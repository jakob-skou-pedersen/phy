/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/PhyloTreeIO.h"

namespace phy {


  PhyloTree readPhyloTree(string const file, string const & treeString)
  {
    ifstream f;
    openInFile(f, file);
    return readPhyloTree(f, treeString);
  }


  PhyloTree readPhyloTree(istream & str, string const & treeStringArg)
  {
    // declare empty input values
    string substModel;  // rate matrix
    vector_t equiFreqs;
    vector_t rateParameters;
    matrix_t rawRateMatrix;
    number_t rateScale = -1;  
    string treeString;

    // first read alphabet
    StateMap alphabet("");
    str >> alphabet;

    // then the rest
    skipWhiteSpaceAndComments(str);
    while ( moreTags(str) ) {  // end of model definition defined by line starting with whitespace
    string tag;
    str >> tag;

    if (tag == "SUBST_MODEL:")
      str >> substModel;

    else if (tag == "BACKGROUND:")
      str >> equiFreqs;

    else if (tag == "RATE_PARAMETERS:")
      str >> rateParameters;

    else if (tag == "RATE_MATRIX:")
      str >> rawRateMatrix;

    else if (tag == "RATE_SCALE:")
      str >> rateScale;

    else if (tag == "TREE:")
      treeString = readNewickStr(str);

    else if (tag == "TREE_FILE:") {
      string fileName;
      str >> fileName;
      treeString = readNewickStr(fileName);
    }

    else
      errorAbort("Unknown tag ('" + tag + "') in phylogenetic model definition.");
    skipLine(str); // skip rest of line
    }

    if (not treeStringArg.empty() ) {
      if ( (not treeString.empty() ) )
	errorAbort(string("from readPhyloTree: Tree specified both as argument and in configuration file. ") + 
		   "Each tree must be specified in only one place. " + 
		   "You probably want to remove the TREE or TREE_FILE tag in your configuration file. " + 
		   "The trees are gives below (tree from configuration file given first):\n\n" + 
		   treeString +"\n\n" + 
		   treeStringArg + "\n\n");
      else 
	treeString = treeStringArg;
    }
    // use rawRateMatrix to derive rateParameters and equiFreqs for GTR models if needed
    if (substModel == "GTR")
      if (rateParameters.size() == 0)
        deriveEquiFreqAndGtrParamsForReversibleRM(rawRateMatrix, equiFreqs, rateParameters);

    // sanity checks
    if (rateScale < 0)
      errorAbort("RATE_SCALE must be defined and have a positive value (>= 0)");
    if ( abs(sum(equiFreqs) - 1) > 0.001 )
      errorAbort("For phyloTree with subtModel of type '" + substModel + "', sum of equiFreqs (" + toString( sum(equiFreqs) ) + ") deviate more than 0.01 from 1.0");

    boost::shared_ptr<BaseTRRateMatrix> rateMatrixPointer = TRRateMatrixDispatcher(substModel, rateParameters, equiFreqs, alphabet);
    TreeInfo treeInfo( newickToTreeInfo(treeString) );
    TrctmcFactorSet factorSet = TrctmcFactorSet(toNumVector(treeInfo.branchMap), rateMatrixPointer, rateScale);
    DFG dfg = mkPhyloDfg(treeInfo, factorSet);
    StateMaskMapSet stateMaskMapSet( alphabet, dfg.variables.size() );

    // debug
    //    writePhyloTree("tmp_phyloTree_" + alphabet.name() + "_" + substModel + ".txt", PhyloTree(alphabet, rateMatrixPointer, treeInfo, stateMaskMapSet, factorSet, dfg) );

    return PhyloTree(alphabet, rateMatrixPointer, treeInfo, stateMaskMapSet, factorSet, dfg);
  }


  void writePhyloTree(string const file, PhyloTree const & phyloTree)
  {
    ofstream f;
    f.open( file.c_str() );
    writePhyloTree(f, phyloTree);
  }


  void writePhyloTree(ostream & str, PhyloTree const & phyloTree)
  {
    str << phyloTree.alphabet;
    str << "SUBST_MODEL:\t" << phyloTree.rmPtr->name() << endl;
    str << "BACKGROUND:\t" << phyloTree.rmPtr->equiFreqs() << endl;
    str << "RATE_PARAMETERS:\t" << phyloTree.rmPtr->rateParameters() << endl;
    str << "RATE_MATRIX:\t" << phyloTree.rmPtr->mkRateMatrix() << endl;
    str << "RATE_SCALE:\t" << phyloTree.factorSet.rateScale() << endl;
    str << "TREE:\t" << writeNewick(phyloTree.treeInfo) << endl;
    str << endl;
  }


} // end namespace phy

