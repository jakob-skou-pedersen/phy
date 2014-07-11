/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/EmitWrapIO.h"

namespace phy {

  // global tree string (used by overloaded version of readEmitWrappers)
  string globalTreeString = ""; 

  // example:
  // KIND:	phylo
  // NAME:	single
  // TYPE:	single
  // FILE:	sngNucPhyloModel.v3.txt

  void readEmitWrappers(istream & str, EmitWrappers & ew, string const & filePrefix)
  {
    skipWhiteSpaceAndComments(str);
    // modelCount
    unsigned modelCount;
    getFeature(str, "EMIT_MODEL_COUNT:", modelCount);

    // reading the models
    skipWhiteSpaceAndComments(str);
    //while ( str.good() ) {
    for (unsigned i = 0; i < modelCount; i++) {
      string kind, name, type, file, dir;
      unsigned minDist = 1; // define default value
      unsigned tagCount = 0;
      while ( moreTags(str) ) {
	string tag;
	str >> tag;
	if (tag == "KIND:")
	  str >> kind;
	else if (tag == "TYPE:")
	  str >> type;
	else if (tag == "NAME:")
	  str >> name;
	else if (tag == "FILE:")
	  str >> file;
	else if (tag == "DIR:")
	  str >> dir;
	else if (tag == "MIN_DIST:")
	  str >> minDist;
	else
	  errorAbort("Unknown tag ('" + tag + "') in emitWrapper (emitModel) with name  '" + name + "'.");
	skipLine(str); // skip rest of line
	tagCount++;
      }
      if (tagCount < 4) 
	errorAbort("emitWrapper (emitModel) with name  '" + name + "' has wrong number of tags.");
	
      // check type
      if (type != "single" and type != "pair")
	errorAbort("from operator>> for EmitWrapper: emit model with name '" + name +"' has type '" + type + "'. the type must be either 'single' or 'pair'.");

      // mk wrapper
      if (kind == "phylo") { 
	PhyloTree phyloTree = readPhyloTree(filePrefix + file, globalTreeString);
	if (type == "single") {
	  EmitWrappers::VectorEmitWrapPtr_t ptr = boost::shared_ptr<PhyloVectorWrap>(new PhyloVectorWrap(name, phyloTree) );
	  ew.addVectorEmitWrap(ptr);
	}
	else if (type == "pair") {
	  EmitWrappers::MatrixEmitWrapPtr_t ptr = boost::shared_ptr<PhyloMatrixWrap>(new PhyloMatrixWrap(name, minDist, phyloTree) );
	  ew.addMatrixEmitWrap(ptr);
	}
      }

      // mk wrapper
      else if (kind == "dfg") { 
	if (dir == "")
	  errorAbort("from operator>> for EmitWrapper: emit model with name '" + name 
		     + "' is of kind '" + kind + "', and must specify a DIR: (directory) for each emission model, " 
		     + "which should hold the following five files:\n stateMaps.txt factorPotentials.txt variables.txt factorGraph.txt seqToVarSymbol.txt");
	
	DfgInfo dfgInfo = readDfgInfo(filePrefix + dir + "/stateMaps.txt", 
				      filePrefix + dir + "/factorPotentials.txt", 
				      filePrefix + dir + "/variables.txt", 
				      filePrefix + dir + "/factorGraph.txt"); 
	vector<SeqToVarSymbol> stvVec = readSeqToVarSymbols(filePrefix + dir + "/seqToVarSymbol.txt");
	
	if (type == "single") {
	  assertSymbolSizeConsistency(stvVec, dfgInfo.varNames, dfgInfo.stateMaskMapSet, false); // a bit of runtime checking
	  EmitWrappers::VectorEmitWrapPtr_t ptr = boost::shared_ptr<DfgVectorWrap>(new DfgVectorWrap(name, dfgInfo, stvVec) );
	  ew.addVectorEmitWrap(ptr);
	}
	else if (type == "pair") {
	  assertSymbolSizeConsistency(stvVec, dfgInfo.varNames, dfgInfo.stateMaskMapSet, true); // a bit of runtime checking
	  EmitWrappers::MatrixEmitWrapPtr_t ptr = boost::shared_ptr<DfgMatrixWrap>(new DfgMatrixWrap(name, minDist, dfgInfo, stvVec) );
	  ew.addMatrixEmitWrap(ptr);
	}
      }

      else // kind not phylo or dfg
	errorAbort("from operator>> for EmitWrapper: emit model with name '" + name +"' is of kind '" + kind + "', which is currently not implemented");
      skipWhiteSpaceAndComments(str);
    }
  }


  void readEmitWrappers(string const & file, EmitWrappers & ew)
  {
    string filePrefix = pathPrefix(file);
    ifstream f;
    openInFile(f, file);
    readEmitWrappers(f, ew, filePrefix);
  }


  EmitWrappers readEmitWrappers(string const & file)
  {
    EmitWrappers ew;
    readEmitWrappers(file, ew);
    return ew;
  }

  void readEmitWrappers(string const & file, EmitWrappers & ew, string const & treeString)
  {
    globalTreeString = treeString;
    readEmitWrappers(file, ew);
    globalTreeString = "";
  }

  EmitWrappers readEmitWrappers(string const & file, string const & treeString)
  {
    EmitWrappers ew;
    readEmitWrappers(file, ew, treeString);
    return ew;
  }


} // Namespace phy 


