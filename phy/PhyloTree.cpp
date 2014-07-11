/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/PhyloTree.h"

namespace phy {

  TrctmcFactorSet mkTrctmcFactorSet(TreeInfo const & treeInfo, boost::shared_ptr<BaseTRRateMatrix> rateMatrixPtr, number_t rateScale)
  {
    vector_t branchLengths = toNumVector(treeInfo.branchMap);
    return TrctmcFactorSet(branchLengths, rateMatrixPtr, rateScale);
  }


  DFG mkPhyloDfg(TreeInfo const & treeInfo, AbstractBaseFactorSet const & factorSet)
  {
    assert( treeInfo.branchMap.size() == factorSet.facCount() - 1); // rootprior must be included in factor set and not in branchMap
    vector<matrix_t> factorVec = factorSet.mkFactorVec();
    unsigned stateCount = factorVec[0].size2(); // exploiting that all factors will have the same second dimension

    // check that all factor potentials are of same size
    for (unsigned i = 0; i < factorVec.size(); i++) {
      assert(factorVec[i].size1() == 1 or factorVec[i].size1() == stateCount);
      assert(factorVec[i].size2() == stateCount);
    }

    // define varDimensions
    unsigned varCount = treeInfo.nodeMap.size();
    vector<unsigned> varDimensions(varCount, stateCount); 

    // define factor neighbors and add neighbors of the rootPrior factor
    // (not defined in phylo tree and therefore not in TreeInfo)
    vector<vector<unsigned> > facNeighbors = treeInfo.branchNeighbors;
    vector<unsigned> rootPriorFacNeighbors(1, treeInfo.root);
    facNeighbors.push_back(rootPriorFacNeighbors); 

    return DFG(varDimensions, factorSet.mkFactorVec(), facNeighbors);
  }


  // Convert input observation symbols to state masks. Overloaded version of function defined in Observations.h. offset defined to symbol length.
  stateMask2DVec_t const mkStateMask2DVec(SeqData const & seqData, PhyloTree const & phyloTree, char missingDataChar)
  {
    return mkStateMask2DVec(seqData, phyloTree, phyloTree.alphabet.symbolSize(), missingDataChar);
  }


  stateMask2DVec_t const mkStateMask2DVec(SeqData const & seqData, PhyloTree const & phyloTree, unsigned offSet, char missingDataChar)
  {
    return mkStateMask2DVec(seqData, phyloTree.treeInfo.nodeMap, phyloTree.stateMaskMapSet, offSet, missingDataChar);
  }


  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, PhyloTree const & phyloTree, char missingDataChar)
  {
    return mkStateMask3DVec(seqDataVec, phyloTree, phyloTree.alphabet.symbolSize(), missingDataChar);
  }


  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, PhyloTree const & phyloTree, unsigned offSet, char missingDataChar)
  {
    return mkStateMask3DVec(seqDataVec, phyloTree.treeInfo.nodeMap, phyloTree.stateMaskMapSet, offSet, missingDataChar);
  }


  xvector_t calcLikelihoodVector(SeqData const & seqdata, PhyloTree & phyloTree, unsigned offSet, char missingDataChar)
  {
    xvector_t result( symbolCount(seqdata.sequences, offSet) );
    calcLikelihoodVector(result, seqdata, phyloTree, offSet, missingDataChar);
    return result;
  }


  void calcLikelihoodVector(xvector_t & result, SeqData const & seqdata, PhyloTree & phyloTree, unsigned offSet, char missingDataChar)
  {
    probMap_t probMap; // not used!
    calcLikelihoodVector(result, seqdata, phyloTree, probMap, false, 1, offSet, missingDataChar);
  }

  // to do: refactor this code into generic functions defined in DiscreteFactorGraph.h 
  void calcLikelihoodVector(xvector_t & result, SeqData const & seqData, PhyloTree & phyloTree, probMap_t & probMap, bool useProbMap, number_t minMapProb, unsigned offSet, char missingDataChar)
  {
    // seqSize and symbolCount
    vector<string> const & strVec = seqData.sequences;
    long  seqSize = seqData.seqSize();
    unsigned long symCount = symbolCount(strVec, offSet);
    assert( symCount <= result.size() );

    // setup symbol vector
    unsigned commonSymSize = phyloTree.alphabet.symbolSize();
    vector<symbol_t> symVec(strVec.size(), symbol_t(commonSymSize, missingDataChar) );

    // map of sequences to variables
    vector<unsigned> seqToVarMap = mkSeqToVarMap(phyloTree.treeInfo.nodeMap, seqData.seqNames);

    // mk stateMaskVecs, call dfg, and fill in result vector
    typedef probMap_t::iterator mapIter_t;  // map iterator
    stateMaskVec_t stateMaskVec( phyloTree.stateMaskMapSet.stateMapCount() );
    bool found = false;
    for (unsigned i = 0; i < symCount; i++) {
      mkSymbolVector(symVec, commonSymSize, strVec, seqSize, i * offSet, missingDataChar);
      phyloTree.stateMaskMapSet.symbols2StateMasks(stateMaskVec, symVec, seqToVarMap);

      if (useProbMap) {
	mapIter_t iter = probMap.find(stateMaskVec);
	if ( iter != probMap.end() ) {
	  result[i] = iter->second;
	  found = true;
	}
	else
	  found = false;
      }	

      if (not found) {
	result[i] = phyloTree.dfg.calcNormConst(stateMaskVec);
	if (useProbMap)
	  if (result[i] > minMapProb)
	    probMap[stateMaskVec] = result[i];
      }
    }
  }



  xmatrix_t calcLikelihoodMatrix(SeqData const & seqData, PhyloTree & phyloTree, unsigned leftSize, unsigned rightSize, char missingDataChar)
  {
    unsigned seqSize = seqData.seqSize();
    xmatrix_t result(seqSize, seqSize);
    reset(result);
    calcLikelihoodMatrix(result, seqData, phyloTree, leftSize, rightSize, missingDataChar);
    return result;
  }


  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, PhyloTree & phyloTree, unsigned leftSize, unsigned rightSize, char missingDataChar)
  {
    probMap_t probMap; // not used!
    calcLikelihoodMatrix(result, seqData, phyloTree, probMap, false, 1, leftSize, rightSize, missingDataChar);
  }


  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, PhyloTree & phyloTree, probMap_t & probMap, bool useProbMap, number_t minMapProb, unsigned leftSize, unsigned rightSize, char missingDataChar)
  {
    // seqSize
    vector<string> const & strVec = seqData.sequences;
    unsigned seqSize = seqData.seqSize();
    assert( seqSize <= result.size1() and seqSize <= result.size2() );

    // setup symbol vector
    unsigned commonSymSize = phyloTree.alphabet.symbolSize();
    assert(commonSymSize == rightSize + leftSize);

    vector<symbol_t> symVec(strVec.size(), symbol_t(commonSymSize, missingDataChar) );

    // map of sequences to variables
    vector<unsigned> seqToVarMap = mkSeqToVarMap(phyloTree.treeInfo.nodeMap, seqData.seqNames);

    // mk stateMaskVecs, call dfg, and fille in result vector
    typedef probMap_t::iterator mapIter_t;  // map iterator
    stateMaskVec_t stateMaskVec( phyloTree.stateMaskMapSet.stateMapCount() );
    bool found = false;
    for (unsigned j = leftSize; j < seqSize; j++) 
      for (unsigned i = 0; i <= j - leftSize; i++) {
	mkDiSymbolVector(symVec, leftSize, rightSize, strVec, seqSize, i, j, missingDataChar); 
	phyloTree.stateMaskMapSet.symbols2StateMasks(stateMaskVec, symVec, seqToVarMap);

      if (useProbMap) {
	mapIter_t iter = probMap.find(stateMaskVec);
	if ( iter != probMap.end() ) {
	  result(i, j) = iter->second;
	  found = true;
	}
	else
	  found = false;
      }	

      if (not found) {
	result(i, j) = phyloTree.dfg.calcNormConst(stateMaskVec);
	if (useProbMap)
	  if (result(i, j) > minMapProb)
	    probMap[stateMaskVec] = result(i, j);
      }
    }
  }

  //  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqdata, PhyloTree const & phyloTree, char missingDataChar)



} // end namespace phy

