/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/EmitWrapDfg.h"

namespace phy {


  xvector_t calcLikelihoodVector(SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar)
  {
    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    vector<long> seqSizes = getSeqSizes(seqData);
    long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);

    // make output data structure
    xvector_t result(symCount);

    calcLikelihoodVector(result, seqData, dfgInfo, stvVec, missingDataChar);
    return result;
  }



  void calcLikelihoodVector(xvector_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar)
  {
    probMap_t probMap; // not used!
    calcLikelihoodVector(result, seqData, dfgInfo, stvVec, probMap, false, 1, missingDataChar);
  }


  void calcLikelihoodVector(xvector_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, probMap_t & probMap, bool useProbMap, number_t minMapProb, char missingDataChar)
  {
    // make vector with maps corresponding to non-sub and sub vars respectively
    //TODO: Might refactor
    vector<SeqToVarSymbol> stvVecNonSub;
    vector<SeqToVarSymbol> stvVecSub;
    for(int i = 0; i < stvVec.size(); ++i){
      if(stvVec[i].subscription)
	stvVecSub.push_back(stvVec[i]);
      else
	stvVecNonSub.push_back(stvVec[i]);
    }

    // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames); //only used in the symCount check
    vector<int> stvToSeqMapNonSub = mkStvToSeqMap(stvVecNonSub, seqData.seqNames);
    vector<int> stvToSeqMapSub = mkStvToSeqMap(stvVecSub, seqData.seqNames);
    onlyOptionalSeqMissingOrAbort(stvVecNonSub, stvToSeqMapNonSub, seqData.entryId);
    onlyOptionalSeqMissingOrAbort(stvVecSub, stvToSeqMapSub, seqData.entryId); //Subscription variables should never be optional. This is also enforced at construction

    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<long> seqSizes = getSeqSizes(seqData);
    unsigned long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);
    assert( symCount <= result.size() );

    // setup symbol vector
    //TODO: change to vector<symbol_t>
    vector<string> symVec = initSymVec(stvVecNonSub, false);
    vector<string> symVecSub = initSymVec(stvVecSub, false);

    // map from the observed random variables / stv to total set of random variables 
    vector<unsigned> stvToVarMap = mkStvToVarMap(stvVecNonSub, dfgInfo.varNames);
    vector<unsigned> stvToSubVarInvMap = mkStvToVarInvMap(stvVecSub, dfgInfo.subNames);

    // mk stateMaskVecs, call dfg, and fill in result vector
    typedef probMap_t::iterator mapIter_t;  // map iterator
    stateMaskVec_t stateMaskVec( dfgInfo.stateMaskMapSet.stateMapCount() );

    bool found = false;
    for (unsigned long i = 0; i < symCount; i++) {
      mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVecNonSub, stvToSeqMapNonSub, i, missingDataChar);
      mkSymbolVector(symVecSub, seqData.sequences, seqSizes, stvVecSub, stvToSeqMapSub, i, missingDataChar);
      dfgInfo.stateMaskMapSet.symbols2StateMasks(stateMaskVec, symVec, stvToVarMap);
      dfgInfo.updateFactors(symVecSub, stvToSubVarInvMap);
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
	result[i] = dfgInfo.dfg.calcNormConst(stateMaskVec);
	if (useProbMap)
	  if (result[i] > minMapProb)
	    probMap[stateMaskVec] = result[i];
      }
    }
  }


  xmatrix_t calcLikelihoodMatrix(SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, bool useProbMap, number_t minMapProb, char missingDataChar)
  {
    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    vector<long> seqSizes = getSeqSizes(seqData);
    long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);

    // make output data structure
    xmatrix_t result(symCount, symCount);
    reset(result);

    calcLikelihoodMatrix(result, seqData, dfgInfo, stvVec, missingDataChar);
    return result;
  }

  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar)
  {
    probMap_t probMap; // not used!
    calcLikelihoodMatrix(result, seqData, dfgInfo, stvVec, probMap, false, 1, missingDataChar);
  }

  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, probMap_t & probMap, bool useProbMap, number_t minMapProb, char missingDataChar)
  {
    // make vector with maps corresponding to non-sub and sub vars respectively
    //TODO: Might refactor
    vector<SeqToVarSymbol> stvVecNonSub;
    vector<SeqToVarSymbol> stvVecSub;
    for(int i = 0; i < stvVec.size(); ++i){
      if(stvVec[i].subscription)
	stvVecSub.push_back(stvVec[i]);
      else
	stvVecNonSub.push_back(stvVec[i]);
    }

    // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    vector<int> stvToSeqMapNonSub = mkStvToSeqMap(stvVecNonSub, seqData.seqNames);
    vector<int> stvToSeqMapSub = mkStvToSeqMap(stvVecSub, seqData.seqNames);
    onlyOptionalSeqMissingOrAbort(stvVecNonSub, stvToSeqMapNonSub, seqData.entryId);
    onlyOptionalSeqMissingOrAbort(stvVecSub, stvToSeqMapSub, seqData.entryId); //Subscription variables should never be optional. This is also enforced at construction

    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<long> seqSizes = getSeqSizes(seqData);
    unsigned long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);
    assert( symCount <= result.size1() and symCount <= result.size2() );

    // setup symbol vector
    vector<string> symVec = initSymVec(stvVec, true);
    vector<string> symVecSub = initSymVec(stvVecSub, true);

    // map from the observed randome variables / stv to total set of random variables 
    vector<unsigned> stvToVarMap = mkStvToVarMap(stvVecNonSub, dfgInfo.varNames);
    vector<unsigned> stvToSubVarInvMap = mkStvToVarInvMap(stvVecSub, dfgInfo.subNames);

    // mk stateMaskVecs, call dfg, and fill in result vector
    typedef probMap_t::iterator mapIter_t;  // map iterator
    stateMaskVec_t stateMaskVec( dfgInfo.stateMaskMapSet.stateMapCount() );
    bool found = false;
    for (unsigned long j = 1; j < symCount; j++) 
      for (unsigned long i = 0; i < j; i++) {  // left and right side cannot overlap
	mkDiSymbolVector(symVec, seqData.sequences, seqSizes, stvVecNonSub, stvToSeqMapNonSub, i, j, missingDataChar);
	mkDiSymbolVector(symVecSub, seqData.sequences, seqSizes, stvVecSub, stvToSeqMapSub, i, j, missingDataChar);
	dfgInfo.stateMaskMapSet.symbols2StateMasks(stateMaskVec, symVec, stvToVarMap);
	dfgInfo.updateFactors(symVecSub, stvToSubVarInvMap);
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
	  result(i, j) = dfgInfo.dfg.calcNormConst(stateMaskVec);
	  if (useProbMap)
	    if (result(i, j) > minMapProb)
	      probMap[stateMaskVec] = result(i, j);
	}
      }
  }


} // Namespace phy 
