/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __EmitWrapDfg_h
#define __EmitWrapDfg_h

#include "phy/EmitWrap.h"
#include "phy/PhyloTree.h"  // needed for probMap_t
#include "phy/DfgIO.h"


namespace phy {


  /** Return vector with normalization constants (likelihoods) of a symbol vector
      given the DFG and the SeqToVarSymbol definitions. Overloaded
      versions return result by reference and allow only unseen sites
      to be evaluated by use of hashing (a map for now).*/
  xvector_t calcLikelihoodVector(SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar = '.');
  void calcLikelihoodVector(xvector_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar = '.');
  void calcLikelihoodVector(xvector_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, probMap_t & probMap, bool useProbMap = true, number_t minMapProb = 1e-7, char missingDataChar = '.');

  /** Return matrix with normalization constants (likelihoods) of all
      di-symbols columns made from a left and a right part. Overloaded
      versions return result by reference and allow only unseen sites
      to be evaluated by use of hashing.*/
  xmatrix_t calcLikelihoodMatrix(SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, bool useProbMap = true, number_t minMapProb = 1e-7, char missingDataChar = '.');
  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, char missingDataChar = '.');
  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, DfgInfo & dfgInfo, vector<SeqToVarSymbol> const & stvVec, probMap_t & probMap, bool useProbMap = true, number_t minMapProb = 1e-7, char missingDataChar = '.');



  /** Wrapper of generic DFG models. */
  class DfgVectorWrap : public VectorEmitWrap {
  public:
    virtual ~DfgVectorWrap() {}
    DfgVectorWrap(string const & name, DfgInfo const & dfgInfo, vector<SeqToVarSymbol> const & stvVec) : VectorEmitWrap(name), dfgInfo_(dfgInfo), stvVec_(stvVec) {};
    virtual void emitVector(xvector_t & emitVec, SeqData const & seqData) {calcLikelihoodVector(emitVec, seqData, dfgInfo_, stvVec_, probMap_, true, 1e-7);}
    virtual void clearHashMap() {probMap_.clear();}
    
  protected:
    probMap_t probMap_;
    DfgInfo dfgInfo_; // to do: consider using a smart pointer
    vector<SeqToVarSymbol> stvVec_;
  };


  /** Wrapper of generic DFG models for pairwise combined symbols. */
  class DfgMatrixWrap : public MatrixEmitWrap {
  public:
    virtual ~DfgMatrixWrap() {}
    DfgMatrixWrap(string const & name, unsigned minDist, DfgInfo const & dfgInfo, vector<SeqToVarSymbol> const & stvVec) : MatrixEmitWrap(name, minDist), dfgInfo_(dfgInfo), stvVec_(stvVec) {};
    virtual void clearHashMap() {probMap_.clear();}
    
  protected:
    virtual void emitMatrixImpl(xmatrix_t & emitVec, SeqData const & seqData) {calcLikelihoodMatrix(emitVec, seqData, dfgInfo_, stvVec_, probMap_, true, 1e-7);}
    probMap_t probMap_;
    DfgInfo dfgInfo_; // to do: consider using a smart pointer
    vector<SeqToVarSymbol> stvVec_;
  };


} // end namespace phy

#endif  // __EmitWrapDfg_h
