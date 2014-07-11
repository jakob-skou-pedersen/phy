/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __EmitWrapPhylo_h
#define __EmitWrapPhylo_h

#include "phy/EmitWrap.h"
#include "phy/PhyloTreeIO.h"

namespace phy {

  /** Wrapper of phylogenetic models. */
  class PhyloVectorWrap : public VectorEmitWrap {
  public:
    virtual ~PhyloVectorWrap() {}
    PhyloVectorWrap(string const & name, PhyloTree const & phyloTree) : VectorEmitWrap(name), phyloTree_(phyloTree) {};
    virtual void emitVector(xvector_t & emitVec, SeqData const & seqData) {calcLikelihoodVector(emitVec, seqData, phyloTree_, probMap_, true, 1e-7);}
    virtual void clearHashMap() {probMap_.clear();}
    
  protected:
    probMap_t probMap_;
    PhyloTree phyloTree_; // to do: consider using a smart pointer
  };


  /** Wrapper of phylogenetic models for pairwise combined symbols. */
  class PhyloMatrixWrap : public MatrixEmitWrap {
  public:
    virtual ~PhyloMatrixWrap() {}
    PhyloMatrixWrap(string const & name, unsigned minDist, PhyloTree const & phyloTree) : MatrixEmitWrap(name, minDist), phyloTree_(phyloTree) {};
    virtual void clearHashMap() {probMap_.clear();}
    
  protected:
    virtual void emitMatrixImpl(xmatrix_t & emitVec, SeqData const & seqData) {calcLikelihoodMatrix(emitVec, seqData, phyloTree_, probMap_, true, 1e-7);}
    probMap_t probMap_;
    PhyloTree phyloTree_; // to do: consider using a smart pointer
  };


} // end namespace phy

#endif  // __EmitWrapPhylo_hh
