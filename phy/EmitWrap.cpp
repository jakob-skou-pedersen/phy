/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/EmitWrap.h"

namespace phy {

  xvector_t VectorEmitWrap::emitVector(SeqData const & seqData)
  {
    xvector_t emitVec( seqData.seqSize() );
    reset(emitVec);
    emitVector(emitVec, seqData);
    return emitVec;
  }


  xmatrix_t MatrixEmitWrap::emitMatrix(SeqData const & seqData)
  {
    unsigned n = seqData.seqSize();
    xmatrix_t emitMat(n, n);
    reset(emitMat);
    emitMatrix(emitMat, seqData);
    return emitMat;
  }


  void MatrixEmitWrap::emitMatrix(xmatrix_t & emitMat, SeqData const & seqData)
  {
    emitMatrixImpl(emitMat, seqData);
    if (minDist_ > 1) {
      unsigned n = emitMat.size1();
      for (unsigned i = 0; i < n; i++) {
	unsigned m = min(i + minDist_, n);
	for (unsigned j = i + 1; j < m; j++)
	  emitMat(i, j) = 0;
      }
    }
  }


  void EmitWrappers::vectorEmissions(vector<xvector_t> & v, SeqData const & seqData)
  {
    assert( v.size() == vv_.size() );
    for (unsigned i = 0; i < vv_.size(); i++)
      vv_[i]->emitVector(v[i], seqData);
  }


  vector<xvector_t> EmitWrappers::vectorEmissions(SeqData const & seqData)
  {
    vector<xvector_t> v(vv_.size(), xvector_t( seqData.seqSize() ) );
    vectorEmissions(v, seqData);
    return v;
  }


  void EmitWrappers::matrixEmissions(vector<xmatrix_t> & v, SeqData const & seqData)
  {
    assert( v.size() == vm_.size() );
    for (unsigned i = 0; i < vm_.size(); i++)
      vm_[i]->emitMatrix(v[i], seqData);
  }


  vector<xmatrix_t> EmitWrappers::matrixEmissions(SeqData const & seqData)
  {
    unsigned n = seqData.seqSize();
    vector<xmatrix_t> v(vm_.size(), xmatrix_t(n, n) );
    BOOST_FOREACH(xmatrix_t & m, v)
      reset(m);
    matrixEmissions(v, seqData);
    return v;
  }


  vector<string> EmitWrappers::vecWrapNames() const
  {
    vector<string> v;
    BOOST_FOREACH(VectorEmitWrapPtr_t const & vPtr, vv_) 
      v.push_back( vPtr->name() );
    return v;
  }


  vector<string> EmitWrappers::matWrapNames() const
  {
    vector<string> v;
    BOOST_FOREACH(MatrixEmitWrapPtr_t const & vPtr, vm_) 
      v.push_back( vPtr->name() );
    return v;
  }


  void EmitWrappers::clearHashMap() {
    BOOST_FOREACH(VectorEmitWrapPtr_t & vPtr, vv_) 
      vPtr->clearHashMap(); 
    BOOST_FOREACH(MatrixEmitWrapPtr_t & mPtr, vm_) 
      mPtr->clearHashMap();
  }

} // namespace phy
