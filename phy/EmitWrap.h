/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __EmitWrap_h
#define __EmitWrap_h

#include "phy/PhyDef.h"
#include "phy/utilsLinAlg.h"
#include "phy/SeqIO.h"
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

namespace phy {


  /** Wrapper class of probility a distribution defined on sequential data. */
  class VectorEmitWrap {
  public:
    VectorEmitWrap(string const & name) : name_(name) {}
    virtual ~VectorEmitWrap() {}
    virtual void emitVector(xvector_t & emitVec, SeqData const & seqData) = 0;
    xvector_t emitVector(SeqData const & seqData);
    virtual void clearHashMap() {errorAbort("clearHashMap is not implemted for vectorEmitWrap.\n");};
    string name() const {return name_;}
  protected:
    string name_;

  };


  /** Wrapper class of probility a distribution defined on sequential data. */
  class MatrixEmitWrap {
  public:
    MatrixEmitWrap(string const & name, unsigned minDist = 1) : name_(name), minDist_(minDist) {}
    virtual ~MatrixEmitWrap() {}
    void emitMatrix(xmatrix_t & emitMat, SeqData const & seqData);
    xmatrix_t emitMatrix(SeqData const & seqData);
    virtual void clearHashMap() {errorAbort("clearHashMap is not implemted for matrixEmitWrap.\n");};
    string name() const {return name_;}
  protected:
    virtual void emitMatrixImpl(xmatrix_t & emitMat, SeqData const & seqData) = 0;
    string name_;
    unsigned minDist_; // minimum distance between left and right part of pair
  };


  /** Class that holds collection of emission wrappers and allows easy generation of vectors of emission probabilities. */
  class EmitWrappers {
  public:
    typedef boost::shared_ptr<VectorEmitWrap> VectorEmitWrapPtr_t;
    typedef boost::shared_ptr<MatrixEmitWrap> MatrixEmitWrapPtr_t;

    /** Constructor */
    EmitWrappers( vector<VectorEmitWrapPtr_t> vv = vector<VectorEmitWrapPtr_t>(), vector<MatrixEmitWrapPtr_t> vm = vector<MatrixEmitWrapPtr_t>() ) : vv_(vv), vm_(vm) {};
    
    virtual ~EmitWrappers() {};

    /** returns vector of vectors of emissions probabilities. */
    void vectorEmissions(vector<xvector_t> & v, SeqData const & seqData);
    vector<xvector_t> vectorEmissions(SeqData const & seqData);

    /** returns vector of matrices of emissions probabilities. */
    void matrixEmissions(vector<xmatrix_t> & v, SeqData const & seqData);
    vector<xmatrix_t> matrixEmissions(SeqData const & seqData);

    /** Count wrappers*/
    unsigned vCount() const {return vv_.size();}
    unsigned mCount() const {return vm_.size();}

    /** return wrapper names */
    vector<string> vecWrapNames() const;
    vector<string> matWrapNames() const;

    /** Clear hashMaps of all constructors */
    virtual void clearHashMap();

    /** Add additional Wrappers */
    void addVectorEmitWrap(VectorEmitWrapPtr_t & vPtr) {vv_.push_back(vPtr);}
    void addMatrixEmitWrap(MatrixEmitWrapPtr_t & mPtr) {vm_.push_back(mPtr);}

  protected:
    
    vector<VectorEmitWrapPtr_t> vv_;
    vector<MatrixEmitWrapPtr_t> vm_;
  };



} // end namespace phy

#endif  // __EmitWrap_h
