/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Emit_h
#define __Emit_h

#include "phy/PhyDef.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace phy {

/** Abstract base class for both single and pair emission distributions. Only exists to allow common base pointers. */
class BaseEmit {
public:
  virtual ~BaseEmit() {}
  simulate() const {errorAbort("Simulate not implemted for this emit class.\n");}
};


/** Abstract base class for both single emission distributions. */
class BaseSingleEmit {
public:
  virtual ~BaseSingleEmit() {}
  virtual emit(unsigned long i) const = 0;
};


/** Abstract base class for both pair emission distributions. */
class BaseSingleEmit {
public:
  virtual ~BasePairEmit() {}
  virtual emit(unsigned long i, unsigned long j) = 0;
};


/** Emit distribution based on precalculated stored emission probabilities. */
public SingleEmitStore {
  SingleEmitStore(vector_t const & probsVec) : probsVec_(probsVec) {}
  virtual ~BaseSingleEmit() {}
  virtual emit(unsigned long i) const {return probsVec[i];}
  reset(vector_t const & probsVec) {probsVec_ = probsVec;}

 protected:
  vector_t probsVec_;
};


/** Emit distribution based on precalculated stored emission probabilities. */
public PairEmitStore {
  PairEmitStore(matrix_t const & probs2DVec) : probs2DVec_(probs2DVec) {}
  virtual ~BasePairEmit() {}
  virtual emit(unsigned long i, unsigned long j) const {return probs2DVec(i, j);}
  reset(matrix_t const & probs2DVec) {probs2DVec_ = probs2DVec;}

 protected:
  matrix_t probs2DVec_;
};



} // end namespace phy

#endif  // __Emit_h
