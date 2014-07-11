/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __RateMatrix_h
#define __RateMatrix_h

#include "phy/PhyDef.h"
#include "phy/utilsLinAlg.h"
#include "phy/Observations.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>


////////////////////////////////////////////////////////////////
// To do: Incorporate JK matrix in hierarchy


/** Classes for specifying continuous time rate matrices. */

namespace phy {

  /** return a normalized jukes cantor rate matrix (no free parameters) */
  matrix_t jukesCantorRM(); 
  

  // base class. The name defines the type of the derived rate matrix.
  class BaseRateMatrix {
  public:
    BaseRateMatrix(string const & name = "") : name_(name) {}
    virtual ~BaseRateMatrix();
    // Returns a rateMatrix based on the input parameters. The input parameters are split into rateParameters and eqiuFrequencies for convenience.
    virtual matrix_t mkRateMatrix() const = 0;  
    virtual string const & name() const {return name_;}
  protected:
    string name_; 
  };


  // base class for time reversible rate matrices
  class BaseTRRateMatrix : public BaseRateMatrix {
  public:
    /** Constructor that provides default parameter values */
    BaseTRRateMatrix (unsigned rateParameterCount, unsigned equiFreqCount, string const & name = "", bool equiFreqsFixed = false);

    /** Constructor that takes explicit parameter values */
    BaseTRRateMatrix (vector_t const & rateParameters, vector_t const & equiFreqs, string const & name = "", bool equiFreqsFixed = false);

    virtual ~BaseTRRateMatrix() {};

    // set parameters
    virtual void resetRateParameters(vector_t const & rateParameters);
    virtual void resetEquiFreqs(vector_t const & equiFreqs);

    // get parameters
    vector_t const & rateParameters() const {return rateParameters_;} 
    vector_t const & equiFreqs() const {return equiFreqs_;} 

    // report parameter set sizes. Convenient for optimization.
    unsigned rateParameterCount() const;
    unsigned equiFreqCount() const;

    // returns true if equiFreqs are fixed (i.e., cannot be optimized)
    bool const equiFreqsFixed() {return equiFreqsFixed_;}

  protected:

    virtual vector_t defaultRateParameters(unsigned rateParameterCount) const;
    virtual vector_t defaultEquiFreqs(unsigned equiFreqCount) const;

    unsigned rateParameterCount_;
    unsigned equiFreqCount_;

    vector_t rateParameters_;
    vector_t equiFreqs_;

    bool const equiFreqsFixed_;

  };


  // base class for time reversible rate matrices
  class GTRRateMatrix : public BaseTRRateMatrix {
  public:

    /** Constructors. */
    GTRRateMatrix(StateMap const & staMap);
    GTRRateMatrix(StateMap const & staMap, vector_t const & rateParameters, vector_t const & equiFreqs);
    virtual ~GTRRateMatrix() {};

    virtual matrix_t mkRateMatrix() const;

  protected:
    ublas::symmetric_matrix<number_t, ublas::upper> gtrRateParametersToMatrix() const;

    // data
    StateMap staMap_;
  };


  // free functions with general RM functionality

  /** set each diagonal entry to minus the sum of the row (excluding the diagonal entry). */
  void setDiagonalRM(matrix_t & Q);
  
  /** Scale the matrix to 'subs' expected substitutions per position. */
  void normalizeRM(matrix_t & Q, StateMap const & staMap, float subs);
  /** Scale the matrix to symbolSize expected substitutions per position. */
  void normalizeRM(matrix_t & Q, StateMap const & staMap);

  /** Return vector of equilibrium frequencies for a time reversible rate ublas::matrix Q */
  vector_t deriveEquiFreqForReversibleRM(matrix_t const & Q);

  /** Return vector of equilibrium frequencies and parameters of the GTR substitution model for a time reversible rate ublas::matrix Q */
  void deriveEquiFreqAndGtrParamsForReversibleRM(matrix_t const & Q, vector_t & equiFreq, vector_t & gtrParams);

  /** Construct rate matrix of requested type and return shared (smart) pointer to base class */
  boost::shared_ptr<BaseTRRateMatrix> TRRateMatrixDispatcher(string const & substModel, vector_t const & rateParameters, vector_t const & equiFreq, StateMap const & alphabet);

  /** Class for defining hamming distance between two symbols. */
  class HammingDistance 
  {
  public:
    HammingDistance(StateMap const & staMap);
    unsigned operator()(state_t i, state_t j);
    unsigned operator()(symbol_t s, symbol_t t);

  protected:
    StateMap staMap_;
  };


} // end namespace phy

#endif  // __RateMatrix_h
