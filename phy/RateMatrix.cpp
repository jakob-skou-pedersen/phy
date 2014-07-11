/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <boost/numeric/ublas/io.hpp>
#include "phy/RateMatrix.h"

namespace phy {


  matrix_t jukesCantorRM()
  {
    matrix_t Q(4, 4);
    for (unsigned i = 0; i < Q.size1(); ++ i)
      for (unsigned j = 0; j < Q.size2(); ++ j)
	Q (i, j) = 1;
    setDiagonalRM(Q);
    StateMap nucMap( mkNucStateMap() );
    normalizeRM(Q, nucMap, 1.0);
    return Q;
  }


  BaseRateMatrix::~BaseRateMatrix() {}


  BaseTRRateMatrix::BaseTRRateMatrix (unsigned rateParameterCount, unsigned equiFreqCount, string const & name, bool equiFreqsFixed) :  
    BaseRateMatrix(name), rateParameterCount_(rateParameterCount), equiFreqCount_(equiFreqCount), rateParameters_( defaultRateParameters(rateParameterCount) ), equiFreqs_( defaultEquiFreqs(equiFreqCount) ), equiFreqsFixed_(equiFreqsFixed)
  {}


  BaseTRRateMatrix::BaseTRRateMatrix (vector_t const & rateParameters, vector_t const & equiFreqs, string const & name, bool equiFreqsFixed) 
    : BaseRateMatrix(name), rateParameterCount_(rateParameters.size() ), equiFreqCount_(equiFreqs.size() ), rateParameters_(rateParameters), equiFreqs_(equiFreqs), equiFreqsFixed_(equiFreqsFixed)
  {}


  void BaseTRRateMatrix::resetRateParameters(vector_t const & rateParameters) 
  {
    assert(rateParameters.size() == rateParameterCount_);
    rateParameters_ = rateParameters;
  }


  void BaseTRRateMatrix::resetEquiFreqs(vector_t const & equiFreqs) 
  {
    assert(not equiFreqsFixed() );
    assert(equiFreqs.size() == equiFreqCount_);
    equiFreqs_ = equiFreqs;
  }

    // report parameter set sizes. Convenient for optimization.
  unsigned BaseTRRateMatrix::rateParameterCount() const 
  {
    return rateParameterCount_;
  }


  unsigned BaseTRRateMatrix::equiFreqCount() const 
  {
    return equiFreqCount_;
  }


  vector_t BaseTRRateMatrix::defaultRateParameters(unsigned rateParameterCount) const 
  {
    return ublas::scalar_vector<number_t>(rateParameterCount, 1);
  }


  vector_t BaseTRRateMatrix::defaultEquiFreqs(unsigned equiFreqCount) const 
  {
    return ublas::scalar_vector<number_t>(equiFreqCount, 1.0 / equiFreqCount);
  }


  GTRRateMatrix::GTRRateMatrix(StateMap const & staMap) 
    : BaseTRRateMatrix( (staMap.stateCount() * ( staMap.stateCount() - 1) ) / 2, staMap.stateCount(), "GTR"), staMap_(staMap) 
  {}


  GTRRateMatrix::GTRRateMatrix(StateMap const & staMap, vector_t const & rateParameters, vector_t const & equiFreqs) 
    : BaseTRRateMatrix(rateParameters, equiFreqs, "GTR"), staMap_(staMap) 
    { 
      assert( (staMap.stateCount() * ( staMap.stateCount() - 1) ) / 2 == rateParameterCount_);
      assert(staMap.stateCount() == equiFreqCount_);
    }


  ublas::symmetric_matrix<number_t, ublas::upper> GTRRateMatrix::gtrRateParametersToMatrix() const
  {
    unsigned n = equiFreqCount_;
    assert( n > 0 and rateParameters_.size() == (n * (n - 1) / 2) );

    ublas::symmetric_matrix<number_t, ublas::upper> q(n, n);
    reset(q); // clear
    unsigned idx = 0;
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = i + 1; j < n; j++) {
	q(i, j) = rateParameters_[idx];
	idx++;
      }
    return q;
  }

  /** Returns general time reversible (GTR) rate matrix based on rateParameters and equiFreqs.
   *  The rateParameters specify the q_ij (= q_ji) parameters.
   *  The entries of the rate matrix are:
   *  m_ij = pi_j * q_ij, where pi is the equiFreq. 
   */
  matrix_t GTRRateMatrix::mkRateMatrix() const
  {
    matrix_t q = gtrRateParametersToMatrix();
    unsigned n = staMap_.stateCount();
    matrix_t Q(n, n);
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < n; j++)
	Q(i, j) = equiFreqs_[j] * q(i, j);
    setDiagonalRM(Q);
    normalizeRM(Q, staMap_);
    return Q;
  }

  // free functions with general RM functionality


  void setDiagonalRM(matrix_t & Q)
  {
    for (unsigned i = 0; i < Q.size1(); ++ i)
      Q(i, i) = 0;
    for (unsigned i = 0; i < Q.size1(); ++ i)
      Q(i, i) = - sum( row (Q, i) );
  }


  void normalizeRM(matrix_t & Q, StateMap const & staMap, float subs)
  {
    vector_t equiFreq = deriveEquiFreqForReversibleRM(Q);
    HammingDistance hammingDistance(staMap);
    number_t normConst = 0;
    for (unsigned i = 0; i < Q.size1(); ++ i)
      for (unsigned j = 0; j < Q.size2(); ++ j)
	normConst += equiFreq(i) * Q(i,j) * hammingDistance(i,j);

    Q = Q * (subs / normConst);
  }


  void normalizeRM(matrix_t & Q, StateMap const & staMap)
  {
    normalizeRM(Q, staMap, staMap.symbolSize() );
  }


  unsigned nonZeroDiagonalEntryOrDie(matrix_t const & Q)
  {
    assert( Q.size1() == Q.size1() ); // Q is quadratic

    for (unsigned i = 0; i < Q.size1(); i++)
      if (Q(i, i) != 0.0)
	return i;

    errorAbort("All zero diagonal entries. Bails out.");
    return 0; // will never reach this, only to quiet compiler
  }


  unsigned countNonZeroDiagonalEntries(matrix_t const & Q)
  {
    assert( Q.size1() == Q.size1() ); // Q is quadratic

    unsigned n = 0;
    for (unsigned i = 0; i < Q.size1(); i++)
      if (Q(i, i) != 0.0)
	n++;
    return n;
  }

  unsigned countNonZeroEntries(vector_t const & v)
  {

    unsigned n = 0;
    for (unsigned i = 0; i < v.size(); i++)
      if ( v(i) )
	n++;
    return n;
  }

  void deriveEquiFreqs(matrix_t const & Q, vector_t & equiFreq)
  {
    unsigned size = Q.size1();

    for (unsigned int j = 0; j < size; j++)
      if (equiFreq(j) == 0 && Q(j, j) != 0.0 )
	for (unsigned int i = 0; i < size; i++) 
	  if (equiFreq(i) != 0.0 && Q(j, i) != 0.0 )
	    equiFreq(j) = equiFreq(i) * Q(i, j) / Q(j, i);
  }


  // Q is a quadratic rate ublas::matrix fulfilling detailed balance. 
  vector_t deriveEquiFreqForReversibleRM(matrix_t const & Q)
  {
    assert( Q.size1() == Q.size1() ); // Q is quadratic
    unsigned int size = Q.size1();

    // define refState
    vector_t equiFreq(size);
    equiFreq.clear();
    unsigned refState = nonZeroDiagonalEntryOrDie(Q);
    equiFreq(refState) = 1;

    // calculate un-normalized equi-freqs 
    // This code could be cleaned up given conditions on the form of the rate ublas::matrix
    unsigned remaining = countNonZeroDiagonalEntries(Q) - countNonZeroEntries(equiFreq);
    while ( remaining ) {
      deriveEquiFreqs(Q, equiFreq);

      // any nondetermined?
      unsigned nonZero = countNonZeroEntries(equiFreq);
      if (nonZero == remaining)
	errorAbort("deriveEquiFreqFromReversibleRM: Bad rateUblas::Matrix (non-reversible?), can't derive equiFreq. Caught in loop");

      remaining = countNonZeroDiagonalEntries(Q) - nonZero;
    }

    equiFreq *= 1 / sum(equiFreq);
    return equiFreq;
  }


  void deriveEquiFreqAndGtrParamsForReversibleRM(matrix_t const & Q, vector_t & equiFreq, vector_t & gtrParams)
  {
    equiFreq = deriveEquiFreqForReversibleRM(Q);
    unsigned n = equiFreq.size();
   
    // check
    for (unsigned i = 0; i < n; i++)
      assert(equiFreq[i] > 0);
    
    matrix_t pm = Q;  // parameter matrix
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < n; j++)
	pm(i, j) = Q(i, j) / equiFreq[j];

    // check grtParams size
    unsigned paramCount = (n * (n - 1) / 2);
    if (gtrParams.size() != paramCount)
      gtrParams.resize(paramCount);

    // reverse of gtrRateParametersToMatrix
    unsigned idx = 0;
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = i + 1; j < n; j++) {
	gtrParams[idx] = pm(i, j);
	idx++;
      }
  }


  boost::shared_ptr<BaseTRRateMatrix> TRRateMatrixDispatcher(string const & substModel, vector_t const & rateParameters, vector_t const & equiFreq, StateMap const & alphabet)
  {
    boost::shared_ptr<BaseTRRateMatrix> rmPtr;

    if (substModel == "GTR")
      rmPtr = boost::shared_ptr<GTRRateMatrix>( new GTRRateMatrix(alphabet) );
    // define other substitution models here ...
    else
      errorAbort("substitution model '" + substModel + "' not defined. If newly implemented, add to TRRateMatrixDispatcher.");


    if (equiFreq.size() > 0) // not defined, stick with default
      rmPtr->resetEquiFreqs(equiFreq);
    if (rateParameters.size() > 0) // not defined, stick with default
      rmPtr->resetRateParameters(rateParameters);

    return rmPtr;
  }



//  matrix_t GTRRateMatrix::mkRateMatrix() const
//  {
//    matrix_t q = gtrRateParametersToMatrix();
//    unsigned n = staMap_.stateCount();
//    matrix_t Q(n, n);
//    for (unsigned i = 0; i < n; i++)
//      for (unsigned j = 0; j < n; j++)
//	Q(i, j) = equiFreqs_[j] * q(i, j);
//    setDiagonalRM(Q);
//    normalizeRM(Q, staMap_);
//    return Q;
//  }
//
//  ublas::symmetric_matrix<number_t, ublas::upper> GTRRateMatrix::gtrRateParametersToMatrix() const
//  {
//    unsigned n = equiFreqCount_;
//    assert( n > 0 and rateParameters_.size() == (n * (n - 1) / 2) );
//
//    ublas::symmetric_matrix<number_t, ublas::upper> q(n, n);
//    q *= 0; // clear
//    unsigned idx = 0;
//    for (unsigned i = 0; i < n; i++)
//      for (unsigned j = i + 1; j < n; j++) {
//	q(i, j) = rateParameters_[idx];
//	idx++;
//      }
//    return q;
//  }


  // HammingDistance class
  HammingDistance::HammingDistance(StateMap const & staMap) : staMap_(staMap) {};

  unsigned HammingDistance::operator()(state_t i, state_t j)
  {
    return hammingDistance(staMap_.state2Symbol(i), staMap_.state2Symbol(j));
  }


  unsigned HammingDistance::operator()(symbol_t s, symbol_t t)
  {
    return hammingDistance(s, t);  // from utils.h
  }


} // end namespace phy

