/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __TRCTMCTree_h
#define __TRCTMCTree_h

#include "phy/PhyDef.h"
#include "phy/utils.h"
#include "phy/utilsLinAlg.h"
#include "phy/Tree.h"
#include "phy/Observations.h"
#include "phy/RateMatrix.h"
#include "phy/Factors.h"
#include "phy/optimizexx.h"
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp> 
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>

using namespace std;

namespace phy {

  //// Some general rate matrix manipulations ////

  /** Diagonalize a time reversible rate Q as V * diag(lambda) * VInv,
      where lambda is a vector of eigenvalues, V is the corresponding
      matrix of eigenvectors (columns), VInv is the inverse of V
      (transpose since the eigenvectors orthonormal).*/
  void diagonalizeReversibleRM(matrix_t const & Q, matrix_t & V, matrix_t & VInv, vector_t & lambda); 

  /** Calculate the matrix exponential of Q (exp(tQ)) using the diagonalized from above (exp(tQ) = V exp( diag(lamda) ) VInv). .*/
  void initTransitionMatrix(matrix_t const & V, matrix_t const & VInv, vector_t const & lambda, matrix_t & T, number_t t);

  /** if multiple transition matrices are made the above version should be preferred */
  void initTransitionMatrix(matrix_t const & Q, matrix_t & T, number_t t);

  /** Update a vector of transition matrices (TM) given branch lengths (branchLengths) and a rate matrix (Q) */
  void updateTransitionMatrices(matrix_t const & Q, vector_t branchLengths, std::vector< matrix_t > & TM);

  /** Update a vector of transition matrices (TM) given branch lengths (branchLengths) and the diagonalized form of a rate matrix (Q) */
  void updateTransitionMatrices(matrix_t & V, matrix_t & VInv, vector_t const & lambda, vector_t branchLengths, std::vector< matrix_t > & TM);


  /**  TrctmcFactorSet */
  class TrctmcFactorSet : public AbstractBaseFactorSet {
  public:
    
    // The factor set includes a (transition matrix) factor for each branch length and an additional factor for the prior distribution of the root node.
    TrctmcFactorSet(vector_t const & branchLengths, boost::shared_ptr<BaseTRRateMatrix> const & rateMatrixPtr, number_t rateScale);
   
    virtual ~TrctmcFactorSet();
   
    virtual void submitCounts(matrix_t const & counts, unsigned idx); // idx is the factor index within the class
    using AbstractBaseFactorSet::submitCounts; // bringing other definitions into this name space (hidden otherwise)

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters();

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const; 
    virtual void mkFactor(matrix_t & m, unsigned idx) const;

    /**     clear all submitted expectation counts */
    virtual void clearCounts();

    /** objective function used in parameter optimization */
    double objFunction(vector_t const & x);

    /** Return value of rateScale. */
    number_t rateScale() const {return rateScale_;}

  protected:
   
    void initCounts();
    void initTransitionMatrices();

    // update all factors based on rateMatrix and branchlengths (updates both transition matrices and rootPrior).
    void updateFactors();

    // set transition matrices based on rateMarices, rateScale and branchLengths
    void updateTMs();

    // set rootPrior to rateMatrix equiFreqs
    void updateRootPrior();
 
    vector_t encodeNonEquiFreqParameters() const;
    void decodeNonEquiFreqParameters(vector_t const & parameters);

    // estimate equiFreqs
    void estimateEquiFreq(vector_t & equiFreq) const;
    vector_t estimateEquiFreq() const;

    void estimateEquiFreqFromRootPrior(vector_t & equiFreq) const;
    vector_t estimateEquiFreqFromRootPrior() const;

    // returns complete logLikelihood based on substitution counts
    number_t countsLogLikelihood();

    // parameters
    vector_t const branchLengths_;
    boost::shared_ptr<BaseTRRateMatrix> rateMatrixPtr_;
    number_t rateScale_;

    // counts
    vector<matrix_t> countsVec_;

    // factors (transition matrices)
    vector<matrix_t> transitionMatrices_;
    matrix_t rootPrior_;
  };










 //
 //  /** Returns a vector with branch lengths indexed by branch index. */
 //  vector_t extractBranchLengths(tree_t const & phyTree);
 //
 //  /** Returns a factor tree corresponding to the phylogenetic
 //      tree. Also returns the mapping of branches and nodes of the
 //      phylogenetic tree to the nodes of factor tree. Finally returns
 //      the index of the factor node corresponding to the initial state
 //      distribution of the phylogentic tree. (All output is done by reference.) */
 //  void mkFactorGraph(tree_t const & phyTree, 
 //		     FG_t & factorGraph, 
 //		     std::vector<unsigned> & nodeMap, 
 //		     std::vector<unsigned> & branchMap, 
 //		     std::vector<FGVertex_t> & factorFirstVarMap, 
 //		     unsigned & initDistrNode);
//
//  vector<FGBaseNode *> mkFactorGraphNodes(std::vector<unsigned> const & nodeMap, 
//					     FG_t const & factorGraph, 
//					     std::vector<unsigned> const & branchMap, 
//					     std::vector<FGVertex_t> const & factorFirstVarMap, 
//					     unsigned initDistrNode, 
//					     std::vector< matrix_t > const & TM, 
//					     vector_t const & equiFreq);
//
//
//  // for debugging, no need to be in header
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v);
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v, FGVertex_t firstVar);
//
  


} // end namespace phy

#endif  //__TRCTMCTree_h
