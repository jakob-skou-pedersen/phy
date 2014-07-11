/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/TRCTMCTree.h"

// DEBUG
#include <boost/numeric/ublas/io.hpp>
#include <iostream>


namespace phy {

  void initPiMatrices(matrix_t const & Q, ublas::banded_matrix<number_t> & PISqrt, ublas::banded_matrix<number_t> & PISqrtInv)
  {
    vector_t pi = deriveEquiFreqForReversibleRM(Q);
    for (unsigned i = 0; i < pi.size(); i++)
      PISqrt(i, i) = sqrt( pi(i) );
    for (unsigned i = 0; i < pi.size(); i++)
      PISqrtInv(i, i) = 1 / PISqrt(i, i);
  }

  // If Q is revesible (i.e. it fulfills detailed balance: PI * Q=
  // PI * Q^t), then Q is similar to a symmetric ublas::matrix S: S =
  // PI_Sqrt * Q * PI_minusSqrt. (all the PI matrices are diagonal
  // in the above.) Q thus shares eigenvalues with the symmetrix
  // ublas::matrix S, and eigenvectors (VQ) is given by: VS * PISqrt).
  void diagonalizeReversibleRM(matrix_t const & Q, matrix_t & V, matrix_t & VInv, vector_t & lambda)
  {
    assert( Q.size1() == Q.size1() ); // Q is quadratic
    
    ublas::banded_matrix<number_t> PISqrt(Q.size1(), Q.size2(), 0, 0);
    ublas::banded_matrix<number_t> PISqrtInv(Q.size1(), Q.size2(), 0, 0);
    initPiMatrices(Q, PISqrt, PISqrtInv);
    ublas::matrix<double, ublas::column_major> S( prod( PISqrt, Q) );  // syev requires a colum_major ublas::matrix. May be able to replace double with number_t
    S = prod(S, PISqrtInv);  // in total: S = PISqrt * Q * PISqrtInv
    boost::numeric::bindings::lapack::syev( 'V', 'L', S, lambda, boost::numeric::bindings::lapack::minimal_workspace() );
    V = prod(PISqrtInv, S);
    VInv = prod(trans(S), PISqrt);
  }


  void initTransitionMatrix(matrix_t const & V, matrix_t const & VInv, vector_t const & lambda, matrix_t & T, number_t t)
  {
    T.clear();
    for (unsigned i = 0; i < lambda.size(); i++)
      T(i, i) = exp(t * lambda(i));
    T = prod(V, T);
    T = prod(T, VInv);  // T = V * exp(t * D_lamda) * VInv
  }


  void initTransitionMatrix(matrix_t const & Q, matrix_t & T, number_t t)
  {
    assert( Q.size1() == Q.size2() );
 
    unsigned dim = Q.size1();
    matrix_t V(dim, dim);
    matrix_t VInv(dim,dim);
    vector_t lambda(dim);

    diagonalizeReversibleRM(Q, V, VInv, lambda);
    initTransitionMatrix (V, VInv, lambda, T, t);
  }


  void updateTransitionMatrices(matrix_t & V, matrix_t & VInv, vector_t const & lambda, vector_t branchLengths, std::vector< matrix_t > & TM) 
  {
    unsigned n = branchLengths.size();
    assert ( TM.size() == n); // branchLengths and TM of same size
    for (unsigned i = 0; i < n; i++)
      initTransitionMatrix(V, VInv, lambda, TM[i], branchLengths[i]);
  }


  void updateTransitionMatrices(matrix_t const & Q, vector_t branchLengths, std::vector< matrix_t > & TM) 
  {
    unsigned dim = Q.size1();
    matrix_t V(dim, dim);
    matrix_t VInv(dim,dim);
    vector_t lambda(dim);

    diagonalizeReversibleRM(Q, V, VInv, lambda);
    updateTransitionMatrices(V, VInv, lambda, branchLengths, TM);
  }


  // implementation of TrctmcFactorSet

  // The factor set includes a (transition matrix) factor for each branch length and an additional factor for the prior distribution of the root node.
  TrctmcFactorSet::TrctmcFactorSet(vector_t const & branchLengths, boost::shared_ptr<BaseTRRateMatrix> const & rateMatrixPtr, number_t rateScale) 
    : AbstractBaseFactorSet(branchLengths.size() + 1), 
      branchLengths_(branchLengths), 
      rateMatrixPtr_(rateMatrixPtr), 
      rateScale_(rateScale), 
      countsVec_(facCount_), 
      transitionMatrices_( branchLengths.size() ), 
      rootPrior_(1, rateMatrixPtr->equiFreqCount())
  {
    initCounts();
    initTransitionMatrices();
    updateFactors();
  }


  TrctmcFactorSet::~TrctmcFactorSet() {}

  void TrctmcFactorSet::submitCounts(matrix_t const & counts, unsigned idx) 
  {
    assert(counts.size1() == rateMatrixPtr_->equiFreqCount() or counts.size1() == 1); 
    assert(counts.size2() == rateMatrixPtr_->equiFreqCount() ); 
    countsVec_[idx] += counts;
  }


  /** return value may optionally report on success or other aspects of optimization. */
  int TrctmcFactorSet::optimizeParameters()
  {
    // GenericFunctionObject<TrctmcFactorSet> f_obj(*this, & TrctmcFactorSet::transformingObjectiveFunction);
    //    vector_t parVec = rateMatrixPtr_->rateParameters();  // debug: inc rateScale
    vector_t parVec = encodeNonEquiFreqParameters();
    transformVec(parVec, 0);

    boost::function< double (vector_t const &) > objFct;
    objFct = boost::bind(& TrctmcFactorSet::objFunction, this, _1);

    parVec = minimizexx(parVec, objFct, "output/optOut.txt");
    deTransformVec(parVec, 0);

    //    double MIN_VALUE = 1e-10, MAX_VALUE = 1e10;
    //    x = forceInRange(x, MIN_VALUE, MAX_VALUE);

    // optimizing equiFreqs
    if ( not rateMatrixPtr_->equiFreqsFixed() ) {
      vector_t equiFreqs = estimateEquiFreq();
      rateMatrixPtr_->resetEquiFreqs(equiFreqs);
      updateFactors();
    }

    cout << "Optimization done" << endl << endl;
    cout << "Final parameter values: " << endl;
    cout << "rateParameters: " << rateMatrixPtr_->rateParameters() << endl;
    cout << "rateScale:      " << rateScale_ << endl;
    cout << "equiFreqs:      " << rateMatrixPtr_->equiFreqs() << endl;
    cout << endl;

    return 0;
  }

  // function called by the numerical optimization procedure
  double TrctmcFactorSet::objFunction(vector_t const & x)
  {
    // to do: assert parameter count
    vector_t parVec = x;
    deTransformVec(parVec, 0);

    // correct extreme values for numerical stability
    //double MIN_VALUE = 1e-10, MAX_VALUE = 1e10;
    //deTrans_x = forceInRange(deTrans_x, MIN_VALUE, MAX_VALUE);

    decodeNonEquiFreqParameters(parVec);
    //    rateMatrixPtr_->resetRateParameters(parVec);
    updateFactors();
    return -countsLogLikelihood();
  }


   /** return or set matrix defining factor idx  */
   matrix_t TrctmcFactorSet::mkFactor(unsigned idx) const
   {
     assert(idx < facCount_);
     if (idx < facCount_ - 1)
       return transitionMatrices_[idx];
     else
       return rootPrior_;
   }


   void TrctmcFactorSet::mkFactor(matrix_t & m, unsigned idx) const 
   {
     m = mkFactor(idx);
   }


   /**     clear all submitted expectation counts */
   void TrctmcFactorSet::clearCounts() 
   {
     for (unsigned i = 0; i < countsVec_.size(); i++) 
       reset(countsVec_[i]);
   }


   void TrctmcFactorSet::initCounts()
   {
     unsigned n = rateMatrixPtr_->equiFreqCount();
     for (unsigned i = 0; i < facCount_ - 1; i++) {
       countsVec_[i].resize(n, n);
       reset(countsVec_[i]);
     }
     (* --countsVec_.end() ).resize(1, n); // root prior distr
     reset(* --countsVec_.end() );
   }
 

  void TrctmcFactorSet::initTransitionMatrices()
  {
    unsigned const n = rateMatrixPtr_->equiFreqCount();
    for (unsigned i = 0; i < transitionMatrices_.size(); i++) {
      transitionMatrices_[i].resize(n, n);
      reset( transitionMatrices_[i] );
    }
  }


  void TrctmcFactorSet::updateFactors()
  {
    updateTMs();
    updateRootPrior();
  }


  void TrctmcFactorSet::updateTMs()
  {
    matrix_t scaledRateMatrix(rateMatrixPtr_->mkRateMatrix() * rateScale_);
    updateTransitionMatrices(scaledRateMatrix, branchLengths_, transitionMatrices_);
  }


  void TrctmcFactorSet::updateRootPrior()
  {
    assert(rootPrior_.size2() == rateMatrixPtr_->equiFreqCount());
    for (unsigned i = 0; i < rootPrior_.size2(); i++) 
      rootPrior_(0, i) = rateMatrixPtr_->equiFreqs()[i];
  }


  // return vector with rateParameters and rateScale
  vector_t TrctmcFactorSet::encodeNonEquiFreqParameters() const
  {
    vector_t parameters = rateMatrixPtr_->rateParameters();
    parameters.resize(parameters.size() + 1);
    ( * --parameters.end() ) = rateScale_;
    //    parameters[ parameters.size() ] = rateScale_; // is this correct? debug
    return parameters;
  }

  
  // extract rateParameters and rateScale from parameters
  void TrctmcFactorSet::decodeNonEquiFreqParameters(vector_t const & parameters)
  {
    vector_t rateParameters = parameters;
    rateParameters.resize(rateParameters.size() - 1);
    rateMatrixPtr_->resetRateParameters(rateParameters);
    rateScale_ = * ( --parameters.end() );
  }


  // the estimation is an approximation based on the assumption that
  // at most one substitution has happened on each branch.  
  //
  // CAUTION:
  // the estimate is blind to the presence of missing data, often
  // better to use estimateEquiFreqFromRootPrior ... I think...

  void TrctmcFactorSet::estimateEquiFreq(vector_t & equiFreq) const
  {
    // Summing state usage over all branches, weighted by branch length 
    unsigned n = rateMatrixPtr_->equiFreqCount();
    matrix_t stateUsage(n, n);
    reset(stateUsage);
    for (unsigned i = 0; i < branchLengths_.size(); i++)
      stateUsage += countsVec_[i] * branchLengths_[i];
    reset(equiFreq);
    for (unsigned k = 0; k < n; k++) 
      for (unsigned l = 0; l < n; l++) {
	equiFreq[k] += stateUsage(k, l); // should be weighted by 0.5, but freqs will be normalized anyway, so no need
	equiFreq[l] += stateUsage(k, l); // should be weighted by 0.5, but freqs will be normalized anyway, so no need
      }
    equiFreq *= ( 1.0 / sum(equiFreq) );
  }


  vector_t TrctmcFactorSet::estimateEquiFreq() const
  {
    unsigned n = rateMatrixPtr_->equiFreqCount();
    vector_t equiFreq(n);
    estimateEquiFreq(equiFreq);
    return equiFreq;
  }


  void TrctmcFactorSet::estimateEquiFreqFromRootPrior(vector_t & equiFreq) const
  {
    // sanity checks
    matrix_t const & c = *( --countsVec_.end() );
    unsigned n = rateMatrixPtr_->equiFreqCount();
    assert(c.size1() == n);
    assert(c.size2() == 1);

    equiFreq *= 0;
    for (unsigned i = 0; i < n; i++) 
      equiFreq[i] = c(0, i);
    
    equiFreq *= 1.0 / sum(equiFreq);
  }


  vector_t TrctmcFactorSet::estimateEquiFreqFromRootPrior() const
  {
    // sanity checks
    unsigned n = rateMatrixPtr_->equiFreqCount();
    vector_t equiFreq(n);
    estimateEquiFreqFromRootPrior(equiFreq);
    return equiFreq;
  }


  // The complete log likelihood as calculated below is derived in:
  //
  // A. Siepel and D. Haussler Phylogenetic estimation of
  // context-dependent substitution rates by maximum
  // likelihood. Mol. Biol. Evol. (2004) vol. 21 (3) pp. 468-488
  number_t TrctmcFactorSet::countsLogLikelihood()
  {
    number_t logLik = 0;
    unsigned n = rateMatrixPtr_->equiFreqCount();
    for (unsigned i = 0; i < transitionMatrices_.size(); i++) { // note that countsVec_ is one longer than transitionMatrices_, due to the rootPrior
      matrix_t const & c  = countsVec_[i];
      matrix_t const & tm = transitionMatrices_[i];
      for (unsigned k = 0; k < n; k++) 
	for (unsigned l = 0; l < n; l++) 
	  logLik += c(k, l) * log( tm(k, l) );
    }
    return logLik;
  }



//  vector_t extractBranchLengths(tree_t const & phyTree)
//  {
//    vector_t v( num_edges( phyTree ) );
//    treeTraits_t::edge_iterator ei, ei_end;
//    tie(ei, ei_end) = edges(phyTree);
//    for (; ei < ei_end; ++ei)
//      v[ phyTree[*ei].id ] = phyTree[*ei].length;
//    return v;
//  }
//
//  
//  vertex_t getParent(tree_t const & phyTree, vertex_t const & v)
//  {
//    treeTraits_t::in_edge_iterator ei, ei_end;
//    tie(ei, ei_end) = in_edges(v, phyTree);
//    if (ei >= ei_end) 
//      return treeTraits_t::null_vertex();
//    else
//      return source(*ei, phyTree);
//  }
//
//
//  int getParentBranchId(tree_t const & phyTree, vertex_t const & v)
//  {
//    treeTraits_t::in_edge_iterator ei, ei_end;
//    tie(ei, ei_end) = in_edges(v, phyTree);
//    if (ei_end - ei != 1) 
//      errorAbort("getParentBranchId: Error: Node should have exactly one parent");
//    return phyTree[*ei].id;
//  }
//
//
//  struct PreVisitFactorGraphConstructor {
//    PreVisitFactorGraphConstructor(tree_t const & phyTree, FG_t & factorGraph, std::vector<unsigned> & nodeMap, std::vector<unsigned> & branchMap, std::vector<FGVertex_t> & factorFirstVarMap) 
//      : phyTree_ (phyTree), factorGraph_(factorGraph), nodeMap_(nodeMap), branchMap_(branchMap), factorFirstVarMap_(factorFirstVarMap)
//    {
//      nodeMap_.resize( num_vertices(phyTree) );
//      branchMap_.resize( num_edges(phyTree) );
//      factorFirstVarMap_.resize( num_edges(phyTree) );
//    };
// 
//    void operator()(vertex_t const & v) 
//    {
//      FGVertex_t u = add_vertex(factorGraph_); // add node
//      nodeMap_[v] = u;
//      vertex_t phyParent = getParent(phyTree_, v);
//      if ( phyParent != treeTraits_t::null_vertex() ) {
//	int branchId = getParentBranchId(phyTree_, v);
//	FGVertex_t f = add_vertex(factorGraph_); // add factor
//	branchMap_[branchId] = f;
//	FGVertex_t FGParent = nodeMap_[phyParent];
//	factorFirstVarMap_[branchId] = FGParent;
//	// add edges
//	add_edge(FGParent, f, factorGraph_);
//	add_edge(f, u, factorGraph_);
//      }
//    }
//
//    tree_t const & phyTree_;
//    FG_t & factorGraph_;
//    std::vector<unsigned> & nodeMap_;
//    std::vector<unsigned> & branchMap_;
//    std::vector<FGVertex_t> & factorFirstVarMap_;
//  };
//
//
//  void mkFactorGraph(tree_t const & phyTree, 
//		     FG_t & factorGraph, 
//		     std::vector<unsigned> & nodeMap, 
//		     std::vector<unsigned> & branchMap, 
//		     std::vector<FGVertex_t> & factorFirstVarMap, 
//		     unsigned & initDistrNode)
//  {
//    PreVisitFactorGraphConstructor preVisitFGCons(phyTree, factorGraph, nodeMap, branchMap, factorFirstVarMap);
//    unsigned root = 0;
//    treeTraversal(phyTree, root, preVisitFGCons);
//
//    // add initial distribution factor
//    initDistrNode = add_vertex(factorGraph); // add factor
//    add_edge(root, initDistrNode, factorGraph);
//  }
//
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v)
//  {
//    std::vector<FGVertex_t> neighbors;
//    FGTraits_t::adjacency_iterator ai, ai_end;
//    tie(ai, ai_end) = adjacent_vertices(v, factorGraph);
//    for (; ai < ai_end; ++ai)
//      neighbors.push_back(*ai);
//    return neighbors;
//  }
//
//  // order neighbors according to firstVar variable
//  std::vector<FGVertex_t> getNeighbors(FG_t const & factorGraph, FGVertex_t v, FGVertex_t firstVar)
//  {
//    std::vector<FGVertex_t> neighbors = getNeighbors(factorGraph, v);
//    unsigned n = neighbors.size();
//    assert( n == 2);
//    if ( neighbors[0] == firstVar )
//      return neighbors;
//    else {
//      assert (neighbors[1] == firstVar);
//      reverse(neighbors.begin(), neighbors.end());
//      return neighbors;
//    }
//  }
//
//  vector<phy::FGBaseNode *> mkFactorGraphNodes(std::vector<unsigned> const & nodeMap, 
//					     FG_t const & factorGraph, 
//					     std::vector<unsigned> const & branchMap, 
//					     std::vector<FGVertex_t> const & factorFirstVarMap, 
//					     unsigned initDistrNode, 
//					     std::vector< matrix_t > const & TM, 
//					     vector_t const & equiFreq)
//  {
//    // to do: some check of consistency of factor node ids...
//    map<unsigned, FGBaseNode *> ndMap;
//
//    // dimension vectors
//    unsigned dim = equiFreq.size(); // all random variables of same size
//    vector<unsigned> dimVec1D(1, dim);
//    vector<unsigned> dimVec2D(2, dim);
//
//    for (unsigned i = 0; i < nodeMap.size(); i++)
//      ndMap[ nodeMap[i] ] = new DiscreteVariableNode(getNeighbors( factorGraph, nodeMap[i] ), dimVec1D);
//
//    for (unsigned i = 0; i < branchMap.size(); i++)
//      ndMap[ branchMap[i] ] = new DiscreteFactor2DNode(getNeighbors( factorGraph, branchMap[i], factorFirstVarMap[i] ), dimVec2D, TM[i] );
//				
//    ndMap[initDistrNode] = new DiscreteFactor1DNode(getNeighbors( factorGraph, initDistrNode ), dimVec2D, equiFreq );
//
//    vector<FGBaseNode *> v;
//    for (unsigned i = 0; i < ndMap.size(); i++)
//      v.push_back( ndMap[i] );
//    return v;
//  }
//
  // initFactorGraphNodes
  // input: nodeMap, branchMap, initDistrNode, TMs, 
  // need to define neigbors for each node


} //namespace phy
