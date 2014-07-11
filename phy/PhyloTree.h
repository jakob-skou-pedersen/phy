/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __PhyloTree_h
#define __PhyloTree_h

#include "phy/DiscreteFactorGraph.h"
#include "phy/TRCTMCTree.h"
#include "phy/ObservationsIO.h"
#include "phy/PhyIO.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/unordered_map.hpp> 

using namespace std;

namespace phy {

  // functions and data structures for input / output of time
  // reversible phylogenetic trees and accompanying data structures

  //typedef map<stateMaskVec_t, xnumber_t> probMap_t;
  typedef boost::unordered_map<stateMaskVec_t, xnumber_t> probMap_t;

  struct PhyloTree {
    /* Constructor **/
    PhyloTree(StateMap const & alphabet, boost::shared_ptr<BaseTRRateMatrix> const & rmPtr, TreeInfo const & treeInfo, StateMaskMapSet const & stateMaskMapSet, TrctmcFactorSet const & factorSet, DFG const & dfg)
      : alphabet(alphabet), rmPtr(rmPtr), treeInfo(treeInfo), stateMaskMapSet(stateMaskMapSet), factorSet(factorSet), dfg(dfg) {}

    StateMap alphabet;                   // alphabet
    boost::shared_ptr<BaseTRRateMatrix> rmPtr;  // rate matrix pointer
    TreeInfo treeInfo;                    // basic tree information
    StateMaskMapSet stateMaskMapSet;     // map to create state masks according to observations
    TrctmcFactorSet factorSet;           // set of all transitin matrices as well as the prior root distribution
    DFG dfg;                             // factor graph / Bayesian network implementation of the phylo tree, based on above data structures
  };

  // helper functions. To do: Could also be defined elsewhere

  /** Wrapper for constructing a TrctmcFactorSet. */ 
  TrctmcFactorSet mkTrctmcFactorSet(TreeInfo const & treeInfo, boost::shared_ptr<BaseTRRateMatrix> rateMatrixPtr, number_t rateScale = 1);

  /** Wrapper for constructing a DFG that represent a phylogenetic tree. */ 
  DFG mkPhyloDfg(TreeInfo const & treeInfo, AbstractBaseFactorSet const & factorSet);

  /** Convert input observation symbols to state masks. Overloaded version of function defined in Observations.h. offset defined to symbol length. */
  stateMask2DVec_t const mkStateMask2DVec(SeqData const & seqData, PhyloTree const & phyloTree, char missingDataChar = '.');

  /** Same as above, but allows offSet to be explicitly defined. */
  stateMask2DVec_t const mkStateMask2DVec(SeqData const & seqData, PhyloTree const & phyloTree, unsigned offSet, char missingDataChar = '.');

  /** Convert input observation symbols to state masks. Overloaded version of function defined in Observations.h. offset defined to symbol length. */
  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, PhyloTree const & phyloTree, char missingDataChar = '.');

  /** Same as above, but allows offSet to be explicitly defined. */
  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, PhyloTree const & phyloTree, unsigned offSet, char missingDataChar = '.');

  /** Return vector with likelihood of symbol column (site) given the
      phylogenetic tree and the substitution model. Overloaded
      versions return result by reference and allow only unseenn sites
      to be evaluated by use of hashing (a map for now).*/
  xvector_t calcLikelihoodVector(SeqData const & seqData, PhyloTree & phyloTree, unsigned offSet = 1, char missingDataChar = '.');
  void calcLikelihoodVector(xvector_t & result, SeqData const & seqdata, PhyloTree & phyloTree, unsigned offSet = 1, char missingDataChar = '.');
  void calcLikelihoodVector(xvector_t & result, SeqData const & seqdata, PhyloTree & phyloTree, probMap_t & probMap, bool useProbMap = true, number_t minMapProb = 1e-7, unsigned offSet = 1, char missingDataChar = '.');

  /** Return matrix with likelihood of all di-symbols columns made
      from a left and a right part. Overloaded versions return result
      by reference and allow only unseenn sites to be evaluated by use
      of hashing.*/
  xmatrix_t calcLikelihoodMatrix(SeqData const & seqData, PhyloTree & phyloTree, unsigned leftSize = 1, unsigned rightSize = 1, char missingDataChar = '.');
  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, PhyloTree & phyloTree,  unsigned leftSize = 1, unsigned rightSize = 1, char missingDataChar = '.');
  void calcLikelihoodMatrix(xmatrix_t & result, SeqData const & seqData, PhyloTree & phyloTree, probMap_t & probMap, bool useProbMap = true, number_t minMapProb = 1e-7, unsigned leftSize = 1, unsigned rightSize = 1, char missingDataChar = '.');



} // end namespace phy

#endif  //__PhyloTree_h
