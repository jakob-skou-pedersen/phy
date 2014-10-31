/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __DfgDataWrap_h
#define __DfgDataWrap_h

#include "phy/Observations.h"
#include "phy/Factors.h"
#include "phy/DiscreteFactorGraph.h"
#include "phy/DfgIO.h"

using namespace std;

namespace phy {

  // General assumptions: 
  //
  // The basic observation is a char, however, some random variables
  // are defined over multiple chars. Observation symbols are
  // therefore defined as strings (of chars) as specified by the
  // symbol_t type (in PhyDef.h). The fundamental input is therefore a
  // vector of symbols, one for each of the observed random variables
  // (RVs) in the discrete factor graph (DFG). The assignment of
  // symbols to observed RVs can be given by a varMap, which
  // indexes the observed variables in the DFG. 
  //

  // In many applications input can be defined as a vector of
  // observationStrings (string of chars) with the number of 
  // symbols (= same number of chars if all symbols of same size). 
  // Symbols are then implicitly given as consecutive and potentially overlapping
  // substrings of these.
  // 
  // The DFG opearates on states represented by unsigned integers. The
  // mappings between states and symbols are given by stateMaps. To
  // allow missing data and degenerate symbols, these are then mapped
  // to stateMasks (see stateMask_t defined in PhyDef.h), which are taken
  // as input by the algorithms acting on the DFG. The mapping from
  // states to stateMasks is done by staMasks. The staMaskSet class
  // holds stateMasks for all (observed) RVs of the DFG.


  /** Set up appropriate data sructures and calculate the DFG normalization constant for each position */
  vector<xnumber_t> strVec2NormConstVector(DFG & dfg, 
					   vector<string> const & strVec, 
					   vector<unsigned> const & varMap, 
					   StateMaskMapSet const & stateMaskMapSet, 
					   char missingDataChar = '.');

  /** Same as above, but intended for the specific case where all random variables are defined over the same state space. The stateMaskMapSet will thus be defined in terms of a single stateMap. */
  vector<xnumber_t> strVec2NormConstVector(DFG & dfg, vector<string> const & strVec, vector<unsigned> const & varMap, StateMap const & staMap, char missingDataChar = '.');

  /** factor set must define (normalized) probability distributions. */
  void dfgEm(AbstractBaseFactorSet & factorSet, vector<stateMask2DVec_t> const & dataVec, DFG & dfg, number_t minDeltaLogLik, unsigned maxIter, string const & logFile = "");

  /** If none-empty, specification of the dfg will be written to the logSpecification files after EM iteration. */
  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDataFile, number_t minDeltaLogLik, unsigned maxIter, string const & logFile);
  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDataFile, number_t minDeltaLogLik, unsigned maxIter, 
	     string const & logStateMapsFile, string const & logFactorPotentialsFile, string const & logVariablesFile, string const & logFactorGraphFile, string const & logFile = "");
  /** Subvar data can now be submitted. The other dfgEm's are wrappers for this one for backward compatibility */
  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDatafile, string const & subVarDataFile, number_t minDeltaLogLik, unsigned maxIter, string const & logFile);
  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDataFile, string const & subVarDataFile, number_t minDeltaLogLik, unsigned maxIter, 
	     string const & logStateMapsFile, string const & logFactorPotentialsFile, string const & logVariablesFile, string const & logFactorGraphFile, string const & logFile);

} // end namespace phy

#endif  //__DfgDataWrap_h
