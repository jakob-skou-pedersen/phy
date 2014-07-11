/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Grammar_h
#define __Grammar_h

#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <stack>
#include <map>
#include <boost/tuple/tuple.hpp>
#include "phy/GrammarDef.h"
#include "phy/utils.h"
#include "phy/utilsLinAlg.h"

namespace phy {

  namespace grammar {

    /** Transition types. Allowed types are:

    W_v -> s W_x    : LEFT
    W_v -> W_x s    : RIGHT
    W_v -> s W_x s  : PAIR
    W_v -> W_x      : SILENT
    W_v -> W_x W_y  : BIFURCATE
    W_v -> s        : TERMINAL
    W_v -> eps      : END  (perhaps not needed)

    Where W denote non-terminals, s denote terminals (with eps being the
    null string), and v,x,y are indices:
    */
    // If ttype_t is changed, then also change definition of "transTypeTagVec" in Grammar.cpp!
    enum ttype_t {LEFT, RIGHT, PAIR, SILENT, BIFURCATE, TERMINAL, END, TYPE_COUNT};  // TYPE_COUNT is simply the number of transition types. 

    /** Global variable defining a mapping of enum values to strings for convenience. */
    extern vector<string> const transTypeTagVec;


    struct Transition {
      /* Default constructor **/
      Transition() : from(-1), to(-1), toL(-1), toR(-1), ei(-1) {};

      // string tags for reference and human interaction
      string name;           // some unique name to identity the transition. May be used to map annotation onto parses
      string from_id;        // left hand side non-terminal.
      string to_id;          // right hand side non-terminal for LEFT, RIGHT, PAIR, and SILENT transtion types.
      string toL_id;         // leftmost right hand side non-terminal for BIFURCATE transtion type.
      string toR_id;         // rightmost right hand side non-terminal for BIFURCATE transtion type.
      string e_id;           // Emission distribution name         

      // raw data
      state_t from;           // left hand side non-terminal.
      state_t to;             // right hand side non-terminal for LEFT, RIGHT, PAIR, and SILENT transtion types.
      state_t toL;            // leftmost right hand side non-terminal for BIFURCATE transtion type.
      state_t toR;            // rightmost right hand side non-terminal for BIFURCATE transtion type.

      ttype_t type;           // transition type (see above)
      ynumber_t p;            // transition prob
      unsigned ei;            // Emission distribution index. Indexing is seperate for single and pair emissions (only relevant for LEFT, RIGHT, PAIR, and TERMINAL).
    };

    // convenience functions for handling Transitions */
    bool isEmitType(ttype_t type);


    /** Holds information on the transition which emitted a given position */
    struct EmitInfo {
      string name;       // name of transition
      bool left;         // true if left part of a pair emission (false if right part) 
      unsigned i;        // begin of substring before emission
      unsigned j;        // end (one past) of substring before emission

      // convenience setters
      /** reset EmitInfo data. leftPart only need to be set when dealing with pair types. */
      void reset(string const & n, unsigned iPos, unsigned jPos, bool leftPart = true) {name = n; i = iPos; j = jPos; left = leftPart;}
    };

    /** Indexing:
     ** 
     ** Intervals are generally defined as half open: [i;j). The start
     ** position is included and the end position excluded.
     */

    class Grammar {
    public:
      /**  convenient typedef's: */
      typedef vector<ynumber_t>  vector1D_t;  // One-dimensional table 
      typedef vector<vector1D_t> vector2D_t;  // Two-dimensional table 
      typedef vector<vector2D_t> vector3D_t;  // Three-dimensional table 
      typedef vector< pair<unsigned, unsigned> >  vector1DPair_t;  // one-dimensional table holding traceback values for cyk. the pair values are: 1 = max transition (u); 2 = max split (k), for BIFURCATE trans
      typedef vector<vector1DPair_t> vector2DPair_t;  // two-dimensional table holding traceback values for cyk. 
      typedef vector<vector2DPair_t> vector3DPair_t;  // Three-dimensional table holding traceback values for cyk.

      Grammar(vector<Transition> const & outTrans);
      Grammar(vector<Transition> const & outTrans, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames = vector<string>() );

      /* Data structures defining SCFG or HMM */
      vector<vector<Transition> > outTransitions;  // 2D vector with outgoing transitions (aka production rules). outTransitions[v] defines all the production rules going from state v to other states.
      vector<vector<Transition> > inTransitions;   // 2D vector with incoming transitions. inTransitions[v] defines all the production rules leading to state v from other states.
      unsigned stateCount;                         // number of states (aka non-terminals)

      /** reset transition probs for both outTransitions and inTransitions. Idexing is according to outProbs */
      void resetTransitionProbs(vector<vector_t> const & transitionProbs);

      /** reset emission probabilities. This need to be done before running any of the algorithms below. */
      void resetEmissions(vector<xvector_t> const & singleEmissions);
      void resetEmissions(vector<xmatrix_t> const & pairEmissions);
      void resetEmissions(vector<xvector_t> const & singleEmissions, vector<xmatrix_t> const & pairEmissions);

      /** Run the inside algorithm and return the likelihood of the data. */
      ynumber_t inside();
      ynumber_t inside(vector3D_t & tbl) const;

      /** Run the outside algorithm and return the likelihood of the data. Precondition: outside has been called. */
      void outside();
      void outside(vector3D_t & tbl, vector3D_t const & inTbl) const;

      /** Run inside as well as outside */
      ynumber_t insideOutside();
      ynumber_t insideOutside(vector3D_t & inTbl, vector3D_t & outTbl);

      /** Precondition for marginal calculations: inside as well as
	  outside has been run. Emissions must not have been reset.*/
      /** Calculate transition marginals. (State marginals can be
	  found as sum of marginals for outgoing transitions.) Useful
	  for EM optimization. Indexing is as for outTransition. The
	  size of sv and tv will be adjusted if necessary. */
      void calcAccumulatedMarginals(vector2D_t & tv) const;
      void calcAccumulatedMarginals(vector2D_t & tv, vector3D_t const & inTbl, vector3D_t const & outTbl) const;

      /** Run the inside algorithm and return the likelihood of the data. */
      ynumber_t cyk();
      ynumber_t cyk(vector3D_t & tbl);
      ynumber_t cyk(vector3D_t & tbl, vector3DPair_t & tbTbl) const;

      /** To do: maximal expected accuracy version of cyk. Assumes
	  emission probs to be filled with post probs values. Ignores
	  transition probs, all set to one. Emission probs should be
	  added instead of multiplied. Since accuracy is per position,
	  pairs should be weighted by two. */

      /** Derive the most likely parse and return the transitions that generated each positin of the sequence. */
      vector<EmitInfo> const maxParse() const;  // assumes that cyk() has just been run (filling out the tb_ table)
      vector<EmitInfo> const maxParse(vector3DPair_t const & tbTbl) const;  // assumes tbTbl to be based on sequence of same length as current emission models

      /** Precondition: inside and outside have been called. */
      /** Lookup posterior probability of each Transition defined by EmitInfo and return in vector. Can be used to find post prob of every emission in maxParse. */ 
      vector_t const postProbParse(vector<EmitInfo> const & emitInfoVec) const;
      vector_t const postProbParse(vector<EmitInfo> const & emitInfoVec, vector3D_t const & inTbl, vector3D_t const & outTbl) const;

      /** Same as above, but also takes an equivalence map between emitting transitions. This can be created based on an AnnoMap and allows calculation of the marginal probability of a given type of annotaion at each position. */
      vector_t const postProbParse(vector<EmitInfo> const & emitInfoVec, map<string, vector<string> > const & equivalenceMap) const;
      vector_t const postProbParse(vector<EmitInfo> const & emitInfoVec, vector3D_t const & inTbl, vector3D_t const & outTbl, map<string, vector<string> > const & equivalenceMap) const;

      /** Initialize dynamic programming table. M denotes stateCout and L sequence length.*/
      void initDpTable(vector3D_t & tbl, unsigned M, unsigned long L) const;
      void initDpTable(vector3D_t & tbl) const; // will only init if current size smaller than needed

      /** Set all entries of dp table to zero..*/
      void clearDpTable(vector3D_t & tbl) const;
      void clearDpTable(vector3D_t & tbl, int M, long L) const;

      /** Initialize backtrace table used with cyk. M denotes stateCout and L sequence length. */
      void initTbTable(vector3DPair_t & tbl, unsigned M, unsigned long L) const;
      void initTbTable(vector3DPair_t & tbl) const; // will only init if current size smaller than needed

      ////////////////////////////////////////////////////////////////
      // helper functions
      ////////////////////////////////////////////////////////////////

      /** Precondition: inside and outside have been called. */
      /** Returns the probability (p) of using Transition t for
	  deriving subsequence [i;j). This probability devided by the
	  likelihood of the data gives the posterior probability of
	  using the given transition (pp( t-> [i;j) ) = p /
	  inside() ).*/
      inline ynumber_t transitionUseProb(Transition const & t, unsigned i, unsigned j, vector3D_t const & inTbl, vector3D_t const & outTbl) const;
      inline ynumber_t transitionUseProb(Transition const & t, unsigned i, unsigned j) const;
   
      /** Precondition: inside and outside have been called. */
      /** Returns the probability (p) of using Transition t emitting
	  position i. This version will sum over all
	  sub-intervals. The left flag should be set for PAIR
	  transitions, specifying if position was left part of pair (left = true) or right part (left = false). */
      inline ynumber_t emitTransitionAccUseProb(Transition const & t, unsigned i, vector3D_t const & inTbl, vector3D_t const & outTbl, bool left = true) const;
      inline ynumber_t emitTransitionAccUseProb(Transition const & t, unsigned i, bool left = true) const;

      /** Precondition: inside has been callled */
      /** Returns the inside probability for state (non-terminal) v and for subsequence [i;j). */
      inline ynumber_t subSeqInsideProb(unsigned v, unsigned i, unsigned j, vector3D_t const & inTbl) const;
      inline ynumber_t subSeqInsideProb(unsigned v, unsigned i, unsigned j) const;

      /** Precondition: cyk has been callled */
      /** Returns the cyk probability for state (non-terminal) v and for subsequence [i;j). */
      inline ynumber_t subSeqCykProb(unsigned v, unsigned i, unsigned j, vector3D_t const & cykTbl) const;
      inline ynumber_t subSeqCykProb(unsigned v, unsigned i, unsigned j) const;

      /** Returns the index of state 'name'. */
      unsigned stateIndex(string const & name) const;

      /** convert outTransition into a map with name keys. */
      map<string, Transition> mkTransitionMap() const;

      /** Returns vector with all emit transitions */
      vector<Transition> emitTransitions() const;

    protected:

      // helper functions
      void init(vector<Transition> transVec, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames); // init data structures
      void setIndices(vector<Transition> & transVec, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames) const; // convert string_ids to indices
      vector<string> mkStateOrder(vector<Transition> const & transVec) const;  // derives order of states from order of transitions (using from_id)
      vector<string> deriveEmissionOrder(vector<Transition> const & transVec, bool pairType) const;  // derives order in which emission distributions are referred to in transVec. 
      void consistencyCheck() const;
      void mkInTransitions();
      vector<vector<Transition> > const orderByState(vector<Transition> const & transVec) const;
      
      // interface to emissions
      ynumber_t e(unsigned idx, unsigned long i) const;  // emission are indexed by their start position (if overlapping)
      ynumber_t e(unsigned idx, unsigned long i, unsigned long j) const; // both i of the left symbol and j refers to the start position of the right symbol.  

      /** position specific emission probabilities */
      vector<xvector_t> singleEmissions_;          // vector of probability of terminals given a production rule P(s|rule).
      vector<xmatrix_t> pairEmissions_;            // vector of probability of pairs of terminals given a production rule P(s_l,s_r|rule).
      unsigned inLength_;                         // length of input emssion tables. Corresponds to length of original input sequence. Could be changed to long
    
      // dynamic programming tables
      vector3D_t cyk_;
      vector3DPair_t bt_;  // backtrace table
      vector3D_t in_;
      vector3D_t out_;
    };


    // helper functions
    ostream & printTable(ostream & str, Grammar::vector3D_t const & tbl);


    ////////////////////////////////////////////////////////////////
    // AnnoMaps
    ////////////////////////////////////////////////////////////////

    /** Structure holding annotation defined in terms of emission distributions */
    struct EmitAnno {
      EmitAnno() {};
      EmitAnno(string const & n, string const & a, string const & al, string const & ar): name(n), anno(a), annoLeft(al), annoRight(ar) {};
      string name;       // name of emission model
      string anno;       // use for anno of single emitting transitions  ("" when not set.)
      string annoLeft;   // use for anno symbol of left part of pair emitting
      string annoRight;  // use for anno symbol of right part of pair emitting

      /** returns true if annotation is copatible with EmitAnno. */
      bool singleCompatible(string const & a, string const & wildCard = "*") const;
      bool pairCompatible(string const & l, string const & r, string const & wildCard = "*") const;
    };

    typedef map<string, EmitAnno> EmitAnnoMap;

    /** sets emission to zero if not compatible with annotation. */
    void maskEmissions(vector<xvector_t> & singleEmissions, vector<string> const & emitModelNames, string const & anno, EmitAnnoMap const & eam);
    void maskEmissions(vector<xvector_t> & singleEmissions, vector<string> const & emitModelNames, vector<string> const & anno, EmitAnnoMap const & eam);
    void maskEmissions(xvector_t & singleEmissions, vector<string> const & anno, EmitAnno ea);
    void maskEmissions(vector<xmatrix_t> & pairEmissions, vector<string> const & emitModelNames, string const & anno, EmitAnnoMap const & eam);
    void maskEmissions(vector<xmatrix_t> & pairEmissions, vector<string> const & emitModelNames, vector<string> const & anno, EmitAnnoMap const & eam);
    void maskEmissions(xmatrix_t & pairEmissions, vector<string> const & anno, EmitAnno ea);

    /** Structure holding annotation defined directly in terms of transitions. Normally derived from the above based on a set of transitions. */
    struct TransAnno {
      TransAnno() {};
      TransAnno(string const & n, string const & a, string const & al, string const & ar): name(n), anno(a), annoLeft(al), annoRight(ar) {};
      string name;       // name of transition
      string anno;       // use for anno of single emitting transitions  ("" when not set.)
      string annoLeft;   // use for anno symbol of left part of pair emitting
      string annoRight;  // use for anno symbol of right part of pair emitting
    };

    class TransAnnoMap {
    public:
      TransAnnoMap(map<string, TransAnno> const & annoMap) : annoMap_(annoMap) {};
      TransAnnoMap(EmitAnnoMap const & emitAnnoMap, vector<Transition> const & emitTransVec);
    
      string convert(EmitInfo const & emitInfo) const;
      string convertToString(vector<EmitInfo> const & emitInfoVec, string const & sep = "") const;
      vector<string> convertToVector(vector<EmitInfo> const & emitInfoVec) const;

      /* Returns names of all transitions with equivalent anno to the given 'transAnno' **/
      vector<string> equivalentTransitions(TransAnno const & transAnno) const;

      /* Same as above, but defined anno in terms of a transition name. **/
      vector<string> equivalentTransitions(string const & transName) const;

    protected:
      map<string, TransAnno> annoMap_;
    };

    /** Return map between transition with same annotation in TransAnnoMap. Transitions are defined by their names. */
    map<string, vector<string> > mkEquivalenceMap(Grammar const & g, TransAnnoMap const & am);

    ////////////////////////////////////////////////////////////////
    // Inline function implementations
    ////////////////////////////////////////////////////////////////


    inline ynumber_t Grammar::e(unsigned idx, unsigned long i) const
    {
      return singleEmissions_[idx][i];
    }


    inline ynumber_t Grammar::e(unsigned idx, unsigned long i, unsigned long j) const
    {
      return pairEmissions_[idx](i, j);
    }


    inline ynumber_t Grammar::transitionUseProb(Transition const & t, unsigned i, unsigned j, vector3D_t const & inTbl, vector3D_t const & outTbl) const
    {
      if (t.type == LEFT and j - i > 0)
	return outTbl[t.from][i][j]* t.p * e(t.ei, i) * inTbl[t.to][i + 1][j];
      else if (t.type == RIGHT and j - i > 0)
	return outTbl[t.from][i][j]* t.p * e(t.ei, j - 1) * inTbl[t.to][i][j - 1];
      else if (t.type == PAIR and j - i > 0)
	return outTbl[t.from][i][j]* t.p * e(t.ei, i, j - 1) * inTbl[t.to][i + 1][j - 1];
      else if (t.type == SILENT)
	return outTbl[t.from][i][j]* t.p * inTbl[t.to][i][j];
      else if (t.type == BIFURCATE) {
	ynumber_t sumSplit = 0;
	for (unsigned k = i; k <= j; k++) 
	  sumSplit += inTbl[t.toL][i][k] * inTbl[t.toR][k][j];
	return outTbl[t.from][i][j] * t.p * sumSplit;
      }
      else if (t.type == TERMINAL and j - i == 1)
	return outTbl[t.from][i][j] * t.p * e(t.ei, i); 
      else if (t.type == END and j == i) 
	return outTbl[t.from][i][j] * t.p;
      else
	return 0; // certain combinations of t.type(), i, j.
    }


    inline ynumber_t Grammar::emitTransitionAccUseProb(Transition const & t, unsigned i, bool left) const
    {
      return emitTransitionAccUseProb(t, i, in_, out_, left);
    }


    inline ynumber_t Grammar::emitTransitionAccUseProb(Transition const & t, unsigned i, vector3D_t const & inTbl, vector3D_t const & outTbl, bool left) const
    {
      if (not isEmitType(t.type) )
	errorAbort("From Grammar::emitTransitionAccUseProb: Transition must be emit type.");

      ynumber_t sum = 0;
      if (t.type == LEFT or (t.type == PAIR and left) )
	for (unsigned j = i + 1; j <= inLength_; j++)
	  sum += transitionUseProb(t, i, j, inTbl, outTbl);
      else if (t.type == RIGHT or (t.type == PAIR and not left) ) {
	unsigned const l = i + 1; // since l denotes the end of sub int [k,l)
	for (unsigned k = 0; k < l; k++)
	  sum += transitionUseProb(t, k, l, inTbl, outTbl);
      }
      else if (t.type == TERMINAL)
	sum += transitionUseProb(t, i, i + 1, inTbl, outTbl);
      return sum;
    }
  

    inline ynumber_t Grammar::subSeqInsideProb(unsigned v, unsigned i, unsigned j, vector3D_t const & inTbl) const
    {
      assert(i <= inLength_ and j <= inLength_);
      return inTbl[v][i][j];
    }
  

    inline ynumber_t Grammar::subSeqInsideProb(unsigned v, unsigned i, unsigned j) const
    {
      return subSeqInsideProb(v, i, j, in_);
    }


    inline ynumber_t Grammar::subSeqCykProb(unsigned v, unsigned i, unsigned j, vector3D_t const & cykTbl) const
    {
      assert(i <= inLength_ and j <= inLength_);
      return cykTbl[v][i][j];
    }
  

    inline ynumber_t Grammar::subSeqCykProb(unsigned v, unsigned i, unsigned j) const
    {
      return subSeqCykProb(v, i, j, cyk_);
    }


  } // end namespace grammar

} // end namespace phy

#endif  // __Grammar_h
