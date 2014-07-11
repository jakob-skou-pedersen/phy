/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/Grammar.h"

namespace phy {

  namespace grammar {

    /** global variable defining a mapping of enum values to strings for convenience. */
    static char const * arr[] = {"LEFT", "RIGHT", "PAIR", "SILENT", "BIFURCATE", "TERMINAL", "END", "TYPE_COUNT"}; 
    vector<string> const transTypeTagVec(arr, arr + sizeof(arr) / sizeof( arr[0] ) );


    ////////////////////////////////////////////////////////////////
    // Transition
    ////////////////////////////////////////////////////////////////

    // convenience functions for dealing with transitions
    bool isEmitType(ttype_t type)
    {
      return (type == LEFT or type == RIGHT or type == TERMINAL or type == PAIR); // emitting 
    }

    ////////////////////////////////////////////////////////////////
    // Grammar
    ////////////////////////////////////////////////////////////////

    Grammar::Grammar(vector<Transition> const & outTrans) //: outTransitions( orderByState(transVec) ), stateCount(outTransitions.size() )
    {
      vector<string> singleEmitModelNames = Grammar::deriveEmissionOrder(outTrans, false);
      vector<string> pairEmitModelNames   = Grammar::deriveEmissionOrder(outTrans, true);
      init(outTrans, singleEmitModelNames, pairEmitModelNames);
    }


    Grammar::Grammar(vector<Transition> const & outTrans, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames)
    {
      init(outTrans, singleEmitModelNames, pairEmitModelNames);
    }


    void Grammar::init(vector<Transition> outTrans, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames)
    {
      setIndices(outTrans, singleEmitModelNames, pairEmitModelNames);

      // init data members
      outTransitions = orderByState(outTrans);
      stateCount = outTransitions.size();
    
      consistencyCheck();
      mkInTransitions();
    }


    void Grammar::setIndices(vector<Transition> & transVec, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames) const
    {
      vector<string> stateOrder = mkStateOrder(transVec);
      for (unsigned i = 0; i < transVec.size(); i++) {
	Transition & t = transVec[i];
	t.from = getIndex(stateOrder, t.from_id);
	if (t.to_id.size() != 0)
	  t.to = getIndex(stateOrder, t.to_id);
	if (t.toL_id.size() != 0)
	  t.toL = getIndex(stateOrder, t.toL_id);
	if (t.toR_id.size() != 0)
	  t.toR = getIndex(stateOrder, t.toR_id);
	if ( isEmitType(t.type) ) {
	  if (t.e_id.size() == 0)
	    errorAbort("from: Grammar::setIndices: Emit types must define emission distributions! Transition with name '" + t.name + "' did not");
	  if (t.type == PAIR)
	    t.ei = getIndex(pairEmitModelNames, t.e_id);
	  else
	    t.ei = getIndex(singleEmitModelNames, t.e_id);
	}
      }
    }


    // derives (from) states in the order they occur in the transition list
    vector<string> Grammar::mkStateOrder(vector<Transition> const & transVec) const
    {
      vector<string> v;
      for (unsigned i = 0; i < transVec.size(); i++) {
	string const & from_id =  transVec[i].from_id;
	if ( not hasElement(v, from_id) )
	  v.push_back(from_id);
      }
      return v;
    }

    // derives (from) states in the order they occur in the transition list. If pairType is true, then pair emission distributions are 
    vector<string> Grammar::deriveEmissionOrder(vector<Transition> const & transVec, bool pairType) const
    {
      vector<string> v;
      for (unsigned i = 0; i < transVec.size(); i++) {
	Transition const & t = transVec[i];
	if ( isEmitType(t.type) )
	  if ( (pairType and t.type == PAIR) or (not pairType and t.type != PAIR) )
	    if ( not hasElement(v, t.e_id) )
	      v.push_back(t.e_id);
      }
      return v;
    }


    vector<vector<Transition> > const Grammar::orderByState(vector<Transition> const & transVec) const
    {
      unsigned stateCount = 0;
      for (unsigned i = 0; i < transVec.size(); i++)
	if (transVec[i].from + 1 > stateCount)
	  stateCount = transVec[i].from + 1;

      vector<vector<Transition> > traVec2D(stateCount, vector<Transition>() ); // due to zero-based enumeration of states
      for (unsigned i = 0; i < transVec.size(); i++) {
	Transition const & t = transVec[i];
	traVec2D[t.from].push_back(t);
      }
      return traVec2D;
    }


    void Grammar::resetTransitionProbs(vector<vector_t> const & transitionProbs)
    {
      assert( outTransitions.size() == transitionProbs.size() );
      for (unsigned i = 0; i < outTransitions.size(); i++) {
	assert( outTransitions[i].size() == transitionProbs[i].size() );
	for (unsigned j = 0; j < outTransitions[i].size(); j++)
	  outTransitions[i][j].p = transitionProbs[i][j];
      }
      mkInTransitions(); // updating transition prob for in transitions
    }

    void Grammar::resetEmissions(vector<xvector_t> const & singleEmissions)
    {
      // unsigned n = singleEmissions.size();
      //      singleEmissions_.resize(n);
      //      for (unsigned i = 0; i < n; i++) 
      //	singleEmissions_[i] = toYNumber(singleEmissions[i]);
      singleEmissions_ = singleEmissions;
      assert(singleEmissions.size() > 0);
      inLength_ = singleEmissions[0].size();
    }

    void Grammar::resetEmissions(vector<xvector_t> const & singleEmissions, vector<xmatrix_t> const & pairEmissions)
    {
      resetEmissions(singleEmissions);
      resetEmissions(pairEmissions);
    }


    void Grammar::resetEmissions(vector<xmatrix_t> const & pairEmissions)
    {
      //      unsigned n = pairEmissions.size();
      //      pairEmissions_.resize(n);
      //      for (unsigned i = 0; i < n; i++) 
      //	pairEmissions_[i] = toYNumber(pairEmissions[i]);
      pairEmissions_ = pairEmissions;
      assert(pairEmissions.size() > 0);
      inLength_ = pairEmissions[0].size1();
    }


    ostream & printTable(ostream & str, Grammar::vector3D_t const & tbl)
    {
      str << "dim 1: " << tbl.size() << endl;
      if ( tbl.size() ) {
	str << "dim 2: " << tbl[0].size() << endl;
	if ( tbl[0].size() ) {
	  str << "dim 3: " << tbl[0][0].size() << endl;
	}
      }
      str << endl;

      for (unsigned v = 0; v < tbl.size(); v++) { 
	str << "v=" << v << endl;
	for (unsigned i = 0; i < tbl[v].size(); i++) { 
	  for (unsigned j = 0; j < tbl[v][i].size(); j++)
	    str << tbl[v][i][j] << "\t";
	  str << endl;
	}
	str << endl << endl;
      }
      return str;
    }


    ostream & printTbTable(ostream & str, Grammar::vector3DPair_t const & tbl)
    {
      str << "dim 1: " << tbl.size() << endl;
      if ( tbl.size() ) {
	str << "dim 2: " << tbl[0].size() << endl;
	if ( tbl[0].size() ) {
	  str << "dim 3: " << tbl[0][0].size() << endl;
	}
      }
      str << endl;

      for (unsigned v = 0; v < tbl.size(); v++) { 
	str << "v=" << v << endl;
	for (unsigned i = 0; i < tbl[v].size(); i++) { 
	  for (unsigned j = 0; j < tbl[v][i].size(); j++)
	    str << toString(tbl[v][i][j]) << "\t";
	  str << endl;
	}
	str << endl << endl;
      }
      return str;
    }


    ynumber_t Grammar::inside(vector3D_t & tbl) const
    {
      const int M = stateCount;
      const int L = inLength_;

      //      clearDpTable(tbl, M, L); // debug -- probably not needed...

      // ensure zero-length-interval entires are set to zero (only needed if tbl is reused)
      for (int j = 0; j < L + 1; j++) 
	for (int v = M - 1; v >= 0; v--) 
	  tbl[v][j][j] = 0;
      
      // initialisation of zero length intervals
      for (int j = 0; j < L + 1; j++) { 
	for (int v = M - 1; v >= 0; v--) { 
	  ynumber_t sum = 0;
	  for (unsigned u = 0; u < outTransitions[v].size(); u++) { 	// for each out transition
	    Transition const & t = outTransitions[v][u];
	    if (t.type == END)
	      sum += t.p;
	    else if (t.type == SILENT)
	      sum += t.p * tbl[t.to][j][j]; // silent self loops are not allowed, so fine to add to tmp sum.
	    else if (t.type == BIFURCATE)
	      sum += t.p * tbl[t.toL][j][j] * tbl[t.toR][j][j];
	    // remaining transition types adds zero prob
	  }
	  tbl[v][j][j] = sum;
	}
      }

      //    printTable(cout, tbl); // debug
      // Transition types
      //
      //      W_v -> s W_x    : LEFT
      //      W_v -> W_x s    : RIGHT
      //      W_v -> s W_x s  : PAIR
      //      W_v -> W_x      : SILENT
      //      W_v -> W_x W_y  : BIFURCATE
      //      W_v -> s        : TERMINAL
      //      W_v -> eps      : END  (perhaps not needed)

      // recursion
      for (int j = 1; j < L + 1; j++) { 
	for (int i = j - 1; i >= 0; i--) {
	  for (int v = M - 1; v >= 0; v--) {
	    ynumber_t sum = 0;
	    for (unsigned u = 0; u < outTransitions[v].size(); u++) { 	// for each out transition
	      Transition const & t = outTransitions[v][u];
	      if (t.type == LEFT)
		sum += t.p * e(t.ei, i) * tbl[t.to][i + 1][j];
	      else if (t.type == RIGHT)
		sum += t.p * e(t.ei, j - 1) * tbl[t.to][i][j - 1]; 
	      else if (t.type == PAIR and j-i >= 2)
		sum += t.p * e(t.ei, i, j - 1) * tbl[t.to][i + 1][j - 1]; 
	      else if (t.type == SILENT)
		sum += t.p * tbl[t.to][i][j];
	      else if (t.type == BIFURCATE) {
		ynumber_t sumSplit = 0;
		for (int k = i; k <= j; k++) 
		  sumSplit += tbl[t.toL][i][k] * tbl[t.toR][k][j];
		sum += t.p * sumSplit;
	      }
	      else if (t.type == TERMINAL and j - i == 1) // TERMINALs emit one position and can therefore only parse one long subsequences. Could be done in init step
		sum += t.p * e(t.ei, i); 
	      // END type transitions adds zero prob for subseqs of positive length
	    }
	    tbl[v][i][j] = sum;
	  }
	}
      }

      //    cout << "inside table: " << endl;
      //    printTable(cout, tbl); //debug
      return tbl[0][0][L];
    }

    ynumber_t Grammar::inside()
    {
      initDpTable(in_);
      return inside(in_);
    }

    void Grammar::outside()
    {
      initDpTable(out_);
      return outside(out_, in_);
    }


    void Grammar::outside(vector3D_t & tbl, vector3D_t const & inTbl) const
    {
      const unsigned M = stateCount;
      const int L = inLength_;

      //      clearDpTable(tbl, M, L); // debug -- probably not needed...
    
      // initialisation
      tbl[0][0][L] = 1;
      for (unsigned v = 1; v < M; v++)
	tbl[v][0][L] = 0;  // this should not be needed, since parents are evaluated before their children for a given [i;j) interval.

      //    printTable(cout, tbl); // debug
      // Transition types
      //
      //      W_v -> s W_x    : LEFT
      //      W_v -> W_x s    : RIGHT
      //      W_v -> s W_x s  : PAIR
      //      W_v -> W_x      : SILENT
      //      W_v -> W_x W_y  : BIFURCATE
      //      W_v -> s        : TERMINAL
      //      W_v -> eps      : END  (perhaps not needed)

      // loop variables
      ynumber_t sum;
      ynumber_t sumSplit = 0;

      // recursion
      for (int i = 0; i <= L; i++) {
	for (int j = L; j >= i; j--) { // need zero length intervals?
	  for (unsigned v = 0; v < M; v++) {
	    sum = 0;
	    for (unsigned u = 0; u < inTransitions[v].size(); u++) { 	// for each *in* transition
	      Transition const & t = inTransitions[v][u];  
	      if (t.type == LEFT and i > 0)
		sum += t.p * e(t.ei, i - 1) * tbl[t.from][i - 1][j];
	      else if (t.type == RIGHT and j < L)
		sum += t.p * e(t.ei, j) * tbl[t.from][i][j + 1]; 
	      else if (t.type == PAIR and i > 0 and j < L and j-i>=0) 
		sum += t.p * e(t.ei, i - 1, j) * tbl[t.from][i - 1][j + 1]; 
	      else if (t.type == SILENT)
		sum += t.p * tbl[t.from][i][j];
	      else if (t.type == BIFURCATE) {
		sumSplit = 0;
		if (t.toL == v) 
		  for (int k = j; k <= L; k++) 
		    sumSplit += tbl[t.from][i][k] * inTbl[t.toR][j][k];
		if (t.toR == v) // note that incoming transitions with toL == toR should only be listed once. If listed twice they will contribute double...
		  for (int k = 0; k <= i; k++) 
		    sumSplit += inTbl[t.toL][k][i] * tbl[t.from][k][j];
		sum += t.p * sumSplit;
	      }
	      else if (t.type == TERMINAL or t.type == END) // these transitions cannot be incoming
		errorAbort("From Grammar::outside. Unexpected incoming transition type to state " + toString(v) + ". \nThis should not be possible. ");
	    }
	    tbl[v][i][j] = sum;
	    // initialization is overirdden for state 0
	    if (i == 0)
	      if (j == L)
		if (v == 0)
		  tbl[v][i][j] = 1;
	  }
	}
      }

      //    cout << "outside table:" << endl;
      //    printTable(cout, tbl); //debug
      return;
    }


    ynumber_t Grammar::insideOutside()
    {
      initDpTable(in_);
      initDpTable(out_);
      return insideOutside(in_, out_);
    }


    ynumber_t Grammar::insideOutside(vector3D_t & inTbl, vector3D_t & outTbl)
    {
      ynumber_t p = inside(inTbl);
      outside(outTbl, inTbl);
      return p;
    }


    void Grammar::calcAccumulatedMarginals(vector2D_t & tv) const
    {
      calcAccumulatedMarginals(tv, in_, out_);
    }


    void Grammar::calcAccumulatedMarginals(vector2D_t & tv, vector3D_t const & inTbl, vector3D_t const & outTbl) const
    {
      // init result data structure if necessary
      if (tv.size() != stateCount)
	tv.resize(stateCount);
      for (unsigned v = 0; v < stateCount; v++) 
	if ( tv[v].size() != outTransitions[v].size() )
	  tv[v].resize( outTransitions[v].size() );

      // check table sizes
      assert( stateCount == inTbl.size() and stateCount == outTbl.size() );
      if (stateCount > 0) {
	assert( inTbl[0].size() > inLength_ and outTbl[0].size() > inLength_ );
	assert( inTbl[0][0].size() > inLength_ and outTbl[0][0].size() > inLength_ );
      }

      unsigned M = stateCount;
      unsigned L = inLength_;
    
      for (unsigned v = 0; v < M; v++) 
	for (unsigned u = 0; u < outTransitions[v].size(); u++) { 	// for each out transition
	  Transition const & t = outTransitions[v][u];
	  ynumber_t sum = 0;
	  for (unsigned j = 0; j <= L; j++) 
	    for (unsigned i = 0; i <= j ; i++) 
	      sum += transitionUseProb(t, i, j, inTbl, outTbl);
	  tv[v][u] = sum / inTbl[0][0][L];
	}
    }


    ynumber_t Grammar::cyk()
    {
      initDpTable(cyk_);
      return cyk(cyk_);
    }


    ynumber_t Grammar::cyk(vector3D_t & tbl)
    {
      initTbTable(bt_);
      return cyk(tbl, bt_);
    }


    ynumber_t Grammar::cyk(vector3D_t & tbl, vector3DPair_t & tbTbl) const
    {
      const int M = stateCount;
      const int L = inLength_;
    
      //loop variables
      ynumber_t tmp;
      ynumber_t m;
      ynumber_t maxSplit;
      unsigned maxK = 0;  // =0 to please compiler
      pair<unsigned, unsigned> bt;
    
      // initialisation of zero length intervals
      for (int j = 0; j < L + 1; j++) { 
	for (int v = M - 1; v >= 0; v--) { 
	  m = 0;
	  for (unsigned u = 0; u < outTransitions[v].size(); u++) { 	// for each out transition
	    Transition const & t = outTransitions[v][u];
	    if (t.type == END) {
	      if (t.p > m) {
		m = t.p;
		bt.first = u;
	      }
	    }
	    else if (t.type == SILENT) {
	      if ( (tmp = t.p * tbl[t.to][j][j]) > m) {
		m = tmp;
		bt.first = u;
	      }
	    }
	    else if (t.type == BIFURCATE) {
	      if ( (tmp = t.p * tbl[t.toL][j][j] * tbl[t.toR][j][j]) > m) {
		m = tmp;
		bt.first = u;
		bt.second = 0;
	      }
	    }
	    // remaining transition types adds zero prob
	  }
	  tbl[v][j][j] = m;
	  tbTbl[v][j][j] = bt;
	}
      }

      //    printTable(cout, tbl); // debug

      // Transition types
      //
      //      W_v -> s W_x    : LEFT
      //      W_v -> W_x s    : RIGHT
      //      W_v -> s W_x s  : PAIR
      //      W_v -> W_x      : SILENT
      //      W_v -> W_x W_y  : BIFURCATE
      //      W_v -> s        : TERMINAL
      //      W_v -> eps      : END  (perhaps not needed)

      // recursion
      for (int j = 1; j < L + 1; j++) { 
	for (int i = j - 1; i >= 0; i--) {
	  for (int v = M - 1; v >= 0; v--) {
	    m = 0;
	    for (unsigned u = 0; u < outTransitions[v].size(); u++) { 	// for each out transition
	      Transition const & t = outTransitions[v][u];
	      if (t.type == LEFT) {
		if ( ( tmp = t.p * e(t.ei, i) * tbl[t.to][i + 1][j] ) > m) {
		  m = tmp;
		  bt.first = u;
		}
	      }
	      else if (t.type == RIGHT) {
		if ( ( tmp = t.p * e(t.ei, j - 1) * tbl[t.to][i][j - 1] ) > m) {
		  m = tmp;
		  bt.first = u;
		}
	      }
	      else if (t.type == PAIR and j-i >= 2) {
		if ( ( tmp = t.p * e(t.ei, i, j - 1) * tbl[t.to][i + 1][j - 1] ) > m) {
		  m = tmp;
		  bt.first = u;
		}
	      }
	      else if (t.type == SILENT) {
		if ( ( tmp = t.p * tbl[t.to][i][j] ) > m) {
		  m = tmp;
		  bt.first = u;
		}
	      }
	      else if (t.type == BIFURCATE) {
		maxSplit = 0;
		for (int k = i; k <= j; k++)
		  if ( ( tmp = tbl[t.toL][i][k] * tbl[t.toR][k][j] ) > maxSplit) {
		    maxSplit = tmp;
		    maxK = k;
		  }
		if ( (tmp = t.p * maxSplit) > m) {
		  m = tmp;
		  bt.first = u;
		  bt.second = maxK;
		}
	      }
	      else if (t.type == TERMINAL and j - i == 1) // TERMINALs emit one position and can therefore only parse one long subsequences. Could be done in init step
		if ( ( tmp = t.p * e(t.ei, i) ) > m) {
		  m = tmp;
		  bt.first = u;
		}
	      // END type transitions adds zero prob for subseqs of positive length
	    }
	    tbl[v][i][j] = m;
	    tbTbl[v][i][j] = bt;
	  }
	}
      }
    
      // printTable(cout, tbl); //debug
      // printTbTable(cout, tbTbl); //debug

      return tbl[0][0][L];
    }


    vector<EmitInfo> const Grammar::maxParse() const{
      return maxParse(bt_);
    }

  
    vector<EmitInfo> const Grammar::maxParse(vector3DPair_t const & tbTbl) const {
      assert(tbTbl.size() > 0 );
      assert(tbTbl[0].size() > inLength_);

      vector<EmitInfo> emitInfoVec(inLength_);  // assuming this is the right size

      // initiate stack
      stack<boost::tuple<unsigned, unsigned, unsigned> > maxStack;
      maxStack.push( boost::make_tuple(0, 0, inLength_) ); // start traceback at state 0, entire sequence

      // traceback
      while (maxStack.size() > 0) {
	unsigned v, i, j; // table indexes
	unsigned u, k;    // transition index (u) and maxSplit value (k)
	boost::tie(v, i, j) = maxStack.top();
	maxStack.pop(); // due to design of stl stack
	boost::tie(u, k) = tbTbl[v][i][j];
	Transition const & t = outTransitions[v][u];

	// set name if emission and update maxStack
	if (t.type == LEFT) {
	  emitInfoVec[i].reset(t.name, i, j);
	  maxStack.push( boost::make_tuple(t.to, i + 1, j) );
	}
	else if (t.type == RIGHT) {
	  emitInfoVec[j - 1].reset(t.name, i, j);
	  maxStack.push( boost::make_tuple(t.to, i, j - 1) );
	}
	else if (t.type == PAIR) {
	  emitInfoVec[i].reset(t.name, i, j, true);       // left part
	  emitInfoVec[j - 1].reset(t.name, i, j, false); 	// right part
	  maxStack.push( boost::make_tuple(t.to, i + 1, j - 1) );
	}
	else if (t.type == TERMINAL) 
	  emitInfoVec[i].reset(t.name, i, j);
	else if (t.type == SILENT)
	  maxStack.push( boost::make_tuple(t.to, i, j) );
	else if (t.type == BIFURCATE) {
	  maxStack.push( boost::make_tuple(t.toR, k, j) );
	  maxStack.push( boost::make_tuple(t.toL, i, k) );
	}
	else if (t.type == END)
	  ; // do nothing
	else 
	  errorAbort("From maxParse. Transition of unkonwn type. Should not happen...");
      }

      return emitInfoVec;
    }


    vector_t const Grammar::postProbParse(vector<EmitInfo> const & emitInfoVec) const {
      return postProbParse(emitInfoVec, in_, out_);
    }


    vector_t const Grammar::postProbParse(vector<EmitInfo> const & emitInfoVec, vector3D_t const & inTbl, vector3D_t const & outTbl) const {

//      cout << "in:" << endl;
//      printTable(cout, inTbl);
//      cout << endl;
//
//      cout << "out:" << endl;
//      printTable(cout, outTbl);
//      cout << endl;

      map<string, Transition> tm = mkTransitionMap();
      unsigned L = emitInfoVec.size();
      vector_t ppv(L);
      reset(ppv);

      if (inTbl.size() < stateCount)
	errorAbort("From Grammar::postProbParse: inTbl size smaller than expected. Has inside and outside been called?");

      ynumber_t inProb = inTbl[0][0][inLength_]; // inside probabilty of input seq

      // cout << "inProb: " << inProb << endl;

      for (unsigned i = 0; i < L; i++) {
	EmitInfo const & ei = emitInfoVec[i];
	Transition const & t = tm[ei.name];
	ppv[i] = toNumber(transitionUseProb(t, ei.i, ei.j, inTbl, outTbl) / inProb);

	//	cout << "ppv[" << i << "]: " << ppv[i] << endl;

      }

      return ppv;
    }


    vector_t const Grammar::postProbParse(vector<EmitInfo> const & emitInfoVec, map<string, vector<string> > const & equivalenceMap) const
    {
      return postProbParse(emitInfoVec, in_, out_, equivalenceMap);
    }
 
    vector_t const Grammar::postProbParse(vector<EmitInfo> const & emitInfoVec, vector3D_t const & inTbl, vector3D_t const & outTbl, map<string, vector<string> > const & equivalenceMap) const
    {
      map<string, Transition> tm = mkTransitionMap();
      unsigned L = emitInfoVec.size();
      vector_t ppv(L);
      reset(ppv);

      if (inTbl.size() < stateCount)
	errorAbort("From Grammar::postProbParse: inTbl size smaller than expected. Has inside and outside been called?");

      ynumber_t inProb = inTbl[0][0][inLength_]; // inside probabilty of input seq
      for (unsigned i = 0; i < L; i++) {
	EmitInfo const & ei = emitInfoVec[i];

	map<string, vector<string> >::const_iterator it = equivalenceMap.find(ei.name);
	//      if (it == equivalenceMap.end() )
	//	errorAbort("equivalenceMap does not contain name: '" + ei.name + "'.");
	vector<string> const & transNames = it->second;
	ynumber_t s = 0;
	for (unsigned j = 0; j < transNames.size(); j++) {
	  Transition const & t = tm[ transNames[j] ];
	  s += transitionUseProb(t, ei.i, ei.j, inTbl, outTbl) / inProb;
	}
	ppv[i] = toNumber(s);
      }
      return ppv;
    }


    void Grammar::initDpTable(vector3D_t & tbl, unsigned M, unsigned long L) const
    {
      tbl = vector3D_t(M, vector2D_t(L + 1, vector1D_t(L + 1) ) );  // +1 to include the end position of half open intercal including the last position and to include emptu sequences at both ends of the original sequence.
    }


    void Grammar::initDpTable(vector3D_t & tbl) const
    {
      // only resize if needed...
      if (tbl.size() >= stateCount)
	if (stateCount > 0)
	  if (tbl[0].size() >= inLength_ + 1)
	    if (tbl[0][0].size() >= inLength_ + 1)
	      return; // no need to init

      initDpTable(tbl, stateCount, inLength_);
    }


    void Grammar::clearDpTable(vector3D_t & tbl) const
    {
      unsigned M = tbl.size();
      if (M > 0) {
	long L = tbl[0].size();
	clearDpTable(tbl, M, L - 1);
      }
    }


    void Grammar::clearDpTable(vector3D_t & tbl, int M, long L) const
    {
      for (int v = 0; v < M; v++) 
	for (long i = 0; i <= L; i++) 
	  for (long j = 0; j <= L; j++) 
	    tbl[v][i][j] = 0;
    }



    void Grammar::initTbTable(vector3DPair_t & tbl, unsigned M, unsigned long L) const
    {
      // +1 to include the end position of half open interval including the last position and to include empty sequences at both ends of the original sequence.
      tbl = vector3DPair_t(M, vector2DPair_t(L + 1, vector1DPair_t(L + 1, pair<unsigned, unsigned>() ) ) );
    }


    void Grammar::initTbTable(vector3DPair_t & tbl) const
    {
      // only resize if needed...
      if (tbl.size() >= stateCount)
	if (stateCount > 0)
	  if (tbl[0].size() >= inLength_ + 1)
	    if (tbl[0][0].size() >= inLength_ + 1)
	      return; // no need to init

      initTbTable(tbl, stateCount, inLength_);
    }


    map<string, Transition> Grammar::mkTransitionMap() const
    {
      map<string, Transition> tm;
      for (unsigned i = 0; i < outTransitions.size(); i++) 
	for (unsigned j = 0; j < outTransitions[i].size(); j++) {
	  Transition const & t = outTransitions[i][j];
	  tm[t.name] = t;
	}
      return tm;
    }


    unsigned Grammar::stateIndex(string const & name) const
    {
      for (unsigned i = 0; i < outTransitions.size(); i++) 
	if (outTransitions[i][0].from_id == name)
	  return i;
      errorAbort("From stateIndex: No state (non-terminal) with name '" + name + "'.");
      return 0; // to please compiler
    }


    vector<Transition> Grammar::emitTransitions() const
    {
      vector<Transition> v;
      for (unsigned i = 0; i < outTransitions.size(); i++) 
	for (unsigned j = 0; j < outTransitions[i].size(); j++) {
	  Transition const & t = outTransitions[i][j];
	  if ( isEmitType(t.type) )
	    v.push_back(t);
	}
      return v;
    }


    void Grammar::mkInTransitions() {
      inTransitions.clear();
      inTransitions.resize(stateCount);
      for (unsigned i = 0; i < outTransitions.size(); i++)
	for (unsigned j = 0; j < outTransitions[i].size(); j++) {
	  Transition const & t = outTransitions[i][j];
	  if (t.type == BIFURCATE) {
	    inTransitions[t.toL].push_back(t);
	    if (t.toR != t.toL) // incoming transitions only listed once for every state -- essential for the outside algorithm
	      inTransitions[t.toR].push_back(t);
	  }
	  else if (not (t.type == END or t.type == TERMINAL) )  // END and TERMINAL do not have destination states
	    inTransitions[t.to].push_back(t);
	}
    }

    // predeclaration
    ostream &operator<<(ostream & str, Grammar const & g);

    void Grammar::consistencyCheck() const
    {
      for (unsigned i = 0; i < outTransitions.size(); i++) {
	ynumber_t sum = 0;
	for (unsigned j = 0; j < outTransitions[i].size(); j++) 
	  sum += outTransitions[i][j].p;
	if (sum < 0.9999 or sum > 1.0001) {
	  stringstream ss;
	  ss << *this;
	  errorAbort("From Grammar::consistencyCheck: (outgoing) transitions from state " + outTransitions[i][0].from_id + " sum to " + toString(sum) + " instead of 1.0. \n If reading from file, check that the number of transitions is specified correctly. \nPrinting grammar below:\n\n" + ss.str() );
	}
      }
    }


    ////////////////////////////////////////////////////////////////
    // Free functions related to Grammar 
    ////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////
    // AnnoMaps
    ////////////////////////////////////////////////////////////////

    bool EmitAnno::singleCompatible(string const & a, string const & wildCard) const 
    {
      return (a == wildCard or a == anno);
    }

    bool EmitAnno::pairCompatible(string const & l, string const & r, string const & wildCard) const 
    {
      return ( (l == annoLeft or l == wildCard) and (r == annoRight or r == wildCard) );
    }


    void maskEmissions(vector<xvector_t> & singleEmissions, vector<string> const & emitModelNames, string const & anno, EmitAnnoMap const & eam)
    {
      maskEmissions(singleEmissions, emitModelNames, stringToVectorOfStrings(anno), eam);
    }


    void maskEmissions(vector<xvector_t> & singleEmissions, vector<string> const & emitModelNames, vector<string> const & anno, EmitAnnoMap const & eam)
    {
      assert(singleEmissions.size() == emitModelNames.size() );
      for (unsigned i = 0; i < singleEmissions.size(); i++) {
	// get emitAnno
	EmitAnnoMap::const_iterator it = eam.find( emitModelNames[i] );
	if ( it == eam.end() ) 	// check 
	  errorAbort("From maskEmissions: emission id '" + emitModelNames[i] + "' not found in emitAnnoMap.");
	EmitAnno const & emitAnno = it->second;

	maskEmissions(singleEmissions[i], anno, emitAnno);
      }
    }


    void maskEmissions(xvector_t & singleEmissions, vector<string> const & anno, EmitAnno ea)
    {
      assert( singleEmissions.size() == anno.size() );
      for (unsigned long i = 0; i < anno.size(); i++) {
	if (not ea.singleCompatible( anno[i] ) )
	    singleEmissions[i] = 0;
      }
    }


    void maskEmissions(vector<xmatrix_t> & pairEmissions, vector<string> const & emitModelNames, string const & anno, EmitAnnoMap const & eam)
    {
      maskEmissions(pairEmissions, emitModelNames, stringToVectorOfStrings(anno), eam);
    }


    void maskEmissions(vector<xmatrix_t> & pairEmissions, vector<string> const & emitModelNames, vector<string> const & anno, EmitAnnoMap const & eam)
    {
      assert(pairEmissions.size() == emitModelNames.size() );
      for (unsigned i = 0; i < pairEmissions.size(); i++) {
	// get emitAnno
	EmitAnnoMap::const_iterator it = eam.find( emitModelNames[i] );
	if ( it == eam.end() ) 	// check 
	  errorAbort("From maskEmissions: emission id '" + emitModelNames[i] + "' not found in emitAnnoMap.");
	EmitAnno const & emitAnno = it->second;

	maskEmissions(pairEmissions[i], anno, emitAnno);
      }
    }


    void maskEmissions(xmatrix_t & pairEmissions, vector<string> const & anno, EmitAnno ea)
    {
      for (unsigned long j = 1; j < anno.size(); j++) {
	for (unsigned long i = 0; i < j; i++) {
	  if (not ea.pairCompatible(anno[i], anno[j]) )
	    pairEmissions(i, j) = 0;
	}
      }
    }


    TransAnnoMap::TransAnnoMap(EmitAnnoMap const & emitAnnoMap, vector<Transition> const & emitTransVec)
    {
      for (unsigned i = 0; i < emitTransVec.size(); i++) {
	Transition const & t = emitTransVec[i];
	EmitAnnoMap::const_iterator it = emitAnnoMap.find(t.e_id);
	if ( it == emitAnnoMap.end() ) 	// check 
	  errorAbort("From TransAnnoMap: emission id '" + t.e_id + "' not found in emitAnnoMap.");
	EmitAnno const & emitAnno = it->second;
	annoMap_[t.name] = TransAnno(t.name, emitAnno.anno, emitAnno.annoLeft, emitAnno.annoRight);
      }
    }


    string TransAnnoMap::convert(EmitInfo const & emitInfo) const
    {

      map<string, TransAnno>::const_iterator it = annoMap_.find(emitInfo.name);
      if (it == annoMap_.end() )
	errorAbort("TransAnnoMap does not hold annotation for transition of name: '" + emitInfo.name + "'.");

      TransAnno const & ta = it->second;
      if (ta.anno.size() > 0)
	return ta.anno;
      else {
	if (emitInfo.left)
	  return ta.annoLeft;
	else
	  return ta.annoRight;
      }
    }


    string TransAnnoMap::convertToString(vector<EmitInfo> const & emitInfoVec, string const & sep) const
    {
      return toSepString( convertToVector(emitInfoVec), sep);
    }

    vector<string> TransAnnoMap::convertToVector(vector<EmitInfo> const & emitInfoVec) const
    {
      unsigned L = emitInfoVec.size();
      vector<string> v(L);
      for (unsigned i = 0; i < L; i++) 
	v[i] = convert(emitInfoVec[i]);
      return v;
    }


    vector<string> TransAnnoMap::equivalentTransitions(TransAnno const & transAnno) const
    {
      vector<string> transNameVec;
      map<string, TransAnno>::const_iterator it;
      for (it = annoMap_.begin(); it != annoMap_.end(); it++) {
	TransAnno const & ta = it->second;
	if (transAnno.anno == ta.anno and transAnno.annoLeft == ta.annoLeft and transAnno.annoRight == ta.annoRight)
	  transNameVec.push_back(ta.name);
      }    
      return transNameVec;
    }


    vector<string> TransAnnoMap::equivalentTransitions(string const & transName) const
    {
      map<string, TransAnno>::const_iterator it = annoMap_.find(transName);
      if (it == annoMap_.end() )
	errorAbort("From equivalentTransitions: TransAnnoMap does not hold annotation for transition of name: '" + transName + "'.");

      return equivalentTransitions(it->second);
    }

 
    map<string, vector<string> > mkEquivalenceMap(Grammar const & g, TransAnnoMap const & am)
    {
      map<string, vector<string> > equiMap;
      for (unsigned i = 0; i < g.outTransitions.size(); i++)
	for (unsigned j = 0; j < g.outTransitions[i].size(); j++) {
	  Transition const & t = g.outTransitions[i][j];
	  if ( isEmitType(t.type) )
	    equiMap[t.name] = am.equivalentTransitions(t.name);
	}
      return equiMap;
    }


  } // end namespace grammar

} // end namespace phy
