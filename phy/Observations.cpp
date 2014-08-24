/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/Observations.h"
#include <sstream>

namespace phy {


  /****************************************
   TODO Split into another header
   ****************************************/

  /** Implementations of StateMaps following the PIMPL idiom **/
  class StateMapImplContinuous : public StateMapImpl {
  public:
    StateMapImplContinuous( vector_t const & breakpoints) : breakpoints_(breakpoints), states_(vector<state_t>(breakpoints.size() + 1)) { init(); }
    
    virtual state_t const & symbol2State(symbol_t const & s) const;
    virtual symbol_t const & state2Symbol(state_t i) const;
    virtual unsigned stateCount() const { return stateCount_; }
    virtual unsigned metaStateCount() const { return stateCount_; }
    virtual state_t symbolSize() const{ return 1;}
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const;
    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask_) const;
  private:
    void init();
    vector_t breakpoints_;
    vector<state_t> states_; //TODO Check if this can be avoided
    vector<stateMask_t> metaState2StateMask_; //TODO make dynamic

    //TODO Remove this:
    symbol_t state2Symbol_;
    unsigned stateCount_;

    //TODO Definitely improve(remove)
    static vector<symbol_t> dummyVecSymbol;
  };
  vector<symbol_t> StateMapImplContinuous::dummyVecSymbol = vector<symbol_t>(1,"");

  class StateMapImplSymbol : public StateMapImpl {
  public:
    StateMapImplSymbol(string const & symbols) : state2Symbol_( stringToVectorOfStrings(symbols) ) { init(); }
    StateMapImplSymbol(vector<symbol_t> const & symbols) : state2Symbol_(symbols) { init(); }
    StateMapImplSymbol(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap) : state2Symbol_(symbols), degeneracyMap_(metaSymbolDegeneracyMap){ init(); }
    StateMapImplSymbol(StateMap const & staMap, unsigned n);

    virtual state_t const & symbol2State(symbol_t const & s) const;
    virtual symbol_t const & state2Symbol(state_t i) const { return state2Symbol_[i]; }
    virtual unsigned stateCount() const { return stateCount_; } //Move to StateMapImpl?
    virtual unsigned metaStateCount() const { return metaStateCount_; }
    virtual state_t symbolSize() const { return symbolSize_; } //TODO Move to StateMapImpl
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const;
    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask_) const;

  private:
    void init();

    vector<symbol_t> state2Symbol_;
    boost::unordered_map<symbol_t, state_t> symbol2State_;
    boost::unordered_map<symbol_t, vector<symbol_t> > degeneracyMap_;
    unsigned stateCount_;
    unsigned metaStateCount_;
    unsigned symbolSize_;
  };

  StateMapImplSymbol::StateMapImplSymbol(StateMap const & staMap, unsigned n){
    state2Symbol_ = mkMultiStateSymbols(staMap, n);
    boost::unordered_map<symbol_t, vector<symbol_t> > degMap;
    for (unsigned i = staMap.stateCount(); i < staMap.metaStateCount(); i++) {
      symbol_t sym = staMap.state2Symbol(i);
      degMap[sym] = staMap.degeneracyVector(sym);
    }
    degeneracyMap_ = mkMultiStateSymbolDegeneracyMap(degMap, staMap, n);

    init();
  }

  void StateMapImplSymbol::init(){
    stateCount_ = state2Symbol_.size();
    // add basic symbols to degeneracy map
    BOOST_FOREACH(symbol_t const & sym,  state2Symbol_)
      degeneracyMap_[sym] = vector<symbol_t>(1, sym);
    // add (only) metaSymbols to state2Symbols_
    symbol_t sym;
    vector<symbol_t> degSymVec;
    BOOST_FOREACH(boost::tie(sym, degSymVec), degeneracyMap_)
      if ( find(state2Symbol_.begin(), state2Symbol_.end(), sym) == state2Symbol_.end() ) // not found
	state2Symbol_.push_back(sym);
    metaStateCount_ = state2Symbol_.size();

    // symbol2State
    for (unsigned i = 0; i < metaStateCount_; i++)
      symbol2State_[ state2Symbol_[i] ] = i;
    //symbolSize
    symbolSize_= (stateCount_ > 0) ? state2Symbol_[0].size() : 0;
  }

  void StateMapImplSymbol::setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const{
    //Checks
    if(metaState2StateMask.size() != metaStateCount_)
      errorAbort("StateMapImplSymbol: metaState2StateMask_.size() != metaStateCount_");

    for (unsigned i = 0; i < metaStateCount_; i++) {
      symbol_t sym = state2Symbol(i);
      vector<symbol_t> const & v = degeneracyVector(sym);
      for (vector<symbol_t>::const_iterator it = v.begin(); it < v.end(); it++) {
	state_t j = symbol2State(*it);
	if (j >= stateCount_)
	  errorAbort("StateMaskMap::StateMaskMap: Degenerate symbol '" + sym + "' is defined in terms of another meta-symbol '" + *it + "'.");
	metaState2StateMask.at(i)(j) = true;
      }
    }
  }

  void StateMapImplContinuous::init(){
    stateCount_ = states_.size();
    for(int i = 0; i < states_.size(); ++i){
      states_.at(i) = i;
    }

    //TODO Check if neccesary i.e. who calls state2Symbol?
    state2Symbol_ = "h";
  }

  void StateMapImplContinuous::setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const {
    if(metaState2StateMask.size() != stateCount_)
      errorAbort("StateMapImplSymbol: metaState2StateMask_.size() != stateCount_");

    for(int i = 0 ; i < stateCount_ ; ++i ){
      metaState2StateMask.at(i)(i) = true;
    }    
  }

  state_t const & StateMapImplSymbol::symbol2State(symbol_t const & s) const {
    boost::unordered_map<symbol_t, state_t >::const_iterator it = symbol2State_.find(s);
    if ( it == symbol2State_.end() )
      errorAbort("StateMap::symbol2State: Symbol '" + s + "' not found in stateMap.");
    return it->second;
  }

  state_t const & StateMapImplContinuous::symbol2State(symbol_t const & s) const {
    double sym;
    try{
      sym = boost::lexical_cast<double>(s);
    }
    catch( boost::bad_lexical_cast &){
      errorAbort("Statemap::symbol2State: Symbol '" + s + "' could not be converted to double");
    }
    //TODO Binary search?
    for(state_t i = 0; i < breakpoints_.size(); ++i){
      if(sym < breakpoints_(i)) return states_.at(i);
    }
    return states_.at(breakpoints_.size());
  }

  symbol_t const & StateMapImplContinuous::state2Symbol(state_t i) const {
    //TODO have a function that takes a reference to a string stream. And append the following
    /*
    std::stringstream s;
    s << '(';
    if(i == 0)
      s << '-Inf';
    else
      s << breakpoints_(i-1) ;
    s << ';';
    if(i < breakpoints_.size())
      s << breakpoints_(i);
    else
      s << 'Inf';
    s << ')';
    return s.str();
    */
    return state2Symbol_;
  }

  vector<symbol_t> const StateMapImplContinuous::degeneracyVector(symbol_t const & s) const{
    //TODO Refactor this function away
    return dummyVecSymbol;
  }

  vector<symbol_t> const StateMapImplSymbol::degeneracyVector(symbol_t const & s) const{
    boost::unordered_map<symbol_t, vector<symbol_t> >::const_iterator it = degeneracyMap_.find(s);
    if ( it == degeneracyMap_.end() )
      errorAbort("StateMap::degeneracyVector: Symbol '" + s + "' not found in stateMap.");
    return it->second;
  }

  /*****************************************
     StateMap stuff
   *****************************************/

  StateMap const & StateMap::operator=(StateMap const &rhs)
  {
    name_ = rhs.name_;
    pImpl = rhs.pImpl;
    return *this;
  }

  /* StateMap constructors */
  /** Constructor */
  StateMap::StateMap(string const & symbols, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols)) { }

  /** Constructor */
  StateMap::StateMap(vector<symbol_t> const & symbols, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols)) { }

  /** Constructor */
  StateMap::StateMap(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols, metaSymbolDegeneracyMap)) { }
 
  /** Constructor for continuous type StateMap */
  StateMap::StateMap( vector_t const & breakpoints, string const & name) : name_(name), pImpl(new StateMapImplContinuous(breakpoints)) { }

  StateMap::StateMap(StateMap const & staMap, unsigned n, string const & explicitName) : name_(explicitName), pImpl( new StateMapImplSymbol(staMap, n) )
  {
    if (name_.size() == 0)
      if ( staMap.name().size() )
	name_ = toString(n) + "-" + staMap.name();
  }

  vector<symbol_t> const StateMap::state2Symbol(vector<state_t> v) const
  {
    vector<symbol_t> u;
    u.reserve( v.size() );
    for (unsigned i = 0; i < v.size(); i++) {
      u.push_back( state2Symbol(v[i]) );
    }
    return u;
  }

  vector<state_t> const StateMap::symbol2State(vector<symbol_t> v) const
  {
    vector<state_t> u( v.size() );
    for (unsigned i = 0; i < v.size(); i++) {
      u.push_back( symbol2State(v[i]) );
    }
    return u;
  }

  string mkMultiStateSymbol(unsigned state, StateMap const & sm, unsigned n)
  {
    symbol_t symbol;
    unsigned base = sm.stateCount();
    for (int j = n - 1; j >= 0; j--) {
      unsigned val = power(base, j);
      state_t subState = state / val;
      symbol += sm.state2Symbol(subState);
      state = state % val;
    }
    return symbol;
  }


  // returns a symbol string of all possible n-long symbols constructed from the StateMap symbols
  vector<symbol_t> mkMultiStateSymbols(StateMap const & sm, unsigned n)
  {
    unsigned stateCount = power(sm.stateCount(), n);
    vector<symbol_t> symbols;
    symbols.reserve(stateCount);
    for (unsigned i = 0; i < stateCount; i++)
      symbols.push_back( mkMultiStateSymbol(i, sm, n) );
    return symbols;
  }
    

  void addBasicSymbolsToDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > & degMap, vector<symbol_t> const & symbols)
  {
    BOOST_FOREACH(symbol_t const & sym,  symbols)
      degMap[sym] = vector<symbol_t>(1, sym);
  }


  void addBasicSymbolsToDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > & degMap, StateMap const & staMap)
  {
    vector<symbol_t> symbols;
    for (unsigned i = 0; i < staMap.stateCount(); i++) 
      symbols.push_back( staMap.state2Symbol(i) );
    addBasicSymbolsToDegeneracyMap(degMap, symbols);
  }


  // take a degeneracy map and the corresponding state map and output a degeneracy map for the corresponding multi state symbol set.
  boost::unordered_map<symbol_t, vector<symbol_t> > mkMultiStateSymbolDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > const & degMap, StateMap const & staMap, unsigned n)
  {
    boost::unordered_map<symbol_t, vector<symbol_t> > extDegMap = degMap;
    addBasicSymbolsToDegeneracyMap(extDegMap, staMap);
    boost::unordered_map<symbol_t, vector<symbol_t> > mulDegMap = extDegMap;

    if (power(static_cast<unsigned int>(extDegMap.size()), static_cast<unsigned int>(n)) > 100000) // 
      errorAbort("Silly number of degenerate symbols in multi-symbol degeneracy map. Reimplement 'mkMultiStateDegeneracyMa' and related functions functions");

    // some loop internal variables
    boost::unordered_map<symbol_t, vector<symbol_t> > tmpMulDegMap;
    symbol_t preSym, appSym;
    vector<symbol_t> preDegSymVec, appDegSymVec;
    for (unsigned i = 1; i < n; i ++) {
      tmpMulDegMap.clear();
      BOOST_FOREACH(boost::tie(preSym, preDegSymVec), mulDegMap) {
	BOOST_FOREACH(boost::tie(appSym, appDegSymVec), extDegMap) {
	  symbol_t sym = preSym + appSym;
	  vector<symbol_t> degSymVec;
	  BOOST_FOREACH(symbol_t const & preDegSym, preDegSymVec) 
	    BOOST_FOREACH(symbol_t const & appDegSym, appDegSymVec) 
	    degSymVec.push_back(preDegSym + appDegSym);
	  tmpMulDegMap[sym] = degSymVec;
	}
      }
      mulDegMap = tmpMulDegMap;
    }
    return mulDegMap;
  }

  StateMap mkNucStateMap(void)
  {
    vector<symbol_t> symbols = stringToVectorOfStrings("ACGT");
    return StateMap(symbols, "nuc");
  }


  StateMap mkMetaNucStateMap(void)
  {
    vector<symbol_t> symbols = stringToVectorOfStrings("ACGT");
    boost::unordered_map<symbol_t, vector<symbol_t> > degMap;
    
    degMap["U"] = stringToVectorOfStrings("T"); 
    degMap["R"] = stringToVectorOfStrings("AG");
    degMap["Y"] = stringToVectorOfStrings("CT");
    degMap["K"] = stringToVectorOfStrings("GT");   
    degMap["M"] = stringToVectorOfStrings("CA");   
    degMap["S"] = stringToVectorOfStrings("CG");   
    degMap["W"] = stringToVectorOfStrings("AT");   
    degMap["B"] = stringToVectorOfStrings("CGT");  
    degMap["D"] = stringToVectorOfStrings("AGT");  
    degMap["H"] = stringToVectorOfStrings("ACT");  
    degMap["V"] = stringToVectorOfStrings("ACG");  
    degMap["N"] = stringToVectorOfStrings("ACGT");
    degMap["-"] = stringToVectorOfStrings("ACGT");
    degMap["."] = stringToVectorOfStrings("ACGT");

    return StateMap(symbols, degMap, "metaNuc");
  }


  // creates stateMap where each symbol consists of n nucleotides.
  StateMap mkMultiMetaNucStateMap(unsigned n)
  {
    StateMap metaNucMap = mkMetaNucStateMap();
    return StateMap(metaNucMap, n);
  }


  // creates stateMap where each symbol consists of n nucleotides.
  StateMap mkMultiNucStateMap(unsigned n)
  {
    StateMap nucMap = mkNucStateMap();
    return StateMap(nucMap, n);
  }


  StateMaskMap::StateMaskMap(StateMap const & staMap)
    : metaState2StateMask_( staMap.metaStateCount() , stateMask_t(staMap.stateCount(), false) )      
  {
    staMap.setMetaState2StateMask(metaState2StateMask_);
  }


  StateMaskMapSet::StateMaskMapSet(StateMap const & staMap, unsigned n) : stateMapCount_(n), stateMapVector_(n), stateMaskMapVector_(n)
  {
    StateMapPtr_t staMapPtr(new StateMap(staMap));
    StateMaskMapPtr_t staMasPtr(new StateMaskMap(staMap));

    for (unsigned i = 0; i < n ;i++) {
      stateMapVector_[i] = staMapPtr;
      stateMaskMapVector_[i] = staMasPtr;
    }
  }


  StateMaskMapSet::StateMaskMapSet(vector<StateMapPtr_t> const & stateMapVector) : stateMapCount_( stateMapVector.size() ), stateMapVector_(stateMapVector), stateMaskMapVector_(stateMapVector.size() ) 
  {
    unsigned n = stateMapCount_;
    for (unsigned i = 0; i < n ;i++) 
      if (stateMapVector_[i] != NULL)
	stateMaskMapVector_[i] = StateMaskMapPtr_t(new StateMaskMap(*stateMapVector_[i]));
      else {
	stateMaskMapVector_[i] = StateMaskMapPtr_t();
	//errorAbort("from StateMaskMapSet::StateMaskMapSet: Null StateMapPtr for index " + toString(i) + ". This should not happen. ");
      }
  }

  void StateMaskMapSet::symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec) const
  {
    assert(symVec.size() == stateMapCount_);

    for (unsigned i = 0; i < stateMapCount_; i++) {
      unsigned state = stateMapVector_[i]->symbol2State( symVec[i] );
      obsMasVec[i] = & stateMaskMapVector_[i]->metaState2StateMask(state);
    }
  }


  stateMaskVec_t StateMaskMapSet::symbols2StateMasks(vector<symbol_t> const & symVec) const
  {
    stateMaskVec_t obsMasVec(stateMapCount_);
    symbols2StateMasks(obsMasVec, symVec);
    return obsMasVec;
  }


  void StateMaskMapSet::symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec, vector<unsigned> const & seqToVarMap) const
  {
    assert( symVec.size() == seqToVarMap.size() );
    assert( obsMasVec.size() == stateMapCount_ );

    for (unsigned i = 0; i < stateMapCount_; i++)
      obsMasVec[i] = NULL;
    for (unsigned i = 0; i < symVec.size(); i++) {
      unsigned const & ranVarIdx = seqToVarMap[i];
      state_t state = stateMapVector_[ranVarIdx]->symbol2State( symVec[i] );
      obsMasVec[ ranVarIdx ] = & stateMaskMapVector_[ranVarIdx]->metaState2StateMask(state);
    }
  }


  stateMaskVec_t StateMaskMapSet::symbols2StateMasks(vector<symbol_t> const & symVec, vector<unsigned> const & seqToVarMap) const
  {
    stateMaskVec_t obsMasVec(stateMapCount_);
    symbols2StateMasks(obsMasVec, symVec, seqToVarMap);
    return obsMasVec;
  }


  SeqToVarSymbol::SeqToVarSymbol(string varName, string seqName, unsigned size, int symbolOffset, int indexOffset, bool optional, unsigned rSize, int rIndexOffset)
    : varName(varName), seqName(seqName), size(size), symbolOffset(symbolOffset), indexOffset(indexOffset), optional(optional), rSize(rSize), rIndexOffset(rIndexOffset)
  {
    assert(varName.size() != 0);
    assert(seqName.size() != 0);
    
    unsigned DEFAULT_VALUE_UNSIGNED = numeric_limits<unsigned>::max();
    int DEFAULT_VALUE_INT = numeric_limits<int>::max();

    if ( symbolOffset == DEFAULT_VALUE_INT )
      SeqToVarSymbol::symbolOffset = size;
    if ( rSize == DEFAULT_VALUE_UNSIGNED )
      SeqToVarSymbol::rSize = size;
    if (rIndexOffset == DEFAULT_VALUE_INT)
      SeqToVarSymbol::rIndexOffset = -1 * indexOffset;
  }


  void mkStvToSeqMap(vector<int> & stvToSeqMap, vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames)
  {
    assert(stvToSeqMap.size() == stvVec.size() );
    for (unsigned i = 0; i < stvVec.size(); i++) {
      string const & seqName = stvVec[i].seqName;
      if ( hasElement(seqNames, seqName) )
	stvToSeqMap[i] = getIndex(seqNames, seqName);
      else
	stvToSeqMap[i] = -1;
    }
  }


  vector<int> mkStvToSeqMap(vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames)
  {
    vector<int> stvToSeqMap( stvVec.size() );
    mkStvToSeqMap(stvToSeqMap, stvVec, seqNames);
    return stvToSeqMap;
  }


  vector<unsigned> mkSeqToVarMap(vector<string> const & varNames, vector<string> const & seqNames)
  {
    return mkSubsetMap(varNames, seqNames);
  }


  // preconditions: 
  // strVec contains all input sequences, which are assumed of equal length
  // stateMask2DVec the same length as the input sequences. 
  void mkStateMask2DVec(vector<string> const & strVec, stateMask2DVec_t & stateMask2DVec, vector<unsigned> const & seqToVarMap, StateMaskMapSet const & stateMaskMapSet, unsigned offSet, char missingDataChar)
  {
    long seqSize = strVec[0].size();
    unsigned long symCount = symbolCount(strVec, offSet);
    assert( symCount == stateMask2DVec.size() );
    unsigned symVecSize = strVec.size();

    // get symbol sizes
    vector<unsigned> symSizeVec(symVecSize);
    for (unsigned i = 0; i < symVecSize; i++)
      symSizeVec[i] = stateMaskMapSet.symbolSize( seqToVarMap[i] );
    
    // setup symbol vector
    vector<symbol_t> symVec(symVecSize);
    for (unsigned i = 0; i < symVecSize; i++)
      symVec[i].resize( symSizeVec[i] );

    // fill in stateMask table
    for (unsigned i = 0; i < symCount; i++) {
      mkSymbolVector(symVec, symSizeVec, strVec, seqSize, i * offSet, missingDataChar);
      stateMaskMapSet.symbols2StateMasks(stateMask2DVec[i], symVec, seqToVarMap);
    }      
  }


  stateMask2DVec_t mkStateMask2DVec(vector<string> const & strVec, vector<unsigned> const & seqToVarMap, unsigned varCount, StateMaskMapSet const & stateMaskMapSet, unsigned offSet, char missingDataChar)
  {
    assert(strVec.size() > 0);
    stateMask2DVec_t stateMask2DVec = initStateMask2DVec(symbolCount(strVec, offSet), varCount);
    mkStateMask2DVec(strVec, stateMask2DVec, seqToVarMap, stateMaskMapSet, offSet, missingDataChar);
    return stateMask2DVec;
  }


  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, unsigned offSet, char missingDataChar)
  {
    vector<unsigned> seqToVarMap = mkSeqToVarMap(varNames, seqData.seqNames);
    unsigned varCount = varNames.size();
    return mkStateMask2DVec(seqData, seqToVarMap, varCount, stateMaskMapSet, offSet, missingDataChar);
  }


  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<unsigned> const & seqToVarMap, unsigned varCount, StateMaskMapSet const & stateMaskMapSet, unsigned offSet, char missingDataChar)
  {
    return mkStateMask2DVec(seqData.sequences, seqToVarMap, varCount, stateMaskMapSet, offSet, missingDataChar);
  }


  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, unsigned offSet, char missingDataChar)
  {
    vector<stateMask2DVec_t> resultVec;
    for (unsigned i = 0; i < seqDataVec.size(); i++)
      resultVec.push_back( mkStateMask2DVec(seqDataVec[i], varNames, stateMaskMapSet, offSet, missingDataChar) );
    return resultVec;
  }


  bool onlyOptionalSeqMissingOrAbort(vector<SeqToVarSymbol> const & stvVec, vector<int> const & stvToSeqMap, string const & seqDataEntryId)
  {
    assert( stvVec.size() == stvToSeqMap.size() );
    for (unsigned i = 0; i < stvVec.size(); i++) {
      SeqToVarSymbol const & stv = stvVec[i];
      if (stvToSeqMap[i] < 0 and (not stv.optional) ) // missing and required
	errorAbort(string("from checkOnlyOptionalSeqMissing: SeqData entry '" + seqDataEntryId 
			  + "' lacks sequence of name, '") + stv.seqName + "' requred to make symbol for random variable, '" + stv.varName + "'.");
    }
    return true;
  }

  long sameSymbolCountOrAbort(vector<SeqToVarSymbol> const & stvVec, vector<int> const & stvToSeqMap, vector<long> const & seqSizes, string const & seqDataEntryId)
  {
    assert( stvVec.size() == stvToSeqMap.size() );

    // first collect names symbol lengths of all sequences.
    vector<string> seqNames;
    vector<long> symCounts;
    for (unsigned i = 0; i < stvVec.size(); i++)
      if (stvToSeqMap[i] >= 0) {
	SeqToVarSymbol const & stv = stvVec[i];
	seqNames.push_back(stv.seqName);
	symCounts.push_back( symbolCount(seqSizes[ stvToSeqMap[i] ], stv.symbolOffset) );
      }
    
    // then check sizes and report any problems
    if (symCounts.size() == 0)
      return 0;
    // else
    long symCount = symCounts[0];
    for (unsigned i = 0; i < symCounts.size(); i++) {
      if (symCounts[i] != symCount) // problem...
	errorAbort( string("From sameSymbolCountOrAbort: Sequences in seqData with entry '") + seqDataEntryId 
		    + "' do not have the same symbol count, based on the given symbol offset.\n"
		    + "Note that sequences may be used multiple times by different variables and can therefore be repeated below.\n"
		    + "Names of evaluated sequences:\n"
		    + toString(seqNames) + "\n"
		    + "symCounts of corresponding sequences:\n"
		    + toString(symCounts) + "\n");
	  }
    return symCount;
  }


  void assertSymbolSizeConsistency(vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, bool isPair)
  {
    // map from the observed randome variables / stv to total set of random variables 
    vector<unsigned> stvToVarMap = mkStvToVarMap(stvVec, varNames);

    // assert symbol sizes
    unsigned n = stvVec.size();
    for (unsigned i = 0; i < n; i++) 
      if (isPair)
	assert( (stvVec[i].size + stvVec[i].rSize) == stateMaskMapSet.symbolSize( stvToVarMap[i] ) );
      else // is single
	assert( stvVec[i].size == stateMaskMapSet.symbolSize( stvToVarMap[i] ) );
  }


  vector<long> getSeqSizes(SeqData const & seqData)
  {
    unsigned seqCount = seqData.sequences.size();
    vector<long> seqSizes(seqCount);
    for (unsigned i = 0; i < seqCount; i++)
      seqSizes[i] = seqData.sequences[i].size();

    return seqSizes;
  }

  // versions that take vector<SeqToVarSymbol> as input
  void  mkStateMask2DVec(stateMask2DVec_t & stateMask2DVec, SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<unsigned> const & stvToVarMap, StateMaskMapSet const & stateMaskMapSet, char missingDataChar)
  {
    // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    onlyOptionalSeqMissingOrAbort(stvVec, stvToSeqMap, seqData.entryId);
    
    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<long> seqSizes = getSeqSizes(seqData);
    long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);

    // assert symbol sizes
    unsigned symVecSize = stvVec.size();
    for (unsigned i = 0; i < symVecSize; i++)
      assert( stvVec[i].size == stateMaskMapSet.symbolSize( stvToVarMap[i] ) );

    // setup symbol vector
    vector<symbol_t> symVec(symVecSize);
    for (unsigned i = 0; i < symVecSize; i++)
      symVec[i].resize( stvVec[i].size );

    // fill in stateMask table
    for (long i = 0; i < symCount; i++) {
      mkSymbolVector(symVec, seqData.sequences, seqSizes, stvVec, stvToSeqMap, i, missingDataChar);
      stateMaskMapSet.symbols2StateMasks(stateMask2DVec[i], symVec, stvToVarMap);
    }
  }

  vector<unsigned> mkStvToVarMap(vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames)
  {
    // count and collect names of observed random variables given by the stv's
    unsigned stvCount = stvVec.size();
    vector<string> stvVarNames(stvCount);
    for (unsigned i = 0; i < stvCount; i++)
      stvVarNames[i] = stvVec[i].varName;

    // map from the observed randome variables / stv to total set of random variables 
    return mkSubsetMap(varNames, stvVarNames);
  }


  void mkStateMask2DVec(stateMask2DVec_t & stateMask2DVec, SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, char missingDataChar)
  {
    // map from the observed randome variables / stv to total set of random variables 
    vector<unsigned> stvToVarMap = mkStvToVarMap(stvVec, varNames);
    mkStateMask2DVec(stateMask2DVec, seqData, stvVec, stvToVarMap, stateMaskMapSet, missingDataChar);
  }


  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, char missingDataChar)
  {
    // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    onlyOptionalSeqMissingOrAbort(stvVec, stvToSeqMap, seqData.entryId);

    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<long> seqSizes = getSeqSizes(seqData);
    long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);
    stateMask2DVec_t stateMask2DVec = initStateMask2DVec( symCount, stateMaskMapSet.stateMapCount() );
    mkStateMask2DVec(stateMask2DVec, seqData, stvVec, varNames, stateMaskMapSet, missingDataChar);
    return stateMask2DVec;
  }


  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<unsigned> const & stvToVarMap, StateMaskMapSet const & stateMaskMapSet, char missingDataChar)
  {
    // map from stv to seq indexes. If sequence is missing for stvVec[i], then stvToSeqMap[i] = -1.
    vector<int> stvToSeqMap = mkStvToSeqMap(stvVec, seqData.seqNames);
    onlyOptionalSeqMissingOrAbort(stvVec, stvToSeqMap, seqData.entryId);

    // get sequence sizes and check all have same symbol-lengths / symbol-offset sizes
    vector<long> seqSizes = getSeqSizes(seqData);
    long symCount = sameSymbolCountOrAbort(stvVec, stvToSeqMap, seqSizes, seqData.entryId);
    stateMask2DVec_t stateMask2DVec = initStateMask2DVec( symCount, stateMaskMapSet.stateMapCount() );
    mkStateMask2DVec(stateMask2DVec, seqData, stvVec, stvToVarMap, stateMaskMapSet, missingDataChar);
    return stateMask2DVec;
  }


  unsigned long symbolCount(vector<string> const & strVec, unsigned offSet)
  {
    assert(offSet > 0);
    long seqSize = (strVec.size() > 0) ? strVec[0].size() : 0;
    unsigned long symCount = (seqSize + (offSet - 1) ) / offSet;  // number of symbols in input sequences. The result should be ceil for integer division.
    return symCount;
  }


  long unsigned symbolCount(long seqSize, unsigned offSet)
  {
    assert(offSet > 0);
    long unsigned symCount = (seqSize + (offSet - 1) ) / offSet;  // number of symbols in input sequences. The result should be ceil for integer division.
    return symCount;
  }

  vector<string> initSymVec(vector<SeqToVarSymbol> const & stvVec, bool isPair)
  {
    unsigned n = stvVec.size();
    vector<symbol_t> symVec(n);
    for (unsigned i = 0; i < n; i++)
      if (isPair)
	symVec[i].resize( stvVec[i].size + stvVec[i].rSize);
      else
	symVec[i].resize( stvVec[i].size );

    return symVec;
  }



  
} // end namespace phy
