/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Observations_h
#define __Observations_h

#include "phy/utils.h"
#include "phy/PhyDef.h"
#include "phy/SeqIO.h"
#include "boost/tuple/tuple.hpp"
#include "boost/foreach.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp> 
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>

using namespace std;

namespace phy {

  ////////////////////////////////////////////////////////////////
  // StateMaps: symbol_t -> state_t
  ////////////////////////////////////////////////////////////////


  /** Map input string (symbol_t) symbols to integer observation
      space. Degenerate observation symbols that map to several states
      are allowed. These symbols are mapped to metaStates with values
      greater than the basic states. Multiple symbols may map to the
      same set of states, however, they are assigned different
      metaStates internally.*/

  class StateMapImpl;

  class StateMapImpl{
  public:
    //    virtual ~StateMapImpl();
    virtual state_t const & symbol2State(symbol_t const & s) const = 0;
    //    virtual symbol_t const & state2Symbol(state_t i) const = 0;
    //virtual const state2Symbol(state_t i, symbol_t & s) const = 0;
  };
  class StateMapImplContinuous;
  class StateMapImplSymbol;
  typedef boost::shared_ptr<StateMapImpl> StateMapImplPtr_t;
  
  class StateMap {
  public:
    
    /** Constructor */
    StateMap(string const & symbols, string const & name = "");

    /** Constructor */
    StateMap(vector<symbol_t> const & symbols, string const & name = "") ;

    /** Constructor */
    StateMap(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap, string const & name = "");
 
    /** Constructor of n-long multiSymbols based on a basic stateMap. If staMap has a name, and no explicit name is given, the name of the constructed StateMap will be n-name.*/
    StateMap(StateMap const & staMap, unsigned n, string const & explicitName = "");

    /** Constructor for continuous type StateMap */
    StateMap( vector_t const & breakpoints, string const & name = "");
 
    /** Copy assignment operator. */
    StateMap const & operator=(StateMap const &rhs);

    /** Returns symbol corresponding to state i */
    symbol_t const & state2Symbol(state_t i) const {return state2Symbol_[i];}
    vector<symbol_t> const state2Symbol(vector<state_t> v) const;

    /** Returns state corresponding to symbol s. Aborts on nonexisting symbols. */
    state_t const & symbol2State(symbol_t const & s) const ;
    vector<state_t> const symbol2State(vector<symbol_t> v) const;

    /** returns degeneracy vector for symbol s */
    vector<symbol_t> const degeneracyVector(symbol_t const & s) const;

    /** returns number of basic states */
    unsigned stateCount() const {return stateCount_;}

    /** returns number of meta states (includes basic states) */
    unsigned metaStateCount() const {return metaStateCount_;}

    /** returns length (in chars) of symbols */
    state_t symbolSize() const  {return symbolSize_;}

    /** Returns StateMap name. Name is empty ("") if not a canonical
	type. Useful in IO-functions. */
    string const & name() const {return name_;}
    bool const & isContinuous() const{ return isCont_;}
  protected:

    /** setup data structures */
    void init();
    void initCont();

    vector<symbol_t> state2Symbol_;
    boost::unordered_map<symbol_t, state_t> symbol2State_;
    boost::unordered_map<symbol_t, vector<symbol_t> > degeneracyMap_;
    unsigned stateCount_;
    unsigned metaStateCount_;
    unsigned symbolSize_;
    string name_;

    //TODO I would rather have an interface class with discrete and continous implementations
    bool isCont_; //TODO initialize to false in existing constructors
    vector_t breakpoints_;
    vector<state_t> states_;
  private:
    StateMapImplPtr_t pImpl;
  };

  /** returns an n-state symbol with state 'state' based on stateMap. */
  symbol_t mkMultiStateSymbol(unsigned state, StateMap const & sm, unsigned n);

  /** returns a vector of all n-state symbols composed from single state symbols. Used to create multi nucleotide stateMaps. */
  vector<symbol_t> mkMultiStateSymbols(StateMap const & sm, unsigned n);

  /** Construct the degeneracy map for n-state symbols baed on the degeneracy map given in the input StateMap (staMap). */
  boost::unordered_map<symbol_t, vector<symbol_t> > mkMultiStateSymbolDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > const & degMap, StateMap const & staMap, unsigned n);


  ////////////////////////////////////////////////////////////////
  // Nucleotide specific StateMaps
  ////////////////////////////////////////////////////////////////
   
  /** Return a symbol map for basic nucleotides */
  StateMap mkNucStateMap();

  /** Return a symbol map for basic nucleotides and IUPAC symbols. Can be used with DNA as well as RNA. */
  StateMap mkMetaNucStateMap();

  /** Returns a stateMap where each symbol consists of n nucleotides. */
  StateMap mkMultiNucStateMap(unsigned n);

  /** Returns a stateMap where each symbol consists of n nucleotides. Include degenerate symbols and should not be used for more than n=3.*/
  StateMap mkMultiMetaNucStateMap(unsigned n);


  ////////////////////////////////////////////////////////////////
  // StateMaskMaps: state_t -> stateMask_t
  ////////////////////////////////////////////////////////////////

  /* StateMaskMap allows conversion from (meta)states to stateMask_t,
     which are vectors of bools denoting if a state is observed or
     not. A stateMask_t is the natural input data structure for, e.g.,
     discrete factor graphs. **/
  class StateMaskMap {
  public:

    //TODO: could allow dynamic calculation of stateMasks to save space Especially for continous factors!
    StateMaskMap(StateMap const & staMap);

    /** Return the bit vector corresponding to state i */
    stateMask_t const & metaState2StateMask(state_t i) const {return metaState2StateMask_[i];}

  protected:
    bool isCont_;
    vector<stateMask_t> metaState2StateMask_;
  };


  ////////////////////////////////////////////////////////////////
  // StateMaskMapSet: vector<symbol_t> -> stateMaskVec_t
  ////////////////////////////////////////////////////////////////

  typedef boost::shared_ptr<StateMaskMap> StateMaskMapPtr_t;
  typedef boost::shared_ptr<StateMap> StateMapPtr_t;

  /** Class that holds all the stateMasks for a given indexed set of
      random variables. Since the member functions return pointers to
      internally stored stateMasks, which are normally used in
      calculations on a factorGraph or similar, this datastructure
      must exist as long as calculations are performed or as long as
      the stateMasks are needed. */
  class StateMaskMapSet {
  public:

    /** Constructor to use when all random variables are over the same state space (i.e, all use the same StateMap). */
    StateMaskMapSet(StateMap const & staMap, unsigned n);

    /** Constructor to use when random variables are defined over different state spaces (i.e., they use different StateMaps). If StateMap is not defined for a random variable, a NULL pointer can be given.*/
    StateMaskMapSet(vector<StateMapPtr_t> const & stateMapVector);

    /** Convert a vector of symbols to a vector of observation
	masks. The indexing of symVec corresponds to the indexing of
	StateMaskMaps in StateMaskMapSet (which should correspond to the
	indexing of random variables in, e.g., a factor graph). The
	returned vector will have a stateMask pointer for each random
	variable. For unobserved random variables the version below
	can be used. */
    void symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec) const;
    stateMaskVec_t symbols2StateMasks(vector<symbol_t> const & symVec) const;

    /** Same as above but allows for unobserved variables. The randomVariableMap defines which random variable (and thus stateMask index) each index of symVec pertains to. obsMasVec[i] is set to null (missing data) if not pointed to by varMap. */
    stateMaskVec_t symbols2StateMasks(vector<symbol_t> const & symVec, vector<unsigned> const & varMap) const;
    void symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec, vector<unsigned> const & varMap) const;

    unsigned stateMapCount(void) const {return stateMapCount_;}
    unsigned symbolSize(unsigned idx) const {return stateMapVector_[idx]->symbolSize();}
    
  protected:
   
    unsigned stateMapCount_;
    vector<StateMapPtr_t> stateMapVector_;
    vector<StateMaskMapPtr_t> stateMaskMapVector_;
   };

  ////////////////////////////////////////////////////////////////
  // SeqToVarSymbol holds mapping information for: sequences -> symbols
  ////////////////////////////////////////////////////////////////

  /** SeqToVarSymbol holds all the information needed to map from sequence data to symbols for a set of (observed) random variables of the DFG. */
  struct SeqToVarSymbol
  {
    /**SeqToVarSymbol constructor. The negative default values will be
       overwritten as follows: symbolOffset=size, rSize=size,
       andrIndexOffset=indexOffset in the definition. Note that
       numeric_limits max values are used as default values in some
       cases to be able to determine when these were not set by the user.*/
    SeqToVarSymbol(string varName, 
		   string seqName, 
		   unsigned size = 1, 
		   int symbolOffset = numeric_limits<int>::max(),
		   int indexOffset = 0, 
		   bool optional = false, 
		   unsigned rSize = numeric_limits<unsigned>::max(), 
		   int rIndexOffset = numeric_limits<int>::max() );

    string varName;   //name of random variable
    string seqName;   // name of input sequence
    unsigned size;    // length of symbol
    int symbolOffset; // offset from previous symbol (symbols overlap if less than size)
    int indexOffset;  // offset relative to current sequence index
    bool optional;         // if false (default), symbol must always be made for the given random variable, i.e., the named sequence must be available
    unsigned rSize;        // length of right part of di-symbol (in which case "size" then defines left part)
    int rIndexOffset; // index offset for right part of di-symbol (in which case "indexOffset" defines offset for left part)
  };

  /** Defines mapping from stvVec to seqNames. stvVec[i].seqName =
      seqNames[ svtToSeqMap[i] ] or if stvVec[i].seqName is missing
      from seqNames, then svtToSeqMap[i]=-1. svtToSeqMap is the
      returned index vector. Precondition: stvToSeqMap.size =
      stvVec.size. */
  void mkStvToSeqMap(vector<int> & stvToSeqMap, vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames);

  /** Same as above, but returns by value.. */
  vector<int> mkStvToSeqMap(vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames);


  ////////////////////////////////////////////////////////////////
  // wrappers for making symbol_t and vector<symbol_t> from various
  // data types.
  // These are made inline for efficiency.
  ////////////////////////////////////////////////////////////////

  /** precondition: the symbol (sym) is of correct size!
      the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize. */
  inline void mkSymbol(symbol_t & sym, unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar = '.');
  inline symbol_t mkSymbol(unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar = '.');
  
  inline void mkSymbolVector(vector<symbol_t> & symVec, vector<unsigned> const & symSizeVec, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData = '.');
  inline void mkSymbolVector(vector<symbol_t> & symVec, unsigned const commonSymSize, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData = '.');

  // A diSymbol is composed of left and a right part.
  // precondition: the symbol (sym) is of correct size!
  // the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize.
  inline void mkDiSymbol(symbol_t & sym, unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar = '.');
  inline symbol_t mkDiSymbol(unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar = '.');

  // A diSymbol is composed of a left and a right part.
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, unsigned const commonSizeLeft, unsigned const commonSizeRight, vector<string> const & strVec, unsigned const seqSize, unsigned const leftStartPos, unsigned const rightStartPos, char missingDataChar = '.');

  /** make missing data symbol. Precondition: sym is of correct size. */
  inline void mkMissingDataSymbol(string & sym, unsigned symSize, char const missingDataChar = '.');

  /** Make symbol vector from seqVec as specified by seqToVarSymbol (information on seq->sym) and stvToSeqMap (mapping from stv to seq index). 
      Preconditions:
      1) symVec is indexed in the same order as vector<seqToVarSymbol> (and same size).
      2) symVec[i] has the size specified by vector<seqToVarSymbol>[i].
      3) symIndex is given in terms of symbols.
      4) seqSizes[i] gives the total char length of symVec[i].
      5) stvToSeqMap[i] indexes seqVec or has value -1 when the sequence corresponding to stvToSeqMap[i] is missing.
  */
  inline void mkSymbolVector(vector<symbol_t> & symVec, 
			     vector<string> const & seqVec, 
			     vector<long> const & seqSizes, 
			     vector<SeqToVarSymbol> const & stvVec, 
			     vector<int> const & stvToSeqMap, 
			     long symIndex, 
			     char const missingData = '.');

  /** Make di-symbol vector from seqVec as specified by seqToVarSymbol (information on seq->sym) and stvToSeqMap (mapping from stv to seq index). 
      Preconditions:
      1) symVec is indexed in the same order as vector<seqToVarSymbol> (and same size).
      2) symVec[i] has the size specified by vector<seqToVarSymbol>[i].
      3) rightSymIndex and leftSymIndex are given in terms of symbols.
      4) seqSizes[i] gives the total char length of symVec[i].
      5) stvToSeqMap[i] indexes seqVec or has value -1 when the sequence corresponding to stvToSeqMap[i] is missing.
  */
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, 
			     vector<string> const & seqVec, 
			     vector<long> const & seqSizes, 
			     vector<SeqToVarSymbol> const & stvVec, 
			     vector<int> const & stvToSeqMap, 
			     long leftSymIndex, 
			     long rightSymIndex, 
			       char const missingData = '.');


  ////////////////////////////////////////////////////////////////
  // wrappers for making stateMask2DVec from various data types.
  ////////////////////////////////////////////////////////////////

  /** Initiates and returns the stateMask2DVec data structure. For convenience. */
  inline stateMask2DVec_t initStateMask2DVec(long seqLength, unsigned varCount) {return stateMask2DVec_t(seqLength, stateMaskVec_t(varCount) );}

  /** Convert input observation-symbol data to 2D vectors of state
      masks. These can be used as input for the Discrete Factor Graph
      framework. 'offSet' defines the offset from one symbol start
      position to the next in the sequences of input symbols. The
      character denoting missing symbols can be defined by the
      missingDataChar variable.

      preconditions: 
      1) All sequences within a data entry (strVec or SeqData) are assumed to be of equal length.
      2) When the state maske table is given as a reference, it is assumed have correct dimension: 
         1st. dim is the number of symbols in the input sequences (taking offSet into account).
	 2nd. dim is the number of variables as defined in the StateMaskMapSet.
  */
  void mkStateMask2DVec(vector<string> const & strVec, stateMask2DVec_t & stateMask2DVec, vector<unsigned> const & seqToVarMap, StateMaskMapSet const & stateMaskMapSet, unsigned offSet = 1, char missingDataChar = '.');
  stateMask2DVec_t mkStateMask2DVec(vector<string> const & strVec, vector<unsigned> const & seqToVarMap, unsigned varCount, StateMaskMapSet const & stateMaskMapSet, unsigned offSet = 1, char missingDataChar = '.');
  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, unsigned offSet = 1, char missingDataChar = '.');
  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<unsigned> const & varMap, unsigned varCount, StateMaskMapSet const & stateMaskMapSet, unsigned offSet = 1, char missingDataChar = '.');

  // equivalent version for vectors of data sets
  vector<stateMask2DVec_t> const mkStateMask3DVec(vector<SeqData> const & seqDataVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, unsigned offSet = 1, char missingDataChar = '.');


  /** Convert input observation-symbol data to 2D vectors of state
      masks. The functions below define mappings from sequences to
      symbols by SeqToVarSymbol structures. Note that the above
      requirement of equal length sequences in the SeqData objects is
      relaxed in this case.
  */
  void  mkStateMask2DVec(stateMask2DVec_t & stateMask2DVec, SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<unsigned> const & stvToVarMap, StateMaskMapSet const & stateMaskMapSet, char missingDataChar = '.');
  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, StateMaskMapSet const & stateMaskMapSet, char missingDataChar = '.');
  stateMask2DVec_t mkStateMask2DVec(SeqData const & seqData, vector<SeqToVarSymbol> const & stvVec, vector<unsigned> const & stvToVarMap, StateMaskMapSet const & stateMaskMapSet, char missingDataChar = '.');


  ////////////////////////////////////////////////////////////////
  // Helper functions
  ////////////////////////////////////////////////////////////////

  /** return the number of symbols in SeqData given a certain fixed offset between each. */
  unsigned long symbolCount(vector<string> const & strVec, unsigned offSet);
  unsigned long symbolCount(long seqSize, unsigned offSet);

  /** Make symbol vector of correct size from vector of SeqToVarSymbol structures. */
  vector<string> initSymVec(vector<SeqToVarSymbol> const & stvVec, bool isPair=false);
  
  /** varNames[i] gives the name of variabe with index i (may be the
      empty string). seqNames[i] gives the name of input sequence i
      (must match a unique varName). The returned varMap maps the
      observed sequences to variable indexes. */
  vector<unsigned> mkSeqToVarMap(vector<string> const & varNames, vector<string> const & seqNames);

  /** Define mapping from vector of seqToVarSymbol data structures to a sequence vector. */
  void mkStvToSeqMap(vector<int> & stvToSeqMap, vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames);
  vector<int> mkStvToSeqMap(vector<SeqToVarSymbol> const & stvVec, vector<string> const & seqNames);

  /** Define mapping from vector of seqToVarSymbol data structures to an indexing of variables given by varNames. */
  vector<unsigned> mkStvToVarMap(vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames);


  /** stvToSeqMap[i] = -1 if the sequence name defined in SeqToVarSymbol[i] is missing from the input. This functions checks that only optional sequences are missing. It will abort otherwise. */
  bool onlyOptionalSeqMissingOrAbort(vector<SeqToVarSymbol> const & stvVec, vector<int> const & stvToSeqMap, string const & seqDataEntryId);

  /** This function checks that all sequence sizes correspond to the same symbol count, given the symbol offsets defined in stv. */
  long sameSymbolCountOrAbort(vector<SeqToVarSymbol> const & stvVec, vector<int> const & stvToSeqMap, vector<long> const & seqSizes, string const & seqDataEntryId);

  /** Check consistecy in the symbol sizes specified by the SeqToVarSymbol data structures and those given by the corresponding state map set. */
  void assertSymbolSizeConsistency(vector<SeqToVarSymbol> const & stvVec, vector<string> const & varNames, StateMaskMapSet const & stateMaskMapSet, bool isPair);

  /** Collect all sequence sizes and return in vector indexed in the same order as the sequences in seqData.*/
  vector<long> getSeqSizes(SeqData const & seqData);


  ////////////////////////////////////////////////////////////////
  // inline function definitions
  ////////////////////////////////////////////////////////////////

  inline void mkSymbol(symbol_t & sym, unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar)
  {
    long pos;
      for (unsigned i = 0; i < symSize; i++) {
	pos = symStartPos + i;
      if (pos < 0)
	sym[i] = missingDataChar;
      else if (pos >= seqSize)
	sym[i] = missingDataChar;
      else 
	sym[i] = seq[pos];
    }
  }


  inline symbol_t mkSymbol(unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar)
  {
    symbol_t sym(symSize, missingDataChar);
    mkSymbol(sym, symSize, seq, seqSize, symStartPos, missingDataChar);
    return sym;
  }
  
  
  inline void mkSymbolVector(vector<symbol_t> & symVec, vector<unsigned> const & symSizeVec, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkSymbol(symVec[i], symSizeVec[i], strVec[i], seqSize, symStartPos, missingData);
  }


  inline void mkSymbolVector(vector<symbol_t> & symVec, unsigned const commonSymSize, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkSymbol(symVec[i], commonSymSize, strVec[i], seqSize, symStartPos, missingData);
  }


  // A diSymbol is composed of left and a right part.
  // precondition: the symbol (sym) is of correct size!
  // the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize.
  inline void mkDiSymbol(symbol_t & sym, unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar)
  {
    // left part
    long pos;
    for (unsigned i = 0; i < sizeLeft; i++) {
      pos = leftStartPos + i;
      if (pos < 0)
	sym[i] = missingDataChar;
      else if (pos >= seqSize)
	sym[i] = missingDataChar;
      else 
	sym[i] = seq[pos];
    }

    // right part
    for (unsigned i = 0; i < sizeRight; i++) {
      pos = rightStartPos + i;
      if (pos < 0)
	sym[sizeLeft + i] = missingDataChar;
      else if (pos >= seqSize)
	sym[sizeLeft + i] = missingDataChar;
      else 
	sym[sizeLeft + i] = seq[pos];
    }
  }


  inline symbol_t mkDiSymbol(unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar)
  {
    symbol_t sym(sizeLeft + sizeRight, missingDataChar);
    mkDiSymbol(sym, sizeLeft, sizeRight, seq, seqSize, leftStartPos, rightStartPos, missingDataChar);
    return sym;
  }


  // A diSymbol is composed of a left and a right part.
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, unsigned const commonSizeLeft, unsigned const commonSizeRight, vector<string> const & strVec, unsigned const seqSize, unsigned const leftStartPos, unsigned const rightStartPos, char missingDataChar)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkDiSymbol(symVec[i], commonSizeLeft, commonSizeRight, strVec[i], seqSize, leftStartPos, rightStartPos, missingDataChar);
  }

  inline void mkMissingDataSymbol(string & sym, unsigned symSize, char const missingDataChar)
  {
    for (unsigned i = 0; i < symSize; i++)
      sym[i] = missingDataChar;
  }


  // Make symbol vector from seqVec as specified by seqToVarSymbol (information on seq->sym) and stvToSeqMap (mapping from stv to seq index). 
  //    Preconditions:
  //    1) symVec is indexed in the same order as vector<seqToVarSymbol> (and same size).
  //    2) symVec[i] has the size specified by vector<seqToVarSymbol>[i].
  //    3) symIndex is given in terms of symbols.
  //    4) seqSizes[i] gives the total char length of symVec[i].
  //    5) stvToSeqMap[i] indexes seqVec or has value -1 when the sequence corresponding to stvToSeqMap[i] is missing.
  inline void mkSymbolVector(vector<symbol_t> & symVec, 
			     vector<string> const & seqVec, 
			     vector<long> const & seqSizes, 
			     vector<SeqToVarSymbol> const & stvVec, 
			     vector<int> const & stvToSeqMap, 
			     long symIndex, 
			     char const missingData)
  {
    assert( symVec.size() == stvVec.size() );
    assert( stvVec.size() == stvToSeqMap.size() );
    for (unsigned i = 0; i < symVec.size(); i++) {
      int seqIdx = stvToSeqMap[i];
      if (seqIdx < 0)  // corresponding sequence is missing 
	if (stvVec[i].optional)
	  mkMissingDataSymbol(symVec[i], stvVec[i].size, missingData);
	else // not optional
	  errorAbort(string("from mkSymbolVector: Sequence of name, '") + stvVec[i].seqName + "' requred to make symbol for random variable, '" 
		     + stvVec[i].varName + "', is missing.");
      else {
	long symStartPos = symIndex * stvVec[i].symbolOffset + stvVec[i].indexOffset;
	mkSymbol(symVec[i], stvVec[i].size, seqVec[seqIdx], seqSizes[seqIdx], symStartPos, missingData);
      }
    }
  }


  //  Make di-symbol vector from seqVec as specified by seqToVarSymbol (information on seq->sym) and stvToSeqMap (mapping from stv to seq index). 
  //  Preconditions:
  //  1) symVec is indexed in the same order as vector<seqToVarSymbol> (and same size).
  //  2) symVec[i] has the size specified by vector<seqToVarSymbol>[i].
  //  3) rightSymIndex and leftSymIndex are given in terms of symbols.
  //  4) seqSizes[i] gives the total char length of symVec[i].
  //  5) stvToSeqMap[i] indexes seqVec or has value -1 when the sequence corresponding to stvToSeqMap[i] is missing.
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, 
			     vector<string> const & seqVec, 
			     vector<long> const & seqSizes, 
			     vector<SeqToVarSymbol> const & stvVec, 
			     vector<int> const & stvToSeqMap, 
			     long leftSymIndex, 
			     long rightSymIndex, 
			     char const missingData)
  {
    assert( symVec.size() == stvVec.size() );
    assert( stvVec.size() == stvToSeqMap.size() );
    for (unsigned i = 0; i < symVec.size(); i++) {
      int seqIdx = stvToSeqMap[i];
      if (seqIdx < 0)  // corresponding sequence is missing 
	if (stvVec[i].optional)
	  mkMissingDataSymbol(symVec[i], stvVec[i].size + stvVec[i].rSize, missingData);
	else // not optional
	  errorAbort(string("from mkSymbolVector: Sequence of name, '") + stvVec[i].seqName + "' requred to make symbol for random variable, '" 
		     + stvVec[i].varName + "', is missing.");
      else {
	long leftSymStartPos  = leftSymIndex  * stvVec[i].symbolOffset + stvVec[i].indexOffset;
	long rightSymStartPos = rightSymIndex * stvVec[i].symbolOffset + stvVec[i].rIndexOffset;
	mkDiSymbol(symVec[i], stvVec[i].size, stvVec[i].rSize, seqVec[seqIdx], seqSizes[seqIdx], leftSymStartPos, rightSymStartPos, missingData);
      }
    }
  }


} // end namespace phy


#endif  //__Observations_h
