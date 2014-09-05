/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/ObservationsIO.h"

namespace phy {

  // convert string with space separated symbols into vector of strings
  vector<symbol_t> symStrToSymVec(string const & symStr)
  {
    return split( strip(symStr) );
  }

  // convert string with semicolon separated meta-symbols into map of meta-symbol to vector of symbols
  boost::unordered_map<symbol_t, vector<symbol_t> > metaSymStrToMetaSymMap(string const & metaSymStr)
  {
    boost::unordered_map<symbol_t, vector<symbol_t> > metaSymMap;
    vector<string> metaSymSpecs = split( strip(metaSymStr, "; \t"), ";");
    BOOST_FOREACH(string const & s, metaSymSpecs) {
      vector<string> v = split(s, "=");
      assert(v.size() == 2);
      metaSymMap[ strip(v[0]) ] = symStrToSymVec( v[1] );
    }
    return metaSymMap;
  }


  const StateMap stateMapDispatcher(string const & name)
  {
    if (name == "nuc")
      return mkNucStateMap();
    else if (name == "metaNuc")
      return mkMetaNucStateMap();
    else
      errorAbort("No StateMap named :'" + name + "'."); 

    return mkNucStateMap(); // to please compiler
  }


  istream & operator>>(istream & str, StateMap & staMap)
  {
    skipWhiteSpaceAndComments(str);

    unsigned mult = 1; // multiplicity
    string name;
    string tag;
    string real;
    bool is_real = false;
    double minv;
    double maxv;
    int no_bp;
    str >> tag;

    if (tag == "ALPHABET_NAME:") {  // use pre-defined state map
      str >> name;

      unsigned mult = 1;
      if ( isdigit(name[0]) ) {  // if numeric then assume n-name form, where n gives the multiplicity and name names a basic StateMap
	vector<string> v = split(name, "-");
	assert(v.size() == 2);
	fromString(v[0], mult);
	name = v[1];
      }
      
      staMap = stateMapDispatcher(name);

      // use multiplicity
      assert(mult > 0 and mult <= 4); // limit on possible sizes of n 
      if (mult > 1)
	staMap = StateMap(staMap, mult);
    }

    else if (tag == "NAME:") {  // state map definition follows
      str >> name;
      skipWhiteSpaceAndComments(str);

      string symbols;
      string metaSymbols;
      while ( moreTags(str) ) {
	str >> tag;
	if (tag == "SYMBOLS:")
	  getline(str, symbols);
	else if (tag == "META_SYMBOLS:")
	  getline(str, metaSymbols);
	else if (tag == "MULTIPLICITY:") {
	  str >> mult;
	  skipLine(str); // skip rest of line
	}
	else if (tag == "REAL:") {
	  str >> real;
	  if( real == "TRUE" or real == "T")
	    is_real = true;
	  else if ( real == "FALSE" or real == "F")
	    is_real = false;
	  else
	    errorAbort("REAL: tag should be either T, TRUE, F, FALSE, note default is FALSE");
	  skipLine(str);
	}
	else if ( tag == "BREAKPOINTS:" ){
	  str >> no_bp;
	  skipLine(str);
	}
	else if ( tag == "MIN:"){
	  str >> minv;
	  skipLine(str);
	}
	else if ( tag == "MAX:"){
	  str >> maxv;
	  skipLine(str);
	}
	else
	  errorAbort("Unknown tag ('" + tag + "') in state map specification with name '" + name + "'.");
      }
      // check
      if (symbols.size() == 0 && !is_real)
	errorAbort("Missing symbol string in state map spscification with name  '" + name + "'.");
	
      if(!is_real)
	staMap = StateMap(symStrToSymVec(symbols), metaSymStrToMetaSymMap(metaSymbols), name);
      if(is_real){
	//TODO Redundant code
	//	vector_t breakpoints(no_bp);
	//	for(int i = 0; i < no_bp; ++i)
	//	  breakpoints(i) = minv + i*(maxv-minv)/(no_bp-1);
	staMap = StateMap( no_bp+1, minv, maxv, name);
      }

      // use multiplicity
      assert(mult > 0 and mult <= 4); // limit on possible sizes of n 
      if (mult > 1)
	staMap = StateMap(staMap, mult, name);
    }

    else
      errorAbort("Unknown initial tag '" + tag + "' used in state map specification");

    return str;
  }


  ostream & operator<<(ostream & str, StateMap const & staMap)
  {
    if (staMap.name().size() == 0) // no defined name
      errorAbort("Currently only named StateMaps can be output");
    str << "ALPHABET_NAME:\t" << staMap.name() << endl;
    return str;
  }


  void writeStateMap2(ostream & str, StateMap const & staMap)
  {
    if (staMap.name().size() == 0) // no defined name
      errorAbort("From writeStateMap2: Only state maps with names can be output");

    //statemap separator
    str << endl;

    // name
    str << "NAME:\t" << staMap.name() << endl;

    // symbols
    unsigned stateCount = staMap.stateCount();
    str << "SYMBOLS:\t";
    for (unsigned i = 0; i < stateCount; i++)
      str << " " << staMap.state2Symbol(i);
    str << endl;

    // meta symbols
    unsigned metaStateCount = staMap.metaStateCount();
    str << "META_SYMBOLS:\t";
    for (unsigned i = stateCount; i < metaStateCount; i++) {
      vector<symbol_t> degVec = staMap.degeneracyVector( staMap.state2Symbol(i) );
      str << " " << staMap.state2Symbol(i) << " =";
      BOOST_FOREACH(string const & s, degVec)
	str << " " << s;
      str << ";";
    }
    str << endl;
  }


  map<string, StateMapPtr_t> readStateMapFile(string const & file)
  {
    ifstream f;
    openInFile(f, file);
    return readStateMapFile(f);
  }


  map<string, StateMapPtr_t> readStateMapFile(istream & str)
  {
    map<string, StateMapPtr_t> smMap;
    skipWhiteSpaceAndComments(str);
    while ( str.good() ) {
      StateMapPtr_t smPtr(new StateMap("") );
      str >> *smPtr;
      smMap[ smPtr->name() ] = smPtr;

      skipWhiteSpaceAndComments(str);
    }
    return smMap;
  }


  void writeStateMapFile(string const & file, map<string, StateMapPtr_t> const & smMap)
  {
    ofstream f( file.c_str() );
    writeStateMapFile(f, smMap);
  }


  void writeStateMapFile(ostream & str, map<string, StateMapPtr_t> const & smMap)
  {
    for (map<string, StateMapPtr_t>::const_iterator iter = smMap.begin();
	 iter != smMap.end();
	 iter++) {
      str << *(iter->second);
      str << endl;
    }
  }


  void writeStateMapFile2(string const & file, map<string, StateMapPtr_t> const & smMap)
  {
    ofstream f( file.c_str() );
    writeStateMapFile2(f, smMap);
  }


  void writeStateMapFile2(ostream & str, map<string, StateMapPtr_t> const & smMap)
  {
    for (map<string, StateMapPtr_t>::const_iterator iter = smMap.begin();
	 iter != smMap.end();
	 iter++)
      writeStateMap2(str, *(iter->second) );
  }



  // Overloaded input / output operator for SeqToVarSymbol entries.
  istream & operator>>(istream & str, SeqToVarSymbol & seqToVarSym)
  {
    // the max values are used as the non-allowed defaults in the SeqToVarSymbol constructor - hence we used it here also 
    // this is messy and error prone! Better if default could be defined in class or similar
    unsigned DEFAULT_VALUE_UNSIGNED = numeric_limits<unsigned>::max();
    int DEFAULT_VALUE_INT = numeric_limits<int>::max();

    string varName;
    string seqName;
    unsigned size = 1;
    int symbolOffset = DEFAULT_VALUE_INT;
    int indexOffset = 0;
    bool optional = false;
    int rSize = DEFAULT_VALUE_UNSIGNED;
    int rIndexOffset = DEFAULT_VALUE_INT;

    skipWhiteSpaceAndComments(str);
    while ( moreTags(str) ) {
      string tag;
      str >> tag;
      if (tag == "VAR_NAME:") {
	str >> varName;
	skipLine(str); // skip rest of line
      }
      else if (tag == "SEQ_NAME:") {
	str >> seqName;
	skipLine(str); // skip rest of line
      }
      else if (tag == "SIZE:") {
	str >> size;
	skipLine(str); // skip rest of line
      }
      else if (tag == "SYMBOL_OFFSET:") {
	str >> symbolOffset;
	skipLine(str); // skip rest of line
      }
      else if (tag == "INDEX_OFFSET:") {
	str >> indexOffset;
	skipLine(str); // skip rest of line
      }
      else if (tag == "OPTIONAL:") {
	optional = true;
	skipLine(str); // skip rest of line
      }
      else if (tag == "R_SIZE:") {
	str >> rSize;
	skipLine(str); // skip rest of line
      }
      else if (tag == "R_INDEX_OFFSET:") {
	str >> rIndexOffset;
	skipLine(str); // skip rest of line
      }
      else
	errorAbort("From operator>> for seqToVarSymbol: use of undefined tag '" + tag + "' for entry with VAR_NAME '" + varName + "'.");
    }
    if (varName.size() == 0 or seqName.size() == 0)
      errorAbort("From operator>> for seqToVarSymbol: Incomplete specification of seqToVarSymbol entry with varName '" + varName + "' and seqName '" + seqName + "'.");
    
    seqToVarSym = SeqToVarSymbol(varName, seqName, size, symbolOffset, indexOffset, optional, rSize, rIndexOffset);
    return str;
  }


  ostream & operator<<(ostream & str, SeqToVarSymbol const & stv)
  {
    str << "VAR_NAME:\t" << stv.varName << endl;
    str << "SEQ_NAME:\t" << stv.seqName << endl;
    if (stv.size != stv.rSize or stv.size != 1)
      str << "SIZE:\t" << stv.size << endl;
    if (stv.symbolOffset != static_cast<int>(stv.size))
      str << "SYMBOL_OFFSET:\t" << stv.symbolOffset << endl;
    if (stv.indexOffset != 0)
      str << "INDEX_OFFSET:\t" << stv.indexOffset << endl;
    if (stv.optional != false)
      str << "OPTIONAL:\t" << 1 << endl;
    if (stv.rSize != stv.size)
      str << "R_SIZE:\t" << stv.rSize << endl;
    if ( (-1 * stv.rIndexOffset) != stv.indexOffset)
      str << "R_INDEX_OFFSET:\t" << stv.rIndexOffset << endl;

    str << endl;
    return str;
  }


  vector<SeqToVarSymbol> readSeqToVarSymbols(string const & file)
  {
    ifstream str;
    openInFile(str, file);
    return readSeqToVarSymbols(str);
  }


  vector<SeqToVarSymbol> readSeqToVarSymbols(istream & str)
  {
    vector<SeqToVarSymbol> v;

    // parse specification
    skipWhiteSpaceAndComments(str);
    while ( str.good() ) {
      SeqToVarSymbol stv(" "," ");
      str >> stv;
      v.push_back(stv);	
      skipWhiteSpaceAndComments(str);
    }
    return v;
  }

  void writeSeqToVarSymbols(string const & file, vector<SeqToVarSymbol> const & v)
  {
    ofstream str;
    openOutFile(str, file);
    writeSeqToVarSymbols(str, v);
  }

  void writeSeqToVarSymbols(ostream & str, vector<SeqToVarSymbol> const & v)
  {
    BOOST_FOREACH(SeqToVarSymbol const & stv, v)
      str << stv;
  }

} // end namespace phy
