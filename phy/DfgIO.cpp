/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/DfgIO.h"

namespace phy {


  // helper functions, not declared in header
  vector<AbsBasFacPtr_t> mkFacVec(vector<string> const & potNames, map<string, AbsBasFacPtr_t> const & facMap)
  {
    vector<AbsBasFacPtr_t> facVec;
    BOOST_FOREACH(string const & pot, potNames) {
      map<string, AbsBasFacPtr_t>::const_iterator iter = facMap.find(pot);
      if ( iter == facMap.end() )
	errorAbort("From mkFacVec: factor with name '" + pot + "' not found in factor map.");
      facVec.push_back(iter->second);
    }
    return facVec;
  }


  vector<StateMapPtr_t> mkStateMapVec(vector<string> const & varNames, map<string, string> const & var2smMap, map<string, StateMapPtr_t> const & smMap)
  {
    vector<StateMapPtr_t> smVec;
    BOOST_FOREACH(string const & var, varNames) {
      map<string, string>::const_iterator iter1 = var2smMap.find(var);
      if ( iter1 == var2smMap.end() )
	errorAbort("From mkStateMapVec: variable with name '" + var + "' not found in variable to stateMap map.");
      string const sm = iter1->second;

      map<string, StateMapPtr_t>::const_iterator iter2 = smMap.find(sm);
      if ( iter2 == smMap.end() )
	errorAbort("From mkStateMapVec: state map with name '" + sm + "' not found in stateMap map.");
      smVec.push_back(iter2->second);
    }
    return smVec;
  }


  vector<unsigned> mkVarDimensions(vector<StateMapPtr_t> const & smVec)
  {
    vector<unsigned> varDim;
    BOOST_FOREACH(StateMapPtr_t const & smPtr, smVec)
      varDim.push_back( smPtr->stateCount() );
    return varDim;
  }


  map<string, StateMapPtr_t> smVecToSmMap(vector<StateMapPtr_t> const & smVec)
  {
    map<string, StateMapPtr_t> smMap;
    BOOST_FOREACH(StateMapPtr_t const & smPtr, smVec) {
      if (smPtr->name().size() == 0)
	errorAbort("from smVecToSmMap: attempt to add stateMap lacking name a to stateMap map");
      smMap[smPtr->name()] = smPtr;
    }
    return smMap;
  }


  vector<string> mkPotNames(vector<AbsBasFacPtr_t> const & facVec)
  {
    vector<string> potVec;
    BOOST_FOREACH(AbsBasFacPtr_t const & facPtr, facVec)
      potVec.push_back( facPtr->name() );
    return potVec;
  }

  // definition of header declared functions
  DfgInfo::DfgInfo(vector<string> const & varNames, 
		   vector<string> const & facNames, 
		   vector<string> const & potNames, 
		   vector<vector<unsigned> > const & facNeighbors, 
		   map<string, StateMapPtr_t> const & smMap, 
		   map<string, AbsBasFacPtr_t> const & facMap, 
		   map<string, string> const & var2smMap) :
    varNames(varNames), 
    facNames(facNames), 
    facVec( mkFacVec(potNames, facMap) ),
    facSet(facVec),
    stateMapVec( mkStateMapVec(varNames, var2smMap, smMap) ), 
    stateMaskMapSet(stateMapVec),
    dfg( mkVarDimensions(stateMapVec), facSet.mkFactorVec(), facNeighbors ) 
  {};


  map<string, string> readVariables(string const & file)
  {
    ifstream f;
    openInFile(f, file);
    return readVariables(f);
  }


  map<string, string> readVariables(istream & str)
  {
    map<string, string> varMap;

    skipWhiteSpaceAndComments(str);
    while ( str.good() ) {
      while ( moreTags(str) ) {
	// get state map name
	string smName;
	getFeatureAndSkipLine(str, "STATE_MAP_NAME:", smName);

	// get variable names
	string tag;
	str >> tag;
	if (tag != "VAR_NAMES:")
	  errorAbort("From readVariables: Unknown tag ('" + tag + "') in specification of variable to stateMap map for state map with name '" + smName + "'.");
	string varStr;
	getline(str, varStr);
	vector<string> varVec = split( strip(varStr) );

	// insert name pairs into map
	BOOST_FOREACH(string const & var, varVec) {
	  map<string, string>::const_iterator iter = varMap.find(smName);
	  if ( iter != varMap.end() )
	    errorAbort("From readVariables: variable with name '" + smName + "' is assigned more than one state map.");
	  varMap[var] = smName;
	}

	skipWhiteSpaceAndComments(str);
      }
    }
    return varMap;
  }


  map<string, vector<string> > mkStateMap2varNamesMap(vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec)
  {
    assert( varNames.size() == stateMapVec.size() );
    map<string, vector<string> > sm2varNames;
    for (unsigned i = 0; i < varNames.size(); i++) {
      string smName  = stateMapVec[i]->name();
      string varName = varNames[i];
      map<string, vector<string> >::const_iterator iter=sm2varNames.find(smName);
      if ( iter == sm2varNames.end() ) // not found
	sm2varNames[smName] = vector<string>();
      sm2varNames[smName].push_back(varName); 
    }
    return sm2varNames;
  }


  void writeVariables(string const & file, vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec)
  {
    ofstream f( file.c_str() );
    writeVariables(f, varNames, stateMapVec);
  }


  void writeVariables(ostream & str, vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec)
  {
    map<string, vector<string> > sm2varNames = mkStateMap2varNamesMap(varNames, stateMapVec);
    typedef map<string, vector<string> >::value_type pair_t;
    BOOST_FOREACH(pair_t const & p, sm2varNames) {
      string smName = p.first;
      vector<string> varVec = p.second;
      str << "STATE_MAP_NAME:\t" << smName << endl;
      str << "VAR_NAMES:\t";
      BOOST_FOREACH(string const & var, varVec)
	str << var << " ";
      str << endl << endl; 
    }
  }


  vector<unsigned> mkFacNbs(vector<string> & varNames, vector<string> facNbsNames)
  // facNbsNames is an (ordered) list of random variables neighboring a specific factor
  // new variable names will be added to varNames and the facNbsNames converted into a list of varNames indices
  {
    vector<unsigned> facNbs;
    // add var names if not present
    BOOST_FOREACH(string const & var, facNbsNames) {
      vector<string>::const_iterator iter = find(varNames.begin(), varNames.end(), var);
      if ( iter == varNames.end() )
	varNames.push_back(var);
      facNbs.push_back( getIndex(varNames, var) );
    }
    return facNbs;
  }


  void readFacInfo(istream & str, string & facName, string & potName, vector<string> & facNbsNames)
  {
    // clear and define variables 
    facNbsNames.clear();
    facName.clear();
    potName.clear();
    string nb1, nb2;

    while ( moreTags(str) ) {
      string tag;
      str >> tag;
      if (tag == "NAME:") {
	str >> facName;
	skipLine(str); // skip rest of line
      }
      else if (tag == "NB1:") {
	str >> nb1;
	skipLine(str); // skip rest of line
      }
      else if (tag == "NB2:") {
	str >> nb2;
	skipLine(str); // skip rest of line
      }
      else if (tag == "POT:") {
	str >> potName;
	skipLine(str); // skip rest of line
      }
      else
	errorAbort("From readFacInfo: use of undefined tag '" + tag + "' for factor with name '" + facName + "'.");
    }
    if (facName.size() == 0 or potName.size() == 0 or nb1.size() == 0)
      errorAbort("From readFacInfo: Incomplete specification of factor with name '" + facName + "' and potential '" + potName + "' and NB1 '" + nb1 + "'.");
    
    facNbsNames.push_back(nb1);
    if (nb2.size() != 0)
      facNbsNames.push_back(nb2);
  }


  void readFactorGraph(string const & file, vector<string> & varNames, vector<string> & facNames, vector<string> & potNames, vector<vector<unsigned> > & facNeighbors)
  {
    ifstream f;
    openInFile(f, file);
    readFactorGraph(f, varNames, facNames, potNames, facNeighbors);
  }


  void readFactorGraph(istream & str, vector<string> & varNames, vector<string> & facNames, vector<string> & potNames, vector<vector<unsigned> > & facNeighbors)
  {
    // clear referenced datas structures
    varNames.clear();
    facNames.clear();
    potNames.clear();
    facNeighbors.clear();

    // parse specification
    skipWhiteSpaceAndComments(str);
    while ( str.good() ) {
      string fac, pot;
      vector<string> facNbsNames;
      readFacInfo(str, fac, pot, facNbsNames);
      facNames.push_back(fac);
      potNames.push_back(pot);
      facNeighbors.push_back( mkFacNbs(varNames, facNbsNames) );	
      skipWhiteSpaceAndComments(str);
    }
  }


  void writeFactorGraph(string const & file, vector<string> const & varNames, vector<string> const & facNames, vector<string> const & potNames, vector<vector<unsigned> > const & facNeighbors)
  {
    ofstream f( file.c_str() );
    writeFactorGraph(f, varNames, facNames, potNames, facNeighbors);
  }


  void writeFactorGraph(ostream & str, vector<string> const & varNames, vector<string> const & facNames, vector<string> const & potNames, vector<vector<unsigned> > const & facNeighbors)
  {
    assert( facNames.size() == facNeighbors.size() );
    assert( facNames.size() == potNames.size() );
    for (unsigned i = 0; i < facNames.size(); i++) {
      // fac
      str << "NAME:" << "\t" << facNames[i] << endl;

      // fac neighbors
      vector<unsigned> const & nbs = facNeighbors[i];
      assert(nbs.size() > 0);
      for (unsigned j = 0; j < nbs.size(); j++) {
	assert( nbs[j] < varNames.size() );
	str << "NB" << j + 1 << ":" << "\t" << varNames[ nbs[j] ] << endl;
      }

      // fac potential
      str << "POT:" << "\t" << potNames[i] << endl;
      str << endl;
    }
  }


  DfgInfo readDfgInfo(string const & stateMapsFile, string const & factorPotentialsFile, string const & variablesFile, string const & factorGraphFile)
  {
    vector<string> varNames;
    vector<string> facNames;
    vector<string> potNames; 
    vector<vector<unsigned> > facNeighbors;
    readFactorGraph(factorGraphFile, varNames, facNames, potNames, facNeighbors);

    return DfgInfo(varNames, 
		   facNames, 
		   potNames, 
		   facNeighbors, 
		   readStateMapFile(stateMapsFile), 
		   readFactorFile(factorPotentialsFile), 
		   readVariables(variablesFile) );
  }


  void writeDfgInfo(DfgInfo const & dfgInfo, string const & stateMapsFile, string const & factorPotentialsFile, string const & variablesFile, string const & factorGraphFile)
  {
    writeStateMapFile2( stateMapsFile, smVecToSmMap(dfgInfo.stateMapVec) );
    writeFactorVec(factorPotentialsFile, dfgInfo.facVec);
    writeVariables(variablesFile, dfgInfo.varNames, dfgInfo.stateMapVec);

    // need to extact facNeigbors from dfg, where they are encode together with variable neighbors 
    // (somewhat hacky solution since it relies in implementation details)
    unsigned facOffset = dfgInfo.dfg.convFacToNode(0);
    vector<vector<unsigned> > facNeighbors(dfgInfo.dfg.neighbors.begin() + facOffset, dfgInfo.dfg.neighbors.end() );
    writeFactorGraph(factorGraphFile, dfgInfo.varNames, dfgInfo.facNames, mkPotNames(dfgInfo.facVec), facNeighbors);
  }

} // end namespace phy
