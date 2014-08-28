/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __DfgIO_h
#define __DfgIO_h

#include "phy/DiscreteFactorGraph.h"
#include "phy/ObservationsIO.h"
#include "phy/FactorsIO.h"
#include "phy/DataIO.h"

namespace phy {

  using namespace std;

  /** Structure holding a discrete factor graph (dfg) and convenient information for using it. */
  struct DfgInfo 
  {
    /* Constructors **/
    DfgInfo(vector<string> const & varNames, vector<string> const & facNames, vector<AbsBasFacPtr_t> const & facVec, vector<StateMapPtr_t> const & stateMapVec, DFG const & dfg) :
      varNames(varNames), facNames(facNames), facVec(facVec), facSet(facVec), stateMapVec(stateMapVec), stateMaskMapSet(stateMapVec), dfg(dfg) {};

    DfgInfo(vector<string> const & varNames, 
	    vector<string> const & facNames, 
	    vector<string> const & potNames, 
	    vector<vector<unsigned> > const & facNeighbors, 
	    map<string, StateMapPtr_t> const & smMap, 
	    map<string, AbsBasFacPtr_t> const & facMap, 
	    map<string, string> const & var2smMap);

    /* Data **/
    vector<string> varNames;           // defines enumeration of random variables
    vector<string> facNames;           // defines enumeration of factors
    vector<AbsBasFacPtr_t> facVec;     // factor potential (pointer) for each factor. Each factor also holds the potential name
    CompositeFactorSet facSet;         // data structure for holding and optimizing factor potentials
    vector<StateMapPtr_t> stateMapVec; // stateMap for each variable
    StateMaskMapSet stateMaskMapSet;   // defines mapping between symbols and stateMasks for all variables
    DFG dfg;                           // Discrete factor graph

    void writeInfo( ostream & str );
  };


  /** Read mapping of variables to stateMaps.
      The returned map is indexed by variable name. */
  map<string, string> readVariables(string const & file);
  map<string, string> readVariables(istream & str);


  /** Write mapping of variables to stateMaps. */
  void writeVariables(string const & file, vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec);
  void writeVariables(ostream & str, vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec);


  /** Read factor graph structure;
      Output is by reference. varNames, facNames, and facNeighbors will be overwritten. */
  void readFactorGraph(string const & file, vector<string> & varNames, vector<string> & facNames, vector<string> & potNames, vector<vector<unsigned> > & facNeighbors);
  void readFactorGraph(istream & str, vector<string> & varNames, vector<string> & facNames, vector<string> & potNames, vector<vector<unsigned> > & facNeighbors);


  /** write factor graph structure to specification file. */
  void writeFactorGraph(string const & file, vector<string> const & varNames, vector<string> const & facNames, vector<string> const & potNames, vector<vector<unsigned> > const & facNeighbors);
  void writeFactorGraph(ostream & str, vector<string> const & varNames, vector<string> const & facNames, vector<string> const & potNames, vector<vector<unsigned> > const & facNeighbors);


  /** read and write dfg info. */
  DfgInfo readDfgInfo(string const & stateMapsFile, string const & factorPotentialsFile, string const & variablesFile, string const & factorGraphFile); 
  void writeDfgInfo(DfgInfo const & dfgInfo, string const & stateMapsFile, string const & factorPotentialsFile, string const & variablesFile, string const & factorGraphFile); 


  /** Class for reading symbolic variable data */
  class VarData : public NamedData<symbol_t> {
  public:

    /** Constructor */
    VarData(string const & file, vector<string> const & allNames) : NamedData<symbol_t>(file, allNames) {};
  };


  /** Class for reading potential matrix data */
  class FacData : public NamedData<xmatrix_t> {
  public:

    /** Constructor */
    FacData(string const & file, vector<string> const & allNames) : NamedData<xmatrix_t>(file, allNames) {};
  };

  /** Helper functions declared here for debug purposes */
  vector<AbsBasFacPtr_t> mkFacVec(vector<string> const & potNames, map<string, AbsBasFacPtr_t> const & facMap);
  vector<StateMapPtr_t> mkStateMapVec(vector<string> const & varNames, map<string, string> const & var2smMap, map<string, StateMapPtr_t> const & smMap);
  vector<unsigned> mkVarDimensions(vector<StateMapPtr_t> const & smVec);
  map<string, StateMapPtr_t> smVecToSmMap(vector<StateMapPtr_t> const & smVec);
  vector<string> mkPotNames(vector<AbsBasFacPtr_t> const & facVec);
} // end namespace phy

#endif  // __DfgIO_h
