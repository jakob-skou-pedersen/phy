/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __ObservationsIO_h
#define __ObservationsIO_h

#include "phy/Observations.h"
#include <iostream>
#include <cctype> 
#include <map> 

namespace phy {

  /** Overloaded input / output operator for StateMaps. */
  istream & operator>>(istream & str, StateMap & staMap);
  ostream & operator<<(ostream & str, StateMap const & staMap);

  /** Writes explicit output specification of StateMaps. */
  void writeStateMap2(ostream & str, StateMap const & staMap);

  /** Read state map specification file and return map between names and stateMap pointers */
  map<string, StateMapPtr_t> readStateMapFile(string const & file);
  map<string, StateMapPtr_t> readStateMapFile(istream & str);

  /** Write state map specification file (only useful with named state maps) */
  void writeStateMapFile(string const & file, map<string, StateMapPtr_t> const & smMap);
  void writeStateMapFile(ostream & str, map<string, StateMapPtr_t> const & smMap);

  /** Write state map specification file (outputs complete state map specifications) */
  void writeStateMapFile2(string const & file, map<string, StateMapPtr_t> const & smMap);
  void writeStateMapFile2(ostream & str, map<string, StateMapPtr_t> const & smMap);

  /** Dispatch named StateMaps. */
  const StateMap stateMapDispatcher(string const & name);

  ////////////////////////////////////////////////////////////////
  // The functions below do io of SeqToVarSymbol data structures,
  // which hold the necessary information form mapping sequences to
  // random variable symbols.
  ////////////////////////////////////////////////////////////////

  /** Overloaded input / output operator for SeqToVarSymbol entries. */
  istream & operator>>(istream & str, SeqToVarSymbol & seqToVarSym);
  ostream & operator<<(ostream & str, SeqToVarSymbol const & seqToVarSym);

  /** Read SeqToVarSymbol entries from file. */
  vector<SeqToVarSymbol> readSeqToVarSymbols(string const & file);
  vector<SeqToVarSymbol> readSeqToVarSymbols(istream & str);

  /** Write SeqToVarSymbol entries to file. */
  void writeSeqToVarSymbols(string const & file, vector<SeqToVarSymbol> const & v);
  void writeSeqToVarSymbols(ostream & str, vector<SeqToVarSymbol> const & v);


} // end namespace phy

#endif  //__ObservationsIO_h
