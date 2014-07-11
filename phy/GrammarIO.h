/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __GrammarIO_h
#define __GrammarIO_h

#include "phy/Grammar.h"

namespace phy {

  namespace grammar {


    ////////////////////////////////////////////////////////////////
    // Transitions
    ////////////////////////////////////////////////////////////////

    /** Overloaded input / output operator for Transitions. */
    istream & operator>>(istream & str, Transition & trans);
    ostream & operator<<(ostream & str, Transition const & trans);

    // helper functions
    ttype_t getTransitionType(istream & str);


    ////////////////////////////////////////////////////////////////
    // Grammar
    ////////////////////////////////////////////////////////////////

    /** Overloaded input operator for vector of Transitions. Basically reads a grammar. */
    istream & operator>>(istream & str, vector<Transition> & transVec);

    /** Overloaded input / output operator for Grammar. */
    istream &operator>>(istream & str, Grammar & g);
    ostream &operator<<(ostream & str, Grammar const & g);

    /** Read a grammar with the order of emission distributions defined
	explicitly by the singleEmitModelNames and pairEmitModelNames
	parameters. The emitModel names may be given by the emitWrappers
	class (defined in emitWrap.h). */
    Grammar readGrammar(string const & grammarFile, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames = vector<string>() );

    /** Read / write grammar from a file. Note that the order of
	emissions is implicitly giving by the order they occur in the
	file. Use above version of readGrammar to avoid ambiguity.*/
    Grammar readGrammar(string const & file);
    void writeGrammar(string const & file, Grammar const & g);

    ////////////////////////////////////////////////////////////////
    // AnnoMap
    ////////////////////////////////////////////////////////////////

    EmitAnnoMap readEmitAnnoMap(string const & file);
    void writeEmitAnnoMap(string const & file, EmitAnnoMap const & annoMap);

    /** Overloaded input / output operator for EmitAnnoMaps. */
    istream & operator>>(istream & str, EmitAnnoMap & annoMap);
    ostream & operator<<(ostream & str, EmitAnnoMap const & annoMap);

    TransAnnoMap mkTransAnnoMap(string const & emitAnnoMapFile, Grammar const & g);


  } // end namespace grammar

} // end namespace phy

#endif  // __GrammarIO_h
