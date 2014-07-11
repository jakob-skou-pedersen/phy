/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/GrammarIO.h"

namespace phy {

  namespace grammar {


    ////////////////////////////////////////////////////////////////
    // Transitions
    ////////////////////////////////////////////////////////////////

    ttype_t getTransitionType(istream & str)
    {
      string typeTag;
      str >> typeTag;

      return static_cast<ttype_t>( getIndex(transTypeTagVec, typeTag) );
    }
  
    // check inputs...
    istream & operator>>(istream & str, Transition & trans)
    {
      skipWhiteSpaceAndComments(str);
    
      unsigned tagCount = 0;

      trans.type = TYPE_COUNT; // check that we set the type

      while ( moreTags(str) ) {
	string tag;
	str >> tag;
	if (tag == "NAME:")
	  str >> trans.name;
	else if (tag == "TYPE:")
	  trans.type = getTransitionType(str);
	else if (tag == "FROM:")
	  str >> trans.from_id;
	else if (tag == "TO:")
	  str >> trans.to_id;
	else if (tag == "TO_LEFT:")
	  str >> trans.toL_id;
	else if (tag == "TO_RIGHT:")
	  str >> trans.toR_id;
	else if (tag == "PROB:")
	  str >> trans.p;
	else if (tag == "EMIT_MODEL:")
	  str >> trans.e_id;
	else
	  errorAbort("Unknown tag ('" + tag + "') in transition with name  '" + trans.name + "'.");
	skipLine(str); // skip rest of line
    
	tagCount++;
      }

      // sanity checks
      if (trans.type == TYPE_COUNT)
	errorAbort("Transition (with name '" + trans.name + "') has no TYPE.");
      if (trans.type == END and tagCount != 4)
	errorAbort("Transition (with name '" + trans.name + "') has wrong number of TAGS. (TYPE END should have exactly four.)" );
      if (trans.type == BIFURCATE and tagCount != 6)
	errorAbort("Transition (with name '" + trans.name + "') has wrong number of TAGS. (TYPE BIFURCATE should have exactly six.)" );
      if ( (trans.type == LEFT or trans.type == RIGHT or trans.type == PAIR) and tagCount != 6)
	errorAbort("Transition (with name '" + trans.name + "') has wrong number of TAGS. (emitting, non-terminal types should have exactly six.)" );
      if ( (trans.type == TERMINAL) and tagCount != 5)
	errorAbort("Transition (with name '" + trans.name + "') has wrong number of TAGS. (TYPE TERMINAL should have exactly five.)" );

      return str;
    }


    ostream & operator<<(ostream & str, Transition const & trans)
    {
      str << "NAME:\t" << trans.name << endl;
      str << "TYPE:\t" << transTypeTagVec[trans.type] << endl;
      str << "FROM:\t" << trans.from_id << endl;
      if (trans.type == BIFURCATE) {
	str << "TO_LEFT:\t" << trans.toL_id << endl;
	str << "TO_RIGHT:\t" << trans.toR_id << endl;
      }
      else if (not (trans.type == END or trans.type == TERMINAL) )  // END and TERMINAL do not have destination states
	str << "TO:\t" << trans.to_id << endl;
      str << "PROB:\t" << trans.p << endl;
      if ( isEmitType(trans.type) )
	str << "EMIT_MODEL:\t" << trans.e_id << endl;
      str << endl; // white space line to denote end of Transition

      return str;
    }


    ////////////////////////////////////////////////////////////////
    // Grammar
    ////////////////////////////////////////////////////////////////


    istream & operator>>(istream & str, vector<Transition> & transVec)
    {
      skipWhiteSpaceAndComments(str);

      // format version
      int ver = 0;
      getFeature(str, "GRAMMAR_VERSION:", ver);
      if (ver != 1)
	errorAbort("Sorry, only grammar version 1 is supported.");

      // transition count
      //      unsigned n = 0;
      //      getFeature(str, "TRANSITION_COUNT:", n);

      // Transitions
      //      for (unsigned i = 0; i < n; i++) {
      //	if ( not str.good() )
      //	  errorAbort("An error occurred while parsing the '" + toString(n) + "' of the grammar specification. \nPlease check number of transitions and the grammar format.");
      
      skipWhiteSpaceAndComments(str);
      while ( str.good() ) {
      	while ( moreTags(str) ) {
	  Transition trans;
	  str >> trans;
	  transVec.push_back(trans);
	  skipWhiteSpaceAndComments(str);
	}
      }
      return str;
    }

    /** Overloaded input / output operator for Transitions. */
    istream & operator>>(istream & str, Grammar & g)
    {
      vector<Transition> transVec;
      str >> transVec;
      g = Grammar(transVec);

      return str;
    }


    Grammar readGrammar(string const & grammarFile, vector<string> const & singleEmitModelNames, vector<string> const & pairEmitModelNames)
    {
      ifstream f;
      openInFile(f, grammarFile);
      vector<Transition> outTransVec;
      f >> outTransVec;

      return Grammar( outTransVec, singleEmitModelNames, pairEmitModelNames );
    }


    Grammar readGrammar(string const & file)
    {
      ifstream f;
      openInFile(f, file);
      vector<Transition> outTransVec;
      f >> outTransVec;

      return Grammar(outTransVec);
    }
  

    ostream & operator<<(ostream & str, Grammar const & g)
    {
      str << "GRAMMAR_VERSION:\t1" << endl;
      str << endl;

      for (unsigned i = 0; i < g.outTransitions.size(); i++)
	for (unsigned j = 0; j < g.outTransitions[i].size(); j++) 
	  str << g.outTransitions[i][j];

      return str;
    }


    void writeGrammar(string const & file, Grammar const & g)
    {
      ofstream f;
      f.open( file.c_str() );
      f << g;
      f.close();
    }


    ////////////////////////////////////////////////////////////////
    // AnnoMap
    ////////////////////////////////////////////////////////////////

    EmitAnnoMap readEmitAnnoMap(string const & file)
    {
      ifstream f;
      f.open( file.c_str() );
      EmitAnnoMap annoMap;
      f >> annoMap;
      return annoMap;
    }

    void writeEmitAnnoMap(string const & file, EmitAnnoMap const & annoMap)
    {
      ofstream f;
      f.open( file.c_str() );
      f << annoMap;
      f.close();
    }

    // Overloaded input / output operator for EmitAnnoMaps. 
    istream & operator>>(istream & str, EmitAnnoMap & annoMap)
    {
      skipWhiteSpaceAndComments(str);
      while ( str.good() ) {
	EmitAnno emitAnno("","","","");
	while ( moreTags(str) ) {
	  string tag;
	  str >> tag;
	  if (tag == "NAME:")
	    str >> emitAnno.name;
	  else if (tag == "ANNO:")
	    str >> emitAnno.anno;
	  else if (tag == "ANNO_LEFT:")
	    str >> emitAnno.annoLeft;
	  else if (tag == "ANNO_RIGHT:")
	    str >> emitAnno.annoRight;
	  else
	    errorAbort("Unknown tag ('" + tag + "') in transition annotation map with name  '" + emitAnno.name + "'.");
	  skipLine(str); // skip rest of line
	}

	// sanity check
	if (not (emitAnno.name.size() > 0) )
	  errorAbort("Missing name in annotation map entry.");
	if (emitAnno.anno.size() > 0 and (emitAnno.annoLeft.size() > 0 or emitAnno.annoRight.size() > 0) )
	  errorAbort("Annotation map with name  '" + emitAnno.name + "' is wrongly defined.");
	if  ( (emitAnno.annoLeft.size() > 0 or emitAnno.annoRight.size() > 0) and not (emitAnno.annoLeft.size() * emitAnno.annoRight.size() > 0) ) // if one, then both ... 
	  errorAbort("Annotation map with name  '" + emitAnno.name + "' is wrongly defined.");

	annoMap[emitAnno.name] = emitAnno;

	skipWhiteSpaceAndComments(str);
      }

      return str;
    }


    ostream & operator<<(ostream & str, EmitAnnoMap const & annoMap)
    {
      for ( EmitAnnoMap::const_iterator it = annoMap.begin(); it != annoMap.end(); it++) {
	EmitAnno const & emitAnno = it->second;
	str << "NAME:\t" << emitAnno.name << endl;
	if (emitAnno.anno.size() > 0)
	  str << "ANNO:\t" << emitAnno.anno << endl;    
	else {
	  str << "ANNO_LEFT:\t" << emitAnno.annoLeft << endl;    
	  str << "ANNO_RIGHT:\t" << emitAnno.annoRight << endl;    
	}
	str << endl;
      }
      return str;
    }

    
    TransAnnoMap mkTransAnnoMap(string const & emitAnnoMapFile, Grammar const & g)
    {
      return TransAnnoMap(readEmitAnnoMap(emitAnnoMapFile), g.emitTransitions() );
    }


  } //namespace grammar

} //namespace phy
