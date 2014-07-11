/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __ama_h
#define __ama_h

#include "phy/utils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <stdlib.h>

namespace phy {

  using namespace std;

  /** Sequence structure holding the information given in the maf format (http://eyebrowse.cit.nih.gov/genome/goldenPath/help/maf.html) */
  struct amaSeq {
  public:

    /** constructor */
    amaSeq(string const & src, long start, long size, char strand, long srcSize, string const & text);

    /**  Parse string l and create amaSeq object */
    amaSeq(string const & l);

    string src;
    long start;
    long size;    
    char strand;  
    long srcSize; 
    string text;    
  };


  /** annotated multiple alignment format (modeled over the maf format:
      http://genome.ucsc.edu/FAQ/FAQformat#format5
      (or originally: http://eyebrowse.cit.nih.gov/genome/goldenPath/help/maf.html.) */
  struct ama {
  public:

    /** pair type. */
    typedef pair<string, string> keyText_t;

    /** vector type. */
    typedef vector<keyText_t> pairVector_t;

    /** sequence vector. */
    typedef vector<amaSeq> seqVector_t;

    /** lAnno text type. */
    typedef vector<string> lAnno_t;

    /** lAnno pair  type . */
    typedef pair<string, lAnno_t> lAnnoPair_t;

    /** lAnno pair vector type . */
    typedef vector<lAnnoPair_t> lAnnoPairVector_t;

    /** Constructor. */
    ama(pairVector_t const & features   = pairVector_t(),
	seqVector_t const & sequences   = seqVector_t(),
	pairVector_t const & anno       = pairVector_t(), 
	lAnnoPairVector_t const & lAnno = lAnnoPairVector_t(), 
	vector<string> const & comments = vector<string>() );
  
    /** Returns true in ama is empty               */     
    bool empty();

    /** Clear all values               */     
    void clear();

    /** vector of tag and value pairs               */     
    pairVector_t features;

    /** vector of name and text pairs		   */
    seqVector_t sequences;

    /** vector of name and text pairs		   */
    pairVector_t anno;

    /** vector of name and long label list pairs	   */
    lAnnoPairVector_t lAnno;

    /** vector of comments                          */
    vector<string> comments;

    bool hasFeature(string const & key) const;
    string getFeature(string const & key) const;

  };

  /** Overloaded output operator for ama objects. */
  ostream &operator<<(ostream &str, ama const & entry);

  /** Overloaded input operator for ama objects. */
  istream &operator>>(istream &str, ama & entry);

  /** Overloaded output operator for vector of ama objects. */
  ostream &operator<<(ostream &str, vector<ama> const & amaVector);

  /** Overloaded input operator for vector of ama objects. */
  istream &operator>>(istream &str, vector<ama> & amaVector);

  /** File input / output */
  vector<ama> readAmaFile(string const & fileName);
  void writeAmaFile(string const & fileName, vector<ama> const & amaVec);

}

#endif // __ama_h
