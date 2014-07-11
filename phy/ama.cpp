/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "ama.h"

namespace phy {

  bool ama::hasFeature(string const & key) const
  {
    for (unsigned i = 0; i < features.size(); i++)
      if (features[i].first == key)
	return true;
    return false;
  }


  string ama::getFeature(string const & key) const
  {
    for (unsigned i = 0; i < features.size(); i++)
      if (features[i].first == key)
	return features[i].second;
    // only reaches here if no match
    errorAbort("From getFeature: No key mathcing: " + key);
    return "";  // to silence compiler
  }


  amaSeq::amaSeq(string const & src, long start, long size, char strand, long srcSize, string const & text)
    : src(src), start(start), size(size), strand(strand), srcSize(srcSize), text(text)
  {}


  amaSeq::amaSeq(string const & l)
  {
    vector<string> v = split(l);
    if (v.size() != 6)
      throw string("amaSeq: parse error: string has wrong number of elements :\n") + l;
    src     = v[0];
    start   = atol(v[1].c_str());
    size    = atol(v[2].c_str());
    strand  = v[3][0];
    srcSize = atol(v[4].c_str());
    text    = v[5];
  }


  ama::ama(pairVector_t const & features, 
	   seqVector_t const & sequences, 
	   pairVector_t const & anno, 
	   lAnnoPairVector_t const & lAnno, 
	   vector<string> const & comments)
    : features(features), sequences(sequences), anno(anno), lAnno(lAnno), comments(comments)
  {}


  void ama::clear()
  {
    features.clear();
    sequences.clear();
    anno.clear();
    lAnno.clear();
    comments.clear();
  }


  bool ama::empty()
  {
    if (features.empty() and  sequences.empty() and anno.empty() and lAnno.empty() and comments.empty())
      return true;
    else
      return false;
  }


  istream &operator>>(istream &str, ama & entry)
  {
    entry.clear();
    string currentLine;
    while ( getline(str, currentLine) ) 
      if (strip(currentLine).size())
	if (currentLine[0] != '#')
	  break;

    if (str.eof())
      return str;

    // parse features
    if (currentLine[0] != 'a')
      throw string("parse error in line at start of ama entry:\n") + currentLine;
    vector<string> splits = split(currentLine);
    for (string::size_type i = 1; i < splits.size(); i++) {
      if (splits[i].find('=') == string::npos)
	throw string("bad feature segment:\n") + splits[i];
      vector<string> tagAndValue = split(splits[i], "=");
      entry.features.push_back(ama::keyText_t(tagAndValue[0],tagAndValue[1]));
    }
    // parse the rest
  
    getline(str, currentLine);
    while( strip(currentLine).size() ) {
      char c = currentLine[0] ;
      if (c == 's') 
	entry.sequences.push_back( amaSeq( currentLine.substr(1, currentLine.size() - 1 ) ) );
      else if  (c == 'l') {
	splits = split(currentLine);
	entry.anno.push_back(ama::keyText_t(splits[1], splits[2]));
      }
      else if  (c == 'L') {
	splits             = split(currentLine);
	ama::lAnno_t text  = split(splits[2], ",");
	pair<string, ama::lAnno_t> lAnnoPair = pair<string, ama::lAnno_t>(splits[1], text);
	entry.lAnno.push_back(lAnnoPair);
      }
      else if  (c == '#') {
	string comment = currentLine.substr(1,currentLine.size()-1);
	entry.comments.push_back(comment);
      }

      if (not getline(str, currentLine) )
	break;
    }	 

    return str;
  }


  ostream &operator<<(ostream &str, ama const & entry)
  {
    str << "a ";
    for (unsigned i  = 0; i < entry.features.size(); i++)
      str << entry.features[i].first << "=" << entry.features[i].second << " ";
    str << endl;

    for (unsigned i  = 0; i < entry.sequences.size(); i++)
      str << "s\t" << 
	entry.sequences[i].src     << "\t"  << 
	entry.sequences[i].start   << "\t"  << 
	entry.sequences[i].size    << "\t"  << 
	entry.sequences[i].strand  << "\t"  << 
	entry.sequences[i].srcSize << "\t"  << 
	entry.sequences[i].text    << endl;

    for (unsigned i  = 0; i < entry.anno.size(); i++)
      str << "l\t" << entry.anno[i].first << "\t\t\t\t" << entry.anno[i].second << endl;

    for (unsigned i  = 0; i < entry.lAnno.size(); i++) {
      str << "L\t" << entry.lAnno[i].first << "\t\t\t\t";
      for (unsigned j  = 0; j < entry.lAnno[i].second.size(); j++)
	str << entry.lAnno[i].second[j] << ",";
      str << endl;
    }

    for (unsigned i  = 0; i < entry.comments.size(); i++)
      str << "# " << entry.comments[i] << endl;

    return str;
  }


  istream &operator>>(istream &str, vector<ama> & amaVector)
  {
    amaVector.clear();
    ama entry;
    for (str >> entry; not entry.empty(); str >> entry)
      amaVector.push_back(entry);
    return str;
  }


  ostream &operator<<(ostream &str, vector<ama> const & amaVector)
  {
    for (unsigned i = 0; i < amaVector.size(); i++)
      str << amaVector[i] << endl;;
    return str;
  }


  vector<ama> readAmaFile(string const & fileName)
  {
    ifstream f;
    f.open(fileName.c_str());
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");

    vector<ama> amaVec;
    f >> amaVec;
    f.close();
    return amaVec;
  }


  void writeAmaFile(string const & fileName, vector<ama> const & amaVec)
  {
    ofstream f;
    f.open(fileName.c_str());
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");

    f << amaVec;
    f.close();
  }



}
