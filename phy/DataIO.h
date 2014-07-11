/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __DataIO_h
#define __DataIO_h

#include "phy/utils.h"


////////////////////////////////////////////////////////////////
// This header defines data input / output rutines 
////////////////////////////////////////////////////////////////


// to do: 
// Make the names data class a general interface to data. 
// Properly make it a base class, with a file interface as one of the derived classes.
// Various functions could then take the base class as argument, giving flexibility in their use.

namespace phy {

  // free functions

  /** Write a line of named data to ofstream. Output will be of this form:

  data_id   val_1   val_2    val_3

  */
  template<typename T>
  void writeNamedData(ostream & str, string const & id, vector<T> const & v, unsigned prec = 0, string const & sep = "\t");

  /** Same as above, but for vectors of corresponding IDs and values. Will add a header line based on the names of the form:

  NAME:   id_val_1   id_val_2    id_val_3

  */
  template<typename T>
  void writeNamedData(ostream & str, vector<string> const & names, vector<string> const & idVec, vector< vector<T> > const & dataVec, unsigned prec = 0, string const & sep = "\t");

  template<typename T> 
  void getAllNamedData(string const & file, vector<string> & allNames, vector<string> & idVec, vector<vector<T> > & dataVec);


  // classes

  /** Class for handling generic types of names data.
      The expected format is:

      NAME:      id_val_1     id_val_2
      id_data1   val_1.1      val_1.2
      id_data2   val_2.1      val_2.2

      The first line defines variable names of input data and all
      subsequent lines the actual input data, starting with an id.
      Note that 'NAME:' is a keyword.
  */
  template <typename T> 
  class NamedData {
  public:

    /** Constructor that takes the input data file as well as the
	complete list of variable names as input. This is useful when
	only a subset of the variables are specified in input
	file. The allNames then defines an enumeration of all
	variables, and the mapping onto these can be retrieved using
	the map() member function. */
    NamedData(string const & file, vector<string> const & allNames);

    /** Constructor that assumes that all variables are defined in the
	input file. In this case, map() will return an identity
	function. */
    NamedData(string const & file);

    /** Reset class*/
    void reset(string const & file, vector<string> const & allNames);
    void reset(string const & file);

    /** Returns next line of data as a pair where first member is the
	data id and the second a vector of symbols. Returns an empty
	id string on end of file. */
    pair< string, vector<T> > const next();
    /** Same as above, but returns output by reference. v must be of correct size. Returns false on failure. */
    bool next(string & id, vector<T> & v);

    /** Returns (ordered) vector of the input variable names (as defined in data file) */
    vector<string> const & names() const {return names_;}

    /** Returns the mapping of the input variables to the allNames vector given in the constructor */
    vector<unsigned> const & map() const {return map_;}

    /** Returns name count */
    unsigned const & count() const {return count_;}


  protected:

    // helper functions
    /** initialize data structures. Called by constructors.*/
    void init(vector<string> const & allNames);
    void init();
    vector<string> const parseElemNameLine();
    bool parseDataLine(string & id, vector<T> & v);

    // data
    vector<string> names_;
    vector<unsigned> map_;
    unsigned count_;

    // input fstream
    ifstream str_;
  };


  // definition of templated functions and members

  template<typename T> 
  void getAllNamedData(string const & file, vector<string> & allNames, vector<string> & idVec, vector<vector<T> > & dataVec)
  {
    NamedData<T> namedData(file);
    allNames = namedData.names();

    string id;
    vector<T> v( namedData.count() );
  
    while ( namedData.next(id, v) ) { 
      idVec.push_back(id);
      dataVec.push_back(v);
    }
  }


  template<typename T>
  void writeNamedData(ostream & str, string const & id, vector<T> const & v, unsigned prec, string const & sep)
  {
    if (prec)
      str.precision(prec);
    str << id;
    for (unsigned i = 0; i < v.size(); i++)
      str << sep << v[i];
    str << endl;
  }


  template<typename T>
  void writeNamedData(ostream & str, vector<string> const & names, vector<string> const & idVec, vector< vector<T> > const & dataVec, unsigned prec, string const & sep)
  {
    if (dataVec.size > 0)
      assert( names.size() == dataVec[0].size() );

    writeNamedData(str, "NAME:", names, prec, sep);
    assert( idVec.size() == dataVec.size() );
    for (unsigned i = 0; i < idVec.size(); i++)
      writeNamedData(str, idVec[i], dataVec[i], prec, sep);
  }


  template <typename T> 
  NamedData<T>::NamedData(string const & file, vector<string> const & allNames) 
  {
    openInFile(str_, file);
    init(allNames);
  }


  template <typename T> 
  NamedData<T>::NamedData(string const & file) 
  {
    openInFile(str_, file);
    init();
  }


  template <typename T> 
  void NamedData<T>::reset(string const & file, vector<string> const & allNames)
  {
    str_.close();
    openInFile(str_, file);
    init(allNames);
  }


  template <typename T> 
  void NamedData<T>::reset(string const & file)
  {
    str_.close();
    openInFile(str_, file);
    init();
  }


  template <typename T> 
  pair< string, vector<T> > const NamedData<T>::next()
  {
    string id;
    vector<T> v(count_);
    parseDataLine(id, v);
    return pair<string, vector<T> >(id, v);
  }


  template <typename T> 
  bool NamedData<T>::next(string & id, vector<T> & v)
  {
    return parseDataLine(id, v);
  }


  template <typename T> 
  void NamedData<T>::init(vector<string> const & allNames)
  {
    names_ = parseElemNameLine();    
    map_   = mkSubsetMap(allNames, names_);
    count_ = names_.size();
  }


  template <typename T> 
  void NamedData<T>::init()
  {
    names_ = parseElemNameLine();    
    map_   = mkSubsetMap(names_, names_);
    count_ = names_.size();
  }


  template <typename T> 
  vector<string> const NamedData<T>::parseElemNameLine()
  {
    skipWhiteSpaceAndComments(str_);

    string tag;
    str_ >> tag;
    if (tag != "NAME:")
      errorAbort("From NamedData::parseElemNameLine: Named data specification lacks 'NAME:' on fist non-comment line. Found tag: '" + tag + "'.");

    // fetch var names
    string namesStr;
    getline(str_, namesStr);
    return split( strip(namesStr) );
  }


  template <typename T> 
  bool NamedData<T>::parseDataLine(string & id, vector<T> & v)
  {
    // id
    skipWhiteSpaceAndComments(str_);
    if ( not str_.good() )
      return false;
    str_ >> id;
    
    // data elements
    string s;
    getline(str_, s);
    vector<string> u = split( strip(s) );

    if (u.size() != count_)
      errorAbort("From NamedData::parseDataLine: Corrupted data line in named data input file. Number of elements (" 
		 + toString( u.size() ) 
		 + ") must match number of names given in header line (" 
		 + toString(count_) 
		 + "). Note that elements are separated by white-spaces, hence no whitespace allowed in elements spcification.\n"
		 + "Corrupted information follows: \nid:\t" 
		 + id + "\nelement data:" + s + "\n");

    fromString(u, v);
    return true;
  }




} // end namespace phy

#endif  //__DataIO_h
