/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "phy/DataIO.h"
#include "phy/StatDist.h"
#include "phy/xdoubleMod.h"

namespace po = boost::program_options;
using namespace phy;
using namespace NTL;

template<typename T>
ublas::matrix<T> toRowMat(vector<T> v)
{
  ublas::matrix<T> m( 1, v.size() );
  for (unsigned i = 0; i < v.size(); i++)
    m(0, i) = v[i];
  return m;  
}

template<typename T>
ublas::matrix<T> toColMat(vector<T> v)
{
  ublas::matrix<T> m( v.size(), 1);
  for (unsigned i = 0; i < v.size(); i++)
    m(i, 0) = v[i];
  return m;  
}


void writeOutputHeader(ostream & outStr, vector<string> const & parIdVec, string const & outputFormat)
{
  if ( outputFormat.size() ) {
    if (outputFormat == "vector" )
      outStr << "NAME:\t" << toString(parIdVec) << endl;
    else if  (outputFormat == "rowMat" )
      outStr << "NAME:\t" << toString( toRowMat(parIdVec) ) << endl;
    else if  (outputFormat == "colMat" )
      outStr << "NAME:\t" << toString( toColMat(parIdVec) ) << endl;
    else
      errorAbort("From main: unknown outputFormat: '" + outputFormat + "'.");
  }
  else
    writeNamedData(outStr, "NAME:", parIdVec);
}


void writeOutput(ostream & outStr, string const & id, vector<xdoubleMod> const & u, string const & outputFormat, unsigned prec)
{
  if ( outputFormat.size() ) {
    if (outputFormat == "vector" )
      outStr << id << "\t" << toString(u, prec) << endl;
    else if  (outputFormat == "rowMat" )
      outStr << id << "\t" << toString(toRowMat(u), prec) << endl;
    else if  (outputFormat == "colMat" )
      outStr << id << "\t" << toString(toColMat(u), prec) << endl;
  }
  else // default
    writeNamedData(outStr, id, u, prec);
}


void checkParameters(vector<string> const & parIdVec, vector< vector<xdoubleMod> > const & parVec)
{
  for (unsigned i = 0; i < parVec.size(); i++) {
    xdoubleMod sum = 0;
    for (unsigned j = 0; j < parVec[i].size(); j++)
      sum += parVec[i][j];
    if ( abs( sum - 1.0 ) > 1e-5)
      errorAbort("form main: parameters for multinomial distribution with name '" + parIdVec[i] + "' sum to " + toString(sum) + " instead of 1.0 (deviation > 1e-5).");
  }
}


int main(int argc, char * argv[])
{
  // option and argument variables
  string parameterFile, countsFile, outputFile, outputFormat;
  unsigned prec;
  bool coefficients, takeLog;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("parameterFile", po::value<string>(& parameterFile), "Input parameters.")
    ("countsFile", po::value<string>(& countsFile), "Input counts.")
    ("outputFile", po::value<string>(& outputFile), "Output probabilities or coefficients.");

  // define help message and options
  po::options_description visible(string("multinomial calculates probabilities of category counts under the multinomial distribution.\n\n")
				  + "  Usage: multinomial [options] <parameters.tab> <counts.tab> <output.tab> \n\n"
				  + "By default, all input and output files use named data format.\n\n"
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("precision,p", po::value<unsigned>(& prec)->default_value(5), "Output precision of real numbers.")
    ("coefficients,c", po::bool_switch(& coefficients)->default_value(false), "Output coefficients instead of probabilities.")
    ("outputFormat,f", po::value<string>(& outputFormat)->default_value(""), "Use alternative output format. Possible values are: vector, rowMat, and colMat, which all use ublas style formatting.")
    ("logarithm,l", po::bool_switch(& takeLog)->default_value(false), "Output natural logarithm of result values.");
  
  // setting up options parser
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  po::positional_options_description p;
  p.add("parameterFile", 1).add("countsFile", 1).add("outputFile", 1);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  // print help message
  if (vm.count("help")) {
    cout << visible << endl;
    return 1;
  }

  // check arguments
  if (vm.count("parameterFile") != 1 or vm.count("countsFile") != 1 or vm.count("outputFile") != 1)
    errorAbort("\nWrong number of arguments. Try -h for help");

  
  // read parameters
  vector<string> parNames;
  vector<string> parIdVec;
  vector< vector<xdoubleMod> > parVec;
  getAllNamedData(parameterFile, parNames, parIdVec, parVec);

  // init input counts
  NamedData<unsigned> countsData(countsFile);
  
  // checks
  if (countsData.names() != parNames)
    errorAbort("From main: Names defined in header line of parameters file and counts file must match");
  checkParameters(parIdVec, parVec);

  // open output stream
  ofstream f;
  if (outputFile != "-")
    openOutFile(f, outputFile);
  ostream & outStr = ( f.is_open() ) ? f : cout;
  xdoubleMod::SetOutputPrecision(prec); // set precision of xdouble

  // write output header
  writeOutputHeader(outStr, parIdVec, outputFormat);

  // loop over data and output
  string id;
  vector<unsigned> v( countsData.count() );
  vector<xdoubleMod> u( parVec.size() );
  while ( countsData.next(id, v) ) { 
    for (unsigned i = 0; i < parVec.size(); i++)

      if (coefficients)
	u[i] = multinomialCoefficient<xdoubleMod>(v);
      else
	u[i] = multinomialDistribution(parVec[i], v);

    if (takeLog)
      for (unsigned i = 0; i < parVec.size(); i++)
	u[i] = log(u[i]);

    // output
    writeOutput(outStr, id, u, outputFormat, prec);
  }
  
  return 0;
}
