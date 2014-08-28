/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <boost/program_options.hpp>
#include "phy/DfgIO.h"
#include "phy/DfgDataWrap.h"

namespace po = boost::program_options;
using namespace phy;

int main(int argc, char * argv[])
{
  // option and argument variables
  string dfgSpecPrefix, stateMapsFile, factorPotentialsFile, variablesFile, factorGraphFile;  // specification files
  string varFile, facFile, logFile, outSpecPrefix, tmpSpecPrefix, dotFile; // input / ouput files
  unsigned prec, maxIter;
  number_t minDeltaLogLik;
  bool emTrain, writeInfo;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("varFile", po::value<string>(& varFile), "Input variable file in named data format.")
    ("facFile", po::value<string>(& facFile), "Input factor file in named data format. Must use same identifiers in same order as varFile.");

  // define help message and options
  po::options_description visible(string("dfgTrain allows implementation of discrete factor graphs and evaluates the probability of data sets under these models.\n\n")
				  + "  Usage: dfgTrain [options] <inputVarData.tab> [inputFacData.tab]\n\n"
				  + "The arguments inputVarData.tab and inputFacData.tab are both in named data format.\n\n"
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("precision,p", po::value<unsigned>(& prec)->default_value(5), "Output precision of real numbers.")
    ("minDeltaLogLik,d", po::value<number_t>(& minDeltaLogLik)->default_value(1e-4), "Defines stopping criteria for the EM training. The training will stop when the difference in log likelihood is below minDeltaLogLik (default is 1e-4).")
    ("maxIter,i", po::value<unsigned>(& maxIter)->default_value(100), "Max numbr if iterations of the EM training (default is 100).")
    ("logFile,l", po::value<string>(& logFile)->default_value("logFile.txt"), "Log file for EM training.")
    ("emTrain,e", po::bool_switch(& emTrain)->default_value(false), "Perform EM training.")
    ("dotFile", po::value<string>(& dotFile)->default_value(""), "Output dfg in dot format to given fileName. (To convert to ps format, e.g. run: \"cat fileName.dot | dot -Tps -Kneato -Gsize=\"7.5,10\" -o dfg.ps\".)")
    ("dfgSpecPrefix,s", po::value<string>(& dfgSpecPrefix)->default_value("./dfgSpec/"), "Prefix of DFG specification files.")
    ("outSpecPrefix,o", po::value<string>(& outSpecPrefix)->default_value("out_"), "Prefix of DFG specification files. Any included prefix directory must already exist.")
    ("tmpSpecPrefix,t", po::value<string>(& tmpSpecPrefix)->default_value(""), "Prefix of DFG specification files written during each iteration of training. Allows state of EM to be saved - especially useful for large datasets. Any included prefix directory must already exist. Not defined and not performed by default. ")
    ("factorGraphFile", po::value<string>(& factorGraphFile)->default_value("factorGraph.txt"), "Specification of the factor graph structure.")
    ("variablesFile", po::value<string>(& variablesFile)->default_value("variables.txt"), "Specification of the state map used by each variable.")
    ("stateMapFile", po::value<string>(& stateMapsFile)->default_value("stateMaps.txt"), "Specification of state maps.")
    ("facPotFile", po::value<string>(& factorPotentialsFile)->default_value("factorPotentials.txt"), "Specification of factor potentials.")
    ("writeInfo", po::bool_switch(& writeInfo)->default_value(false), "Print factor graph info. Useful for debugging factor graph specification.");

  // setting up options parser
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  po::positional_options_description p;
  p.add("varFile", 1).add("facFile", -1);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  // print help message
  if (vm.count("help")) {
    cout << visible << endl;
    return 1;
  }

  // check arguments
  if (vm.count("varFile") != 1)
    errorAbort("\nWrong number of arguments. Try -h for help");


  // read dfg
  DfgInfo dfgInfo = readDfgInfo(dfgSpecPrefix + stateMapsFile, 
				dfgSpecPrefix + factorPotentialsFile, 
				dfgSpecPrefix + variablesFile, 
				dfgSpecPrefix + factorGraphFile);

  // evaluate and output according to options
  if ( emTrain ) {
    if ( tmpSpecPrefix.size() )
	dfgEm(dfgInfo, varFile, facFile, minDeltaLogLik, maxIter, 
	      tmpSpecPrefix + stateMapsFile, 
	      tmpSpecPrefix + factorPotentialsFile, 
	      tmpSpecPrefix + variablesFile, 
	      tmpSpecPrefix + factorGraphFile, 
	      logFile);
    else
      dfgEm(dfgInfo, varFile, facFile, minDeltaLogLik, maxIter, logFile);

    // write final dfg version0
    writeDfgInfo(dfgInfo,
		 outSpecPrefix + stateMapsFile, 
		 outSpecPrefix + factorPotentialsFile, 
		 outSpecPrefix + variablesFile, 
		 outSpecPrefix + factorGraphFile);
  }

  writeDfgInfo(dfgInfo,
	       outSpecPrefix + stateMapsFile, 
	       outSpecPrefix + factorPotentialsFile, 
	       outSpecPrefix + variablesFile, 
	       outSpecPrefix + factorGraphFile);

  if ( dotFile.size() )
    dfgInfo.dfg.writeDot(dotFile);

  if (writeInfo){
    dfgInfo.dfg.writeInfo(cerr, dfgInfo.varNames, dfgInfo.facNames);
    dfgInfo.writeInfo(cerr);
  }
  return 0;
}
