/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <boost/program_options.hpp>
#include "phy/DfgIO.h"

namespace po = boost::program_options;
using namespace phy;

void takeMinusLog(xvector_t & v, string const & id)
{
  for (unsigned i = 0; i < v.size(); i++) {
    if (v[i] < 0)
      errorAbort("From takeMinusLog: Trying to logarithm of negative number resulting from data line with id: '" + id + "'.");
    v[i] = - log(v[i]);
  }
}


// useful with post probs close to one, where precision is lost  
void ppSumOther(xvector_t & v)
{
  xvector_t u( v.size(), 0);
  for (unsigned i = 0; i < v.size(); i++) {
    u[i] = 0;
    for (unsigned j = 0; j < v.size(); j++)
      if (i != j)
	u[i] += v[j];
  }
  v = u;
}


void transformByOptions(xvector_t & v, bool minusLogarithm, bool sumOther, string const & id)
{
  if (sumOther)
    ppSumOther(v);
  if (minusLogarithm)
    takeMinusLog(v, id);
}


// check identity of ids
void checkIds(string const & idVar, string const & idFac, unsigned lineCount)
{ 
  if (idVar != idFac)
    errorAbort("From main: problem with input number " + toString(lineCount) + ": idFac '" + idFac + "' and idVar '" + idVar + "' differ.");
}


// if facDataPtr not NULL, then read next line of potentials and reset dfg
void resetFactorPotential(FacData * facDataPtr, string const id, unsigned const lineCount, DFG & dfg)
{
  if (facDataPtr != NULL) {
    vector<xmatrix_t> facVec( facDataPtr->count() );
    string idFac;
    facDataPtr->next(idFac, facVec);
    checkIds(id, idFac, lineCount);
    dfg.resetFactorPotentials( facVec, facDataPtr->map() );
    dfg.consistencyCheck();
  }
}

void updateFactorPotentials(VarData * subVarDataPtr, string const id, unsigned const lineCount, DfgInfo & dfgInfo){
  if ( subVarDataPtr != NULL){
    vector<symbol_t> subVarVec( subVarDataPtr->count() );
    string idSubVar;
    subVarDataPtr->next(idSubVar, subVarVec);
    checkIds(id, idSubVar, lineCount);
    dfgInfo.updateFactors( subVarVec, subVarDataPtr->invMap() );
  }
}

// generate 2D vector of state symbols according to vector of stateMaps
vector< vector<string> > mkStateSymbolTable(vector<StateMapPtr_t> stateMapVec)
{ 
  vector< vector<string> > ssTable( stateMapVec.size() );
  for (unsigned i = 0; i < stateMapVec.size(); i++) {
    StateMap const & sm = *stateMapVec[i];
    for (unsigned j = 0; j < sm.stateCount(); j++)
      ssTable[i].push_back( sm.state2Symbol(j) );
  }
  return ssTable;
}
  

// convert vector of states to vector of symbols
void stateVecToSymbolVec(vector<StateMapPtr_t> const & stateMapVec, vector<state_t> const & maxVarStates, vector<symbol_t> & maxVarSymbols)
{
  assert( maxVarStates.size() == maxVarSymbols.size() );
  assert( maxVarStates.size() == stateMapVec.size() );
  for (unsigned i = 0; i < maxVarStates.size(); i++)
    maxVarSymbols[i] = stateMapVec[i]->state2Symbol( maxVarStates[i] );
}


// parse ppVarVecStr, which is of the form "X = a b c; Y = a b", and
// return variables (X and Y in this case) and states (a, b, and c in
// this case) as commonly indexed vectors.
void mkVarAndStateSymbolList(string const & varSpecStr,  vector<string> & varNames, vector< vector<symbol_t> > & varStates)
{
  varNames.clear();
  varStates.clear();
  vector<string> specs = split( strip(varSpecStr, "; \t"), ";");
  BOOST_FOREACH(string & s, specs) {
    // check for empty specs
    if (s.size() == 0 or strip(s).size() == 0)
      errorAbort("From mkVarAndStateSymbolList: Empty variable specification '" + s + "' found in this line:\n" + varSpecStr + "\n");

    vector<string> v = split(s, "=");
  
    //check
    if (not (v.size() == 1 or v.size() == 2) )
      errorAbort("From mkVarAndStateSymbolList: Error in specification of variable and state list.:\n" + varSpecStr + "\n");

    // add var name
    varNames.push_back( strip(v[0]) );

    // add var states (or none)
    vector<symbol_t> states;
    if (v.size() == 2) // states defined
      states = split( strip( v[1] ) );
    varStates.push_back(states);
  }
}


void writePostProbLegend(ostream & str, vector<string> const & varNames, vector< vector<symbol_t> > const & varStates)
{
  str << "#Definition of state order for each variable:" << endl;
  str << "#" << "NAME\t" << "ranVar" << "\tstate order ..." << endl;
  for (unsigned i = 0; i < varNames.size(); i++) {
    str << "#";
    writeNamedData(str, "name\t" + varNames[i], varStates[i]);
  }
}


int main(int argc, char * argv[])
{
  // option and argument variables
  string dfgSpecPrefix, stateMapsFile, factorPotentialsFile, variablesFile, factorGraphFile;  // specification files
  string varFile, facFile, subVarFile, postProbFile, normConstFile, maxProbStateFile; // input / ouput files
  string mpsVarVecStr, ppVarVecStr; // other options
  unsigned prec;
  bool minusLogarithm, ppSumOther;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("varFile", po::value<string>(& varFile), "Input variable file in named data format.")
    ("facFile", po::value<string>(& facFile), "Input factor file in named data format. Must use same identifiers in same order as varFile.");

  // define help message and options
  po::options_description visible(string("dfgEval allows implementation of discrete factor graphs and evaluates the probability of data sets under these models.\n\n")
				  + "  Usage: dfgEval [options] <inputVarData.tab> [inputFacData.tab]\n\n"
				  + "The arguments inputVarData.tab and inputFacData.tab are both in named data format.\n"
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("ppFile,o", po::value<string>(& postProbFile)->default_value(""), "Calculate posterior probabilities for each state of each random variable and output to file.")
    ("ncFile,n", po::value<string>(& normConstFile)->default_value(""), "Calculate normalization constant output to file.")
    ("mpsFile,m", po::value<string>(& maxProbStateFile)->default_value(""), "Calculate most probable state for each random variable and output to file.")
    ("precision,p", po::value<unsigned>(& prec)->default_value(5), "Output precision of real numbers.")
    ("ppSumOther", po::bool_switch(& ppSumOther)->default_value(false), "For post probs, for each state output sum of post probs for all the other states for that variable. This retains precision for post probs very close to one.")
    ("minusLogarithm,l", po::bool_switch(& minusLogarithm)->default_value(false), "Output minus the natural logarithm of result values (program will terminate on negative results...).")
    ("mpsVars", po::value<string>(& mpsVarVecStr)->default_value(""), "Define the random variables for which the most probable state (mps) should be output. Default is to output the mps for all random variables. The specification string must be enclosed in citation marks and whitespace separated if it includes more than one random variable, e.g.: \"X Y\".")
    ("ppVars", po::value<string>(& ppVarVecStr)->default_value(""), "Define the random variables for which the posterior state probabilities (pp) should be calculated. Default is to output the pp for all states of all random variables (may generate much output!). Random variables are specified similar to mpsVars, but must be semicolon (';') separated. It is possible to only output pp's for certain states, in which case the following specification format is used: \"X=a b c; Y=a b\".")
    ("dfgSpecPrefix,s", po::value<string>(& dfgSpecPrefix)->default_value("./dfgSpec/"), "Prefix of DFG specification files..")
    ("factorGraphFile", po::value<string>(& factorGraphFile)->default_value("factorGraph.txt"), "Specification of the factor graph structure.")
    ("variablesFile", po::value<string>(& variablesFile)->default_value("variables.txt"), "Specification of the state map used by each variable.")
    ("stateMapFile", po::value<string>(& stateMapsFile)->default_value("stateMaps.txt"), "Specification of state maps.")
    ("facPotFile", po::value<string>(& factorPotentialsFile)->default_value("factorPotentials.txt"), "Specification of factor potentials.")
    ("subVarFile", po::value<string>(& subVarFile), "Input subscribed variables file in named data format. Must use same identifiers in same order as varFile");
  
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

  // set output precision at this point in case of xdouble type
#ifdef XNUMBER_IS_XDOUBLE
  xnumber_t::SetOutputPrecision(prec);
#endif

  // read dfg
  DfgInfo dfgInfo = readDfgInfo(dfgSpecPrefix + stateMapsFile, 
				dfgSpecPrefix + factorPotentialsFile, 
				dfgSpecPrefix + variablesFile, 
				dfgSpecPrefix + factorGraphFile);

  // define convenient bools
  bool calcMps =  maxProbStateFile.size() ? true : false;
  bool calcNc  =  normConstFile.size() ? true : false;
  bool calcPp  =  postProbFile.size() ? true : false;
  bool calcExpect = true;
//  // debug
//  cout << endl << endl;
//  cout << "calcMps\t" << calcMps << endl;
//  cout << "calcNc\t" << calcNc << endl;
//  cout << "calcPp\t" << calcPp << endl;
//  cout << "dfgSpecPrefix:\t\t" << dfgSpecPrefix  << endl;
//  cout << "stateMapsFile:\t\t" << stateMapsFile  << endl;
//  cout << "factorPotentialsFile:\t\t" << factorPotentialsFile  << endl;
//  cout << "variablesFile:\t\t" << variablesFile  << endl;
//  cout << "factorGraphFile:\t\t" << factorGraphFile  << endl;
//  cout << "varFile:\t\t" << varFile  << endl;
//  cout << "facFile:\t\t" << facFile  << endl;
//  cout << "postProbFile:\t\t" << postProbFile  << endl;
//  cout << "normConstFile:\t\t" << normConstFile  << endl;
//  cout << "maxProbStateFile:\t\t" << maxProbStateFile  << endl;

  // setup input data structures
  VarData varData(varFile, dfgInfo.varNames);
  FacData * facDataPtr = NULL;
  VarData * subVarDataPtr = NULL;
  if (facFile.size() != 0) {
    facDataPtr = new FacData(facFile, dfgInfo.facNames);
  }
  if (subVarFile.size() != 0){
    subVarDataPtr = new VarData(subVarFile, dfgInfo.subNames);
  }
  else if( dfgInfo.subNames.size() > 0 ){
    errorAbort("dfgEval.cpp::main There are subscribed factors but no subscribed variable file were provided");
  }

  // open output streams
  // max prob state
  ofstream mpsF;
  if (calcMps and maxProbStateFile != "-")
    openOutFile(mpsF, maxProbStateFile);
  ostream & mpsStr = ( mpsF.is_open() ) ? mpsF : cout;

  // norm const
  ofstream ncF;
  if (calcNc and normConstFile != "-")
    openOutFile(ncF, normConstFile);
  ostream & ncStr = ( ncF.is_open() ) ? ncF : cout;

  // post prob
  ofstream ppF;
  if (calcPp and postProbFile != "-")
    openOutFile(ppF, postProbFile);
  ostream & ppStr = ( ppF.is_open() ) ? ppF : cout;

  // variables needed in data loop
  string idVar;
  vector<symbol_t> varVec( varData.count() );
  stateMaskVec_t stateMasks( dfgInfo.varNames.size() );
  
  // mps init data structures 
  vector<unsigned> maxVarStates( dfgInfo.varNames.size() );
  vector<symbol_t> maxVarSymbols( dfgInfo.varNames.size() );
  // mps output random variable subset
  vector<string> mpsVarNames = split(strip(mpsVarVecStr) );
  if (mpsVarNames.size() == 0)
    mpsVarNames = dfgInfo.varNames;
  vector<unsigned> mpsVarMap = mkSubsetMap(dfgInfo.varNames, mpsVarNames);
  vector<symbol_t> mpsSymVec( mpsVarNames.size() );
 
  // pp init data structures 
  vector<string>  ppVarNames;
  vector< vector<symbol_t> >  ppVarStates;
  mkVarAndStateSymbolList(ppVarVecStr,  ppVarNames, ppVarStates); 
  // pp output data structures
  vector<xvector_t> variableMarginals;
  initGenericVariableMarginals(variableMarginals, dfgInfo.dfg);
  // define which random variables and which states to output post probs for
  assert( ppVarNames.size() == ppVarStates.size() );
  if (ppVarNames.size() == 0) {
    ppVarNames = dfgInfo.varNames;
    ppVarStates = vector< vector<symbol_t> >( ppVarNames.size() ); // reset varStates as well
  }
  vector<unsigned> ppVarMap = mkSubsetMap(dfgInfo.varNames, ppVarNames);
  // state maps
  vector< vector<string> > ssTable = mkStateSymbolTable(dfgInfo.stateMapVec); 
  vector< vector<unsigned> > ppVarStateMap;
  for (unsigned i = 0; i < ppVarStates.size(); i++) {
    if (ppVarStates[i].size() == 0)
      ppVarStates[i] = ssTable[ ppVarMap[i] ];
    ppVarStateMap.push_back( mkSubsetMap( ssTable[ ppVarMap[i] ], ppVarStates[i] ) );
  }

  // write headers
  // mps
  if (calcMps)
    writeNamedData(mpsStr, "NAME:\tprob", mpsVarNames);

  // nc
  if (calcNc)
    ncStr << "NAME:\t" << "normConst" << endl;

  // pp
  if (calcPp) {
    if (ppVarNames.size() == 1) 
      writeNamedData(ppStr, "NAME:\tranVar", ppVarStates[0]);
    else {
      writePostProbLegend(ppStr, ppVarNames, ppVarStates);
      ppStr << "NAME:\t" << "ranVar\t" << "..." << endl;
    }
  }

  // data loop: calc and output
  unsigned lineCount = 1;
  while ( varData.next(idVar, varVec) ) {
    // read and reset factor potentials
    resetFactorPotential(facDataPtr, idVar, lineCount, dfgInfo.dfg);
    updateFactorPotentials(subVarDataPtr, idVar, lineCount, dfgInfo);
    dfgInfo.stateMaskMapSet.symbols2StateMasks(stateMasks, varVec, varData.map());

    if (calcMps) {
      xnumber_t p = dfgInfo.dfg.runMaxSum(stateMasks, maxVarStates);
      if (minusLogarithm)
	p = - log(p);
      stateVecToSymbolVec(dfgInfo.stateMapVec, maxVarStates, maxVarSymbols);
      mkSubset(maxVarSymbols, mpsVarMap, mpsSymVec); // map maxVarStates to the requested subset of symbols
      writeNamedData(mpsStr, idVar + "\t" + toString(p, prec), mpsSymVec);
    }

    if (calcNc or calcPp) 
      dfgInfo.dfg.runSumProduct(stateMasks);  

    if (calcNc) {
      xnumber_t p = dfgInfo.dfg.calcNormConst2(stateMasks);
      if (minusLogarithm)
	p = - log(p);
      ncStr << idVar << "\t" << toString(p, prec) << endl;
    }

    if (calcPp) {
      dfgInfo.dfg.calcVariableMarginals(variableMarginals, stateMasks);
      for (unsigned i = 0; i < ppVarNames.size(); i++) {
	xvector_t ppVec = variableMarginals[ ppVarMap[i] ];
	transformByOptions(ppVec, minusLogarithm, ppSumOther, idVar);
	writeNamedData(ppStr, idVar + "\t" + ppVarNames[i], mkSubset(toStdVector(ppVec), ppVarStateMap[i]), prec);
      }
    }

    if ( calcExpect ){
      //Start with function_a equal to potentials
      dfgInfo.dfg.runSumProduct(stateMasks);
      xnumber_t expect = dfgInfo.dfg.calcExpect();
      xnumber_t expect2 = dfgInfo.dfg.calcExpect2(stateMasks);
      std::cout << idVar << "\tExpect:\t" << toString(expect, prec) << "\tExpect2:\t" << toString(expect2, prec) << std::endl;
      //Use counts at each node
    }

    lineCount++;
  }

  // clean up
  if (facDataPtr != NULL)
    delete facDataPtr;

  return 0;
}
