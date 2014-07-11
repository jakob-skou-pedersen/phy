/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <fstream>
#include <boost/program_options.hpp>
#include "phy/utils.h"
#include "phy/utilsLinAlg.h"
#include "phy/GrammarIO.h"
#include "phy/GrammarWrap.h"
#include "phy/EmitWrapIO.h"

namespace po = boost::program_options;

using namespace phy;
using namespace phy::grammar;

////////////////////////////////////////////////////////////////
//                    Global constants
// Names of 'begin' non-terminals defining the structural and
// non-structural grammars

string const STR_BEGIN = "MS";
string const NON_STR_BEGIN = "E_NS";

struct PredRecord
{
  string id;                      // id
  unsigned begin;                 // begin of structure (0) / substructure (i)
  unsigned end;                   // end of structure (L) / substructure (j)
  unsigned bpCount;               // base pair count
  ynumber_t cykProbStr;           // CYK probabilty under str grammar
  ynumber_t cykProbNoStr;         // CYK probabilty under no-str str grammar
  ynumber_t insideProbStr;        // inside probabilty under str grammar
  ynumber_t insideProbNoStr;      // inside probabilty under no-str grammar
  number_t cykScore;              // log-odds score betweek cyk prob of str and no-str grammars
  number_t score;                 // log-odds score between str and no-str grammars
  number_t strPostProb;           // posterior probability of enitire cyk structure pred
  string fold;                    // secondary structure in parenthesis format
  vector_t posScore;              // position specific scores (posterior proababilities of predicted annotation)

  // string representation
  vector<string> toStrVec(unsigned prec1, unsigned prec2) const;
  string toString(unsigned prec1, unsigned prec2) const;
};

vector<string> PredRecord::toStrVec(unsigned prec1, unsigned prec2) const {
  vector<string> v;
  v.push_back(id);
  v.push_back( ::toString(begin) );
  v.push_back( ::toString(end) );
  v.push_back( ::toString(bpCount) );
  v.push_back( ::toString(cykProbStr, prec1) );
  v.push_back( ::toString(cykProbNoStr, prec1) );
  v.push_back( ::toString(insideProbStr, prec1) );
  v.push_back( ::toString(insideProbNoStr, prec1) );
  v.push_back( ::toString(cykScore, prec1) );
  v.push_back( ::toString(score, prec1) );
  v.push_back( ::toString(strPostProb, prec1) );
  v.push_back(fold);
  v.push_back( toSepString( toStdVector(posScore), ",", prec2) );
  return v;
}

string PredRecord::toString(unsigned prec1, unsigned prec2) const {
  return ::toSepString( toStrVec(prec1, prec2), "\t", 0);
}

vector<SeqData> readAlignments(istream & algFile)
{
  vector<ama> amaVec;
  algFile >> amaVec;
  return amaToSeqData(amaVec);
}

/** Return number of base pairs in fold */
unsigned countBasePairs(string const & fold, char leftPart = '(')
{
  return count(fold.begin(), fold.end(), leftPart);
}


void defineSubStructureIntervals(string const & fold, vector<unsigned> & begin, vector<unsigned> & end)
{
  unsigned stck = 0;
  for (unsigned i = 0; i < fold.size(); i++) {
    char const & c = fold[i];
    if (c == '(') {
      if (stck == 0)
	begin.push_back(i);
      stck++;
    }
    else if (c == ')') {
      stck--;
      if (stck == 0)
	end.push_back(i + 1);
    }
  }
  assert (begin.size() == end.size() );
}

PredRecord mkStrPredRecord(Grammar const & g, string const & entryId, unsigned begin, unsigned end, string const & fold, vector_t const & posScore)
{
    PredRecord pr;
    pr.id                = entryId;
    pr.begin             = begin;
    pr.end               = end;
    pr.cykProbStr        = g.subSeqCykProb(g.stateIndex(STR_BEGIN), begin, end);
    pr.cykProbNoStr      = g.subSeqCykProb(g.stateIndex(NON_STR_BEGIN), begin, end);
    pr.insideProbStr     = g.subSeqInsideProb(g.stateIndex(STR_BEGIN), begin, end);
    pr.insideProbNoStr   = g.subSeqInsideProb(g.stateIndex(NON_STR_BEGIN), begin, end);
    pr.cykScore          = (pr.cykProbNoStr != 0) ? log(pr.cykProbStr / pr.cykProbNoStr) : INFINITY;
    pr.score             = (pr.insideProbNoStr != 0) ? log(pr.insideProbStr / pr.insideProbNoStr) : INFINITY;
    if (pr.cykProbStr !=0)
      if (pr.insideProbStr != 0)
	pr.strPostProb       = toNumber(pr.cykProbStr / pr.insideProbStr);
      else // denom = 0
	pr.strPostProb = INFINITY;
    else
      pr.strPostProb = 0;
    pr.fold              = string(fold, begin, end - begin);
    pr.bpCount           = countBasePairs(fold);
    pr.posScore.resize(end - begin);
    for (unsigned i = 0; i < end - begin; i++)
      pr.posScore[i] = posScore[i + begin];
    return pr;
}

void predStructure(vector<PredRecord> & subStrPredVec, vector<PredRecord> & completeStrPredVec, bool calcComplStr, vector<SeqData> const & seqDataVec, Grammar & g, EmitWrappers & ew, EmitAnnoMap & eam, string const & annoName)
{
  TransAnnoMap am(eam, g.emitTransitions() );
  map<string, vector<string> > annoEquiMap = mkEquivalenceMap(g, am);

  for (unsigned i = 0; i < seqDataVec.size(); i++) {
    SeqData const & seqData = seqDataVec[i];

    // setup data structures for emissions probs (could be done outside loop if grammar accepted seqLength argument!)
    calcAndResetEmissions(seqData, g, ew, eam, annoName);

    // call main algorithms
    g.cyk();
    vector<EmitInfo> vei = g.maxParse();
    g.insideOutside();
    string fold = am.convertToString(vei);
    vector_t posScore = toNumber( g.postProbParse(vei, annoEquiMap) );

    // substructures
    vector<unsigned> begin;
    vector<unsigned> end;
    defineSubStructureIntervals(fold, begin, end);
    for (unsigned i = 0; i < begin.size(); i++) {
      string entryId = seqData.entryId + "." + toString(i);
      subStrPredVec.push_back( mkStrPredRecord(g, entryId, begin[i], end[i], fold, posScore) );
    }
    
    // complete structures
    if (calcComplStr) {
      completeStrPredVec.push_back( mkStrPredRecord(g, seqData.entryId, 0, seqData.seqSize(), fold, posScore) );
    }
  }
}


void printPredRecordVec(vector<PredRecord> & strPredVec, ostream & str, unsigned prec)
{
  str << "# seqId\tbeginPos\tendPos\tbasePairCount\tstrCykProb\tbgCykProb\tstrProb\tbgProb\tcykScore\tscore\tstrPostProb\tfold\tposScore" << endl;
  for (unsigned i = 0; i < strPredVec.size(); i++) {
    PredRecord const & pr = strPredVec[i];
    str << pr.toString(prec, 2) << endl; 
  }
}


int main(int argc, char * argv[])
{
  // option and argument variables
  string alignmentFile, treeFile, configFilePath, annoName, completeFile, outputFile; 
  unsigned prec;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("alignmentFile", po::value<string>(& alignmentFile), "Alignment file in ama format.")
    ("treeFile", po::value<string>(& treeFile), "Tree in newick format.");

  // define help message and options
  po::options_description visible(string("EvoFold predicts and scores secondary structures in mulitple sequence alignments.\n\n")
				  + "  Usage: EvoFold [options] <alg.ama> <tree.neiwck>\n\n"
				  + "The tabular output contains the following columns:\n\n"
				  + "  seqId beginPos endPos basePairCount strCykProb bgCykProb strProb bgProb cykScore score strPostProb fold posScore\n\n"
				  + "Setting alg.ama to '-' reads from stdin.\n" 
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("configFilePath,c", po::value<string>(& configFilePath)->default_value("./"), "Path to EvoFold configuration files.")
    ("completeFile,f", po::value<string>(& completeFile)->default_value(""), "Output complete structure predictions for each input element in addition to the sub-structures.")
    ("annoName,n", po::value<string>(& annoName)->default_value(""), "Name of annotation to use (see annoMap file for definition of annotation symbols. Note that '*' can be used as wildcard. Specifying annotation is useful for adding constraints on the predicted structure.")
    ("decimals", po::value<unsigned>(& prec)->default_value(5), "Output precision of score.")
    ("outputFile,o", po::value<string>(& outputFile)->default_value("-"), "Output file (default is stdout).");
  
  // setting up options parser
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  po::positional_options_description p;
  p.add("alignmentFile", 1);
  p.add("treeFile", 1);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  // print help message
  if (vm.count("help")) {
    cout << visible << endl;
    return 1;
  }

  // check arguments
  if (vm.count("alignmentFile") != 1 or vm.count("treeFile") != 1)
    errorAbort("\nWrong number of arguments. Try -h for help");

//  // debug
//  cerr << "configFilePath:\t" << configFilePath << endl
//       << "completeFile:\t" << subStrFile << endl
//       << "outputFile:\t" << outputFile << endl
//       << "alignmentFile:\t" << alignmentFile << endl
//       << "treeFile:\t" << treeFile << endl;

  // read alignment 
  ifstream algF;
  if (alignmentFile != "-")
    openInFile(algF, alignmentFile);
  istream & algFile = ( algF.is_open() ) ? algF : cin;
  vector<SeqData> seqData = readAlignments(algFile);

  // read tree
  string const treeString = readNewickStr(treeFile);

  // open output file
  ofstream outF;
  if (outputFile != "-")
    openOutFile(outF, outputFile);
  ostream & outFile = ( outF.is_open() ) ? outF : cout;


  ////////////////////////////////////////////////////////////////
  // Default configuration file names
  configFilePath += "/";
  string grammarStrFile       = configFilePath + "EvoFoldGrammarStr.txt";
  string emitModels           = configFilePath + "EvoFoldEmitModels.txt";
  string annoMapFile          = configFilePath + "EvoFoldAnnoMap.txt";

  EmitWrappers ew;
  readEmitWrappers(emitModels, ew, treeString);

  Grammar g       = readGrammar( grammarStrFile, ew.vecWrapNames(), ew.matWrapNames() );
  EmitAnnoMap eam = readEmitAnnoMap(annoMapFile);

  vector<PredRecord> subStrPredVec;
  vector<PredRecord> completeStrPredVec;
  bool calcComplStr = (completeFile != "") ? true : false;
  predStructure(subStrPredVec, completeStrPredVec, calcComplStr, seqData, g, ew, eam, annoName);

  // output
  printPredRecordVec(subStrPredVec, outFile, prec);
  if (calcComplStr) {
    ofstream f;
    openOutFile(f, completeFile);
    printPredRecordVec(completeStrPredVec, f, prec);
    f.close();
  }

//   if (vm.count("treeFile")) {
//     cout << "Tree file was defined as: " 
// 	 << vm["treeFile"].as<string>() << ".\n";
//   } 
// 
//   if (vm.count("subStrFile")) {
//     cout << "Sub-structure file was defined as: " 
// 	 << vm["subStrFile"].as<string>() << ".\n";
//   } 

  return 0;
}
