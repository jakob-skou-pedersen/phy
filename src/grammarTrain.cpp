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
#include "phy/GrammarTrain.h"
#include "phy/EmitWrapIO.h"

namespace po = boost::program_options;

using namespace phy;
using namespace phy::grammar;

vector<SeqData> readAlignments(istream & algFile)
{
  vector<ama> amaVec;
  algFile >> amaVec;
  return amaToSeqData(amaVec);
}



int main(int argc, char * argv[])
{
  // option and argument variables
  string grammarFile, emitModelsFile, alignmentFile, treeFile, annoMapFile, annoName, outputFile, logFile, tmpGrammar; 
  number_t pseudoCounts, minDeltaLogLik;
  unsigned maxIter;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("grammarFile", po::value<string>(& grammarFile), "Grammar in version 1.0 format.")
    ("emitModelsFile", po::value<string>(& emitModelsFile), "File defining emit models.")
    ("alignmentFile", po::value<string>(& alignmentFile), "Alignment file in ama format.");

  // define help message and options
  po::options_description visible(string("grammarTrins estimate transition probabilities (production rule probabilities) of stochastic context free grammars. Currently ama type input data is supported.\n\n")
				  + "  Usage: grammarTrain [options] <grammar.txt> <emitModels.txt> <alg.ama>\n\n"
				  + "Setting alg.ama to '-' reads from stdin.\n" 
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("treeFile,t", po::value<string>(& treeFile)->default_value(""), "File with Newick tree used with phylo grammars.")
    ("annoMapFile,a", po::value<string>(& annoMapFile)->default_value(""), "Anno map file.")
    ("annoName,n", po::value<string>(& annoName)->default_value(""), "Name of annotation to use.")
    ("pseudoCounts,p", po::value<number_t>(& pseudoCounts)->default_value(0), "Defines total number of pseudocounts used for each transition distribution (For each transition, the number of pseudocounts is defined as the initital transition probs in the input file times the given pseudoCounts value).")
    ("minDeltaLogLik,d", po::value<number_t>(& minDeltaLogLik)->default_value(1e-4), "Defines stopping criteria for the EM training. The training will stop when the difference in log likelihood is below minDeltaLogLik (default is 1e-4).")
    ("maxIter,i", po::value<unsigned>(& maxIter)->default_value(100), "Max numbr if iterations of the EM training (default is 100).")
    ("logFile,l", po::value<string>(& logFile)->default_value("grammarLogFile.txt"), "Log file for EM training (default is ./grammarLogFile.txt).")
    ("outputGrammar,o", po::value<string>(& outputFile)->default_value("-"), "Output file for trained grammar (default is stdout).")
    ("tmpGrammar,o", po::value<string>(& tmpGrammar)->default_value("tmpGrammar.txt"), "Output file for partly trained grammar, printed in each interation (default is tmpGrammar.txt).");
  
  // setting up options parser
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  po::positional_options_description p;
  p.add("grammarFile", 1);
  p.add("emitModelsFile", 1);
  p.add("alignmentFile", 1);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  // print help message
  if ( vm.count("help") ) {
    cout << visible << endl;
    return 1;
  }

  // check arguments
  if (vm.count("grammarFile") != 1 or vm.count("emitModelsFile") != 1 or vm.count("alignmentFile") != 1)
    errorAbort("\nWrong number of arguments. Try -h for help");

  // input
  // read tree
  if ( not treeFile.size() )
    errorAbort("TreeFile not specified. Currently only phyloGrammars are supported, which need a tree to be defined.");
  string const treeString = readNewickStr(treeFile);

  EmitWrappers ew;
  readEmitWrappers(emitModelsFile, ew, treeString);
  Grammar g       = readGrammar( grammarFile, ew.vecWrapNames(), ew.matWrapNames() );
  EmitAnnoMap eam = readEmitAnnoMap(annoMapFile);

  // read alignment 
  ifstream algF;
  if (alignmentFile != "-")
    openInFile(algF, alignmentFile);
  istream & algFile = ( algF.is_open() ) ? algF : cin;
  vector<SeqData> seqDataVec = readAlignments(algFile);

  // train grammar
  grammarEm(seqDataVec, g, ew, eam, annoName, pseudoCounts, minDeltaLogLik, maxIter, logFile, tmpGrammar);

  // open output file
  ofstream outF;
  if (outputFile != "-")
    openOutFile(outF, outputFile);
  ostream & outFile = ( outF.is_open() ) ? outF : cout;

  // output
  outFile << g;

  return 0;
}
