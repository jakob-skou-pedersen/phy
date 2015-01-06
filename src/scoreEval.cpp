/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <boost/program_options.hpp>
#include "phy/DfgIO.h"

#include <iterator>
#include <algorithm>
#include <utility>
#include <cmath>

namespace po = boost::program_options;
using namespace phy;

typedef std::pair< number_t, number_t> wSample;

int main(int argc, char * argv[])
{
  unsigned no_samples;
  number_t alpha; //Importance sampling tuning variable
  string dfgSpecPrefix, stateMapsFile, factorPotentialsFile1, factorPotentialsFile2, variablesFile, factorGraphFile;  // specification files

  po::options_description visible(string("dfgEval allows implementation of discrete factor graphs and evaluates the probability of data sets under these models.\n\n")
				  + "  Usage: dfgEval [options] <inputVarData.tab> [inputFacData.tab]\n\n"
				  + "The arguments inputVarData.tab and inputFacData.tab are both in named data format.\n"
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("samples,n", po::value<unsigned>(& no_samples)->default_value(10000), "No of samples used for methods naive sampling and importance sampling")
    ("alpha,a", po::value<number_t>(& alpha)->default_value(0.5), "Importance sampling tuning factor. Alpha=0 corresponds to naive sampling. The higher alpha the more the distribution is skewed towards high scores.")
    ("dfgSpecPrefix,s", po::value<string>(& dfgSpecPrefix)->default_value("./dfgSpec/"), "Prefix of DFG specification files..")
    ("factorGraphFile", po::value<string>(& factorGraphFile)->default_value("factorGraph.txt"), "Specification of the factor graph structure.")
    ("variablesFile", po::value<string>(& variablesFile)->default_value("variables.txt"), "Specification of the state map used by each variable.")
    ("stateMapFile", po::value<string>(& stateMapsFile)->default_value("stateMaps.txt"), "Specification of state maps.")
    ("facPotFile1", po::value<string>(& factorPotentialsFile1)->default_value("factorPotentials1.txt"), "Specification of factor potentials in Null model.")
    ("facPotFile2", po::value<string>(& factorPotentialsFile2)->default_value("factorPotentials2.txt"), "Specification of factor potentials in Positive model.")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(visible);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  po::notify(vm);    

  // print help message
  if (vm.count("help")) {
    cout << visible << endl;
    return 1;
  }


  DfgInfo dfgInfo1 = readDfgInfo(dfgSpecPrefix + stateMapsFile, 
				dfgSpecPrefix + factorPotentialsFile1, 
				dfgSpecPrefix + variablesFile, 
				dfgSpecPrefix + factorGraphFile);
  DfgInfo dfgInfo2 = readDfgInfo(dfgSpecPrefix + stateMapsFile, 
				dfgSpecPrefix + factorPotentialsFile2, 
				dfgSpecPrefix + variablesFile, 
				dfgSpecPrefix + factorGraphFile);

  //Generate marginal distributions used generating IS distribution and calculating score
  stateMaskVec_t stateMasks( dfgInfo1.dfg.variables.size() );
  dfgInfo1.dfg.runSumProduct(stateMasks);
  dfgInfo2.dfg.runSumProduct(stateMasks);
  number_t normConst1 = toNumber(dfgInfo1.dfg.calcNormConst2(stateMasks));
  number_t normConst2 = toNumber(dfgInfo2.dfg.calcNormConst2(stateMasks));

  //Naive sampling

  //Importance sampling
  {
    DfgInfo dfgInfoIS = dfgInfo1;

    //Calculate IS distribution
    vector<xmatrix_t> ISfacPotentials(dfgInfoIS.dfg.factors.size());
    for(int f = 0; f < dfgInfoIS.dfg.factors.size(); ++f){
      
      xmatrix_t const & pot1 = dfgInfo1.dfg.nodes[ dfgInfo1.dfg.convFacToNode(f) ].potential;
      xmatrix_t const & pot2 = dfgInfo2.dfg.nodes[ dfgInfo2.dfg.convFacToNode(f) ].potential;
      
      xmatrix_t potIS(pot1.size1(), pot1.size2());
      for(int i = 0; i < potIS.size1(); ++i){
	for(int j = 0; j < potIS.size2(); ++j){
	  potIS(i,j) = power(pot1(i,j), 1-alpha)*power(pot2(i,j), alpha);
	}
      }
      ISfacPotentials.at(f) = potIS;
    }

    dfgInfoIS.dfg.resetFactorPotentials( ISfacPotentials);
    dfgInfoIS.dfg.runSumProduct(stateMasks);
    dfgInfoIS.dfg.calcFactorMarginals();
    dfgInfoIS.dfg.calcVariableMarginals(stateMasks);
    dfgInfo1.dfg.calcFactorMarginals();
    dfgInfo1.dfg.calcVariableMarginals(stateMasks);


    //Setup
    boost::mt19937 gen(std::time(0));
    vector<wSample> samples(no_samples);
    const vector<vector_t> & varMarBg = dfgInfo1.dfg.getVariableMarginals();
    const vector<matrix_t> & facMarBg = dfgInfo1.dfg.getFactorMarginals();
    const vector<vector_t> & varMarIS = dfgInfoIS.dfg.getVariableMarginals();
    const vector<matrix_t> & facMarIS = dfgInfoIS.dfg.getFactorMarginals();

    std::cout << "Generating IS samples:" << endl;
    for(int i = 0; i < no_samples; ++i){
      //Sample
      vector<unsigned> sampleIS;
      number_t weight = 0.0;

      dfgInfo1.dfg.sampleIS(gen,
			    varMarBg, 
			    facMarBg, 
			    varMarIS, 
			    facMarIS, 
			    sampleIS, 
			    weight);

      //Calculate score
      number_t score1 = dfgInfo1.dfg.calcFullLikelihood(sampleIS)/normConst1;
      number_t score2 = dfgInfo2.dfg.calcFullLikelihood(sampleIS)/normConst2;
      number_t score = log(score2/score1);
      
      samples.at(i).first = score;
      samples.at(i).second = weight;
    }
    
    std::cout << "Sorting samples according to score:" << endl;
    std::sort(samples.begin(), samples.end() );
    
    //Report some quantiles
    //TODO make an option for which quantiles to report
    number_t normalized_score_sum = 0;
    number_t unnormalized_score_sum = 0;
    number_t sum_of_weights = 0;
    for(vector<wSample>::reverse_iterator it = samples.rbegin(); it != samples.rend(); ++it){
      if(it == samples.rbegin()){
	std::cout << "MAX:\t" << (*it).first << endl;
      }
      sum_of_weights += (*it).second;
      unnormalized_score_sum += (*it).first;
      normalized_score_sum += (*it).first*(*it).second;
      static bool q9999 = true;
      static bool q999  = true;
      static bool q99   = true;
      static bool q95   = true;

      if(sum_of_weights/no_samples > 0.0001 and q9999){
	q9999 = false;
	std::cout << "99.99%:\t" << (*it).first << endl;
      }
      if(sum_of_weights/no_samples > 0.001 and q999){
	q999 = false;
	std::cout << "99.9%:\t" << (*it).first << endl;
      }
      if(sum_of_weights/no_samples > 0.01 and q99){
	q99 = false;
	std::cout << "99.0%:\t" << (*it).first << endl;
      }
      if(sum_of_weights/no_samples > 0.05 and q95){
	q95 = false;
	std::cout << "95.0%:\t" << (*it).first << endl;
      }
    }
    
    std::cout << endl << "IS Diagnostics" << endl;
    std::cout << "Unnormalized score mean:\t" << unnormalized_score_sum/no_samples << endl;
    std::cout << "Normalized score mean:\t" << normalized_score_sum/sum_of_weights << endl;
  }
  

  return 0;
}
