/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/GrammarTrain.h"

namespace phy {

  namespace grammar {

    // helper functions

    /** Creates a compositeFactorSet with a factor for each
	transitionDistribution of the grammar. pseudoCounts defines
	the total number of pseudoCounts for each transition
	distribution. */
    CompositeFactorSet mkTransFactorSet(Grammar & g, number_t totalPseudoCounts)
    {
      unsigned n = g.outTransitions.size();
      std::vector<AbsBasFacPtr_t> factorPtrs(n);

      for (unsigned i = 0; i < n; i++) {
	unsigned m = g.outTransitions[i].size();
	matrix_t mat(m, 1);
	for (unsigned j = 0; j < m; j++)
	  mat(j, 0) = toNumber(g.outTransitions[i][j].p);
	matrix_t pseudoCounts = mat * totalPseudoCounts;
	factorPtrs[i] = AbsBasFacPtr_t(new ColumnNormFactor(g.outTransitions[i][0].from_id + ".dist", mat, pseudoCounts) );
      }

      return CompositeFactorSet(factorPtrs);
    }


    /** extract column number "col" from each matrix in vm and return vector of vectors */
    vector<vector_t>  matVecToVecVec(vector<matrix_t> const & vm, unsigned col) 
    {
      vector<vector_t> vv;
      for (unsigned i = 0; i < vm.size(); i++) {
	assert(vm[i].size2() > col);
	vv.push_back( ublas::column(vm[i], col) );
      }
      return vv;
    }


    /** convert vector of vectors to vector of single-column matrices */
    vector<matrix_t>  vecVecToVecMat(Grammar::vector2D_t const & vv) 
    {
      vector<matrix_t> vm;
      for (unsigned i = 0; i < vv.size(); i++) {
	unsigned n = vv[i].size();
	ymatrix_t m(n, 1);
	vectorToMatrix(m, vv[i]);
	vm.push_back( toNumber(m) );
      }
      return vm;
    }


    void grammarEm(vector<SeqData> const & seqDataVec, Grammar & g, EmitWrappers & ew, EmitAnnoMap & eam, string const & annoName, number_t totalPseudoCounts, number_t minDeltaLogLik, unsigned maxIter, string const & logFile, string const & tmpOutputFile)
    {
      unsigned iter = 0;
      number_t logLik = 0;
      number_t prevLogLik = 0;
      number_t deltaLogLik = 0;

      CompositeFactorSet transitionDist = mkTransFactorSet(g, totalPseudoCounts);
      vector<matrix_t> expCountsMat     = transitionDist.mkFactorVec();     // Expectation counts
      Grammar::vector2D_t expCountsVec;  // size adjusted on first call of calcAccumulatedMarginals

      // open log file
      ofstream f;
      if ( logFile.size() ) {
	f.open(logFile.c_str(), ios::out);
	if (!f)
	  errorAbort("Cannot open file: " + logFile + "\n");
      
	f << "Annotation use:" << endl;;
	f << "Using annotation: " << ( ( eam.size() and annoName.size() ) ? "yes" : "no" ) << endl;
	f << "annoMap size:     " << eam.size() << endl;
	f << "annoName:         " << annoName << endl;
	f << endl;

	f << "Stopping criteria:" << endl;
	f << "minDeltaLogLik:   " << minDeltaLogLik << endl;
	f << "maxIter:          " << maxIter << endl;
	f << endl;

	f << "iter" << "\t" << "logLik" << "\t" << "deltaLogLik" << endl;
      }

      while ( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter) {
	logLik = 0;
	reset(expCountsMat);
	for (unsigned i = 0; i < seqDataVec.size(); i++) {
	  SeqData const & seqData = seqDataVec[i];
	  
	  calcAndResetEmissions(seqData, g, ew, eam, annoName);
	  ynumber_t likelihood = g.insideOutside();
	  g.calcAccumulatedMarginals(expCountsVec); // note that counts are NOT added, but overwritten
	  transitionDist.submitCounts( vecVecToVecMat(expCountsVec) );
	  
	  logLik += - log(likelihood);
	}
	transitionDist.optimizeParameters();
	g.resetTransitionProbs( matVecToVecVec(transitionDist.mkFactorVec(), 0) );  // convert matrices to vectors
	transitionDist.clearCounts();

	deltaLogLik = abs(prevLogLik - logLik); // set at abs(-logLik) in iter 0

	if ( logFile.size() )
	  f << iter << "\t" << logLik << "\t" << deltaLogLik << endl;

	// print grammar in each iter to tmp output file
	ofstream tmpOutF;
	openOutFile(tmpOutF, tmpOutputFile);
	tmpOutF << "# iter no. = " << iter << endl;
	tmpOutF << g;
	tmpOutF.close();

	prevLogLik = logLik;
	iter++;
      }
    }


  } // end namespace grammar

} // end namespace phy
