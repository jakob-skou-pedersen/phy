/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/DfgDataWrap.h"

// debug
#include <boost/numeric/ublas/io.hpp>

namespace phy {

  // Set up appropriate data structures and calculate the DFG normalization constant for each position 
  vector<xnumber_t> strVec2NormConstVector(DFG & dfg, 
					   vector<string> const & strVec, 
					   vector<unsigned> const & varMap, 
					   StateMaskMapSet const & stateMaskMapSet, 
					   char missingDataChar)
  {
    long seqLength = strVec[0].size();
    unsigned ranVarCount = dfg.variables.size();
    stateMask2DVec_t stateMask2DVec = initStateMask2DVec(seqLength, ranVarCount);
    mkStateMask2DVec(strVec, stateMask2DVec, varMap, stateMaskMapSet, missingDataChar);
    return calcNormConsMultObs(stateMask2DVec, dfg);
  }

  // Same as above, but intended for the specific case where all random variables are defined over the same state space. The stateMaskMapSet will thus be defined in terms of a single stateMap. 
  vector<xnumber_t> strVec2NormConstVector(DFG & dfg, vector<string> const & strVec, vector<unsigned> const & varMap, StateMap const & staMap, char missingDataChar)
  {
    unsigned ranVarCount = dfg.variables.size();
    StateMaskMapSet stateMaskMapSet(staMap, ranVarCount);
    return strVec2NormConstVector(dfg, strVec, varMap, stateMaskMapSet, missingDataChar);
  }


  void dfgEm(AbstractBaseFactorSet & factorSet, vector<stateMask2DVec_t> const & dataVec, DFG & dfg, number_t minDeltaLogLik, unsigned maxIter, string const & logFile)
  {
    unsigned iter = 0;
    number_t logLik = 0;
    number_t prevLogLik = 0;
    number_t deltaLogLik = 0;

    vector<matrix_t> facExpCounts;     // Expectation counts
    initAccFactorMarginals(facExpCounts, dfg);
    xvector_t normConstVec;  // will be resized by calcFacAccMarAndNormConstMultObs

    // open log file
    ofstream f;
    if ( logFile.size() ) {
      f.open(logFile.c_str(), ios::out);
      if (!f)
	errorAbort("Cannot open file: " + logFile + "\n");
      
      f << "minDeltaLogLik: " << minDeltaLogLik << endl;
      f << "maxIter:        " << maxIter << endl;
      f << endl;
      f << "iter" << "\t" << "logLik" << "\t" << "deltaLogLik" << endl;
    }

    while ( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter) {
      logLik = 0;
      reset(facExpCounts);
      for (unsigned i = 0; i < dataVec.size(); i++) {
	stateMask2DVec_t const & stateMask2DVec = dataVec[i];
	calcFacAccMarAndNormConstMultObs(facExpCounts, normConstVec, stateMask2DVec, dfg);
	logLik +=  - accLog(normConstVec);
      }
      factorSet.submitCounts(facExpCounts);
      factorSet.optimizeParameters();
      dfg.resetFactorPotentials( factorSet.mkFactorVec() );
      factorSet.clearCounts();

      deltaLogLik = abs(prevLogLik - logLik); // set at abs(-logLik) in iter 0

      if ( logFile.size() )
	f << iter << "\t" << logLik << "\t" << deltaLogLik << endl;
      //
      //      // debug
      //      cout << "from dfgEm:" << endl;
      //      cout << "iter:        " << iter << endl;
      //      cout << "logLik:      " << logLik << endl;
      //      cout << "prevLogLik:  " << prevLogLik << endl;
      //      cout << "deltaLogLik: " << deltaLogLik << endl;
      //      cout << "minDeltaLogLik: " << minDeltaLogLik << endl;
      //      cout << "maxIter: " << maxIter << endl;
      //      cout << "( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter): " << ( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter) << endl;
      //
      prevLogLik = logLik;
      iter++;
    }
  }

  // to do: refactorize such that the functions below are defined in a header and reused by dfgEval
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


  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDataFile, number_t minDeltaLogLik, unsigned maxIter, string const & logFile)
  {
    dfgEm(dfgInfo, varDataFile, facDataFile, minDeltaLogLik, maxIter, "", "", "", "", logFile);
  }

  void dfgEm(DfgInfo & dfgInfo, string const & varDataFile, string const & facDataFile, number_t minDeltaLogLik, unsigned maxIter, 
	     string const & logStateMapsFile, string const & logFactorPotentialsFile, string const & logVariablesFile, string const & logFactorGraphFile, string const & logFile)
  {
    DFG & dfg = dfgInfo.dfg;     // convenient
    AbstractBaseFactorSet & factorSet = dfgInfo.facSet;

    unsigned iter = 0;
    number_t logLik = 0;
    number_t prevLogLik = 0;
    number_t deltaLogLik = 0;

    // setup input data structures
    VarData varData(varDataFile, dfgInfo.varNames);
    FacData * facDataPtr = NULL;
    if (facDataFile.size() != 0) {
      facDataPtr = new FacData(facDataFile, dfgInfo.facNames);
    }

    // init variables
    // input data variables
    string idVar;
    vector<symbol_t> varVec( varData.count() ); // for input data
    // dfg variables    
    stateMaskVec_t stateMaskVec( dfgInfo.varNames.size() );
    vector<matrix_t> facExpCounts; // Expectation counts
    initAccFactorMarginals(facExpCounts, dfg);
    vector<xmatrix_t> tmpFacMar;  // workspace
    dfg.initFactorMarginals(tmpFacMar);
    xnumber_t normConst;

    // open log file
    ofstream f;
    if ( logFile.size() ) {
      f.open(logFile.c_str(), ios::out);
      if (!f)
	errorAbort("Cannot open file: " + logFile + "\n");
      
      f << "minDeltaLogLik: " << minDeltaLogLik << endl;
      f << "maxIter:        " << maxIter << endl;
      f << endl;
      f << "iter" << "\t" << "logLik" << "\t" << "deltaLogLik" << endl;
    }

    while ( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter) {
      logLik = 0;
      reset(facExpCounts);

      // reset varData and facData
      varData.reset(varDataFile, dfgInfo.varNames);
      if (facDataFile.size() != 0) {
	facDataPtr->reset(facDataFile, dfgInfo.facNames);
      }

      unsigned lineCount = 0;
      while ( varData.next(idVar, varVec) ) {
	lineCount++;
	resetFactorPotential(facDataPtr, idVar, lineCount, dfgInfo.dfg);
	dfgInfo.stateMaskMapSet.symbols2StateMasks( stateMaskVec, varVec, varData.map() );

	//	//debug begin
	//	string facName("CpG_P_1.likelihood");
	//	unsigned facIdx = find(dfgInfo.facNames.begin(), dfgInfo.facNames.end(), facName) - dfgInfo.facNames.begin();
	//	cout << "factor name: " << facName << "; facIdx: " << facIdx << "; potential " << dfg.getFactor(facIdx).potential <<  endl;
	//	//debug end

	calcFacAccMarAndNormConst(facExpCounts, tmpFacMar, normConst, stateMaskVec, dfg);

	logLik += - log(normConst);
      }

      factorSet.submitCounts(facExpCounts);
      factorSet.optimizeParameters();
      dfg.resetFactorPotentials( factorSet.mkFactorVec() );
      factorSet.clearCounts();

      deltaLogLik = abs(prevLogLik - logLik); // set at abs(-logLik) in iter 0

      if ( logFile.size() )
	f << iter << "\t" << logLik << "\t" << deltaLogLik << endl;


      if ( logStateMapsFile.size() )
	writeDfgInfo(dfgInfo, logStateMapsFile, logFactorPotentialsFile, logVariablesFile, logFactorGraphFile);

      //
      //      // debug
      //      cout << "from dfgEm:" << endl;
      //      cout << "iter:        " << iter << endl;
      //      cout << "logLik:      " << logLik << endl;
      //      cout << "prevLogLik:  " << prevLogLik << endl;
      //      cout << "deltaLogLik: " << deltaLogLik << endl;
      //      cout << "minDeltaLogLik: " << minDeltaLogLik << endl;
      //      cout << "maxIter: " << maxIter << endl;
      //      cout << "( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter): " << ( ( (deltaLogLik > minDeltaLogLik) or iter == 0) and iter < maxIter) << endl;
      //

      //      //debug begin
      //      cout << endl << "iteration #: " << iter << endl << endl;
      //      //debug end

      prevLogLik = logLik;
      iter++;
    }

    // clean up
    if (facDataPtr != NULL)
      delete facDataPtr;
  }



} // end namespace phy
