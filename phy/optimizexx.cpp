/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "optimizexx.h"

# ifndef NO_OPTPP

//global pointers used by initWrapper and fctWrapper
boost::function< double (phy::vector_t const &) > const * glbObjFct;
phy::vector_t * glbInitParams;


void initWrapper(int ndim, NEWMAT::ColumnVector & x)
{
  assert( ndim == x.Nrows() );
  assert( ndim == static_cast<signed>( glbInitParams->size() ) );
  x = toColumnVector( *glbInitParams );
}


// wrap boost function in interface required by opt++
void fctWrapper(int ndim, const NEWMAT::ColumnVector & x, double & fx, int & result)
{
  fx = (*glbObjFct)( toNumber(x) );
  result = NLPFunction;
}


phy::vector_t const minimizexx(phy::vector_t x, boost::function< double (phy::vector_t const &) > const & objFct, std::string const & statusFileName, unsigned maxFEval, unsigned maxIter)
{
  glbObjFct = & objFct;
  glbInitParams = & x;

  FDNLF1 nlp(x.size(), fctWrapper, initWrapper);

  OptQNewton objfcn(&nlp);
  //OptFDNewton objfcn(&nlp);
  //OptCG objfcn(&nlp);

  objfcn.setSearchStrategy(LineSearch); //alt: TrustRegion
  // set stopping criteria: see http://csmr.ca.sandia.gov/opt%2B%2B/opt++2.4_doc/html/ControlParameters.html
  objfcn.setMaxFeval(maxFEval); 
  objfcn.setMaxIter(maxIter);

  // The "0" in the second argument says to create a new file.  A "1"
  // would signify appending to an existing file.
  if (statusFileName != "")
    if (!objfcn.setOutputFile(statusFileName.c_str(), 1) )
      cerr << "main: output file open failed" << endl;

  objfcn.optimize();

  char status[] = "Solution from quasi-newton"; // hack to avoid compiler warnings
  objfcn.printStatus(status);
  objfcn.cleanup();

  return toNumber( objfcn.getXPrev() );
}

NEWMAT::ColumnVector const toColumnVector(phy::vector_t const & v)
{
  unsigned n = v.size();
  NEWMAT::ColumnVector u(n);
  for (unsigned i = 0; i < n; i++)
    u[i] = v[i];
  return u;
}


phy::vector_t const toNumber(NEWMAT::ColumnVector const & v)
{
  unsigned n = v.Nrows();
  phy::vector_t u(n);
  for (unsigned i = 0; i < n; i++)
    u[i] = v[i];
  return u;
}

#else /* NO_OPTPP */

phy::vector_t const minimizexx(phy::vector_t x, boost::function< double (phy::vector_t const &) > const & objFct, std::string const & statusFileName, unsigned maxFEval, unsigned maxIter)
{
  std::cerr << "minimizezz cannot be called as compiled with NO_OPTPP preprocessing directive, which excludes the opt++ numeric optimization libraries. To include opt++ functionality recompile without the NO_OPTPP flag. Program termninates." << std::endl;
  exit(0);
}


#endif /* NO_OPTPP */


// map ]l;oo[ into ]-oo;oo[
double transform(double x, double l)
{
  return log(10*(x-l) );
}


// map ]-oo;oo[ back into ]l;oo[  
double deTransform(double y, double l)
{
  return l + exp(y)/10;
}


void transformVec(phy::vector_t & v, double l)
{
  for (unsigned i = 0; i < v.size(); i++)
    v[i] = transform(v[i], l);
}


void deTransformVec(phy::vector_t & v, double l)
{
  for (unsigned i = 0; i < v.size(); i++)
    v[i] = deTransform(v[i], l);
}

// 
// void InitObject::operator()(phy::vector_t & x) 
// {
//   assert( x.size() == x_.size() ); 
//   x = x_;
// }
// 
// 
// double FctObject::operator()(phy::vector_t const & x)
// {
//   return objFct_(x);
// }
// 


phy::vector_t const forceInRange(phy::vector_t const &x, double min_value, double max_value)
{
  int n     = x.size();
  phy::vector_t y(x); 

  for (int i = 0; i < n; i++) {
    if (x[i] < min_value)
      y[i] = min_value;
    if (x[i] > max_value)
      y[i] = max_value;
  }
  return y;
}
