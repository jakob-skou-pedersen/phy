/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __optimizexx_h
#define __optimizexx_h

#include <fstream>

#include "boost/function.hpp"
#include "phy/PhyDef.h"

# ifndef NO_OPTPP // allow compilation without installation of OPTPP:
#define SETUP_C_SUBSCRIPTS // allow zero-based []-indexation of newmat vectors and matrices
#include "newmat/newmatio.h"
#include "optpp/OptQNewton.h"
#include "optpp/Constraint.h"
#include "optpp/BoundConstraint.h"
#include "optpp/OptBCQNewton.h"
#include "optpp/OptBaQNewton.h"

using namespace OPTPP;
# endif /* NO_OPTPP */

/** global optimization */
phy::vector_t const minimizexx(phy::vector_t x, boost::function< double (phy::vector_t const &) > const & objFct, std::string const & statusFileName = "", unsigned maxFEval = 10000, unsigned maxIter = 1000);


////////////////////////////////////////////////////////////////
// useful transformations and conversions
////////////////////////////////////////////////////////////////

/** transform maps ]l;oo[ into ]-oo;oo[ */
double transform(double x, double l);
/** deTransforms maps ]-oo;oo[ into ]l;oo[ (inverse of transform for given l) */
double deTransform(double y, double l);

/** same as transform and deTransform, but acts on whole vectors. Mapping done in place. */
void transformVec(phy::vector_t & v, double l);
void deTransformVec(phy::vector_t & v, double l);

/** force x[i] to be in [min_value; max_value]. */
phy::vector_t const forceInRange(phy::vector_t const &x, double min_value, double max_value);

# ifndef NO_OPTPP
/** Conversion from NEWMAT::ColumnVector to phy::vector_t and visa versa. */
NEWMAT::ColumnVector const toColumnVector(phy::vector_t const & v);
phy::vector_t const toNumber(NEWMAT::ColumnVector const & v);
# endif /* NO_OPTPP */

#endif /* __optimizexx_h */
