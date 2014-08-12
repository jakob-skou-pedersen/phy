/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/Factors.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/numeric/ublas/io.hpp> /* Used in print function mostly used for debugging could be removed later */

namespace phy {

  AbstractFullyParameterizedFactor::AbstractFullyParameterizedFactor(string const & type, string const & id, matrix_t const & m, matrix_t const & pseudoCounts) 
    : AbstractBaseFactor(type, id, m.size1(), m.size2()), m_(m), pseudoCounts_(pseudoCounts)
  {
    if (pseudoCounts.size1() != 0) {
      assert(pseudoCounts.size1() == size1_); 
      assert(pseudoCounts.size2() == size2_); 
    }
  }


  int GlobalNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    number_t n = sumMatrix(m_);
    if (n == 0.0)
      errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible.");
    m_ *= 1.0 / n;
    return 1;
  }


  int ColumnNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    for (unsigned i = 0; i < size2_; i++) {
      // normalize each column
      number_t n = sum( ublas::column(m_, i) );
      if (n == 0.0)
	errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible for column i=" + toString(i) + ".");
      column(m_, i) *= 1.0 / n;
    }      
    return 1;
  }
	

  int RowNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    for (unsigned i = 0; i < size1_; i++) {
      // normalize each row
      number_t n = sum( row(m_, i) );
      if (n == 0.0)
	errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible for row i=" + toString(i) + ".");
      row(m_, i) *= 1.0 / n;
    }      
    return 1;
  }

  void NormalFactor::init(){
    //TODO improve on decision of midpoints for outermost bins
    midpoints_ = vector_t(breakpoints_.size()+1);
    midpoints_(0) = breakpoints_(0);
    midpoints_(breakpoints_.size()) = breakpoints_(breakpoints_.size()-1);

    for(unsigned i = 1; i < breakpoints_.size(); ++i)
      midpoints_(i) = (breakpoints_(i-1)+breakpoints_(i))/2;

    calcPotentials();
  }

  int NormalFactor::optimizeParametersImpl()
  {
    //TODO Stabilitity of singel pass variance calculation
    //TODO Check if maximum likelihood based estimates are better
    number_t S = inner_prod(row(counts_,0), midpoints_);
    number_t USS = inner_prod(row(counts_,0), element_prod(midpoints_, midpoints_));
    number_t n = sum(row(counts_,0));
    mean_ = S/n;
    var_ = (USS-S*S/n)/(n-1); 
    return 1;
  }

  void NormalFactor::calcPotentials() 
  {
    //Loop over breakpoints and set potentials
    boost::math::normal norm(mean_, sqrt(var_) );
    number_t prevCdf = 0;
    for(int i = 0; i < breakpoints_.size(); ++i){
      m_(0,i) = cdf(norm, breakpoints_(i)) - prevCdf;
      prevCdf = cdf(norm, breakpoints_(i));
    }
    m_(0,breakpoints_.size()) = 1-prevCdf;   std::cout << "DEBUG l36 Factors.cpp" << std::endl;
  }

  void NormalFactor::print()
  {
    std::cout << "Mean:\t" << mean_ << std::endl;
    std::cout << "Variance:\t" << var_ << std::endl;
    std::cout << "Matrix m:\t" << m_ << std::endl;
    std::cout << "Breakpoints:\t" << breakpoints_ << std::endl;
    std::cout << "Midpoints:\t" << midpoints_ << std::endl;
    std::cout << "Matrix counts:\t" << counts_ << std::endl;
  }


  // return matrix defining factor idx 
  vector<matrix_t> AbstractBaseFactorSet::mkFactorVec() const
  {
    vector<matrix_t> v; 
    for (unsigned i = 0; i < facCount_; i++) 
      v.push_back(mkFactor(i) );
    return v;
  } 


  // set matrix defining factor idx 
  void AbstractBaseFactorSet::mkFactorVec(vector<matrix_t> & v) const 
  {
    for (unsigned i = 0; i < facCount_; i++) 
      mkFactor(v[i], i);
  }
  

} // end namespace phy
