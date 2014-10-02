/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/Factors.h"
#include "phy/utils.h"
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

  void AbstractFullyParameterizedFactor::serialize(ostream& os) const{
    os << "POT_MAT:\t" << m_ << endl;
    if ( pseudoCounts_.size1() != 0 )
      os << "PC_MAT:\t"  << pseudoCounts_ << endl;
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

  /** Constructor without parameter values, set to reasonable defaults */
  DiscContFactor::DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins ) : AbstractBaseFactor("discCont", name, states, bins), minv_(minv), maxv_(maxv), bins_(bins),states_(states) { 
    mixDist_ = MixPtr_t(new NormalMixture(states, minv, maxv, bins));
  }
  
  /** Constructor taking Mixture object */
  DiscContFactor::DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins, MixPtr_t mixDist ) : AbstractBaseFactor("discCont", name, states, bins), minv_(minv), maxv_(maxv), bins_(bins), states_(states), mixDist_(mixDist) {

  }

  void DiscContFactor::serialize(ostream & os) const
  {
    os << "BINS:\t" << bins_ << endl;
    os << "MIN:\t" << minv_ << endl;
    os << "MAX:\t" << maxv_ << endl;
    os << "STATES:\t" << states_ << endl;
    mixDist_->serialize(os);
  }

  int DiscContFactor::optimizeParametersImpl()
  {
    mixDist_->optimizeParameters(counts_);
  }

  void DiscContFactor::mkFactor(matrix_t &m) const
  {
    mixDist_->mkFactor(m);
  }

  void ContContFactor::serialize(ostream & os) const {
    ConvertIndexNumber is1(minv1_, maxv1_, bins1_);
    ConvertIndexNumber is2(minv2_, maxv2_, bins2_);

    os << "ALPHA:\t" << alpha_* is2.getToNumberAlpha()/is1.getToNumberAlpha() << endl;
    os << "BETA:\t" << beta_*is2.getToNumberAlpha()+is2.getToNumberBeta()-alpha_*is2.getToNumberAlpha()*is1.getToNumberBeta()/is1.getToNumberAlpha() << endl;
    os << "VAR:\t" << var_*is2.getToNumberAlpha()*is2.getToNumberAlpha() << endl << endl;
  }

  void ContContFactor::mkFactor(matrix_t &m) const{
    for(matrix_t::iterator1 it1 = m.begin1(); it1 != m.end1(); ++it1){
      boost::math::normal norm( it1.index1()*alpha_+beta_, sqrt(var_) );   
      for(matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	//TODO Lazy evaluation?
	*it2 = cdf(norm, it2.index2()) - cdf(norm, it2.index2()-1);
      }
    }
  }

  int ContContFactor::optimizeParametersImpl() {
    number_t SP_xy = 0, S_x = 0, S_y = 0, USS_x = 0, USS_y = 0, N = 0;
    for(matrix_t::iterator1 it1 = counts_.begin1(); it1 != counts_.end1(); ++it1){
      for(matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	SP_xy += it2.index1()*it2.index2()  * (*it2);
	S_x += it2.index1() * (*it2);
	S_y += it2.index2() * (*it2);
	USS_x += it2.index1()*it2.index1() * (*it2);
	USS_y += it2.index2()*it2.index2() * (*it2);
	N += *it2;
      }
    }
    
    number_t SSD_x = USS_x-S_x*S_x/N;
    number_t SSD_y = USS_y-S_y*S_y/N;
    number_t SPD_xy = SP_xy - S_x*S_y/N;
    alpha_ = SPD_xy/SSD_x;
    beta_ = (S_y-S_x*alpha_)/N;
    var_ = (SSD_y - SPD_xy*SPD_xy/SSD_x)/N;

    return 1;
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
