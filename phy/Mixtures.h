/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Mixtures_h
#define __Mixtures_h

#include <phy/PhyDef.h>
#include <phy/utils.h>
#include <vector>
#include <stdio.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/shared_ptr.hpp>

/** 
    Mixtures contains distributions that can be used to model a finite mixture
    of some class of distriubution. Making it possible to store parametersets
    in a convenient format for each distribution.
    Notice cdf calculations are done in "binspace".
 */

namespace phy {

  class Mixture; // forward declaration

  /** Pointer to Mixture */
  typedef boost::shared_ptr<Mixture> MixPtr_t;

  class Mixture {
  public:
    Mixture(number_t minv, number_t maxv, unsigned bins, unsigned states) : minv_(minv), maxv_(maxv), bins_(bins), states_(states) {}

    /** In derived classes implement constructor taking inputstream, to be used in FactorsIO for discCont */

    virtual ~Mixture() {};

    /** <Returns P( x in [i; i+1[ ) for X|S=s */
    virtual number_t binProb(int i, int s) const = 0; 

    virtual void serialize( std::ostream & os) const = 0;

    virtual void mkFactor(matrix_t &m) const = 0;

    virtual int optimizeParameters(matrix_t & counts) = 0;

  protected:
    number_t minv_;
    number_t maxv_;
    unsigned bins_;
    unsigned states_;
  };

  class NormalMixture : public Mixture {
  public:
    /** Constructor */
    NormalMixture(vector_t const & means, vector_t const & vars, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Constructor */
    NormalMixture( std::istream & str, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Constructor with defaults? */
    NormalMixture(int states, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Dtor */
    virtual ~NormalMixture() {} ;

    virtual number_t binProb(int i, int s) const;
    
    virtual void serialize( std::ostream & os) const;

    virtual void mkFactor(matrix_t &m) const;

    virtual int optimizeParameters(matrix_t & counts) ;

  private:
    void setDists();
    vector_t means_;
    vector_t vars_;

    vector<boost::math::normal> dists_;
  };


  class GammaMixture : public Mixture {
  public:
    /** Constructor */
    GammaMixture(vector_t const & alphas, vector_t const & betas, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Constructor with inputstream */
    GammaMixture(std::istream & str, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Constructor with defaults? */
    GammaMixture(unsigned states, number_t const & minv, number_t const & maxv, unsigned const & bins);

    /** Dtor */
    virtual ~GammaMixture() {} ;

    virtual number_t binProb(int i, int s) const ;

    virtual void serialize( std::ostream & os) const;    

    virtual void mkFactor(matrix_t &m) const;
 
    virtual int optimizeParameters(matrix_t & counts) ;


  private:
    void setDists();
    vector_t alphas_;
    vector_t betas_;
    vector< boost::math::gamma_distribution<> > dists_;
  };
  
  //IO functions
  MixPtr_t readMixture(std::istream & str, string dist, number_t const & minv, number_t const & maxv, unsigned const & bins, unsigned const & states);
}

#endif // __Mixtures_h
