/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EOS_HPP_
#define EOS_HPP_

#include "simd.hpp"

class EquationOfState
{
  public:
    struct params_t
    {
      double omega;
      double mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    EquationOfState(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &pressure, data_t &omega,
		    	 data_t &dpdrho, data_t &dpdenergy,
                           const vars_t<data_t> &vars) const
   {
               /*
	       pressure = m_params.omega * vars.density * ( 1  + vars.energy);
	       enthalpy = 1 + vars.energy + pressure / vars.density;
	       dpdrho = m_params.omega * ( 1  + vars.energy);
	       dpdenergy = m_params.omega * vars.density ;
               */
               omega = m_params.omega;

   }
};

#endif /* EquationOfState_HPP_ */
