/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(GRHYDRO_HPP_)
#error "This file should only be included through GRHydro.hpp"
#endif

#ifndef GRHYDRO_IMPL_HPP_
#define GRHYDRO_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class eos_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> PerfectFluid<eos_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;
    Tensor<1, data_t> V_i;   // 4-velocity with lower indices

    data_t my_D =  simd_max(vars.D,  -vars.E + 1e-20);
    data_t  fluidT = my_D + vars.E + vars.pressure;

    FOR1(i)
    {
      V_i[i] = vars.Z[i] / fluidT;
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
          fluidT  * V_i[i] * V_i[j] +
          vars.pressure * vars.h[i][j]/vars.chi;
    }

    // S_i (note lower index) = - n^a T_ai
    FOR1(i) { out.Si[i] =  vars.Z[i]; }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // rho = n^a n^b T_ab
    out.rho =  my_D + vars.E;

    return out;
}

// Adds in the RHS for the matter vars
template <class eos_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void PerfectFluid<eos_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    data_t V_dot_dchi = 0;
    FOR1(i){
      V_dot_dchi += vars.V[i] * d1.chi[i];
      total_rhs.Z[i] =0;
     }

    data_t Z_dot_dchi = 0;
    FOR2(i,j){ Z_dot_dchi += vars.Z[i] * d1.chi[j] * h_UU[i][j]; }


    data_t my_D =  simd_max(vars.D,  -vars.E +1e-20);

    Tensor<2, data_t> covdtildeZ;
    Tensor<2, data_t> covdZ;     // D_k Z_l
    FOR2(k, l)
    {
        covdtildeZ[k][l] = d1.Z[k][l];
        FOR1(m) { covdtildeZ[k][l] -= chris.ULL[m][k][l] * vars.Z[m]; }
        covdZ[k][l] =
            covdtildeZ[k][l] +
            (0.5/vars.chi) * (vars.Z[k] * d1.chi[l] + d1.chi[k] * vars.Z[l] -
                   vars.h[k][l] * Z_dot_dchi);
    }

    data_t covdV = 0;    // D_m V^m
    FOR1(m) { covdV += - 3/(2*vars.chi)* d1.chi[m]*vars.V[m]  +  d1.V[m][m]; }


	  total_rhs.D = 0;
    total_rhs.E = 0;

    Tensor<2, data_t> K_tensor;
    FOR2(i, j)
    {
      K_tensor[i][j] = (vars.A[i][j] + vars.h[i][j] * vars.K / 3.) / vars.chi;
    }

    total_rhs.D += advec.D + vars.lapse * vars.K * my_D;
    total_rhs.E += advec.E + vars.lapse * vars.K *
                            (vars.pressure + vars.E);


    FOR1(i)
    {
        total_rhs.D += - vars.lapse * (d1.D[i] * vars.V[i]
                                        + my_D * covdV/GR_SPACEDIM)
                       - d1.lapse[i] * my_D * vars.V[i];


        total_rhs.E += - vars.lapse * (d1.E[i] * vars.V[i]
                                    + vars.E * covdV/GR_SPACEDIM)
                       - d1.lapse[i] * vars.E * vars.V[i]
                       - vars.lapse * (d1.pressure[i] * vars.V[i]
                                    + vars.pressure * covdV/GR_SPACEDIM)
                       - d1.lapse[i] * vars.pressure * vars.V[i]
                       - (my_D + vars.E + vars.pressure) *
                                    vars.V[i] * d1.lapse[i];



        total_rhs.Z[i] += advec.Z[i]
                       - vars.lapse * d1.pressure[i]
                       - d1.lapse[i] * vars.pressure
                       - (vars.E + my_D) * d1.lapse[i]
                       + vars.lapse * vars.K * vars.Z[i];


    }

    FOR2(i, j)
    {
        total_rhs.Z[i] += - vars.lapse * ( covdV/GR_SPACEDIM * vars.Z[i] +
                                           covdZ[i][j] * vars.V[j])
                          - d1.lapse[j] * vars.V[j] * vars.Z[i];


        total_rhs.E +=  (my_D + vars.E + vars.pressure) *
                            vars.lapse * vars.V[i] * vars.V[j] * K_tensor[i][j];

    }
}



template <class eos_t>
template <class data_t>
void PerfectFluid<eos_t>::compute(
  Cell<data_t> current_cell) const
{

    const auto vars = current_cell.template load_vars<Vars>();
    const auto geo_vars = current_cell.template load_vars<GeoVars>();
    auto up_vars = current_cell.template load_vars<Vars>();

    Tensor<1, data_t> V_i; // with lower indices: V_i
    up_vars.D =  simd_max(vars.D,  -vars.E + 1e-20);  
    data_t V2_max = 1 - 1e-15;

    // Inverse metric
    const auto h_UU = TensorAlgebra::compute_inverse_sym(geo_vars.h);
    data_t pressure = 0.0;
    data_t Lorentz, fl_dens;

    // Calculate S_iS^i
    data_t S2 = 0.0;
    FOR2(i, j)
    {
      S2 += vars.Z[i] * vars.Z[j] * h_UU[i][j] * geo_vars.chi;
    }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    recover_primvars_bartropic(current_cell, fl_dens, pressure, Lorentz, S2);
    // Redefine variables
    up_vars.W = 1./Lorentz;
    data_t density = up_vars.D / up_vars.W;
    data_t energy = fl_dens/density - 1;
    up_vars.pressure = pressure;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    FOR1(i) {
       up_vars.V[i] = 0;
     }
    FOR2(i,j) {
      up_vars.V[i] += vars.Z[j] * h_UU[i][j] * geo_vars.chi / (vars.E + up_vars.D + pressure);
    }

    // Overwrite new values for fluid variables
    current_cell.store_vars(up_vars.D, c_D);
    current_cell.store_vars(up_vars.pressure, c_pressure);
    current_cell.store_vars(up_vars.V, GRInterval<c_V1, c_V3>());
    current_cell.store_vars(up_vars.W, c_W);
}




template <class eos_t>
template <class data_t>
void PerfectFluid<eos_t>::recover_primvars_bartropic(Cell<data_t> current_cell,
                           data_t &fl_dens,
                           data_t &pressure,
                           data_t &Lorentz,
                           const data_t &S2) const
{
  const auto vars = current_cell.template load_vars<Vars>();
  auto up_vars = current_cell.template load_vars<Vars>();
  data_t omega = 0;
  data_t  enthalpy, dpdrho, dpdenergy ;
  enthalpy = dpdrho = dpdenergy = 0.0;
  data_t in_sqrt = 0;

  my_eos.compute_eos(pressure, omega, dpdrho, dpdenergy, up_vars);

  up_vars.D =   simd_max(vars.D,  -vars.E +1e-20); 

  /*  This part of the code is for omega < 0, needs testing. 
   * 

  if (omega == -1){
    fl_dens = vars.E + up_vars.D;
  }
  if (omega == 0){
    fl_dens = vars.E + up_vars.D - S2/(vars.E + up_vars.D);
  }

  if (omega > 0 && omega <= 1 ){
      in_sqrt = (omega-1)*(omega-1)*(vars.E + up_vars.D)*(vars.E + up_vars.D)
              - 4*omega*S2 + 4*omega*(vars.E + up_vars.D)*(vars.E + up_vars.D);
      in_sqrt = (in_sqrt > 0) ? in_sqrt : 0;
      fl_dens = ((omega -1)*(vars.E + up_vars.D) + sqrt(in_sqrt))/(2*omega);
  }

  if (omega > 1){
      data_t eplus, eminus;
      in_sqrt = (omega-1)*(omega-1)*(vars.E + up_vars.D)*(vars.E + up_vars.D)
              - 4*omega*S2 + 4*omega*(vars.E + up_vars.D)*(vars.E + up_vars.D);
      in_sqrt = (in_sqrt > 0) ? in_sqrt : 0;

      eplus = ((omega -1)*(vars.E + up_vars.D) +  sqrt(in_sqrt))/(2*omega);
      eminus = ((omega -1)*(vars.E + up_vars.D) - sqrt(in_sqrt))/(2*omega);

      fl_dens = (eplus > eminus) ? eplus : eminus;
  }
  */

  //Only   EoS =  (0, 1]
  in_sqrt = (omega-1)*(omega-1)*(vars.E + up_vars.D)*(vars.E + up_vars.D)
          - 4*omega*S2 + 4*omega*(vars.E + up_vars.D)*(vars.E + up_vars.D);
  in_sqrt =  simd_max(in_sqrt , 1e-50);
  fl_dens = ((omega -1)*(vars.E + up_vars.D) + sqrt(in_sqrt))/(2*omega);


  fl_dens = simd_min(fl_dens , vars.D + vars.E);
  pressure = fl_dens*omega;
  Lorentz = sqrt( (fl_dens + pressure)/(vars.E + up_vars.D + pressure));

}


#endif /* GRHYDRO_IMPL_HPP_ */
