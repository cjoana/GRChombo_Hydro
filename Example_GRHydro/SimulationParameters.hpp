/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "CCZ4.hpp"
#include "EquationOfState.hpp"

// For BH at initio
// #include "KerrBH.hpp"


// #include "ScalarGauss.hpp" FIXME: needed?

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        Real W_def = 10000;
        // for regridding
        pp.load("regrid_threshold_K", regrid_threshold_K);
        pp.load("regrid_threshold_W", regrid_threshold_W, W_def);
        // pp.load("regrid_threshold_phi", regrid_threshold_phi); 


        // Cosmological parameters
        pp.load("G_Newton", G_Newton, 1.0);

        // Predefine ideal Fluids
        pp.load("omega", eos_params.omega, 0.0);
        pp.load("mass", eos_params.mass, 1.0);

        // Relaxation params  
        pp.load("relaxtime", relaxtime);
        pp.load("relaxspeed", relaxspeed);
        pp.load("omega", omega, 0.333333333);



#ifdef USE_AHFINDER
        double AH_guess = 0.5;
        pp.load("AH_1_initial_guess", AH_1_initial_guess, AH_guess);
        pp.load("AH_2_initial_guess", AH_2_initial_guess, AH_guess);
#endif

    }

    // Regrid parameters
    Real regrid_threshold_K;
    Real regrid_threshold_W;

    // Relaxation params
    Real relaxtime, relaxspeed, omega;


    // Initial data for matter and potential
    double G_Newton;
    EquationOfState::params_t eos_params;
    // KerrBH::params_t kerr_params;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    double AH_2_initial_guess;
#endif



};

#endif /* SIMULATIONPARAMETERS_HPP_ */
