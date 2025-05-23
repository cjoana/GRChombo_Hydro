/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 * 
 * Evoltion of cosmological perfect fluid, e.g. PBH formation studies. 
 * Cristian Joana
 */

// General includes common to most GR problems
#include "PerfectFluidLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "AMRReductions.hpp"
#include "MyGauge.hpp"
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
//#include "MatterConstraints.hpp"
// #include "ExtendedMatterConstraints.hpp"
#include "MyNewMatterConstraints.hpp"

// Problem specific includes
#include "FluidInitData.hpp"
// #include "ChiRelaxation.hpp"   // $GRCHOMBO/Source/Matter/
//#include "KTaggingCriterion.hpp"
#include "KandWTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"
#include "GRHydro.hpp"     // in dir ./GRHlib

#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"


// Things to do at each advance step, after the RK4 is calculated
void PerfectFluidLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(),
			    PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());


}

// Initial data for field and metric variables
void PerfectFluidLevel::initialData()
{
    CH_TIME("PerfectFluidLevel::initialData");
    if (m_verbosity)
        pout() << "PerfectFluidLevel::initialData " << m_level << endl;

}


// Things to do before outputting a checkpoint file
void PerfectFluidLevel::prePlotLevel()
{
    fillAllGhosts();

    EquationOfState eos(m_p.eos_params);                                          
    PerfectFluidWithEOS perfect_fluid(eos);
    BoxLoops::loop(MatterConstraints<PerfectFluidWithEOS>(
                  // perfect_fluid, m_dx, m_p.G_Newton),
                  // m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
              perfect_fluid, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom),
              c_Ham_abs_terms,  Interval(c_Mom_abs_terms, c_Mom_abs_terms),
              c_ricci_scalar, c_trA2,
              c_rho, c_S, c_ricci_scalar_tilde),
          m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    fillAllGhosts();

    BoxLoops::loop(
                      MatterWeyl4<PerfectFluidWithEOS>(perfect_fluid, m_p.extraction_params.center, m_dx),
                      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);




    AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
    // AMRReductions<VariableType::diagnostic> camr_reductions(m_gr_amr);
    double Ham = amr_reductions.norm(c_Ham_abs_terms, 2, true);
    double rho = amr_reductions.norm(c_rho, 2, true);

    pout() << "Ham_Abs_tems " << Ham <<  "  rho  " <<   rho   << endl;
}




// Things to do in RHS update, at each RK4 step
void PerfectFluidLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Relaxation function for chi - this will eventually be done separately
    // with hdf5 as input
    if (m_time <= m_p.relaxtime)
    {
        SetValue set_other_values_zero(0.0, Interval(c_chi, NUM_VARS - 1));
        BoxLoops::loop(
           make_compute_pack(TraceARemoval(), set_other_values_zero),
           a_soln,  a_rhs, INCLUDE_GHOST_CELLS);


    }
    else
    {
        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(),
		    PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)), a_soln,
            a_soln, INCLUDE_GHOST_CELLS);


        //AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        AMRReductions<VariableType::evolution> amr_reductions(m_gr_amr);
        double abs_K_mean = amr_reductions.norm(c_K, 2, true);
        // auto Gauge = MovingPunctureGauge(m_p.ccz4_params);
        m_p.ccz4_params.K_mean = - abs_K_mean;

        // Calculate MatterCCZ4 right hand side with matter_t = PerfectFluid
        // We don't want undefined values floating around in the constraints so
        // zero these
	
        const int ratio = pow(2, 5 * (m_level - m_p.max_level));
        const double sigma = m_p.sigma * double(ratio);
	
        EquationOfState eos(m_p.eos_params);
        PerfectFluidWithEOS perfect_fluid(eos);
        MatterCCZ4RHS<PerfectFluidWithEOS> my_ccz4_matter(
            perfect_fluid, m_p.ccz4_params, m_dx, sigma, m_p.formulation,
            m_p.G_Newton);


        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);


    }
}

// Things to do at ODE update, after soln + rhs
void PerfectFluidLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);


    // Update constraints also for plotfile
    EquationOfState eos(m_p.eos_params);
    PerfectFluidWithEOS perfect_fluid(eos);


    if (m_time <= m_p.relaxtime)
    {

      //EquationOfState eos(m_p.eos_params);
      //PerfectFluidWithEOS perfect_fluid(eos);
      FluidInitData<PerfectFluidWithEOS> initialisation(
          perfect_fluid, m_dx, m_p.relaxspeed, m_p.G_Newton);

      BoxLoops::loop(initialisation, a_soln, a_soln, EXCLUDE_GHOST_CELLS);
    }


    // Evaluate fluid vars (density, energy, etc)  CJ !!!
    BoxLoops::loop(perfect_fluid, a_soln, a_soln,
                     EXCLUDE_GHOST_CELLS); 


    fillAllGhosts();

}

void PerfectFluidLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(
        KandWTaggingCriterion(m_dx, m_p.regrid_threshold_K, m_p.regrid_threshold_W),
        current_state, tagging_criterion);
}

void PerfectFluidLevel::specificPostTimeStep()
{
#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
}


