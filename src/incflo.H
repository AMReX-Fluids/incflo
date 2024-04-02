#ifndef INCFLO_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Math.H>

#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <DiffusionTensorOp.H>
#include <DiffusionScalarOp.H>

class incflo : public amrex::AmrCore
{
public:

    friend DiffusionTensorOp;
    friend DiffusionScalarOp;

    enum struct FluidModel {
        Newtonian, powerlaw, Bingham, HerschelBulkley, deSouzaMendesDutra
    };

    incflo ();
    ~incflo () override;

    // Declare a default move constructor so we ensure the destructor is
    // not called when we return an object of this class by value
    incflo(incflo&&)  noexcept = default;

    // Declare a default move assignment operator
    incflo& operator=(incflo&& other)  noexcept = default;

    // Delete the copy constructor and copy assignment operators because
    // the integrator allocates internal memory that is best initialized
    // from scratch when needed instead of making a copy.

    // Delete the copy constructor
    incflo(const incflo& other) = delete;
    //
    // Delete the copy assignment operator
    incflo& operator=(const incflo& other) = delete;

    // Initialize multilevel AMR data
    void InitData ();

#ifdef AMREX_USE_EB
    void WriteMyEBSurface ();
#endif

    // Evolve solution to final time through repeated calls to Advance()
    void Evolve ();

    // Tag cells for refinement
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping
    // Only used during initialization
    void MakeNewLevelFromScratch (int lev, amrex::Real time,
                                  const amrex::BoxArray& new_grids,
                                  const amrex::DistributionMapping& new_dmap) override;

    // Make a new level using provided BoxArray and DistributionMapping,
    // and fill with interpolated coarse level data
    void MakeNewLevelFromCoarse (int lev, amrex::Real time,
                                 const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping,
    // and fill with existing fine and coarse data
    void RemakeLevel (int lev, amrex::Real time,
                      const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    // Delete level data
    void ClearLevel (int lev) override;

// for cuda

    void ComputeDt (int initialization, bool explicit_diffusion);

    amrex::Real vol_wgt_sum (amrex::Vector<amrex::MultiFab*> const& mf, int icomp);

    void compute_tra_forces (amrex::Vector<amrex::MultiFab*> const& tra_forces,
                             amrex::Vector<amrex::MultiFab const*> const& density);
    void compute_vel_forces (amrex::Vector<amrex::MultiFab*> const& vel_forces,
                             amrex::Vector<amrex::MultiFab const*> const& velocity,
                             amrex::Vector<amrex::MultiFab const*> const& density,
                             amrex::Vector<amrex::MultiFab const*> const& tracer_old,
                             amrex::Vector<amrex::MultiFab const*> const& tracer_new,
                             bool include_pressure_gradient = true);
    void compute_vel_forces_on_level ( int lev,
                                       amrex::MultiFab& vel_forces,
                                       const amrex::MultiFab& velocity,
                                       const amrex::MultiFab& density,
                                       const amrex::MultiFab& tracer_old,
                                       const amrex::MultiFab& tracer_new,
                                       bool include_pressure_gradient = true);


    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////

    std::unique_ptr<amrex::iMultiFab> make_BC_MF (int lev,
                                                  amrex::Gpu::DeviceVector<amrex::BCRec> const& bcs,
                                                  std::string const& field);
    amrex::iMultiFab make_nodalBC_mask (int lev);
    amrex::Vector<amrex::MultiFab> make_robinBC_MFs(int lev, amrex::MultiFab* state = nullptr);
    // void make_ccBC_mask (int lev, const amrex::BoxArray& ba,
    //                      const amrex::DistributionMapping& dm);
    // void make_nodalBC_mask (int lev, const amrex::BoxArray& ba,
    //                         const amrex::DistributionMapping& dm);

#ifdef AMREX_USE_EB
    void set_eb_velocity (int lev, amrex::Real time, amrex::MultiFab& eb_vel, int nghost);
    void set_eb_density (int lev, amrex::Real time, amrex::MultiFab& eb_density, int nghost);
    void set_eb_tracer (int lev, amrex::Real time, amrex::MultiFab& eb_tracer, int nghost);
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // convection
    //
    ///////////////////////////////////////////////////////////////////////////

    void ApplyPredictor(bool incremental_projection = false);
    void ApplyCorrector();
    void compute_convective_term (amrex::Vector<amrex::MultiFab*> const& conv_u,
                                  amrex::Vector<amrex::MultiFab*> const& conv_r,
                                  amrex::Vector<amrex::MultiFab*> const& conv_t,
                                  amrex::Vector<amrex::MultiFab const*> const& vel,
                                  amrex::Vector<amrex::MultiFab const*> const& density,
                                  amrex::Vector<amrex::MultiFab const*> const& tracer,
                                  AMREX_D_DECL(amrex::Vector<amrex::MultiFab*> const& u_mac,
                                               amrex::Vector<amrex::MultiFab*> const& v_mac,
                                               amrex::Vector<amrex::MultiFab*> const& w_mac),
                                  amrex::Vector<amrex::MultiFab      *> const& vel_forces,
                                  amrex::Vector<amrex::MultiFab      *> const& tra_forces,
                                  amrex::Real time);

    void compute_MAC_projected_velocities (
                                  amrex::Vector<amrex::MultiFab const*> const& vel,
                                  amrex::Vector<amrex::MultiFab const*> const& density,
                                  AMREX_D_DECL(amrex::Vector<amrex::MultiFab*> const& u_mac,
                                               amrex::Vector<amrex::MultiFab*> const& v_mac,
                                               amrex::Vector<amrex::MultiFab*> const& w_mac),
                                  amrex::Vector<amrex::MultiFab*> const& vel_forces,
                                  amrex::Real time);

    ///////////////////////////////////////////////////////////////////////////
    //
    // derive
    //
    ///////////////////////////////////////////////////////////////////////////

    void ComputeMagVel (int lev, amrex::Real time, amrex::MultiFab& magvel,
                        amrex::MultiFab const& vel);
    void ComputeVorticity (int lev, amrex::Real time, amrex::MultiFab& vort,
                           amrex::MultiFab const& vel);
    void ComputeDivU (amrex::Real time);
    void ComputeDrag ();
    [[nodiscard]] static amrex::Real ComputeKineticEnergy () ;

    void DiffFromExact (int lev, amrex::Geometry& lev_geom, amrex::Real time, amrex::Real dt,
                        amrex::MultiFab& error, int soln_comp, int err_comp) const;

    virtual void compute_strainrate_at_level (int lev,
                                              amrex::MultiFab* strainrate,
                                              amrex::MultiFab* vel,
                                              amrex::Geometry& lev_geom,
                                              amrex::Real time, int nghost);

    ///////////////////////////////////////////////////////////////////////////
    //
    // diffusion
    //
    ///////////////////////////////////////////////////////////////////////////

    [[nodiscard]] amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
    average_velocity_eta_to_faces (int lev, amrex::MultiFab const& cc_eta) const;

    void compute_divtau  (amrex::Vector<amrex::MultiFab      *> const& divtau,
                          amrex::Vector<amrex::MultiFab const*> const& velocity,
                          amrex::Vector<amrex::MultiFab const*> const& density,
                          amrex::Vector<amrex::MultiFab const*> const& eta);

    void diffuse_velocity (amrex::Vector<amrex::MultiFab      *> const&  velocity,
                           amrex::Vector<amrex::MultiFab      *> const&  density,
                           amrex::Vector<amrex::MultiFab const*> const& eta,
                           amrex::Real dt_diff);

    [[nodiscard]] amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
    average_scalar_eta_to_faces (int lev, int comp, amrex::MultiFab const& cc_eta) const;

    void compute_laps (amrex::Vector<amrex::MultiFab      *> const& laps,
                       amrex::Vector<amrex::MultiFab const*> const& scalar,
                       amrex::Vector<amrex::MultiFab const*> const& density,
                       amrex::Vector<amrex::MultiFab const*> const& eta);

    void diffuse_scalar (amrex::Vector<amrex::MultiFab      *> const& scalar,
                         amrex::Vector<amrex::MultiFab      *> const& density,
                         amrex::Vector<amrex::MultiFab const*> const& eta,
                         amrex::Real dt_diff);

    void fixup_eta_on_domain_faces (int lev,
                                    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fc,
                                    amrex::MultiFab const& cc) const;

    ///////////////////////////////////////////////////////////////////////////
    //
    // prob
    //
    ///////////////////////////////////////////////////////////////////////////

    void prob_init_fluid (int lev);
    void prob_set_BC_MF (amrex::Orientation const& ori, amrex::Box const& bx,
                         amrex::Array4<int> const& mask, int lev,
                         int inflow_val, int outflow_val, std::string const& field);
    void prob_set_MAC_robinBCs (amrex::Orientation const& ori, amrex::Box const& bx,
                                amrex::Array4<amrex::Real> const& robin_a,
                                amrex::Array4<amrex::Real> const& robin_b,
                                amrex::Array4<amrex::Real> const& robin_f,
                                int lev);
    void prob_set_diffusion_robinBCs (amrex::Orientation const& ori, amrex::Box const& bx,
                                      amrex::Array4<amrex::Real> const& robin_a,
                                      amrex::Array4<amrex::Real> const& robin_b,
                                      amrex::Array4<amrex::Real> const& robin_f,
                                      amrex::Array4<amrex::Real const> const& bcval,
                                      int lev);

#include "incflo_prob_I.H"
#include "incflo_prob_usr_I.H"

    ///////////////////////////////////////////////////////////////////////////
    //
    // projection
    //
    ///////////////////////////////////////////////////////////////////////////

    void ApplyProjection(amrex::Vector<amrex::MultiFab const*> density,
                         amrex::Real time, amrex::Real scaling_factor, bool incremental);
    void ApplyNodalProjection(amrex::Vector<amrex::MultiFab const*> density,
                              amrex::Real time, amrex::Real scaling_factor, bool incremental);

    ///////////////////////////////////////////////////////////////////////////
    //
    // rheology
    //
    ///////////////////////////////////////////////////////////////////////////

    void compute_viscosity (amrex::Vector<amrex::MultiFab*> const& eta,
                            amrex::Vector<amrex::MultiFab*> const& rho,
                            amrex::Vector<amrex::MultiFab*> const& vel,
                            amrex::Real time, int nghost);
    void compute_viscosity_at_level (int lev,
                                     amrex::MultiFab* eta,
                                     amrex::MultiFab* rho,
                                     amrex::MultiFab* vel,
                                     amrex::Geometry& lev_geom,
                                     amrex::Real time, int nghost);
    void compute_tracer_diff_coeff (amrex::Vector<amrex::MultiFab*> const& tra_eta, int nghost);

#ifdef AMREX_USE_EB
    ///////////////////////////////////////////////////////////////////////////
    //
    // embedded_boundaries
    //
    ///////////////////////////////////////////////////////////////////////////
    void incflo_correct_small_cells (amrex::Vector<amrex::MultiFab*      > const& vel_in,
                                     AMREX_D_DECL(amrex::Vector<amrex::MultiFab const*> const& u_mac,
                                                  amrex::Vector<amrex::MultiFab const*> const& v_mac,
                                                  amrex::Vector<amrex::MultiFab const*> const& w_mac));

    void redistribute_term ( amrex::MultiFab& result, amrex::MultiFab& result_tmp,
                 amrex::MultiFab const& state,
                 amrex::BCRec const* bc,
                 int lev);
    void redistribute_term ( amrex::MFIter const& mfi,
                 amrex::MultiFab& result, amrex::MultiFab& result_tmp,
                 amrex::MultiFab const& state,
                 amrex::BCRec const* bc,
                 int lev);
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // utilities
    //
    ///////////////////////////////////////////////////////////////////////////

    static bool SteadyStateReached ();


private:

    //
    // member variables
    //

    // Be verbose?
    int m_verbose = 0;

    // Member variables for initial conditions
    int m_probtype = 0;
    amrex::Real m_ic_u = amrex::Real(0.0);
    amrex::Real m_ic_v = amrex::Real(0.0);
    amrex::Real m_ic_w = amrex::Real(0.0);
    amrex::Real m_ic_p = amrex::Real(0.0);
    amrex::Real m_ic_t = amrex::Real(0.0);

    amrex::Vector<amrex::Real> m_t_old;
    amrex::Vector<amrex::Real> m_t_new;

    // Times
    amrex::Real m_cur_time     = amrex::Real( 0.0);
    amrex::Real m_dt           = amrex::Real(-1.0);
    amrex::Real m_prev_dt      = amrex::Real(-1.0);
    amrex::Real m_prev_prev_dt = amrex::Real(-1.0);

    // Time step counter
    int m_nstep         = -1;

    // Stop simulation if cur_time reaches stop_time OR nstep reaches max_step
    // OR steady_state = true AND steady_state_tol is reached
    amrex::Real m_stop_time = amrex::Real(-1.0);
    int m_max_step = -1;
    bool m_steady_state = false;
    amrex::Real m_steady_state_tol = amrex::Real(1.0e-5);

    // Options to control time stepping
    amrex::Real m_cfl = amrex::Real(0.5);
    amrex::Real m_fixed_dt = amrex::Real(-1.);
    amrex::Real m_init_shrink = amrex::Real(0.1);

    // Initial projection / iterations
    bool m_do_initial_proj    = true;
    int  m_initial_iterations = 3;

    // Use Boussinesq approximation for buoyancy?
    bool m_use_boussinesq = false;

    // Use tensor solve (as opposed to separate solves for each velocity component)
    bool use_tensor_solve      = true;

    // Use separate solves for each velocity component with lagged tensor terms
    bool use_tensor_correction = false;

    // AMR / refinement settings
    int m_refine_cutcells = 1;
    int m_regrid_int = -1;

    // ***************************************************************
    // MAC projection
    // ***************************************************************
    // Verbosity and MLMG parameters are now ParmParse with "mac_proj" in the
    // inputs file
    // Examples: mac_proj.verbose = 1
    //           mac_proj.bottom_verbose = 1
    //           mac_proj.maxiter
    //           mac_proj.bottom_maxiter
    //           mac_proj.bottom_rtol
    //           mac_proj.bottom_atol
    //           mac_proj.bottom_solver
    // More info in "AMReX-Hydro/Projections/hydro_MacProjector.cpp"

    std::unique_ptr<Hydro::MacProjector> macproj;

    int m_mac_mg_max_coarsening_level = 100;

#ifdef AMREX_USE_FLOAT
    amrex::Real m_mac_mg_rtol = amrex::Real(1.0e-4);
    amrex::Real m_mac_mg_atol = amrex::Real(1.0e-7);
#else
    amrex::Real m_mac_mg_rtol = amrex::Real(1.0e-11);
    amrex::Real m_mac_mg_atol = amrex::Real(1.0e-14);
#endif

    // ***************************************************************
    // Nodal solve
    // ***************************************************************

    // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
    // inputs file
    // Examples: nodal_proj.verbose = 1
    //           nodal_proj.bottom_verbose = 1
    //           nodal_proj.maxiter
    //           nodal_proj.bottom_maxiter
    //           nodal_proj.bottom_rtol
    //           nodal_proj.bottom_atol
    //           nodal_proj.bottom_solver
    // More info in "AMReX-Hydro/Projections/hydro_NodalProjector.cpp"
    // Verbosity

    // Control MLMG behavior
#ifdef AMREX_USE_FLOAT
    amrex::Real m_nodal_mg_rtol = amrex::Real(1.0e-4);
    amrex::Real m_nodal_mg_atol = amrex::Real(1.0e-7);
#else
    amrex::Real m_nodal_mg_rtol = amrex::Real(1.0e-11);
    amrex::Real m_nodal_mg_atol = amrex::Real(1.0e-14);
#endif

    // Max coarsening level
    int m_nodal_mg_max_coarsening_level = 100;

    // ***************************************************************
    // ***************************************************************

    amrex::Real m_p000 = 0.0; // background pressure at node (0,0,0).
    amrex::Vector<amrex::Real> m_gp0{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_delp{{0.0, 0.0, 0.0}};

    // Scalar diffusive coefficient
    amrex::Vector<amrex::Real> m_mu_s;

    // Density (if constant)
    amrex::Real m_ro_0 = 1.0;

    bool m_constant_density = true;
    bool m_advect_tracer    = false;
    bool m_test_tracer_conservation = false;

    // If false then we update velocity using the convective   form, u dot grad u
    // If  true then we update velocity using the conservative form, del dot (u u)
    bool m_advect_momentum = false;

    std::string m_advection_type = "Godunov";

#ifdef AMREX_USE_EB
    std::string m_redistribution_type = "StateRedist";

    // If using Godunov with EB, default to PLM
    bool m_godunov_ppm         = false;

    // EB flow settings
    struct EBFlow_t {
       bool enabled{false};
       bool has_normal{false};
       bool is_mag{false};
       amrex::Real vel_mag{0.};
       amrex::Vector<amrex::Real> velocity;

       // Normal to exclude EB faces. By default all faces within a specified
       // region will
       amrex::Vector<amrex::Real> normal;
       amrex::Real normal_tol{0.};

       // scalars
       amrex::Real density{1.};
       amrex::Vector<amrex::Real> tracer;

       EBFlow_t ()

       = default;
    };
    EBFlow_t m_eb_flow;

    bool m_write_geom_chk = false;
#else
    // If using Godunov with no EB, default to PPM
    bool m_godunov_ppm         = true;
#endif

    // If using Godunov, default to including viscous/diffusive
    //    terms in the forcing terms
    bool m_godunov_include_diff_in_forcing = true;

    // If using Godunov, when to add the forcing terms in
    //    the construction of the "trans" velocities
    bool m_godunov_use_forces_in_trans = false;

    // If using Godunov, if true then don't include the pressure gradient
    //    in the forcing term passed into the godunov routine; instead
    //    use gradient of mac phi which contains the full pressure
    bool m_use_mac_phi_in_godunov = false;

    enum struct DiffusionType {
        Invalid, Explicit, Crank_Nicolson, Implicit
    };
    DiffusionType m_diff_type = DiffusionType::Implicit;

    // Fluid properties
    FluidModel m_fluid_model;
    amrex::Real m_mu = 1.0;
    amrex::Real m_n_0 = 0.0;
    amrex::Real m_tau_0 = 0.0;
    amrex::Real m_papa_reg = 0.0;
    amrex::Real m_eta_0 = 0.0;

    int m_plot_int = -1;

    // Dump plotfiles at as close as possible to the designated period *without* changing dt
    amrex::Real m_plot_per_approx = -1.0;

    // Dump plotfiles at exactcly the designated period by changing dt
    amrex::Real m_plot_per_exact  = -1.0;

    int m_last_plt = -1;
    std::string m_plot_file{"plt"};

    int m_check_int = -1;
    int m_last_chk = -1;
    int m_KE_int = -1;
    std::string m_check_file{"chk"};
    std::string m_restart_file;
    std::string m_tag_file;

    bool m_plotfile_on_restart = false;

    amrex::Vector<amrex::Real> tag_region_lo;
    amrex::Vector<amrex::Real> tag_region_hi;

    // Flags for saving fluid data in plot files
    int m_plt_velx        = 1;
    int m_plt_vely        = 1;
    int m_plt_velz        = 1;
    int m_plt_gpx         = 1;
    int m_plt_gpy         = 1;
    int m_plt_gpz         = 1;
    int m_plt_rho         = 1;
    int m_plt_tracer      = 1;
    int m_plt_p_nd        = 0;
    int m_plt_macphi      = 0;
    int m_plt_eta         = 0;
    int m_plt_vort        = 1;
    int m_plt_magvel      = 1;
    int m_plt_forcing     = 0;
    int m_plt_strainrate  = 0;
    int m_plt_divu        = 0;
    int m_plt_vfrac       = 1;
    int m_plt_error_u     = 0;
    int m_plt_error_v     = 0;
    int m_plt_error_w     = 0;
    int m_plt_error_p     = 0;
    int m_plt_error_mac_p = 0;

    struct LevelData {
        LevelData () = default;
        LevelData (amrex::BoxArray const& ba,
                   amrex::DistributionMapping const& dm,
                   amrex::FabFactory<amrex::FArrayBox> const& fact,
                   int ntrac, int ng_state,
                   const std::string& advection_type, bool implicit_diffusion,
                   bool use_tensor_correction, bool advect_tracer);
        // cell-centered multifabs
        amrex::MultiFab velocity;
        amrex::MultiFab velocity_o;
        amrex::MultiFab velocity_eb;
        amrex::MultiFab density;
        amrex::MultiFab density_eb;
        amrex::MultiFab density_o;
        amrex::MultiFab tracer;
        amrex::MultiFab tracer_eb;
        amrex::MultiFab tracer_o;

        amrex::MultiFab mac_phi; // cell-centered pressure used in MAC projection

        // nodal pressure multifab
        amrex::MultiFab p_nd;

        // cell-centered pressure gradient
        amrex::MultiFab gp;

        amrex::MultiFab conv_velocity;
        amrex::MultiFab conv_velocity_o;
        amrex::MultiFab conv_density;
        amrex::MultiFab conv_density_o;
        amrex::MultiFab conv_tracer;
        amrex::MultiFab conv_tracer_o;

        amrex::MultiFab divtau;
        amrex::MultiFab divtau_o;
        amrex::MultiFab laps;
        amrex::MultiFab laps_o;
    };

    amrex::Vector<std::unique_ptr<LevelData> > m_leveldata;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

    enum struct BC {
        pressure_inflow, pressure_outflow, mass_inflow, no_slip_wall, slip_wall,
        periodic, mixed, undefined
    };

    bool m_has_mixedBC = false;

    amrex::GpuArray<BC                         , AMREX_SPACEDIM*2> m_bc_type;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_pressure;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_density;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                               , AMREX_SPACEDIM*2> m_bc_velocity;
    amrex::GpuArray<amrex::Vector<amrex::Real> , AMREX_SPACEDIM*2> m_bc_tracer;
    amrex::Vector<amrex::Real>                                     m_bc_eb_velocity;

    // amrex::Vector cannot be used on gpu, so ...
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> m_bc_tracer_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_tracer_raii;
    //
    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Vector<amrex::BCRec> m_bcrec_density;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_density_d;
    amrex::Vector<amrex::BCRec> m_bcrec_tracer;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_tracer_d;
    amrex::Vector<amrex::BCRec> m_bcrec_force;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_force_d;

    amrex::Vector<int> m_iconserv_velocity;
    amrex::Gpu::DeviceVector<int> m_iconserv_velocity_d;
    amrex::Vector<int> m_iconserv_density;
    amrex::Gpu::DeviceVector<int> m_iconserv_density_d;
    amrex::Vector<int> m_iconserv_tracer;
    amrex::Gpu::DeviceVector<int> m_iconserv_tracer_d;

    int m_ntrac = 1;

    std::unique_ptr<DiffusionTensorOp> m_diffusion_tensor_op;
    std::unique_ptr<DiffusionScalarOp> m_diffusion_scalar_op;

    //
    // end of member variables
    //

    [[nodiscard]] amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }
#ifdef AMREX_USE_EB
    [[nodiscard]] amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }
#endif

    // Number of ghost cells for field arrays.
    [[nodiscard]] int nghost_state () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular())
        {
            return 4;
        }
#endif
        {
            return (m_advection_type != "MOL") ? 3 : 2;
        }
    }

    // For Godunov, we need 1 ghost cell in addition to the Box we are filling
    // For MOL    , we need 0 ghost cells
    [[nodiscard]] int nghost_force () const
    {
       if (m_advection_type == "MOL")
           return 0;
       else
           return 1;
    }

    [[nodiscard]] int nghost_mac () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular()) return (m_advection_type == "MOL") ? 3 : 4;
#endif
        return (m_advection_type == "MOL") ? 0 : 1;
    }

#ifdef AMREX_USE_EB
    [[nodiscard]] static int nghost_eb_basic ()
    { return 5; }

    // We need 5 for doing StateRedistribution; otherwise 4 would be enough
    [[nodiscard]] static int nghost_eb_volume ()
    { return 5; }

    [[nodiscard]] static int nghost_eb_full ()
    { return 4; }
#endif

    [[nodiscard]] bool need_divtau () const {
        return ( m_godunov_include_diff_in_forcing || DiffusionType::Implicit != m_diff_type );
    }

    [[nodiscard]] bool AdvectMomentum () const {
        return m_advect_momentum;
    }

    [[nodiscard]] bool useTensorSolve () const {
        return use_tensor_solve;
    }

    [[nodiscard]] bool useTensorCorrection () const {
        return use_tensor_correction;
    }

#ifdef AMREX_USE_EB
    [[nodiscard]] bool hasEBFlow () const {
       return m_eb_flow.enabled;
    }
#endif

    DiffusionTensorOp* get_diffusion_tensor_op ();
    DiffusionScalarOp* get_diffusion_scalar_op ();

    amrex::Vector<amrex::MultiFab*> get_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_velocity_eb () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_eb () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_eb () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_mac_phi () noexcept;
    amrex::Vector<amrex::MultiFab*> get_vel_forces () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tra_forces () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_divtau_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_divtau_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_laps_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_laps_new () noexcept;
    //
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_velocity_old_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_velocity_new_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_velocity_eb () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_density_eb () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_tracer_eb () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_density_old_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_density_new_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_tracer_old_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_tracer_new_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_vel_forces_const () const noexcept;
    [[nodiscard]] amrex::Vector<amrex::MultiFab const*> get_tra_forces_const () const noexcept;

    [[nodiscard]] amrex::Vector<int> const& get_velocity_iconserv () const noexcept { return m_iconserv_velocity; }
    [[nodiscard]] amrex::Vector<int> const& get_density_iconserv () const noexcept { return m_iconserv_density; }
    [[nodiscard]] amrex::Vector<int> const& get_tracer_iconserv () const noexcept { return m_iconserv_tracer; }
    //
    [[nodiscard]] int const* get_velocity_iconserv_device_ptr () const noexcept {
        return m_iconserv_velocity_d.data(); }
    [[nodiscard]] int const* get_density_iconserv_device_ptr () const noexcept {
        return m_iconserv_density_d.data(); }
    [[nodiscard]] int const* get_tracer_iconserv_device_ptr () const noexcept {
        return m_iconserv_tracer_d.data(); }

    [[nodiscard]] amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
    [[nodiscard]] amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    [[nodiscard]] amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }
    [[nodiscard]] amrex::Vector<amrex::BCRec> const& get_force_bcrec () const noexcept { return m_bcrec_force; }
    //
    [[nodiscard]] amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept {
        return m_bcrec_velocity_d.data(); }
    [[nodiscard]] amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept {
        return m_bcrec_density_d.data(); }
    [[nodiscard]] amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept {
        return m_bcrec_tracer_d.data(); }
    [[nodiscard]] amrex::BCRec const* get_force_bcrec_device_ptr () const noexcept {
        return m_bcrec_force_d.data(); }

    [[nodiscard]] amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_nodal_projection_bc (amrex::Orientation::Side side) const noexcept;
    [[nodiscard]] amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_mac_projection_bc (amrex::Orientation::Side side) const noexcept;

    [[nodiscard]] amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_diffuse_velocity_bc (amrex::Orientation::Side side, int comp) const noexcept;

    [[nodiscard]] amrex::Vector<amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> >
    get_diffuse_tensor_bc (amrex::Orientation::Side side) const noexcept;

    [[nodiscard]] amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_diffuse_scalar_bc (amrex::Orientation::Side side) const noexcept;

    void fillpatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillpatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillpatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);
    void fillpatch_gradp (int lev, amrex::Real time, amrex::MultiFab& gp, int ng);
    void fillpatch_force (amrex::Real time, amrex::Vector<amrex::MultiFab*> const& force, int ng);

    void fillcoarsepatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillcoarsepatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillcoarsepatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);
    void fillcoarsepatch_gradp (int lev, amrex::Real time, amrex::MultiFab& gp, int ng);

    void fillphysbc_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillphysbc_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillphysbc_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);

    void copy_from_new_to_old_velocity (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_velocity (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_density  (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_density  (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_tracer   (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_tracer   (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    //
    void copy_from_old_to_new_velocity (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_velocity (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_density  (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_density  (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_tracer   (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_tracer   (int lev, amrex::IntVect const& ng = amrex::IntVect{0});

    void Advance ();
    bool writeNow ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////

    void init_bcs ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // convection
    //
    ///////////////////////////////////////////////////////////////////////////

    void init_advection ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // embedded_boundaries
    //
    ///////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_EB
    void MakeEBGeometry ();
    void make_eb_annulus ();
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_twocylinders ();
    void make_eb_regular ();
    void make_eb_sphere ();
    void make_eb_spherecube ();
    void make_eb_cyl_tuscan ();
    void make_eb_tuscan ();
    void make_eb_chkptfile ();
#ifdef CSG_EB
    void make_eb_csg (const std::string& csg_file);
#endif
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // rheology
    //
    ///////////////////////////////////////////////////////////////////////////

    void ReadRheologyParameters ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // setup
    //
    ///////////////////////////////////////////////////////////////////////////

    void set_background_pressure ();
    void ReadParameters ();
    void ReadIOParameters ();
    void ResizeArrays (); // Resize arrays to fit (up to) max_level + 1 AMR levels
    void InitialProjection ();
    void InitialIterations ();
#ifdef AMREX_USE_EB
    void InitialRedistribution ();
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // utilities
    //
    ///////////////////////////////////////////////////////////////////////////

    void WriteHeader (const std::string& name, bool is_checkpoint) const;
    void WriteJobInfo (const std::string& path) const;
    void WriteCheckPointFile () const;
    void WritePlotFile ();
    void ReadCheckpointFile ();
};

#endif
