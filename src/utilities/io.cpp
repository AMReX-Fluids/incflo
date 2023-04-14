#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include <incflo.H>

#if (AMREX_SPACEDIM)== 2
#include <AMReX_EB_LeastSquares_2D_K.H>
#else
#include <AMReX_EB_LeastSquares_3D_K.H>
#endif

using namespace amrex;

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void incflo::WriteHeader(const std::string& name, bool is_checkpoint) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        std::string HeaderFileName(name + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(),
                        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if(!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        if(is_checkpoint) {
            HeaderFile << "Checkpoint version: 1\n";
        } else {
            HeaderFile << "HyperCLaw-V1.1\n";
        }

        HeaderFile << finest_level << "\n";

        // Time stepping controls
        HeaderFile << m_nstep << "\n";
        HeaderFile << m_cur_time << "\n";
        HeaderFile << m_dt << "\n";
        HeaderFile << m_prev_dt << "\n";
        HeaderFile << m_prev_prev_dt << "\n";

        // Geometry
        for(int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geom(0).ProbLo(i) << ' ';
        }
        HeaderFile << '\n';

        for(int i = 0; i < BL_SPACEDIM; ++i)
            HeaderFile << Geom(0).ProbHi(i) << ' ';
        HeaderFile << '\n';

        // BoxArray
        for(int lev = 0; lev <= finest_level; ++lev)
        {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }
}

void incflo::WriteCheckPointFile() const
{
    BL_PROFILE("incflo::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate(m_check_file, m_nstep);

    amrex::Print() << "\n\t Writing checkpoint " << checkpointname << std::endl;

    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

    bool is_checkpoint = true;
    WriteHeader(checkpointname, is_checkpoint);
    WriteJobInfo(checkpointname);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Write(m_leveldata[lev]->velocity,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "velocity"));

        VisMF::Write(m_leveldata[lev]->density,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "density"));

        if (m_ntrac > 0) {
            VisMF::Write(m_leveldata[lev]->tracer,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "tracer"));
        }

        VisMF::Write(m_leveldata[lev]->gp,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gradp"));

        VisMF::Write(m_leveldata[lev]->p_nd,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p_nd"));
    }
}

void incflo::ReadCheckpointFile()
{
    BL_PROFILE("incflo::ReadCheckpointFile()");

    amrex::Print() << "Restarting from checkpoint " << m_restart_file << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(m_restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file

    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    GotoNextLine(is);

    // Step count
    is >> m_nstep;
    GotoNextLine(is);

    // Current time
    is >> m_cur_time;
    GotoNextLine(is);

    // Time step size
    is >> m_dt;
    GotoNextLine(is);

    is >> m_prev_dt;
    GotoNextLine(is);

    is >> m_prev_prev_dt;
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up problem domain
    RealBox rb(prob_lo, prob_hi);
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(lev, Geometry(Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                                  Geom(lev).isPeriodic()));
    }

#ifdef AMREX_USE_MOVING_EB
    if ( m_verbose > 1 ) { Print()<<"Making new geometry at time "<<m_cur_time<<std::endl; }
    // Create new EB geometry at the current time.
    // MakeNewLevelFromScratch() will build EBFabFactories
    MakeNewEBGeometry(m_cur_time);
#endif

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // Create distribution mapping
        DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, m_cur_time, ba, dm);
    }

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

    // Load the field data
    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Read(m_leveldata[lev]->velocity,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "velocity"));

        VisMF::Read(m_leveldata[lev]->density,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "density"));

        if (m_ntrac > 0) {
            VisMF::Read(m_leveldata[lev]->tracer,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "tracer"));
        }

        VisMF::Read(m_leveldata[lev]->gp,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "gradp"));

        VisMF::Read(m_leveldata[lev]->p_nd,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "p_nd"));
    }

    amrex::Print() << "Restart complete" << std::endl;
}

void incflo::WriteJobInfo(const std::string& path) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = path;
        std::string PrettyLine =
            "===============================================================================\n";

        FullPathJobInfoFile += "/incflo_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " incflo Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

        jobInfoFile << "\n\n";

        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
        jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
        jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        if(std::strlen(githash1) > 0)
        {
            jobInfoFile << "incflo git hash: " << githash1 << "\n";
        }
        if(std::strlen(githash2) > 0)
        {
            jobInfoFile << "AMReX git hash: " << githash2 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for(int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for(int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                jobInfoFile << geom[i].Domain().length(dir) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}

void incflo::WritePlotFile()
{
    BL_PROFILE("incflo::WritePlotFile()");

    if (m_plt_vort || m_plt_divu || m_plt_forcing || m_plt_eta || m_plt_strainrate) {
        for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_EB
            const int ng = (EBFactory(0).isAllRegular()) ? 1 : 2;
#else
            const int ng = 1;
#endif
            fillpatch_velocity(lev, m_cur_time, m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_cur_time, m_leveldata[lev]->density, ng);
            fillpatch_tracer(lev, m_cur_time, m_leveldata[lev]->tracer, ng);
        }
    }

    const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename << " at time " << m_cur_time << std::endl;

    int ncomp = 0;

    // Velocity components
    if (m_plt_velx) ++ncomp;
    if (m_plt_vely) ++ncomp;
#if (AMREX_SPACEDIM == 3)
    if (m_plt_velz) ++ncomp;
#endif

    // Pressure gradient components
    if (m_plt_gpx) ++ncomp;
    if (m_plt_gpy) ++ncomp;
#if (AMREX_SPACEDIM == 3)
    if (m_plt_gpz) ++ncomp;
#endif

    // Density
    if (m_plt_rho) ++ncomp;

    // Tracers
    if (m_plt_tracer) ncomp += m_ntrac;

    // Pressure
    if(m_plt_p_nd) ++ncomp;

    // MAC phi
    if(m_plt_macphi) ncomp += 1;

    // Error in u (computed vs exact)
    if(m_plt_error_u) ncomp += 1;

    // Error in v (computed vs exact)
    if(m_plt_error_v) ncomp += 1;

#if (AMREX_SPACEDIM == 3)
    // Error in w (computed vs exact)
    if(m_plt_error_w) ncomp += 1;
#endif

    // Error in nodal pressure (computed vs exact)
    if(m_plt_error_p) ncomp += 1;

    // Error in MAC pressure (computed vs exact)
    if(m_plt_error_mac_p) ncomp += 1;

    // Apparent viscosity
    if(m_plt_eta) ++ncomp;

    // Vorticity
    if(m_plt_vort) ++ncomp;

    // Forcing terms in velocity update
    if(m_plt_forcing) ncomp += 3;

    // Magnitude of the rate-of-strain tensor
    if(m_plt_strainrate) ++ncomp;

    // Divergence of velocity field
    if(m_plt_divu) ++ncomp;

#ifdef AMREX_USE_EB
    // Cut cell volume fraction
    if(m_plt_vfrac) ++ncomp;
#endif

    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
    }

    Vector<std::string> pltscaVarsName;
    int icomp = 0;
    if (m_plt_velx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velx");
        ++icomp;
    }
    if (m_plt_vely) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vely");
        ++icomp;
    }
#if (AMREX_SPACEDIM == 3)
    if (m_plt_velz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velz");
        ++icomp;
    }
#endif
    if (m_plt_gpx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpx");
        ++icomp;
    }
    if (m_plt_gpy) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpy");
        ++icomp;
    }
#if (AMREX_SPACEDIM == 3)
    if (m_plt_gpz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpz");
        ++icomp;
    }
#endif
    if (m_plt_rho) {
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(mf[lev], m_leveldata[lev]->density, 0, icomp, 1, 0);
        pltscaVarsName.push_back("density");
        ++icomp;
    }
    if (m_plt_tracer) {
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(mf[lev], m_leveldata[lev]->tracer, 0, icomp, m_ntrac, 0);
        for (int i = 0; i < m_ntrac; ++i) {
            pltscaVarsName.push_back("tracer"+std::to_string(i));
        }
        icomp += m_ntrac;
    }
    if (m_plt_p_nd) {
        for (int lev = 0; lev <= finest_level; ++lev)
            amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p_nd, 0, 1);
        pltscaVarsName.push_back("p_nd");
        ++icomp;
    }

    if (m_plt_macphi) {
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(mf[lev], m_leveldata[lev]->mac_phi, 0, icomp, 1, 0);
        pltscaVarsName.push_back("mac_phi");
        ++icomp;
    }

    if (m_plt_error_u) {
        int icomp_err_u = 0;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 0, icomp, 1, 0);
            DiffFromExact(lev, Geom(lev), m_cur_time, m_dt, mf[lev], icomp, icomp_err_u);
            amrex::Print() << "Norm0 / Norm2 of u error " <<
                mf[lev].norm0(icomp) << " " << mf[lev].norm2(icomp) / std::sqrt(mf[lev].boxArray().numPts()) << std::endl;
        }
        pltscaVarsName.push_back("error_u");
        ++icomp;
    }

    if (m_plt_error_v) {
        int icomp_err_v = 1;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 1, icomp, 1, 0);
            DiffFromExact(lev, Geom(lev), m_cur_time, m_dt, mf[lev], icomp, icomp_err_v);
            amrex::Print() << "Norm0 / Norm2 of v error " <<
                mf[lev].norm0(icomp) << " " << mf[lev].norm2(icomp) / std::sqrt(mf[lev].boxArray().numPts()) << std::endl;
        }
        pltscaVarsName.push_back("error_v");
        ++icomp;
    }

#if (AMREX_SPACEDIM == 3)
    if (m_plt_error_w) {
        int icomp_err_w = 2;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 2, icomp, 1, 0);
            DiffFromExact(lev, Geom(lev), m_cur_time, m_dt, mf[lev], icomp, icomp_err_w);
            amrex::Print() << "Norm0 / Norm2 of w error " <<
                mf[lev].norm0(icomp) << " " << mf[lev].norm2(icomp) / std::sqrt(mf[lev].boxArray().numPts()) << std::endl;
        }
        pltscaVarsName.push_back("error_w");
        ++icomp;
    }
#endif

    if (m_plt_error_p) {
        int icomp_err_p = AMREX_SPACEDIM;
        for (int lev = 0; lev <= finest_level; ++lev)
            amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p_nd, 0, 1);

        Real offset = mf[0].sum(icomp,true);
        ParallelDescriptor::ReduceRealSum(offset);
        offset *= Real(1.0)/grids[0].numPts();

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            mf[lev].plus(-offset, icomp, 1);
            DiffFromExact(lev, Geom(lev), m_cur_time, m_dt, mf[lev], icomp, icomp_err_p);
            amrex::Print() << "Norm0 / Norm2 of p error " <<
                mf[lev].norm0(icomp) << " " << mf[lev].norm2(icomp) / std::sqrt(mf[lev].boxArray().numPts()) << std::endl;
        }
        pltscaVarsName.push_back("error_p");
        ++icomp;
    }

    if (m_plt_error_mac_p) {
        int icomp_err_mac_p = AMREX_SPACEDIM+1;
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(mf[lev], m_leveldata[lev]->mac_phi, 0, icomp, 1, 0);

        Real offset = mf[0].sum(icomp,true);
        ParallelDescriptor::ReduceRealSum(offset);
        offset *= Real(1.0)/grids[0].numPts();

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            mf[lev].plus(-offset, icomp, 1);
            DiffFromExact(lev, Geom(lev), m_cur_time, m_dt, mf[lev], icomp, icomp_err_mac_p);
            amrex::Print() << "Norm0 / Norm2 of mac_p error " <<
                mf[lev].norm0(icomp) << " " << mf[lev].norm2(icomp) / std::sqrt(mf[lev].boxArray().numPts()) << std::endl;
        }
        pltscaVarsName.push_back("error_mac_p");
        ++icomp;
    }

    if (m_plt_eta) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab vel_eta(mf[lev], amrex::make_alias, icomp, 1);
            compute_viscosity_at_level(lev,
                                       &vel_eta,
                                       &m_leveldata[lev]->density,
                                       &m_leveldata[lev]->velocity,
                                       Geom(lev),
                                       m_cur_time, 0);
        }
        pltscaVarsName.push_back("eta");
        ++icomp;
    }
    if (m_plt_vort) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            (m_leveldata[lev]->velocity).FillBoundary(geom[lev].periodicity());
            MultiFab vort(mf[lev], amrex::make_alias, icomp, 1);
            ComputeVorticity(lev, m_cur_time, vort, m_leveldata[lev]->velocity);
        }
        pltscaVarsName.push_back("vort");
        ++icomp;
    }
    if (m_plt_forcing) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab forcing(mf[lev], amrex::make_alias, icomp, 3);
            compute_vel_forces_on_level(lev, forcing,
                                        m_leveldata[lev]->velocity,
                                        m_leveldata[lev]->density,
                                        m_leveldata[lev]->tracer,
                                        m_leveldata[lev]->tracer);
        }
        pltscaVarsName.push_back("forcing_x");
        pltscaVarsName.push_back("forcing_y");
        pltscaVarsName.push_back("forcing_z");
        icomp += 3;
    }
    if (m_plt_strainrate) {

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab strainrate(mf[lev], amrex::make_alias, icomp, 1);
            compute_strainrate_at_level(lev,
                                        &strainrate,
                                        &m_leveldata[lev]->velocity,
                                        Geom(lev),
                                        m_cur_time, 0);
        }
        pltscaVarsName.push_back("strainrate");
        ++icomp;
    }
    if (m_plt_divu) {
        amrex::Abort("plt_divu: xxxxx TODO");
        pltscaVarsName.push_back("divu");
        ++icomp;
    }
#ifdef AMREX_USE_EB
    if (m_plt_vfrac) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], EBFactory(lev).getVolFrac(), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vfrac");
        ++icomp;
    }
#endif

#ifdef AMREX_USE_EB
    for (int lev = 0; lev <= finest_level; ++lev) {
        EB_set_covered(mf[lev], 0.0);
    }
#endif

    AMREX_ALWAYS_ASSERT(ncomp == static_cast<int>(pltscaVarsName.size()));

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile,
    // but will never change unless we use subcycling.
    // If we do use subcycling, this should be a incflo class member.
    Vector<int> istep(finest_level + 1, m_nstep);

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf),
                                   pltscaVarsName, Geom(), m_cur_time, istep, refRatio());
    WriteJobInfo(plotfilename);
}


// Copy/Paste drag computation from from drangara
void incflo::PrintDragForce(std::ofstream &drag_file) {
   Print() << "PrintDragForce\n";

   Vector<MultiFab> gradx(finest_level + 1);
   Vector<MultiFab> grady(finest_level + 1);
   Vector<MultiFab> gradeb(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
       gradx[lev].define(grids[lev], dmap[lev], 2, 0, MFInfo(), Factory(lev));
       grady[lev].define(grids[lev], dmap[lev], 2, 0, MFInfo(), Factory(lev));
       gradeb[lev].define(grids[lev], dmap[lev], 2, 0, MFInfo(), Factory(lev));
       gradx[lev].setVal(0.);
       grady[lev].setVal(0.);
       gradeb[lev].setVal(0.);
   }

   Vector<Real> drag(finest_level + 1, 0.);
   Gpu::DeviceVector<Real> dv_drag(finest_level + 1, 0.);
   auto* p_dv_drag = dv_drag.data();

   for (int lev = 0; lev <= finest_level; lev++)
   {
      auto& ld = *m_leveldata[lev];
      auto const& dx = geom[lev].CellSizeArray();

      for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& fact = EBFactory(lev);
         int ncomp = ld.velocity.nComp();

         Array4<Real> const& gradx_arr         = gradx[lev].array(mfi);
         Array4<Real> const& grady_arr         = grady[lev].array(mfi);
         Array4<Real> const& gradeb_arr        = gradeb[lev].array(mfi);
         Array4<const Real> const& phi_arr     = ld.velocity.array(mfi);
         Array4<const Real> const& phi_eb_arr  = get_velocity_eb()[lev]->array(mfi);

         Array4<Real const> const& fcx   = (fact.getFaceCent())[0]->const_array(mfi);
         Array4<Real const> const& fcy   = (fact.getFaceCent())[1]->const_array(mfi);
 #if (AMREX_SPACEDIM == 3)
         Array4<Real const> const& fcz   = (fact.getFaceCent())[2]->const_array(mfi);
 #endif
         Array4<Real const> const& ccent = (fact.getCentroid()).array(mfi);
         Array4<Real const> const& bcent = (fact.getBndryCent()).array(mfi);
         Array4<Real const> const& bnorm = (fact.getBndryNormal()).array(mfi);
         Array4<Real const> const& barea = (fact.getBndryArea()).array(mfi);
         const FabArray<EBCellFlagFab>* flags = &(fact.getMultiEBCellFlagFab());
         Array4<EBCellFlag const> const& flag = flags->const_array(mfi);


         bool is_eb_dirichlet = true;
         bool is_eb_inhomog  = true;


         amrex::ParallelFor(bx, ncomp,
         [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            if (flag(i,j,k).isSingleValued()) {
               Real yloc_on_xface = fcx(i,j,k,0);
               Real xloc_on_yface = fcy(i,j,k,0);
               Real nx = bnorm(i,j,k,0);
               Real ny = bnorm(i,j,k,1);
#if (AMREX_SPACEDIM == 2)
               gradx_arr(i,j,k,n) = amrex::grad_x_of_phi_on_centroids(i,j,k,n,
                                                                      phi_arr,phi_eb_arr,
                                                                      flag,ccent,bcent,
                                                                      yloc_on_xface,
                                                                      is_eb_dirichlet,is_eb_inhomog);

               grady_arr(i,j,k,n) = amrex::grad_y_of_phi_on_centroids(i,j,k,n,
                                                                      phi_arr,phi_eb_arr,
                                                                      flag,ccent,bcent,
                                                                      xloc_on_yface,
                                                                      is_eb_dirichlet,is_eb_inhomog);

               gradeb_arr(i,j,k,n) = amrex::grad_eb_of_phi_on_centroids(i,j,k,n,
                                                                        phi_arr,phi_eb_arr,
                                                                        flag,ccent,bcent,
                                                                        nx, ny, is_eb_inhomog);


               Real cell_drag = barea(i,j,k)*dx[0]*dx[1]*(2*nx*gradx_arr(i,j,k,0) + ny*(grady_arr(i,j,k,0) + gradx_arr(i,j,k,1)));

               Gpu::Atomic::Add(&p_dv_drag[lev], cell_drag);

#endif

               //drag_sum += 2*nx*gradx_arr(i,j,k,0) + ny*(grady_arr(i,j,k,0) + gradx_arr(i,j,k,1));
            }
         });

         Gpu::synchronize();

      }
   }

   Gpu::copy(Gpu::deviceToHost, dv_drag.begin(), dv_drag.end(), drag.begin());

   // Print to drag history file
   drag_file << m_cur_time << ", ";
   for (int lev = 0; lev <= finest_level; lev++)
   {
      if (lev == finest_level){
          drag_file << drag[lev] << std::endl;
      } else {
          drag_file << drag[lev] << ", ";
      }
   }


   Vector<int> istep(finest_level + 1, m_nstep);
   Vector<MultiFab> plotmf(finest_level + 1);

   for (int lev = 0; lev <= finest_level; lev++)
   {
      auto const& fact = EBFactory(lev);
#if (AMREX_SPACEDIM == 2)
      plotmf[lev].define(grids[lev],dmap[lev],9,0);
      MultiFab::Copy(plotmf[lev], gradx[lev], 0, 0, 2, 0);
      MultiFab::Copy(plotmf[lev], grady[lev], 0, 2, 2, 0);
      MultiFab::Copy(plotmf[lev], gradeb[lev], 0, 4, 2, 0);
      {
         MultiFab temp_mf = fact.getBndryArea().ToMultiFab(0.0, 0.0);
         MultiFab::Copy(plotmf[lev], temp_mf, 0, 6, 1, 0);
      }
      {
         MultiFab temp_mf = fact.getBndryNormal().ToMultiFab(0.0, 0.0);
         MultiFab::Copy(plotmf[lev], temp_mf, 0, 7, 2, 0);
      }
#endif
   }

   const std::string& pfname = amrex::Concatenate("drag_terms", m_nstep);

   // Write the plotfile
   if (m_eb_flow_plt_drag && ((m_nstep % m_eb_flow_plt_drag_int) == 0)){
       amrex::WriteMultiLevelPlotfile(pfname, finest_level + 1, GetVecOfConstPtrs(plotmf),
               {"dudx", "dvdx", "dudy", "dvdy", "dudn", "dvdn",
               "barea", "bnrmx", "bnrmy"},
               Geom(), m_cur_time, istep, refRatio());
   }
}

