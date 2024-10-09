#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include <incflo.H>

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

        if (m_use_cc_proj) {
            VisMF::Write(m_leveldata[lev]->p_nd,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p_nd"));
        } else {
            VisMF::Write(m_leveldata[lev]->p_cc,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p_cc"));
        }
    }

#ifdef INCFLO_USE_PARTICLES
   particleData.Checkpoint(checkpointname);
#endif
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

    if ( m_regrid_on_restart ) {
        MakeNewGrids(m_cur_time);
    } else {
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

        if (m_use_cc_proj) {
            VisMF::Read(m_leveldata[lev]->p_nd,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "p_nd"));
        } else {
            VisMF::Read(m_leveldata[lev]->p_cc,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "p_cc"));
        }
    }

#ifdef INCFLO_USE_PARTICLES
   particleData.Restart((ParGDBBase*)GetParGDB(),m_restart_file);
#endif

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

    const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename << " at time " << m_cur_time << std::endl;

    // Write the plotfile
    WritePlotVariables(m_plotVars, plotfilename);

    WriteJobInfo(plotfilename);

#ifdef INCFLO_USE_PARTICLES
    particleData.Checkpoint(plotfilename);
#endif
}

void incflo::WriteSmallPlotFile()
{
    BL_PROFILE("incflo::WriteSmallPlotFile()");

    const std::string& plotfilename = amrex::Concatenate(m_smallplot_file, m_nstep);

    amrex::Print() << "  Writing smallplotfile " << plotfilename << " at time " << m_cur_time << std::endl;

    // Write the plotfile
    WritePlotVariables(m_smallplotVars, plotfilename);
}

void incflo::WritePlotVariables(Vector<std::string> vars, const std::string& plotfilename)
{
    BL_PROFILE("incflo::WritePlotVariables()");

    if ( vars.empty() ) {
        Abort("incflo::WritePlotVariables : Must specify at least one variable to plot");
    }

    for(int n = 0; n < vars.size(); n++) {
        if (vars[n] == "vort" || vars[n] == "divu" || vars[n] == "forcing" ||
            vars[n] == "eta"  || vars[n] == "strainrate") {
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
            break;
        }
    }

    // Determine how many components are needed for plotfile MF
    int ncomp = 0;

    for(int n = 0; n < vars.size(); n++) {
        // Check for options with more than one component
        if ( vars[n] == "tracer" ) { ncomp += m_ntrac; }
        else if ( vars[n] == "forcing" ) { ncomp += AMREX_SPACEDIM; }

        // Only include these options if relevant/enabled
        else if ( vars[n] == "velz" || vars[n] == "gpz") {
#if (AMREX_SPACEDIM == 3)
            ++ncomp;
#endif
        }
        else if ( vars[n] == "vfrac" ) {
#ifdef AMREX_USE_EB
            // Cut cell volume fraction
            ++ncomp;
#endif
        }
        else if( vars[n] == "particle_count" ) {
#ifdef INCFLO_USE_PARTICLES
            // Number of particles per cell
            ++ncomp;
#endif
        }
        else { ncomp++; }
    }


    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
    }

    Vector<std::string> pltscaVarsName;
    pltscaVarsName.reserve(ncomp);
    int icomp = 0;
    for(int n = 0; n < vars.size(); n++) {
        if ( vars[n] == "velx" ) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 0, icomp, 1, 0);
            }
            pltscaVarsName.push_back("velx");
            ++icomp;
        }
        else if (vars[n] == "vely") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 1, icomp, 1, 0);
            }
            pltscaVarsName.push_back("vely");
            ++icomp;
        }
        else if (vars[n] == "velz") {
#if (AMREX_SPACEDIM == 3)
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 2, icomp, 1, 0);
            }
            pltscaVarsName.push_back("velz");
            ++icomp;
#endif
        }
        else if (vars[n] == "gpx") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 0, icomp, 1, 0);
                mf[lev].plus(m_gp0[0],icomp,1,0);
            }
            pltscaVarsName.push_back("gpx");
            ++icomp;
        }
        else if (vars[n] == "gpy") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 1, icomp, 1, 0);
                mf[lev].plus(m_gp0[1],icomp,1,0);
            }
            pltscaVarsName.push_back("gpy");
            ++icomp;
        }
        else if (vars[n] == "gpz") {
#if (AMREX_SPACEDIM == 3)
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 2, icomp, 1, 0);
                mf[lev].plus(m_gp0[2],icomp,1,0);
            }
            pltscaVarsName.push_back("gpz");
            ++icomp;
#endif
        }
        else if (vars[n] == "rho" || vars[n] == "density") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->density, 0, icomp, 1, 0);
            }
            pltscaVarsName.push_back("density");
            ++icomp;
        }
        else if (vars[n] == "tracer") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->tracer, 0, icomp, m_ntrac, 0);
            }
            for (int i = 0; i < m_ntrac; ++i) {
                pltscaVarsName.push_back("tracer"+std::to_string(i));
            }
            icomp += m_ntrac;
        }
        else if (vars[n] == "p") {
            if (m_use_cc_proj) {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    MultiFab::Copy(mf[lev], m_leveldata[lev]->p_cc, 0, icomp, 1, 0);
                }
                pltscaVarsName.push_back("p_cc");
            } else {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p_nd, 0, 1);
                }
                pltscaVarsName.push_back("p_nd");
            }
            ++icomp;
        }
        else if (vars[n] == "macphi") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], m_leveldata[lev]->mac_phi, 0, icomp, 1, 0);
            }
            pltscaVarsName.push_back("mac_phi");
            ++icomp;
        }
        else if (vars[n] == "error_u") {
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
        else if (vars[n] == "error_v") {
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
        else if (vars[n] == "error_w") {
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
        else if (vars[n] == "error_p") {
            int icomp_err_p = AMREX_SPACEDIM;
            for (int lev = 0; lev <= finest_level; ++lev)
                amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p_nd, 0, 1);

            Real offset = mf[0].sum(icomp,true);
            ParallelDescriptor::ReduceRealSum(offset);
            offset *= Real(1.0)/Real(grids[0].numPts());

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
        else if (vars[n] == "error_mac_p") {
            int icomp_err_mac_p = AMREX_SPACEDIM+1;
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Copy(mf[lev], m_leveldata[lev]->mac_phi, 0, icomp, 1, 0);

            Real offset = mf[0].sum(icomp,true);
            ParallelDescriptor::ReduceRealSum(offset);
            offset *= Real(1.0)/Real(grids[0].numPts());

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
        else if (vars[n] == "eta") {
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
        else if (vars[n] == "magvel") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab magvel(mf[lev], amrex::make_alias, icomp, 1);
                ComputeMagVel(lev, m_cur_time, magvel, m_leveldata[lev]->velocity);
            }
            pltscaVarsName.push_back("magvel");
            ++icomp;
        }

        else if (vars[n] == "vort") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                (m_leveldata[lev]->velocity).FillBoundary(geom[lev].periodicity());
                MultiFab vort(mf[lev], amrex::make_alias, icomp, 1);
                ComputeVorticity(lev, m_cur_time, vort, m_leveldata[lev]->velocity);
            }
            pltscaVarsName.push_back("vort");
            ++icomp;
        }
        else if (vars[n] == "forcing") {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab forcing(mf[lev], amrex::make_alias, icomp, AMREX_SPACEDIM);
                compute_vel_forces_on_level(lev, forcing,
                                            m_leveldata[lev]->velocity,
                                            m_leveldata[lev]->density,
                                            m_leveldata[lev]->tracer,
                                            m_leveldata[lev]->tracer);
            }
            AMREX_D_TERM(pltscaVarsName.push_back("forcing_x");,
                         pltscaVarsName.push_back("forcing_y");,
                         pltscaVarsName.push_back("forcing_z"););
            icomp += AMREX_SPACEDIM;
        }
        else if (vars[n] == "strainrate") {
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
        else if (vars[n] == "divu") {
            amrex::Abort("plt_divu: xxxxx TODO");
            pltscaVarsName.push_back("divu");
            ++icomp;
        }
        else if (vars[n] == "particle_count") {
#ifdef INCFLO_USE_PARTICLES
            const auto& particles_namelist( particleData.getNames() );
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab temp_dat(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 0);
                temp_dat.setVal(0);
                particleData[particles_namelist[0]]->Increment(temp_dat, lev);
                MultiFab::Copy(mf[lev], temp_dat, 0, icomp, 1, 0);
            }
            pltscaVarsName.push_back("particle_count");
            ++icomp;
#endif
        }
        else if (vars[n] == "vfrac") {
#ifdef AMREX_USE_EB
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf[lev], EBFactory(lev).getVolFrac(), 0, icomp, 1, 0);
            }
            pltscaVarsName.push_back("vfrac");
            ++icomp;
#endif
        }
        else {
            amrex::Abort("incflo::WritePlotfileVariables : plotfile variable '"+vars[n]+"' not found");
        }
    }

#ifdef AMREX_USE_EB
    for (int lev = 0; lev <= finest_level; ++lev) {
        EB_set_covered(mf[lev], 0.0);
    }
#endif

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile,
    // but will never change unless we use subcycling.
    // If we do use subcycling, this should be a incflo class member.
    Vector<int> istep(finest_level + 1, m_nstep);

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf),
                                   pltscaVarsName, Geom(), m_cur_time, istep, refRatio());
}
