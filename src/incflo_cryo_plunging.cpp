#include <incflo.H>
#include <prob_bc.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void incflo::cryo_set_zero_vel ()
{
    BL_PROFILE("incflo::cryo_set_zero_vel");

    // Grid
    if (2976 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real> const& tracer = ld.tracer.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real rz = 3/2.;
                    Real ry = 0.025;
                    Real disp = m_cur_time*1.;
                    Real d = disp<(1.+2.*rz)?disp:(1.+2.*rz);
                    // if ((x*x+(z+rz)*(z+rz)<=rz*rz)) {
                    // if ((x*x+(z+2*rz)*(z+2*rz)<=rz*rz) && (y>=-ry && y<=ry)) {
                    // if ((x*x+(z-rz+dx[2]+d)*(z-rz+dx[2]+d)<=rz*rz) && (y>=-ry && y<=ry)) {
                    if ((x*x+(z+dx[2]+rz)*(z+dx[2]+rz)<=rz*rz) && (y>=-ry && y<=ry)) {
                        vel(i,j,k,0) = 0.;
                        vel(i,j,k,1) = 0.;
                        vel(i,j,k,2) = disp<(1.+2.*rz)?-1.:0.;
                    }
                    // if ((x*x+(z+rz)*(z+rz)<=rz*rz) && (y>=-ry && y<=ry)) {
                    // // if ((x*x+(z-rz+dx[2]+d)*(z-rz+dx[2]+d)<=rz*rz) && (y>=-ry && y<=ry)) {
                    // // if ((x*x+(z-rz+dx[2]+d)*(z-rz+dx[2]+d)<=rz*rz) && (y>=-ry && y<=ry)) {
                    //     if (z>-1.*dx[2] && disp<=2*rz) tracer(i,j,k) = 239.16524;
                    // }
                });
            }
        }
    }

    // Tcouple (sphere)
    else if (2977 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real> const& tracer = ld.tracer.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real r = 0.050;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);
                    if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<=r*r) {
                        vel(i,j,k,0) = 0.;
                        vel(i,j,k,1) = 0.;
                        vel(i,j,k,2) = v;
                        // vel(i,j,k,2) = m_cur_time<4.2?v:0.;
                        // if (z>-2.*dx[2]) tracer(i,j,k) = 239.16524;
                        if (z>-2.*dx[2]) tracer(i,j,k) = 219.42506;
                    }
                    // Set top row back to ethane temperature
                    if (disp>2.*r+dx[2] && z>-1.*dx[2]) tracer(i,j,k) = 88.53434;
                });
            }
        }
    }

    // Tcouple (sphere) with energy
    else if (2978 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.velocity); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real> const& tracer = ld.tracer.array(mfi);
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real r = 0.050;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);
                    // Measured temperature with 1 m/s
                    Real temp_eth = 88.53434;
                    // Real temp_tcp = 239.16524;
                    Real temp_tcp = 219.42506;
                    Real cv_eth = 1.55;
                    Real cv_tcp = 5.4406;

                    // if (x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2])<r*r) {
                    // if (x*x+y*y+(z+1.5*r+dx[2])*(z+1.5*r+dx[2])<r*r) {
                    if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<r*r) {
                        vel(i,j,k,0) = 0.;
                        vel(i,j,k,1) = 0.;
                        vel(i,j,k,2) = v;
                        if (z>-2.*dx[2]) {
                            energy(i,j,k) = temp_tcp * cv_tcp;
                            tracer(i,j,k) = temp_tcp;// * cv_tcp;
                        }
                    }

                    // // Smoothing
                    // Real phi = x*x+y*y+(z-r+d+dx[2])*(z-r+d+dx[2]) - r*r;
                    // // Real phi = x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2]) - r*r;
                    // Real eps = 1.0*dx[0];
                    // if (phi<=0)
                    // {
                    //     vel(i,j,k,0) = 0.;
                    //     vel(i,j,k,1) = 0.;
                    //     vel(i,j,k,2) = v;
                    // }

                    // if (phi<=eps)
                    // {   
                    //     if (phi<-eps) 
                    //     {   // Tcouple
                    //         if (z>-2.*dx[2]) energy(i,j,k) = temp_tcp * cv_tcp;
                    //     } else {
                    //         // 2020 JFM Heaviside function
                    //         Real xx=phi/eps;
                    //         Real energy_tcp = temp_tcp * cv_tcp;
                    //         Real energy_eth = temp_eth * cv_eth;
                    //         if (z>-2.*dx[2]) energy(i,j,k) = temp_tcp * cv_tcp;
                    //         // if (z>-2.*dx[2]) energy(i,j,k) = energy_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(energy_tcp-energy_eth);
                    //         // tracer(i,j,k) = temp_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(temp_tcp-temp_eth);
                    //     }
                    // }

                    // Set top row back to ethane temperature
                    if (disp>2.*r+dx[2] && z>-1.*dx[2]) {
                        energy(i,j,k) = temp_eth * cv_eth;
                        tracer(i,j,k) = temp_eth;// * cv_eth;
                    }
                });
            }
        }
    }

    // Disk with energy
    else if (2979 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.velocity); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real> const& tracer = ld.tracer.array(mfi);
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real rz = 3/2.;
                    Real ry = 0.025;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);
                    // Measured temperature with 1 m/s
                    Real temp_eth = 88.53434;
                    // Real temp_tcp = 239.16524;
                    Real temp_tcp = 219.42506;
                    Real cv_eth = 1.55;
                    Real cv_tcp = 5.4406;

                    if ((x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d)<rz*rz) && (y>-ry && y<ry)) {
                        vel(i,j,k,0) = 0.;
                        vel(i,j,k,1) = 0.;
                        vel(i,j,k,2) = v;
                        if (z>-2.*dx[2]) {
                            energy(i,j,k) = temp_tcp * cv_tcp;
                            tracer(i,j,k) = temp_tcp;// * cv_tcp;
                        }
                    }

                    // Set top row back to ethane temperature
                    if (disp>2.*rz+dx[2] && z>-1.*dx[2]) {
                        energy(i,j,k) = temp_eth * cv_eth;
                        tracer(i,j,k) = temp_eth;// * cv_eth;
                    }
                });
            }
        }
    }

    // Grid with energy
    else if (2980 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.velocity); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real> const& tracer = ld.tracer.array(mfi);
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real rz = 3.05/2.;
                    Real rd = 0.24/2.;
                    Real ry = 0.025/2.;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);
                    // Measured temperature with 1 m/s
                    Real temp_eth = 88.53434;
                    // Real temp_tcp = 239.16524;
                    Real temp_tcp = 219.42506;
                    Real temp_spl = temp_tcp;
                    Real cv_eth = 1.55;
                    Real cv_tcp = 5.4406;
                    Real cv_spl = 6.5421875;

                    Real outer_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                    Real inner_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                    Real grid_gap = 0.06;
                    Real zz=z-rz+2*dx[2]+d;

                    if (outer_grid<rz*rz && y*y<ry*ry) {
                        if (inner_grid>(rz-rd)*(rz-rd) ||
                            (x<ry && x>-ry) ||
                            (x<ry-1*grid_gap && x>-ry-1*grid_gap) ||
                            (x<ry-2*grid_gap && x>-ry-2*grid_gap) ||
                            (x<ry-3*grid_gap && x>-ry-3*grid_gap) ||
                            (x<ry-4*grid_gap && x>-ry-4*grid_gap) ||
                            (x<ry-5*grid_gap && x>-ry-5*grid_gap) ||
                            (x<ry-6*grid_gap && x>-ry-6*grid_gap) ||
                            (x<ry-7*grid_gap && x>-ry-7*grid_gap) ||
                            (x<ry-8*grid_gap && x>-ry-8*grid_gap) ||
                            (x<ry-9*grid_gap && x>-ry-9*grid_gap) ||
                            (x<ry-10*grid_gap && x>-ry-10*grid_gap) ||
                            (x<ry-11*grid_gap && x>-ry-11*grid_gap) ||
                            (x<ry-12*grid_gap && x>-ry-12*grid_gap) ||
                            (x<ry-13*grid_gap && x>-ry-13*grid_gap) ||
                            (x<ry-14*grid_gap && x>-ry-14*grid_gap) ||
                            (x<ry-15*grid_gap && x>-ry-15*grid_gap) ||
                            (x<ry-16*grid_gap && x>-ry-16*grid_gap) ||
                            (x<ry-17*grid_gap && x>-ry-17*grid_gap) ||
                            (x<ry-18*grid_gap && x>-ry-18*grid_gap) ||
                            (x<ry-19*grid_gap && x>-ry-19*grid_gap) ||
                            (x<ry-20*grid_gap && x>-ry-20*grid_gap) ||
                            (x<ry-21*grid_gap && x>-ry-21*grid_gap) ||
                            (x<ry-22*grid_gap && x>-ry-22*grid_gap) ||
                            (x<ry-23*grid_gap && x>-ry-23*grid_gap) ||
                            (x<ry+1*grid_gap && x>-ry+1*grid_gap) ||
                            (x<ry+2*grid_gap && x>-ry+2*grid_gap) ||
                            (x<ry+3*grid_gap && x>-ry+3*grid_gap) ||
                            (x<ry+4*grid_gap && x>-ry+4*grid_gap) ||
                            (x<ry+5*grid_gap && x>-ry+5*grid_gap) ||
                            (x<ry+6*grid_gap && x>-ry+6*grid_gap) ||
                            (x<ry+7*grid_gap && x>-ry+7*grid_gap) ||
                            (x<ry+8*grid_gap && x>-ry+8*grid_gap) ||
                            (x<ry+9*grid_gap && x>-ry+9*grid_gap) ||
                            (x<ry+10*grid_gap && x>-ry+10*grid_gap) ||
                            (x<ry+11*grid_gap && x>-ry+11*grid_gap) ||
                            (x<ry+12*grid_gap && x>-ry+12*grid_gap) ||
                            (x<ry+13*grid_gap && x>-ry+13*grid_gap) ||
                            (x<ry+14*grid_gap && x>-ry+14*grid_gap) ||
                            (x<ry+15*grid_gap && x>-ry+15*grid_gap) ||
                            (x<ry+16*grid_gap && x>-ry+16*grid_gap) ||
                            (x<ry+17*grid_gap && x>-ry+17*grid_gap) ||
                            (x<ry+18*grid_gap && x>-ry+18*grid_gap) ||
                            (x<ry+19*grid_gap && x>-ry+19*grid_gap) ||
                            (x<ry+20*grid_gap && x>-ry+20*grid_gap) ||
                            (x<ry+21*grid_gap && x>-ry+21*grid_gap) ||
                            (x<ry+22*grid_gap && x>-ry+22*grid_gap) ||
                            (x<ry+23*grid_gap && x>-ry+23*grid_gap) ||
                            (zz<ry && zz>-ry) ||
                            (zz<ry-1*grid_gap && zz>-ry-1*grid_gap) ||
                            (zz<ry-2*grid_gap && zz>-ry-2*grid_gap) ||
                            (zz<ry-3*grid_gap && zz>-ry-3*grid_gap) ||
                            (zz<ry-4*grid_gap && zz>-ry-4*grid_gap) ||
                            (zz<ry-5*grid_gap && zz>-ry-5*grid_gap) ||
                            (zz<ry-6*grid_gap && zz>-ry-6*grid_gap) ||
                            (zz<ry-7*grid_gap && zz>-ry-7*grid_gap) ||
                            (zz<ry-8*grid_gap && zz>-ry-8*grid_gap) ||
                            (zz<ry-9*grid_gap && zz>-ry-9*grid_gap) ||
                            (zz<ry-10*grid_gap && zz>-ry-10*grid_gap) ||
                            (zz<ry-11*grid_gap && zz>-ry-11*grid_gap) ||
                            (zz<ry-12*grid_gap && zz>-ry-12*grid_gap) ||
                            (zz<ry-13*grid_gap && zz>-ry-13*grid_gap) ||
                            (zz<ry-14*grid_gap && zz>-ry-14*grid_gap) ||
                            (zz<ry-15*grid_gap && zz>-ry-15*grid_gap) ||
                            (zz<ry-16*grid_gap && zz>-ry-16*grid_gap) ||
                            (zz<ry-17*grid_gap && zz>-ry-17*grid_gap) ||
                            (zz<ry-18*grid_gap && zz>-ry-18*grid_gap) ||
                            (zz<ry-19*grid_gap && zz>-ry-19*grid_gap) ||
                            (zz<ry-20*grid_gap && zz>-ry-20*grid_gap) ||
                            (zz<ry-21*grid_gap && zz>-ry-21*grid_gap) ||
                            (zz<ry-22*grid_gap && zz>-ry-22*grid_gap) ||
                            (zz<ry-23*grid_gap && zz>-ry-23*grid_gap) ||
                            (zz<ry+1*grid_gap && zz>-ry+1*grid_gap) ||
                            (zz<ry+2*grid_gap && zz>-ry+2*grid_gap) ||
                            (zz<ry+3*grid_gap && zz>-ry+3*grid_gap) ||
                            (zz<ry+4*grid_gap && zz>-ry+4*grid_gap) ||
                            (zz<ry+5*grid_gap && zz>-ry+5*grid_gap) ||
                            (zz<ry+6*grid_gap && zz>-ry+6*grid_gap) ||
                            (zz<ry+7*grid_gap && zz>-ry+7*grid_gap) ||
                            (zz<ry+8*grid_gap && zz>-ry+8*grid_gap) ||
                            (zz<ry+9*grid_gap && zz>-ry+9*grid_gap) ||
                            (zz<ry+10*grid_gap && zz>-ry+10*grid_gap) ||
                            (zz<ry+11*grid_gap && zz>-ry+11*grid_gap) ||
                            (zz<ry+12*grid_gap && zz>-ry+12*grid_gap) ||
                            (zz<ry+13*grid_gap && zz>-ry+13*grid_gap) ||
                            (zz<ry+14*grid_gap && zz>-ry+14*grid_gap) ||
                            (zz<ry+15*grid_gap && zz>-ry+15*grid_gap) ||
                            (zz<ry+16*grid_gap && zz>-ry+16*grid_gap) ||
                            (zz<ry+17*grid_gap && zz>-ry+17*grid_gap) ||
                            (zz<ry+18*grid_gap && zz>-ry+18*grid_gap) ||
                            (zz<ry+19*grid_gap && zz>-ry+19*grid_gap) ||
                            (zz<ry+20*grid_gap && zz>-ry+20*grid_gap) ||
                            (zz<ry+21*grid_gap && zz>-ry+21*grid_gap) ||
                            (zz<ry+22*grid_gap && zz>-ry+22*grid_gap) ||
                            (zz<ry+23*grid_gap && zz>-ry+23*grid_gap) ) {
                            vel(i,j,k,0) = 0.;
                            vel(i,j,k,1) = 0.;
                            vel(i,j,k,2) = v;
                            if (z>-2.*dx[2]) {
                                energy(i,j,k) = temp_tcp * cv_tcp;
                                tracer(i,j,k) = temp_tcp;// * cv_tcp;
                            }
                        }
                    }
                    // Add biological samples
                    if (inner_grid<0.25*0.25 && (y>ry && y<ry+0.025)) {
                    // if (inner_grid<pow(0.04/2,2) && (y>ry && y<ry+0.01)) {
                        vel(i,j,k,0) = 0.;
                        vel(i,j,k,1) = 0.;
                        vel(i,j,k,2) = v;
                        if (z>-2.*dx[2]) {
                            energy(i,j,k) = temp_spl * cv_spl;
                            tracer(i,j,k) = temp_spl;// * cv_tcp;
                        }
                    }

                    // Set top row back to ethane temperature
                    if (disp>2.*rz+dx[2] && z>-1.*dx[2]) {
                        energy(i,j,k) = temp_eth * cv_eth;
                        tracer(i,j,k) = temp_eth;// * cv_eth;
                    }
                });
            }
        }
    }
}

void incflo::cryo_update_thermal () {
    BL_PROFILE("incflo::cryo_update_thermal");

    // TODO: Make this generic
    if (2978 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& density = ld.density.array(mfi);
                Array4<Real> const& vel = ld.density.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real r = 0.050;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;

                    // Thermal conductivity (converted)
                    Real cv_eth = 1.55;
                    Real cv_tcp = 5.4406;
                    Real k_eth = 0.0003875;
                    Real k_tcp = 0.0340625;

                    // Smoothing
                    // Real phi = x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2]) - r*r;
                    Real phi = x*x+y*y+(z-r+d+dx[2])*(z-r+d+dx[2]) - r*r;
                    Real eps = 2.0*dx[0];
                    if (phi>eps)
                    {   // Ethane
                        vol_heat_cap(i,j,k) = cv_eth;
                        ther_conduct(i,j,k) = k_eth;
                    // } else
                    // {   // Tcouple
                    //     vol_heat_cap(i,j,k) = cv_tcp;
                    //     ther_conduct(i,j,k) = k_tcp;
                    // }
                    } else if (phi<-eps) 
                    {   // Tcouple
                        vol_heat_cap(i,j,k) = cv_tcp;
                        ther_conduct(i,j,k) = k_tcp;
                    } else {
                        // 2020 JFM Heaviside function
                        Real xx=phi/eps;
                        vol_heat_cap(i,j,k) = cv_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(cv_tcp-cv_eth);
                        ther_conduct(i,j,k) = k_tcp-0.5*(1.+xx+1./M_PI*sin(M_PI*xx))*(k_tcp-k_eth);
                    }

                    // if (x*x+y*y+(z+1.5*r+dx[2])*(z+1.5*r+dx[2])<r*r) {
                    // // if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<(r+2*dx[0])*(r+2*dx[0])) {
                    //     // tcouple
                    //     vol_heat_cap(i,j,k) = cv_tcp;
                    //     ther_conduct(i,j,k) = k_tcp;
                    // } else {
                    //     // ethane
                    //     vol_heat_cap(i,j,k) = cv_eth;
                    //     ther_conduct(i,j,k) = k_eth;
                    // }
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype

    // Disk with energy
    else if (2979 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& density = ld.density.array(mfi);
                Array4<Real> const& vel = ld.density.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real rz = 3/2.;
                    Real ry = 0.025;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);

                    // Thermal conductivity (converted)
                    Real cv_eth = 1.55;
                    Real cv_tcp = 3.76953125; //plunger
                    Real k_eth = 0.0003875;
                    Real k_tcp = 0.0005390625; // plunger

                    if ((x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d)<rz*rz) && (y>-ry && y<ry)) {
                        // tcouple
                        vol_heat_cap(i,j,k) = cv_tcp;
                        ther_conduct(i,j,k) = k_tcp;
                    } else {
                        // ethane
                        vol_heat_cap(i,j,k) = cv_eth;
                        ther_conduct(i,j,k) = k_eth;
                    }
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype

    // Grid with energy
    else if (2980 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& density = ld.density.array(mfi);
                Array4<Real> const& vel = ld.density.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                Array4<Real> const& ther_conduct = ld.ther_conduct.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
                    Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
                    Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

                    // Set velocity field inside the plunged objects to 0
                    Real rz = 3.05/2.;
                    Real rd = 0.24/2.;
                    Real ry = 0.025/2.;
                    // // 1 m/s
                    // Real v0 = 0.7408041682457519;
                    // Real acc = -0.04623555291981721;
                    // 500 mm/s
                    Real v0 = 0.485047915352808;
                    Real acc = -0.03131555980546892;
                    Real disp = v0*m_cur_time+0.5*acc*m_cur_time*m_cur_time;
                    Real d = disp;
                    Real v = -(v0+acc*m_cur_time);

                    // Thermal conductivity (converted)
                    Real cv_eth = 1.55;
                    Real cv_tcp = 3.76953125; //plunger
                    Real cv_spl = 6.5421875; // sample
                    Real k_eth = 0.0003875;
                    Real k_tcp = 0.0005390625; // plunger
                    Real k_spl = 0.0009344219; // sample

                    Real outer_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                    Real inner_grid = x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d);
                    Real grid_gap = 0.06;
                    Real zz=z-rz+2*dx[2]+d;

                    if (outer_grid<rz*rz && y*y<ry*ry) {
                        if (inner_grid>(rz-rd)*(rz-rd) ||
                            (x<ry && x>-ry) ||
                            (x<ry-1*grid_gap && x>-ry-1*grid_gap) ||
                            (x<ry-2*grid_gap && x>-ry-2*grid_gap) ||
                            (x<ry-3*grid_gap && x>-ry-3*grid_gap) ||
                            (x<ry-4*grid_gap && x>-ry-4*grid_gap) ||
                            (x<ry-5*grid_gap && x>-ry-5*grid_gap) ||
                            (x<ry-6*grid_gap && x>-ry-6*grid_gap) ||
                            (x<ry-7*grid_gap && x>-ry-7*grid_gap) ||
                            (x<ry-8*grid_gap && x>-ry-8*grid_gap) ||
                            (x<ry-9*grid_gap && x>-ry-9*grid_gap) ||
                            (x<ry-10*grid_gap && x>-ry-10*grid_gap) ||
                            (x<ry-11*grid_gap && x>-ry-11*grid_gap) ||
                            (x<ry-12*grid_gap && x>-ry-12*grid_gap) ||
                            (x<ry-13*grid_gap && x>-ry-13*grid_gap) ||
                            (x<ry-14*grid_gap && x>-ry-14*grid_gap) ||
                            (x<ry-15*grid_gap && x>-ry-15*grid_gap) ||
                            (x<ry-16*grid_gap && x>-ry-16*grid_gap) ||
                            (x<ry-17*grid_gap && x>-ry-17*grid_gap) ||
                            (x<ry-18*grid_gap && x>-ry-18*grid_gap) ||
                            (x<ry-19*grid_gap && x>-ry-19*grid_gap) ||
                            (x<ry-20*grid_gap && x>-ry-20*grid_gap) ||
                            (x<ry-21*grid_gap && x>-ry-21*grid_gap) ||
                            (x<ry-22*grid_gap && x>-ry-22*grid_gap) ||
                            (x<ry-23*grid_gap && x>-ry-23*grid_gap) ||
                            (x<ry+1*grid_gap && x>-ry+1*grid_gap) ||
                            (x<ry+2*grid_gap && x>-ry+2*grid_gap) ||
                            (x<ry+3*grid_gap && x>-ry+3*grid_gap) ||
                            (x<ry+4*grid_gap && x>-ry+4*grid_gap) ||
                            (x<ry+5*grid_gap && x>-ry+5*grid_gap) ||
                            (x<ry+6*grid_gap && x>-ry+6*grid_gap) ||
                            (x<ry+7*grid_gap && x>-ry+7*grid_gap) ||
                            (x<ry+8*grid_gap && x>-ry+8*grid_gap) ||
                            (x<ry+9*grid_gap && x>-ry+9*grid_gap) ||
                            (x<ry+10*grid_gap && x>-ry+10*grid_gap) ||
                            (x<ry+11*grid_gap && x>-ry+11*grid_gap) ||
                            (x<ry+12*grid_gap && x>-ry+12*grid_gap) ||
                            (x<ry+13*grid_gap && x>-ry+13*grid_gap) ||
                            (x<ry+14*grid_gap && x>-ry+14*grid_gap) ||
                            (x<ry+15*grid_gap && x>-ry+15*grid_gap) ||
                            (x<ry+16*grid_gap && x>-ry+16*grid_gap) ||
                            (x<ry+17*grid_gap && x>-ry+17*grid_gap) ||
                            (x<ry+18*grid_gap && x>-ry+18*grid_gap) ||
                            (x<ry+19*grid_gap && x>-ry+19*grid_gap) ||
                            (x<ry+20*grid_gap && x>-ry+20*grid_gap) ||
                            (x<ry+21*grid_gap && x>-ry+21*grid_gap) ||
                            (x<ry+22*grid_gap && x>-ry+22*grid_gap) ||
                            (x<ry+23*grid_gap && x>-ry+23*grid_gap) ||
                            (zz<ry && zz>-ry) ||
                            (zz<ry-1*grid_gap && zz>-ry-1*grid_gap) ||
                            (zz<ry-2*grid_gap && zz>-ry-2*grid_gap) ||
                            (zz<ry-3*grid_gap && zz>-ry-3*grid_gap) ||
                            (zz<ry-4*grid_gap && zz>-ry-4*grid_gap) ||
                            (zz<ry-5*grid_gap && zz>-ry-5*grid_gap) ||
                            (zz<ry-6*grid_gap && zz>-ry-6*grid_gap) ||
                            (zz<ry-7*grid_gap && zz>-ry-7*grid_gap) ||
                            (zz<ry-8*grid_gap && zz>-ry-8*grid_gap) ||
                            (zz<ry-9*grid_gap && zz>-ry-9*grid_gap) ||
                            (zz<ry-10*grid_gap && zz>-ry-10*grid_gap) ||
                            (zz<ry-11*grid_gap && zz>-ry-11*grid_gap) ||
                            (zz<ry-12*grid_gap && zz>-ry-12*grid_gap) ||
                            (zz<ry-13*grid_gap && zz>-ry-13*grid_gap) ||
                            (zz<ry-14*grid_gap && zz>-ry-14*grid_gap) ||
                            (zz<ry-15*grid_gap && zz>-ry-15*grid_gap) ||
                            (zz<ry-16*grid_gap && zz>-ry-16*grid_gap) ||
                            (zz<ry-17*grid_gap && zz>-ry-17*grid_gap) ||
                            (zz<ry-18*grid_gap && zz>-ry-18*grid_gap) ||
                            (zz<ry-19*grid_gap && zz>-ry-19*grid_gap) ||
                            (zz<ry-20*grid_gap && zz>-ry-20*grid_gap) ||
                            (zz<ry-21*grid_gap && zz>-ry-21*grid_gap) ||
                            (zz<ry-22*grid_gap && zz>-ry-22*grid_gap) ||
                            (zz<ry-23*grid_gap && zz>-ry-23*grid_gap) ||
                            (zz<ry+1*grid_gap && zz>-ry+1*grid_gap) ||
                            (zz<ry+2*grid_gap && zz>-ry+2*grid_gap) ||
                            (zz<ry+3*grid_gap && zz>-ry+3*grid_gap) ||
                            (zz<ry+4*grid_gap && zz>-ry+4*grid_gap) ||
                            (zz<ry+5*grid_gap && zz>-ry+5*grid_gap) ||
                            (zz<ry+6*grid_gap && zz>-ry+6*grid_gap) ||
                            (zz<ry+7*grid_gap && zz>-ry+7*grid_gap) ||
                            (zz<ry+8*grid_gap && zz>-ry+8*grid_gap) ||
                            (zz<ry+9*grid_gap && zz>-ry+9*grid_gap) ||
                            (zz<ry+10*grid_gap && zz>-ry+10*grid_gap) ||
                            (zz<ry+11*grid_gap && zz>-ry+11*grid_gap) ||
                            (zz<ry+12*grid_gap && zz>-ry+12*grid_gap) ||
                            (zz<ry+13*grid_gap && zz>-ry+13*grid_gap) ||
                            (zz<ry+14*grid_gap && zz>-ry+14*grid_gap) ||
                            (zz<ry+15*grid_gap && zz>-ry+15*grid_gap) ||
                            (zz<ry+16*grid_gap && zz>-ry+16*grid_gap) ||
                            (zz<ry+17*grid_gap && zz>-ry+17*grid_gap) ||
                            (zz<ry+18*grid_gap && zz>-ry+18*grid_gap) ||
                            (zz<ry+19*grid_gap && zz>-ry+19*grid_gap) ||
                            (zz<ry+20*grid_gap && zz>-ry+20*grid_gap) ||
                            (zz<ry+21*grid_gap && zz>-ry+21*grid_gap) ||
                            (zz<ry+22*grid_gap && zz>-ry+22*grid_gap) ||
                            (zz<ry+23*grid_gap && zz>-ry+23*grid_gap) ) {
                            // tcouple
                            vol_heat_cap(i,j,k) = cv_tcp;
                            ther_conduct(i,j,k) = k_tcp;
                        }
                    } else {
                        // ethane
                        vol_heat_cap(i,j,k) = cv_eth;
                        ther_conduct(i,j,k) = k_eth;
                    }
                    if (inner_grid<0.25*0.25 && (y>ry && y<ry+0.025)) {
                    // if (inner_grid<pow(0.04/2,2) && (y>ry && y<ry+0.01)) {
                        // samples
                        vol_heat_cap(i,j,k) = cv_spl;
                        ther_conduct(i,j,k) = k_spl;
                    }
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype
}

void incflo::cryo_compute_temp ()
{
    BL_PROFILE("incflo::cryo_compute_temp");

    // TODO: Make this generic
    if (2978 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.energy); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // temp = energy / vol_heat_cap(x)
                    temp(i,j,k) = energy(i,j,k) / vol_heat_cap(i,j,k);
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype

    else if (2979 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.energy); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // temp = energy / vol_heat_cap(x)
                    temp(i,j,k) = energy(i,j,k) / vol_heat_cap(i,j,k);
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype

    else if (2980 == m_probtype) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
            Box const& domain = geom[lev].Domain();
            auto const& dx = geom[lev].CellSizeArray();
            auto const& problo = geom[lev].ProbLoArray();
            auto const& probhi = geom[lev].ProbHiArray();
            for (MFIter mfi(ld.energy); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                Array4<Real> const& energy = ld.energy.array(mfi);
                Array4<Real> const& temp = ld.temp.array(mfi);
                Array4<Real> const& vol_heat_cap = ld.vol_heat_cap.array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // temp = energy / vol_heat_cap(x)
                    temp(i,j,k) = energy(i,j,k) / vol_heat_cap(i,j,k);
                }); // amrex::ParallelFor
            } // mfi
        } // lev
    } // m_probtype
}

// void incflo::cryo_convect_energy ()
// {
//     BL_PROFILE("incflo::cryo_convect_energy");

//     // TODO: Make this generic
//     if (2978 != m_probtype) return;

//     // Follow incflo_apply_predictor.cpp advect tracer parts

//     Real l_dt = m_dt;
//     constexpr Real m_half = Real(0.5);
//     // We use the new time value for things computed on the "*" state
//     Real new_time = m_cur_time + m_dt;
    
//     // Forcing terms (?)
//     Vector<MultiFab> energy_eta;
//     // Allocate space for the forcing terms
//     for (int lev = 0; lev <= finest_level; ++lev) {
//         energy_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), Factory(lev));
//     }
    
//     // Compute diffusive coefficients
//     compute_tracer_diff_coeff(GetVecOfPtrs(energy_eta), 1);
//     // Compute explicit diffusive term
//     // Note: Here we need grad T
//     compute_laps(get_laps_old(), get_temp_old_const(), get_one_const(),
//                  GetVecOfConstPtrs(energy_eta));
//     // Compute convective term for energy
//     // Done in incflo_apply_predictor.cpp compute_convective_term_energy();

//     // Update the energy field
//     for (int lev = 0; lev <= finest_level; lev++)
//     {
//         auto& ld = *m_leveldata[lev];
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//         for (MFIter mfi(ld.energy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             Box const& bx = mfi.tilebox();
//             Array4<Real const> const& energy_o = ld.energy_o.const_array(mfi);
//             Array4<Real const> const& temp_o   = ld.temp_o.const_array(mfi);
//             Array4<Real const> const& dedt_o   = ld.conv_energy_o.const_array(mfi);
//             Array4<Real> const& energy         = ld.energy.array(mfi);

//             if (m_diff_type == DiffusionType::Explicit)
//             {
//                 // TODO: Check if laps is operator or already on vel
//                 Array4<Real const> const& laps_o = ld.laps_o.const_array(mfi);
//                 amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                 {
//                     // (energy)^new = (energy)^old + dt * ( div(energy vel) + div(k grad temp) )
//                     Real energy_new = energy_o(i,j,k) + l_dt *
//                         ( dedt_o(i,j,k) + laps_o(i,j,k) );
//                     energy(i,j,k) = energy_new;
//                 }); // amrex::ParallelFor
//             } // if (m_diff_type)
//         } // mfi
//     } // lev

//     // *************************************************************************************
//     // Solve diffusion equation for energy
//     // *************************************************************************************
//     if (m_advect_energy)
//     {
//         const int ng_diffusion = 1;
//         for (int lev = 0; lev <= finest_level; ++lev)
//             fillphysbc_density(lev, new_time, m_leveldata[lev]->energy, ng_diffusion);

//         Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : m_half*m_dt;
//         diffuse_scalar(get_energy_new(), get_one_new(), GetVecOfConstPtrs(energy_eta), dt_diff);

//     } // if (m_advect_energy)
// }

// TODO: Decide if need to move to incflo_fillpatch.cpp
void incflo::fillpatch_energy (int lev, Real time, MultiFab& energy, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloEneFill> > physbc(geom[lev], get_energy_bcrec(),
                                                            IncfloEneFill{m_probtype, m_bc_energy});
        FillPatchSingleLevel(energy, IntVect(ng), time,
                             {&(m_leveldata[lev]->energy_o),
                              &(m_leveldata[lev]->energy)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_energy_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloEneFill> > cphysbc
            (geom[lev-1], bcrec, IncfloEneFill{m_probtype, m_bc_energy});
        PhysBCFunct<GpuBndryFuncFab<IncfloEneFill> > fphysbc
            (geom[lev], bcrec, IncfloEneFill{m_probtype, m_bc_energy});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(energy, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->energy_o),
                            &(m_leveldata[lev-1]->energy)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->energy_o),
                            &(m_leveldata[lev]->energy)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, 1, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

void incflo::fillpatch_temp (int lev, Real time, MultiFab& temp, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > physbc(geom[lev], get_density_bcrec(),
                                                            IncfloDenFill{m_probtype, m_bc_density});
        FillPatchSingleLevel(temp, IntVect(ng), time,
                             {&(m_leveldata[lev]->temp_o),
                              &(m_leveldata[lev]->temp)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_density_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > cphysbc
            (geom[lev-1], bcrec, IncfloDenFill{m_probtype, m_bc_density});
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > fphysbc
            (geom[lev], bcrec, IncfloDenFill{m_probtype, m_bc_density});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(temp, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->temp_o),
                            &(m_leveldata[lev-1]->temp)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->temp_o),
                            &(m_leveldata[lev]->temp)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, 1, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

void incflo::copy_from_new_to_old_energy (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_new_to_old_energy(lev, ng);
    }
}

void incflo::copy_from_new_to_old_energy (int lev, IntVect const& ng)
{
    if (m_ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->energy_o,
                       m_leveldata[lev]->energy, 0, 0, 1, ng);
    }
}

void incflo::copy_from_old_to_new_energy (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_old_to_new_energy(lev, ng);
    }
}

void incflo::copy_from_old_to_new_energy (int lev, IntVect const& ng)
{
    if (m_ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->energy,
                       m_leveldata[lev]->energy_o, 0, 0, 1, ng);
    }
}

void incflo::cryo_copy_from_temp_to_tracer (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        cryo_copy_from_temp_to_tracer(lev, ng);
    }
}

void incflo::cryo_copy_from_temp_to_tracer (int lev, IntVect const& ng)
{
    if (m_ntrac == 1) {
        MultiFab::Copy(m_leveldata[lev]->tracer,
                       m_leveldata[lev]->temp, 0, 0, 1, ng);
    }
}
