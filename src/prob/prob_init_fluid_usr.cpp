#include <incflo.H>

using namespace amrex;

// Initialize fluid field for 3D cryo-plunging
// m_probtype = 2976, and the if-statement check is defined in
// prob_init_fluid() in prob_init_fluid.cpp
void incflo::init_cryo_plunging ( Box const& vbx, Box const& /*gbx*/,
                                  Array4<Real> const& vel,
                                  Array4<Real> const& density,
                                  Array4<Real> const& tracer,
                                  Box const& /*domain*/,
                                  GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                  GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                  GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);

        // Set initial velocity and density fields
        // TODO: double check if need incflo.ro_0 here
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        density(i,j,k) = 1.;

        // Set initial temperature field
        // TODO: also need to set vhc and k here (check how to do +2 tracer fields)
        // HACK: pre-set the plunged solid location (not moving already immersed)
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 90.3;
        if ((x>=-0.05 && x<=0.05) && (y<=0 && y>=-0.3))
            tracer(i,j,k) = 293.;
    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
        Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

        // Set initial velocity and density fields
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;
        density(i,j,k) = 1.;

        // Set initial temperature field
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 88.53434;
        // TODO: better way to initialize disk geometry
        // Real r=3./2;
        // if ((x*x+(z+r+0.1)*(z+r+0.1)<=r*r) && (y>=-0.025 && y<=0.025))
        //     tracer(i,j,k) = 239.16524;
    });

#endif
}

// Initialize fluid field for 3D thermo couple
// m_probtype = 2977, and the if-statement check is defined in
// prob_init_fluid() in prob_init_fluid.cpp
void incflo::init_thermo_couple ( Box const& vbx, Box const& /*gbx*/,
                                  Array4<Real> const& vel,
                                  Array4<Real> const& density,
                                  Array4<Real> const& tracer,
                                //   Array4<Real> const& diff_coeff,
                                  Box const& /*domain*/,
                                  GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                  GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                  GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);

        // Set initial velocity and density fields
        // TODO: double check if need incflo.ro_0 here
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        density(i,j,k) = 1.;

        // Set initial temperature field
        // TODO: also need to set vhc and k here (check how to do +2 tracer fields)
        // HACK: pre-set the plunged solid location (not moving already immersed)
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 90.3;
        Real r = 0.025;
        if (x*x+(y+r+0.1)*(y+r+0.1)<=r*r)
            // tracer(i,j,k) = 239.16524;
            tracer(i,j,k) = 219.42506;
    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
        Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

        // Set initial velocity and density fields
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;
        density(i,j,k) = 1.;
        // diff_coeff(i,j,k) = 0.00016; // TODO: make this not a hack

        // Set initial temperature field
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 88.53434;
        // // TODO: better way to initialize tcouple geometry
        // Real r = 0.050;
        // if (x*x+y*y+(z+r+0.1)*(z+r+0.1)<=r*r) {
        //     tracer(i,j,k) = 239.16524;
        //     // diff_coeff(i,j,k) = 0.05547074;
        // }

        // // HACK: set the fluid above the thermocouple to be 1 m/s
        // if (((x*x+y*y<=r*r) && (z>-r-0.1)) && (x*x+y*y+(z+r+0.1)*(z+r+0.1)>r*r))
        //     vel(i,j,k,2) = -1.;
    });

#endif
}

// Initialize fluid field for 3D thermo couple
// m_probtype = 2978, and the if-statement check is defined in
// prob_init_fluid() in prob_init_fluid.cpp
void incflo::init_tcouple_energy ( Box const& vbx, Box const& /*gbx*/,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& density,
                                   Array4<Real> const& tracer,
                                   Array4<Real> const& energy,
                                   Array4<Real> const& temp,
                                   Array4<Real> const& vol_heat_cap,
                                   Array4<Real> const& ther_conduct,
                                   Array4<Real> const& one,
                                   Box const& /*domain*/,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);

        // Set initial velocity and density fields
        // TODO: double check if need incflo.ro_0 here
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        density(i,j,k) = 1.;

        // Set initial energy field
        // TODO: also need to set vhc and k here (check how to do +2 tracer fields)
        // HACK: pre-set the plunged solid location (not moving already immersed)
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 90.3;
        Real r = 0.025;
        if (x*x+(y+r+0.1)*(y+r+0.1)<=r*r)
            tracer(i,j,k) = 239.16524;
    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
        Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

        // Set initial velocity and density fields
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;
        density(i,j,k) = 1.;
        one(i,j,k) = 1.;

        // TODO: Parse physical properties in input file
        Real cv_eth = 1.55;
        Real k_eth = 0.0003875;
        Real cv_tcp = 5.4406;
        Real k_tcp = 0.0340625;
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        Real temp_eth = 88.53434; // final measured temperature 1 m/s
        // Real temp_tcp = 239.16524;
        Real temp_tcp = 219.42506;

        // Set initial energy and temperature field (no solid)
        // Energy = c_v * T
        // energy(i,j,k) = temp_eth;
        energy(i,j,k) = cv_eth * temp_eth;
        temp(i,j,k) = temp_eth;
        tracer(i,j,k) = temp_eth;
        // tracer(i,j,k) = cv_eth * temp_eth;
        vol_heat_cap(i,j,k) = cv_eth;
        ther_conduct(i,j,k) = k_eth;

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

        // if (x*x+y*y+(z+1.5*r+d+dx[2])*(z+1.5*r+d+dx[2])<r*r) {
        // if (x*x+y*y+(z+1.5*r+dx[2])*(z+1.5*r+dx[2])<r*r) {
        if (x*x+y*y+(z-r+dx[2]+d)*(z-r+dx[2]+d)<r*r) {
            vel(i,j,k,0) = 0.;
            vel(i,j,k,1) = 0.;
            vel(i,j,k,2) = v;
            if (z>-2.*dx[2]) {
                energy(i,j,k) = temp_tcp * cv_tcp;
                temp(i,j,k) = temp_tcp;
                tracer(i,j,k) = temp_tcp;// * cv_tcp;
                vol_heat_cap(i,j,k) = cv_tcp;
                ther_conduct(i,j,k) = k_tcp;
            }
        }
    });

#endif
}

// Initialize fluid field for 3D thermo couple
// m_probtype = 2979, and the if-statement check is defined in
// prob_init_fluid() in prob_init_fluid.cpp
void incflo::init_disk_energy ( Box const& vbx, Box const& /*gbx*/,
                                Array4<Real> const& vel,
                                Array4<Real> const& density,
                                Array4<Real> const& tracer,
                                Array4<Real> const& energy,
                                Array4<Real> const& temp,
                                Array4<Real> const& vol_heat_cap,
                                Array4<Real> const& ther_conduct,
                                Array4<Real> const& one,
                                Box const& /*domain*/,
                                GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);

        // Set initial velocity and density fields
        // TODO: double check if need incflo.ro_0 here
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        density(i,j,k) = 1.;

        // Set initial energy field
        // TODO: also need to set vhc and k here (check how to do +2 tracer fields)
        // HACK: pre-set the plunged solid location (not moving already immersed)
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 90.3;
        Real r = 0.025;
        if (x*x+(y+r+0.1)*(y+r+0.1)<=r*r)
            tracer(i,j,k) = 239.16524;
    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
        Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

        // Set initial velocity and density fields
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;
        density(i,j,k) = 1.;
        one(i,j,k) = 1.;

        // TODO: Parse physical properties in input file
        Real cv_eth = 1.55;
        Real cv_tcp = 3.76953125; //plunger
        Real k_eth = 0.0003875;
        Real k_tcp = 0.0005390625; // plunger
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        Real temp_eth = 88.53434; // final measured temperature 1 m/s
        // Real temp_tcp = 239.16524;
        Real temp_tcp = 219.42506;


        // Set initial energy and temperature field (no solid)
        // Energy = c_v * T
        // energy(i,j,k) = temp_eth;
        // Ethane
        energy(i,j,k) = cv_eth * temp_eth;
        temp(i,j,k) = temp_eth;
        tracer(i,j,k) = temp_eth;
        // tracer(i,j,k) = cv_eth * temp_eth;
        vol_heat_cap(i,j,k) = cv_eth;
        ther_conduct(i,j,k) = k_eth;

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

        if ((x*x+(z-rz+2*dx[2]+d)*(z-rz+2*dx[2]+d)<rz*rz) && (y>-ry && y<ry)) {
            vel(i,j,k,0) = 0.;
            vel(i,j,k,1) = 0.;
            vel(i,j,k,2) = v;
            if (z>-2.*dx[2]) {
                // Disk
                energy(i,j,k) = cv_tcp * temp_tcp;
                temp(i,j,k) = temp_tcp;
                tracer(i,j,k) = temp_tcp;
                vol_heat_cap(i,j,k) = cv_tcp;
                ther_conduct(i,j,k) = k_tcp;
            }
        }
    });

#endif
}

// Initialize fluid field for 3D thermo couple
// m_probtype = 2980, and the if-statement check is defined in
// prob_init_fluid() in prob_init_fluid.cpp
void incflo::init_grid_energy ( Box const& vbx, Box const& /*gbx*/,
                                Array4<Real> const& vel,
                                Array4<Real> const& density,
                                Array4<Real> const& tracer,
                                Array4<Real> const& energy,
                                Array4<Real> const& temp,
                                Array4<Real> const& vol_heat_cap,
                                Array4<Real> const& ther_conduct,
                                Array4<Real> const& one,
                                Box const& /*domain*/,
                                GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);

        // Set initial velocity and density fields
        // TODO: double check if need incflo.ro_0 here
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        density(i,j,k) = 1.;

        // Set initial energy field
        // TODO: also need to set vhc and k here (check how to do +2 tracer fields)
        // HACK: pre-set the plunged solid location (not moving already immersed)
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        tracer(i,j,k) = 90.3;
        Real r = 0.025;
        if (x*x+(y+r+0.1)*(y+r+0.1)<=r*r)
            tracer(i,j,k) = 239.16524;
    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0] - 0.5*(probhi[0]-problo[0]);
        Real y = (j+0.5)*dx[1] - 0.5*(probhi[1]-problo[1]);
        Real z = (k+0.5)*dx[2] - (probhi[2]-problo[2]);

        // Set initial velocity and density fields
        vel(i,j,k,0) = 0.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;
        density(i,j,k) = 1.;
        one(i,j,k) = 1.;

        // TODO: Parse physical properties in input file
        Real cv_eth = 1.55;
        Real cv_tcp = 3.76953125; //plunger
        Real cv_spl = 6.5421875; // sample
        Real k_eth = 0.0003875;
        Real k_tcp = 0.0005390625; // plunger
        Real k_spl = 0.0009344219; // sample
        // Liquid ethane: 90.3K, plunged object: 293K (to be determined by the equation)
        Real temp_eth = 88.53434; // final measured temperature 1 m/s
        // Real temp_tcp = 239.16524;
        Real temp_tcp = 219.42506;
        Real temp_spl = temp_tcp;

        // Set initial energy and temperature field (no solid)
        // Energy = c_v * T
        // energy(i,j,k) = temp_eth;
        // Ethane
        energy(i,j,k) = cv_eth * temp_eth;
        temp(i,j,k) = temp_eth;
        tracer(i,j,k) = temp_eth;
        // tracer(i,j,k) = cv_eth * temp_eth;
        vol_heat_cap(i,j,k) = cv_eth;
        ther_conduct(i,j,k) = k_eth;

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
                    // Disk
                    energy(i,j,k) = cv_tcp * temp_tcp;
                    temp(i,j,k) = temp_tcp;
                    tracer(i,j,k) = temp_tcp;
                    vol_heat_cap(i,j,k) = cv_tcp;
                    ther_conduct(i,j,k) = k_tcp;
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
                // sample
                energy(i,j,k) = cv_spl * temp_spl;
                temp(i,j,k) = temp_spl;
                tracer(i,j,k) = temp_spl;
                vol_heat_cap(i,j,k) = cv_spl;
                ther_conduct(i,j,k) = k_spl;
            }
        }
    });

#endif
}
