#include <incflo.H>

using namespace amrex;

void incflo::ReadRheologyParameters()
{
     amrex::ParmParse pp0("incflo");
     pp0.query("do_vof", m_do_vof);
     
     // Initialize Rheology Parameters for Single Fluid
     if (!m_do_vof) 
     {
         amrex::ParmParse pp("incflo");

         std::string fluid_model_s = "newtonian";
         pp.query("fluid_model", fluid_model_s);
         if(fluid_model_s == "newtonian")
         {
             m_fluid_model = FluidModel::Newtonian;
             amrex::Print() << "Newtonian fluid with"
                            << " mu = " << m_mu << std::endl;
         }
         else if(fluid_model_s == "powerlaw")
         {
             m_fluid_model = FluidModel::powerlaw;
             pp.query("n", m_n_0);
             AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                     "No point in using power-law rheology with n = 1");

             amrex::Print() << "Power-law fluid with"
                            << " mu = " << m_mu
                            << ", n = " << m_n_0 <<  std::endl;
         }
         else if(fluid_model_s == "bingham")
         {
             m_fluid_model = FluidModel::Bingham;
             pp.query("tau_0", m_tau_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                     "No point in using Bingham rheology with tau_0 = 0");

             pp.query("papa_reg", m_papa_reg);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                        "Papanastasiou regularisation parameter must be positive");

             amrex::Print() << "Bingham fluid with"
                            << " mu = " << m_mu
                            << ", tau_0 = " << m_tau_0
                            << ", papa_reg = " << m_papa_reg << std::endl;
         }
         else if(fluid_model_s == "hb")
         {
             m_fluid_model = FluidModel::HerschelBulkley;
             pp.query("n", m_n_0);
             AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                     "No point in using Herschel-Bulkley rheology with n = 1");

             pp.query("tau_0", m_tau_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_0 = 0");

             pp.query("papa_reg", m_papa_reg);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             amrex::Print() << "Herschel-Bulkley fluid with"
                            << " mu = " << m_mu
                            << ", n = " << m_n_0
                            << ", tau_0 = " << m_tau_0
                            << ", papa_reg = " << m_papa_reg << std::endl;
         }
         else if(fluid_model_s == "smd")
         {
             m_fluid_model = FluidModel::deSouzaMendesDutra;
             pp.query("n", m_n_0);
             AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);

             pp.query("tau_0", m_tau_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");

             pp.query("eta_0", m_eta_0);
             AMREX_ALWAYS_ASSERT(m_eta_0 > 0.0);

             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << m_mu
                            << ", n = " << m_n_0
                            << ", tau_0 = " << m_tau_0
                            << ", eta_0 = " << m_eta_0 << std::endl;
         }
         else
         {
             amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd");
         }
     }

     // Initialize Rheology Parameters for VOF
     if (m_do_vof) 
     {
         FluidVOF_t fluid0;
         FluidVOF_t fluid1;
     
         amrex::ParmParse ppVOF("incflo.vof");

         amrex::Vector<std::string> fluid_model{{"newtonian", "newtonian"}};;
         ppVOF.queryarr("fluid_model", fluid_model);
         amrex::Vector<amrex::Real> rho;
         ppVOF.queryarr("rho", rho);
         amrex::Vector<amrex::Real> mu;
         ppVOF.queryarr("mu", mu);
         amrex::Vector<amrex::Real> n_0;
         ppVOF.queryarr("n_0", n_0);
         amrex::Vector<amrex::Real> tau_0;
         ppVOF.queryarr("tau_0", tau_0);
         amrex::Vector<amrex::Real> papa_reg;
         ppVOF.queryarr("papa_reg", papa_reg);
         amrex::Vector<amrex::Real> eta_0;
         ppVOF.queryarr("eta_0", eta_0);

         if (fluid_model.size() != 2) {
            amrex::Abort("Need 2 incflo.vof.fluid_model");
         }

         if (rho.size() != 2) {
            amrex::Abort("Need 2 incflo.vof.rho");
         }
         else {
            fluid0.rho = rho[0];
            fluid1.rho = rho[1];
         }

         // fluid 0
         if(fluid_model[0] == "newtonian")
         {
             fluid0.fluid_model = FluidModel::Newtonian;
             fluid0.mu = mu[0];
             amrex::Print() << "Newtonian fluid0 with"
                            << " mu = " << fluid0.mu << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "powerlaw")
         {
             fluid0.fluid_model = FluidModel::powerlaw;
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[0] != 1.0,
                     "No point in using power-law rheology with n = 1");
             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             amrex::Print() << "Power-law fluid0 with"
                            << " mu = " << fluid0.mu
                            << ", n = " << fluid0.n_0 << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "bingham")
         {
             fluid0.fluid_model = FluidModel::Bingham;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using Bingham rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[0] > 0.0,
                        "Papanastasiou regularisation parameter must be positive");
             
             fluid0.mu = mu[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.papa_reg = papa_reg[0];
             amrex::Print() << "Bingham fluid0 with"
                            << " mu = " << fluid0.mu
                            << ", tau_0 = " << fluid0.tau_0
                            << ", papa_reg = " << fluid0.papa_reg << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "hb")
         {
             fluid0.fluid_model = FluidModel::HerschelBulkley;
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[0] != 1.0,
                     "No point in using Herschel-Bulkley rheology with n = 1");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.papa_reg = papa_reg[0];
             amrex::Print() << "Herschel-Bulkley fluid0 with"
                            << " mu = " << fluid0.mu
                            << ", n = " << fluid0.n_0
                            << ", tau_0 = " << fluid0.tau_0
                            << ", papa_reg = " << fluid0.papa_reg << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "smd")
         {
             fluid0.fluid_model = FluidModel::deSouzaMendesDutra;
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT(eta_0[0] > 0.0);

             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.eta_0 = eta_0[0];
             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << fluid0.mu
                            << ", n = " << fluid0.n_0
                            << ", tau_0 = " << fluid0.tau_0
                            << ", eta_0 = " << fluid0.eta_0 << " and rho = " << fluid0.rho << std::endl;
         }
         else
         {
             amrex::Abort("Unknown fluid_model for fluid0! Choose either newtonian, powerlaw, bingham, hb, smd");
         }
         m_fluid_vof.push_back(fluid0);

         //  fluid 1
         if(fluid_model[1] == "newtonian")
         {
             fluid1.fluid_model = FluidModel::Newtonian;
             fluid1.mu = mu[1];
             amrex::Print() << "Newtonian fluid1 with"
                            << " mu = " << fluid1.mu << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "powerlaw")
         {
             fluid1.fluid_model = FluidModel::powerlaw;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[1] != 1.0,
                     "No point in using power-law rheology with n = 1");
             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             amrex::Print() << "Power-law fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0 << " and rho = " << fluid1.rho <<  std::endl;
         }
         else if(fluid_model[1] == "bingham")
         {
             fluid1.fluid_model = FluidModel::Bingham;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using Bingham rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[1] > 0.0,
                        "Papanastasiou regularisation parameter must be positive");
             
             fluid1.mu = mu[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.papa_reg = papa_reg[1];
             amrex::Print() << "Bingham fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", tau_0 = " << fluid1.tau_0
                            << ", papa_reg = " << fluid1.papa_reg << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "hb")
         {
             fluid1.fluid_model = FluidModel::HerschelBulkley;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[1] != 1.0,
                     "No point in using Herschel-Bulkley rheology with n = 1");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.papa_reg = papa_reg[1];
             amrex::Print() << "Herschel-Bulkley fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0
                            << ", tau_0 = " << fluid1.tau_0
                            << ", papa_reg = " << fluid1.papa_reg << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "smd")
         {
             fluid1.fluid_model = FluidModel::deSouzaMendesDutra;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT(eta_0[1] > 0.0);

             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.eta_0 = eta_0[1];
             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0
                            << ", tau_0 = " << fluid1.tau_0
                            << ", eta_0 = " << fluid1.eta_0 << " and rho = " << fluid1.rho << std::endl;
         }
         else
         {
             amrex::Abort("Unknown fluid_model for fluid1! Choose either newtonian, powerlaw, bingham, hb, smd");
         }
         m_fluid_vof.push_back(fluid1);

     }

}
