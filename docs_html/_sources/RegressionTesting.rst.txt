Regression testing
==================

Generating local regression test data
-------------------------------------

Developers are encouraged to create local benchmark data for regression
testing prior to committing code the GitLab repository.

1. Create a location to store benchmark data and clone the incflo and
   AMReX repositories.

   .. code:: shell

       mkdir /home/user/incflo-rt
       mkdir /home/user/incflo-rt/benchmark
       cd /home/user/incflo-rt
       git clone https://github.com/AMReX-Codes/amrex.git
       git clone https://github.com/AMReX-Codes/incflo.git

2. Create a local copy the regression test setup file from the incflo
   repository.

   .. code:: shell

       cp incflo/RegressionTesting/incflo-tests.ini incflo-local.ini

3. Edit the local setup file. Specify the top level directories for test
   and web output under the ``[main]`` heading.

   .. code:: shell

       [main]
       testTopDir = /home/user/incflo-rt/benchmark
       webTopDir =  /home/user/incflo-rt/web

   Specify the AMReX source directory location under the ``[AMReX]``
   heading.

   .. code:: shell

       [AMReX]
       dir = /home/user/incflo-rt/amrex
       branch = development

   Specify the incflo-Exa source directory location under the ``[source]``
   heading.

   .. code:: shell

       [source]
       dir = /home/user/incflo-rt/incflo
       branch = develop

4. Run the AMReX regression test tool. The second argument is a user
   supplied comment.

   .. code:: shell

       cd /home/user/incflo-rt
       amrex/Tools/RegressionTesting/regtest.py --make_benchmarks "incflo" incflo-local.ini

+------------------------------------------------------------------------------------------------------------------------+
| ## Prerequisite: Environment Dependencies on Joule (Joule specific)                                                    |
+------------------------------------------------------------------------------------------------------------------------+
| For the Joule environment, load the gnu module and set environment variables first. If not on Joule, skip this step.   |
+------------------------------------------------------------------------------------------------------------------------+
| \`\`\`shell                                                                                                            |
+------------------------------------------------------------------------------------------------------------------------+
| > module load gnu/6.1.0                                                                                                |
+------------------------------------------------------------------------------------------------------------------------+
| > export CC=/nfs/apps/Compilers/GNU/6.1.0/bin/gcc                                                                      |
+------------------------------------------------------------------------------------------------------------------------+
| > export CXX=/nfs/apps/Compilers/GNU/6.1.0/bin/g++                                                                     |
+------------------------------------------------------------------------------------------------------------------------+
| > export F77=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran                                                                |
+------------------------------------------------------------------------------------------------------------------------+
| > export FC=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran                                                                 |
+------------------------------------------------------------------------------------------------------------------------+
| \`\`\`                                                                                                                 |
+------------------------------------------------------------------------------------------------------------------------+
