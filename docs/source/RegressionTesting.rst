Regression testing
==================

Generating local regression test data
-------------------------------------

Developers are encouraged to create local benchmark data for regression
testing prior to committing code the GitLab repository.

1. Create a location to store benchmark data and clone the MFiX and
   AMReX repositories.

   .. code:: shell

       mkdir /home/user/exa-rt
       mkdir /home/user/exa-rt/benchmark
       cd /home/user/exa-rt
       git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
       git clone https://github.com/AMReX-Codes/amrex.git

2. Create a local copy the regression test setup file from the MFiX
   repository.

   .. code:: shell

       cp mfix/RegressionTesting/MFiX-tests.ini MFiX-local.ini

3. Edit the local setup file. Specify the top level directories for test
   and web output under the ``[main]`` heading.

   .. code:: shell

       [main]
       testTopDir = /home/user/exa-rt/benchmark
       webTopDir =  /home/user/exa-rt/web

   Specify the AMReX source directory location under the ``[AMReX]``
   heading.

   .. code:: shell

       [AMReX]
       dir = /home/user/exa-rt/amrex
       branch = development

   Specify the MFiX-Exa source directory location under the ``[source]``
   heading.

   .. code:: shell

       [source]
       dir = /home/user/exa-rt/mfix
       branch = develop

4. Run the AMReX regression test tool. The second argument is a user
   supplied comment.

   .. code:: shell

       cd /home/user/exa-rt
       amrex/Tools/RegressionTesting/regtest.py --make_benchmarks "MFiX" MFiX-local.ini

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
