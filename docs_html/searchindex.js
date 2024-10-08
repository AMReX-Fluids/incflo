Search.setIndex({"alltitles": {"Algorithm Options": [[5, null]], "Building AMReX": [[4, "building-amrex"]], "Building incflo": [[4, "building-incflo"]], "Building incflo for Summit (OLCF)": [[4, "building-incflo-for-summit-olcf"]], "Building with CMake": [[4, "building-with-cmake"]], "Building with GNU Make": [[4, "building-with-gnu-make"]], "Checkpoint/Restart": [[6, null]], "Contents:": [[18, null]], "Debugging": [[0, null]], "Downloading the code": [[4, "downloading-the-code"]], "Embedded Boundaries": [[1, null]], "Few more notes on building incflo": [[4, "few-more-notes-on-building-incflo"]], "Fluid Equations": [[2, "fluid-equations"]], "Fluid Variables": [[2, null]], "Getting Started": [[4, null]], "Godunov Methods": [[17, "godunov-methods"]], "Gridding and Load Balancing": [[8, null], [16, null]], "Initialization": [[7, null]], "Introduction": [[15, null]], "MOL": [[17, "mol"]], "Multigrid Inputs": [[9, null]], "Plotfiles and Other Output": [[10, null]], "Problem Definition": [[11, null]], "Run-time Inputs": [[14, null]], "Running the code": [[4, "running-the-code"]], "STANDALONE instructions": [[4, "standalone-instructions"]], "SUPERBUILD Instructions (recommended)": [[4, "superbuild-instructions-recommended"]], "Setting the Time Step": [[12, "setting-the-time-step"]], "Solving the Fluid Equations": [[3, null]], "Time Step": [[17, null]], "Time Stepping": [[12, null]], "Verbosity": [[13, null]], "Visualizing the Results": [[4, "visualizing-the-results"]], "Welcome to incflo\u2019s documentation!": [[18, null]]}, "docnames": ["Debugging", "EB", "FluidEquations", "Fluids_Chapter", "GettingStarted", "InputsAlgorithm", "InputsCheckpoint", "InputsInitialization", "InputsLoadBalancing", "InputsMultigrid", "InputsPlotFiles", "InputsProblemDefinition", "InputsTimeStepping", "InputsVerbosity", "Inputs_Chapter", "Introduction", "ManagingGridHierarchy_Chapter", "TimeStep", "index"], "envversion": {"sphinx": 63, "sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx.ext.viewcode": 1}, "filenames": ["Debugging.rst", "EB.rst", "FluidEquations.rst", "Fluids_Chapter.rst", "GettingStarted.rst", "InputsAlgorithm.rst", "InputsCheckpoint.rst", "InputsInitialization.rst", "InputsLoadBalancing.rst", "InputsMultigrid.rst", "InputsPlotFiles.rst", "InputsProblemDefinition.rst", "InputsTimeStepping.rst", "InputsVerbosity.rst", "Inputs_Chapter.rst", "Introduction.rst", "ManagingGridHierarchy_Chapter.rst", "TimeStep.rst", "index.rst"], "indexentries": {}, "objects": {}, "objnames": {}, "objtypes": {}, "terms": {"": [0, 1, 2, 4, 5, 9, 11, 16], "0": [2, 4, 5, 6, 8, 9, 10, 11, 12, 13], "1": [0, 4, 5, 6, 8, 9, 10, 11, 12, 17], "10": 4, "100": 9, "1024000": 8, "105": 4, "11": 9, "14": [4, 9], "2": [4, 5, 15, 17], "200": 9, "3": [4, 7, 15], "32": 8, "3d": 8, "4": [0, 9], "5": 12, "7": [9, 12], "8": 8, "81": 11, "9": 11, "As": 1, "By": [1, 2], "For": [0, 4, 11, 17, 18], "If": [0, 2, 4, 5, 6, 7, 9, 10, 12], "In": [0, 4, 5, 11, 12, 17], "It": [1, 11, 15, 18], "Of": 15, "One": 0, "The": [0, 1, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15], "Then": 4, "There": [4, 16], "These": 4, "To": [1, 4, 9, 11], "_u": [2, 17], "about": [1, 4, 16], "abov": 4, "absolut": [0, 4, 9], "acceler": 2, "accuraci": 17, "achiev": 17, "activ": 18, "actual": 4, "adapt": [15, 18], "add": 4, "addit": [4, 8, 9, 11], "addition": 11, "advect": [5, 11, 12, 15, 17], "advect_momentum": 5, "advect_trac": 5, "advection_typ": 5, "after": 10, "against": 4, "algorithm": [14, 15], "all": [0, 4, 10], "allow": [9, 11], "alon": 4, "also": [0, 1, 4, 11, 15, 16], "altern": [1, 4], "alwai": [4, 12], "amr": [6, 7, 8, 10, 11, 15, 18], "amrex": [0, 1, 8, 9, 11, 15, 16, 17, 18], "amrex_hom": 4, "amrex_hydro_hom": 4, "amrex_root": 4, "amrvi": [0, 15], "an": [0, 4, 15, 17, 18], "ani": 4, "anoth": [0, 15], "answer": 1, "applic": [4, 15], "approach": [1, 15, 17], "appropri": 12, "approxim": [5, 10, 17], "ar": [0, 4, 8, 9, 11, 12, 15, 16, 17], "area": 11, "argument": 4, "art": 0, "assign": 16, "assum": 4, "ast": 17, "atol": 9, "avail": [4, 11, 15], "awar": 4, "backtrac": 0, "balanc": [14, 18], "base": [15, 18], "bash": 4, "bashrc": 4, "basic": [11, 17], "bc": 11, "bd": [5, 17], "befor": [4, 6, 7], "behavior": 0, "believ": 0, "below": [4, 9, 12], "benchmark": 4, "between": [0, 12], "bf": [2, 17], "bicg": 9, "bicgcg": 9, "bicgstab": 9, "binari": 4, "block": 15, "blocking_factor_i": 8, "blocking_factor_x": 8, "blocking_factor_z": 8, "bool": [5, 6, 7, 10, 12], "both": [4, 12, 15], "bottom": 9, "bottom_maxit": 9, "bottom_solv": 9, "bottom_verbos": 9, "boundari": [4, 11, 15, 18], "boxarrai": 16, "branch": 18, "brief": 4, "bug": 0, "build": [1, 15], "builddir": 4, "buildingamrex": 4, "built": [0, 1, 4, 15, 16], "c": 4, "calcul": 7, "call": [4, 9], "can": [0, 1, 4, 11, 12, 15, 16], "cartesian": [11, 15], "case": [0, 12], "cc": 4, "cd": 4, "cdot": [2, 5, 17], "cell": [0, 4, 5, 8, 11, 15], "center": [5, 15, 17], "centroid": 15, "cfl": 12, "cg": 9, "cgbicg": 9, "chang": 5, "chap": 8, "characterist": 4, "check": [0, 4], "check_fil": 6, "check_int": 6, "checkpoint": [4, 14], "chk": 6, "chk00000": 4, "chk00010": 4, "choos": [1, 4], "chosen": 4, "ci": 4, "clone": 4, "cluster": 4, "cmake": 0, "cmake_build_typ": 4, "cmake_cuda_flag": 4, "cmake_cxx_flag": 4, "co": 15, "coars": 1, "coarsen": 1, "coarser": 9, "code": [11, 12, 15, 18], "coexist": 4, "com": [4, 15], "come": 16, "comm": 4, "comm_profil": 4, "command": 4, "comment": 11, "common": 4, "commonli": [8, 9], "comp": 4, "compar": 0, "compil": [0, 4], "complex": [15, 18], "compon": 5, "comput": [12, 16], "condit": [11, 12], "configur": 4, "conserv": [2, 5], "consist": 4, "const": 0, "constant": [5, 11], "constant_dens": 5, "constraint": [2, 12], "construct": [1, 15, 17], "contain": [0, 4, 11, 16], "continu": [4, 5], "control": [6, 9, 10, 16], "convect": [5, 15, 17], "coord_si": 11, "coordin": [11, 15], "corner": 11, "corrector": 17, "correspond": 4, "cpp": 11, "cpu": [8, 15], "crank": 5, "creat": [0, 1, 4, 6, 8, 11, 16], "creation": 16, "criteria": 12, "criterion": 12, "csg": 1, "csgeb_hom": 1, "cuda": [4, 15], "current": [4, 9, 11, 17], "custom": 4, "cut": [4, 15], "cxxflag": 4, "d": [4, 15], "damrex_hom": 4, "damrex_hydro_hom": 4, "damrex_root": 4, "darshan": 4, "data": [0, 4], "databas": 1, "dcmake_build_typ": 4, "dcmake_cxx_compil": 4, "dcmake_install_prefix": 4, "dcp": 15, "ddebug": 0, "debug": [4, 18], "decid": 4, "decompos": 16, "default": [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17], "defin": [1, 4, 11, 15, 16, 17], "definit": [2, 4, 14], "delp": 11, "delta": 17, "densiti": [2, 5, 10, 11, 15, 17, 18], "depend": 1, "deriv": 4, "describ": [1, 3, 4], "descript": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "design": 15, "desir": 4, "detail": [1, 11, 16, 17], "determin": [1, 7, 8, 12], "develop": 18, "diff": 0, "differ": [0, 1, 4, 12, 16], "differenc": 5, "diffus": [2, 5, 9, 11, 17], "diffusion_typ": 5, "diffvar": 0, "dimension": 4, "dincflo_cuda": 4, "dirchlet": 11, "direct": [8, 11], "directori": 4, "discret": [3, 15], "discuss": [1, 4, 17], "disk": 0, "distribut": [4, 16], "distributionmap": 16, "diverg": 10, "divis": 8, "do": [4, 7, 12], "do_initial_proj": 7, "docmument": 16, "docs_html": 4, "document": [0, 1, 4, 8, 9, 11, 17], "doe": [4, 10, 11, 12], "domain": [11, 16], "down": 0, "dp": 10, "drop": 11, "dt": [10, 12], "dt_change_max": 12, "dure": 2, "duse_xsdk_default": 4, "dx": [10, 12], "dy": 10, "dz": 10, "e": [0, 4, 9, 11, 16, 17], "each": [0, 4, 5, 8, 11, 16, 17], "earlier": 4, "eb": [1, 4, 10, 11], "eb2": 1, "ecp": 15, "edit": 4, "either": [4, 12, 17], "embed": [4, 15, 18], "embedded_boundari": 1, "enabl": 4, "encod": 4, "encount": 0, "encourag": 4, "end": 15, "ensur": 4, "entiti": 4, "entri": 5, "environ": 4, "environment": 4, "equat": [5, 15, 18], "error": 4, "etc": 4, "everyon": 0, "evolut": 3, "evolv": 5, "ex": [0, 4], "exact": 10, "examin": 0, "exampl": [0, 4, 12], "except": [0, 4], "execut": [0, 4], "exist": 4, "expand": 4, "explan": 11, "explicit": [5, 15], "explicitli": 4, "export": 4, "extern": 2, "fabarrai": [0, 16], "fabarray_mfit": 8, "face": [11, 15, 17], "factor": 12, "fals": [4, 5, 6, 10, 11, 12], "farraybox": 0, "favorit": 0, "fcompar": 0, "featur": 15, "few": 0, "field": 4, "file": [0, 1, 4, 6, 7, 10, 11, 12], "first": [4, 7, 12], "fix": 12, "fixed_dt": 12, "flag": [0, 1, 4, 10], "flexibl": 16, "float": [0, 4, 9], "flow": [12, 15, 18], "fluid": [1, 4, 12, 15, 18], "folder": 4, "follow": [1, 4, 5, 6, 7, 8, 10, 11, 13], "forc": [2, 10], "format": [1, 10, 15], "found": [0, 4, 15], "fpe_trap_invalid": 0, "frac": [2, 17], "fraction": 10, "framework": [1, 4, 15], "frequenc": [6, 10], "from": [4, 6, 7], "fu": 16, "full": 4, "function": 11, "further": 0, "g": [0, 2, 4, 11], "gcc": 4, "gener": [4, 10], "geom": 1, "geometr": 1, "geometri": [1, 10, 11, 15, 18], "get": 18, "git": 4, "github": [0, 4, 15], "give": 1, "gmake": 0, "gnumakefil": 4, "go": 0, "godunov": [3, 5, 15], "gov": 4, "govern": 3, "gpu": [4, 15], "gravit": 2, "graviti": [2, 11], "green": 4, "grid": [6, 14, 18], "gshop": 1, "guid": 4, "h": [2, 17], "h_": 2, "h_x": 2, "h_y": 2, "h_z": 2, "ha": [0, 17], "have": [0, 4, 12], "hdf5": 15, "heavili": 15, "help": 0, "here": [0, 3, 4, 8, 17], "high": 11, "higher": 4, "hip": 15, "hit": 12, "how": [1, 7, 8, 16, 17], "html": 4, "http": [4, 15], "hybrid": 15, "hydro": [4, 9, 11, 15, 17], "hypr": [4, 9], "i": [0, 1, 2, 4, 5, 9, 11, 12, 15, 16, 17, 18], "iamr": [15, 18], "ideal": 4, "immedi": 1, "implement": 11, "implicit": [5, 15, 17], "impos": 15, "incflo": [0, 1, 5, 7, 11, 12, 13, 15, 16, 17], "incflo2d": 4, "incflo_cuda": 4, "incflo_dim": 4, "incflo_eb": 4, "incflo_fp": 4, "incflo_hypr": 4, "incflo_mpi": 4, "incflo_omp": 4, "includ": [0, 4, 6, 11, 15, 18], "incompress": [2, 15, 18], "incorrect": 0, "increas": 12, "index": 11, "indic": 4, "indirectli": 1, "individu": 16, "infile1": 0, "infile2": 0, "inflow": 11, "info": [4, 11], "inform": 4, "inher": 4, "init_shrink": 12, "initi": [0, 1, 10, 12, 14], "initial_iter": 7, "input": [0, 1, 4, 5, 6, 7, 8, 10, 11, 12, 13, 18], "inputsloadbalanc": 8, "instal": [1, 4], "installdir": 4, "instanc": 11, "instead": 4, "int": [0, 5, 6, 7, 8, 9, 10, 11, 12, 13], "integr": 4, "intel": 4, "intermedi": 15, "introduct": [4, 18], "intvect": [0, 8], "invok": 4, "involv": 17, "io": [4, 10, 15], "is_period": 11, "issu": 0, "iter": [7, 9, 12], "its": 4, "j": 4, "kei": 15, "later": 4, "leak": 0, "left": 17, "let": [4, 12], "lev": 1, "level": [1, 4, 6, 8, 9, 11, 16], "leverag": [15, 16], "librari": [1, 4], "like": 4, "line": [4, 15, 17], "link": 4, "list": [4, 8, 9], "load": [4, 14, 18], "locat": 4, "log": 0, "logic": 8, "longer": 17, "look": [0, 4], "low": [11, 15], "mac": [9, 11, 15, 17], "mac_proj": 9, "mach": 15, "machin": [4, 15, 16], "made": 4, "mai": 4, "make": 0, "makeebgeometri": 1, "makefil": 4, "mani": [0, 1, 4, 7], "mass": [2, 11], "mass_inflow": 11, "massiv": 15, "max": 12, "max_crse_level": 1, "max_grid_size_i": 8, "max_grid_size_x": 8, "max_grid_size_z": 8, "max_level": 11, "max_step": 12, "maximum": [0, 1, 8, 9, 11, 12], "maxit": 9, "mdkir": 4, "meet": 4, "mesh": [15, 18], "method": [0, 3, 4, 15], "methodologi": 17, "mf": 0, "mfiter": 16, "mfiter_tile_s": 16, "mfix": 1, "mg_max_coarsening_level": 9, "mi": 11, "might": 1, "minsizerel": 4, "mix": 11, "mkdir": 4, "mlmg": 1, "mode": [0, 4], "model": [15, 18], "modifi": 10, "modul": 4, "mol": [3, 5, 15], "momentum": [2, 5], "more": [0, 1, 16], "most": [8, 9], "mpi": [0, 4, 15, 16], "mpiexec": 0, "mu": 5, "mu_": [2, 5], "multicompon": 5, "multicor": [15, 16], "multifab": [0, 16], "multigrid": [1, 14], "must": [1, 4, 5, 6, 7, 8, 10, 11, 12, 13], "n": [0, 17], "n_cell": 11, "nabla": [2, 5, 17], "name": [0, 4, 6, 10], "namespac": 1, "nan": 0, "nativ": 15, "navier": [15, 18], "need": [4, 11, 17], "neither": 12, "neumann": 11, "new": [0, 6, 11, 17], "next": 4, "nicholson": 5, "no_slip_wal": 11, "nodal": [5, 9, 11], "nodal_proj": 9, "node": 15, "non": [2, 11], "none": [4, 6, 7, 11], "nor": 12, "note": [0, 1, 12], "now": 4, "nsw": 11, "ntrac": 5, "number": [5, 8, 9, 11, 12, 15], "o": [4, 15], "obviou": 1, "occur": [0, 8], "offer": 0, "often": 8, "old": 0, "omit": 4, "onc": 10, "one": [0, 4, 5, 11, 16], "ongo": 18, "onli": [5, 8, 10, 11], "openmp": [4, 15, 16], "openscad": 1, "option": [1, 4, 8, 9, 11, 12, 14, 15, 17], "order": [4, 17], "origin": 0, "ornl": 4, "other": [4, 14], "otherwis": 5, "out": [0, 10], "outflow": 11, "output": [0, 6, 14], "overrid": [4, 10, 12], "overwritten": 4, "own": [0, 1], "p": [0, 2, 17], "pa": 11, "page": 0, "parallel": 15, "paramet": [0, 1, 4], "paraview": 15, "part": [4, 15], "partial": 2, "pass": 4, "path": 4, "per": [5, 16], "percentag": 12, "perform": 4, "period": [4, 10], "phi": 17, "physic": 11, "pi": 11, "platform": 4, "pleas": [0, 18], "plot": 10, "plot_fil": 10, "plot_int": 10, "plot_per_approx": 10, "plot_per_exact": 10, "plotfil": [0, 4, 6, 14, 15], "plotfile_on_restart": [6, 10], "plt": 10, "plt00000": 4, "plt00000_run1": 0, "plt00000_run2": 0, "plt00010": 4, "plt_ccse_regtest": 10, "plt_divu": 10, "plt_eta": 10, "plt_forc": 10, "plt_gpx": 10, "plt_gpy": 10, "plt_gpz": 10, "plt_p": 10, "plt_rho": 10, "plt_strainrat": 10, "plt_tracer": 10, "plt_veli": 10, "plt_velx": 10, "plt_velz": 10, "plt_vfrac": 10, "plt_vort": 10, "po": 11, "point": [0, 4], "possibl": [4, 15], "preced": [5, 6, 7, 8, 9, 10, 11, 12, 13], "pred": 17, "predictor": 17, "prefix": [6, 10, 12], "presenc": 18, "present": 6, "pressur": [5, 7, 10, 11, 15], "pressure_inflow": 11, "pressure_outflow": 11, "print": 0, "print_stat": 0, "prob": 11, "prob_bc": 11, "prob_hi": 11, "prob_lo": 11, "problem": [4, 12, 14], "procedur": 1, "process": [0, 4], "profil": 4, "project": [5, 7, 9, 11, 15, 17], "provid": [1, 4, 11], "publicli": 15, "pull": 4, "put": 4, "r": 15, "rank": 16, "rate": 10, "rather": [5, 7], "reach": 12, "real": [5, 9, 10, 11, 12], "rectangular": 16, "reduc": 12, "refin": [11, 15, 18], "region": 11, "regress": 4, "regrid": 8, "regrid_int": 8, "regrid_on_restart": 6, "rel": [0, 4, 9], "releas": 4, "relwithdebinfo": 4, "report": 0, "repositori": [1, 4], "represent": 15, "requir": [4, 5, 11, 17], "required_crse_lev": 1, "restart": [4, 7, 10, 14], "restrict": 1, "result": 0, "rho": [2, 17], "rho_0": 5, "right": 17, "routin": [0, 13, 15], "rther": 16, "rtol": 9, "run": [0, 16, 18], "runtim": 2, "safeti": 12, "sai": 4, "same": [4, 11], "save": 10, "scalar": [2, 5, 11], "scheme": [5, 17], "scratch": 7, "search": 4, "second": 17, "section": [0, 4], "see": [0, 1, 4, 8, 9, 11, 15, 16, 17, 18], "select": 1, "separ": [4, 5, 11], "set": [1, 2, 4, 7, 9, 11, 15, 16, 17], "setenv": 4, "setup": 11, "sever": [0, 1, 4], "sh": 4, "shell": 4, "should": [7, 10], "signal": 0, "simpli": [4, 12], "simul": [12, 18], "simultan": 4, "singl": [0, 4, 11], "size": 16, "slightli": 1, "small": 1, "smoother": 9, "snippet": 4, "so": [4, 12], "softwar": 15, "solid": 1, "solv": [1, 5, 15, 17, 18], "solver": [1, 9, 11], "sourc": [2, 4, 15], "space": 11, "specifi": [1, 5, 11, 12, 16], "src": 1, "stagger": 17, "stand": 4, "start": [7, 18], "state": [12, 17], "std": 0, "steadi": 12, "steady_st": 12, "step": [3, 4, 6, 8, 14, 15], "still": [1, 12], "stoke": [15, 18], "stop": 12, "stop_tim": 12, "strain": 10, "strategi": 16, "stress": 2, "string": [0, 5, 6, 7, 9, 10, 11], "strongli": 4, "structur": 15, "subcycl": [15, 18], "subfold": 4, "submission_script": 4, "subsequ": 12, "suggest": [0, 4], "summar": 4, "summit_script": 4, "support": [4, 15, 16], "sure": 4, "system": [4, 15], "t": [2, 17], "tabl": 4, "take": [4, 6, 12], "target": 4, "tau": [2, 5, 17], "taylor": 4, "taylor_green_vortic": 4, "tcsh": 4, "tell": 0, "tensor": 2, "term": [10, 15, 17], "test": 4, "test_2d": 4, "test_3d": 4, "test_no_eb": 4, "test_no_eb_2d": 4, "test_xxx": 4, "than": [0, 5, 7, 16], "therefor": 1, "thi": [0, 1, 4, 5, 7, 10, 11, 12], "those": [0, 16], "thread": 16, "three": 12, "through": [4, 15], "thu": 4, "tile": [8, 16], "tile_s": 8, "time": [3, 4, 6, 10, 15, 18], "timestep": 7, "tini": 4, "tiny_profil": 4, "tip": 0, "toler": [9, 12], "too": 1, "tool": 0, "top": [4, 15], "trac_is_conserv": 5, "trace": 4, "trace_profil": 4, "tracer": [2, 5, 10, 11, 15], "track": 0, "treat": [11, 17], "true": [0, 1, 4, 5, 7, 10, 11], "turn": 0, "two": [0, 4, 12], "type": [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "u": [2, 12, 17], "u_g": 0, "unless": 4, "unload": 4, "unsteadi": 12, "up": 4, "updat": 4, "us": [0, 1, 4, 5, 6, 8, 9, 10, 11, 12, 15, 17], "use_cc_proj": 5, "use_csg": 1, "use_cuda": 4, "use_mpi": 4, "use_omp": 4, "use_tensor_solv": [5, 11], "user": [0, 4], "util": [4, 15], "v": 4, "valgrind": 0, "valid": 4, "vallog": 0, "valu": [1, 4, 11, 12], "variabl": [0, 3, 4, 10, 15, 18], "ve": [0, 4], "vector": 11, "vel": 12, "veloc": [2, 5, 10, 11, 15, 17], "verbos": [9, 14], "version": 4, "via": 15, "viscos": [5, 10, 11, 15], "viscou": [2, 15], "visit": 15, "vismf": 0, "void": 0, "volum": 10, "vortic": [4, 10], "vtp": 10, "wa": 1, "wai": 0, "walk": 4, "want": [0, 4, 12], "we": [0, 1, 3, 4, 5, 7, 8, 10, 11, 16, 17], "well": [4, 10], "were": 4, "what": 10, "when": [0, 1, 4, 10, 11, 12, 16, 17], "where": [1, 4, 15], "whether": [4, 10, 12], "which": [0, 1, 5, 9, 11, 15, 16, 17], "within": [1, 4], "without": [4, 15, 18], "work": 16, "workflow": 17, "write": [0, 6, 10], "write_eb_surfac": 10, "written": [6, 10], "www": 4, "x": [8, 10, 15], "xalt": 4, "xhi": 11, "xlo": 11, "xxx": [1, 12], "y": [8, 10], "ye": [0, 4], "yhi": 11, "yield": 4, "ylo": 11, "you": [0, 4, 12], "your": [1, 4], "yt": 15, "z": [8, 10, 15], "zhi": 11, "zlo": 11}, "titles": ["Debugging", "Embedded Boundaries", "Fluid Variables", "Solving the Fluid Equations", "Getting Started", "Algorithm Options", "Checkpoint/Restart", "Initialization", "Gridding and Load Balancing", "Multigrid Inputs", "Plotfiles and Other Output", "Problem Definition", "Time Stepping", "Verbosity", "Run-time Inputs", "Introduction", "Gridding and Load Balancing", "Time Step", "Welcome to incflo\u2019s documentation!"], "titleterms": {"": 18, "algorithm": 5, "amrex": 4, "balanc": [8, 16], "boundari": 1, "build": 4, "checkpoint": 6, "cmake": 4, "code": 4, "content": 18, "debug": 0, "definit": 11, "document": 18, "download": 4, "embed": 1, "equat": [2, 3], "few": 4, "fluid": [2, 3], "get": 4, "gnu": 4, "godunov": 17, "grid": [8, 16], "incflo": [4, 18], "initi": 7, "input": [9, 14], "instruct": 4, "introduct": 15, "load": [8, 16], "make": 4, "method": 17, "mol": 17, "more": 4, "multigrid": 9, "note": 4, "olcf": 4, "option": 5, "other": 10, "output": 10, "plotfil": 10, "problem": 11, "recommend": 4, "restart": 6, "result": 4, "run": [4, 14], "set": 12, "solv": 3, "standalon": 4, "start": 4, "step": [12, 17], "summit": 4, "superbuild": 4, "time": [12, 14, 17], "variabl": 2, "verbos": 13, "visual": 4, "welcom": 18}})