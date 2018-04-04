# define additional build options
include(CMakeDependentOption)

# enable/disable FlexibleSUSY components
option(ENABLE_COLORS "Colored output (default: no)" OFF)
option(ENABLE_COMPILE "Compile the source code (default: yes)" ON)
option(ENABLE_COMPILER_WARNINGS "Enable compiler warnings (default: no)" OFF)
option(ENABLE_DEBUG "Enable debug output (default: no)" OFF)
option(ENABLE_GM2Calc "Enable use of the GM2Calc addon (default: yes)" ON)
# option(ENABLE_ILP64MKL_WORKAROUND "Enable use ILP64MKL workaround (default: yes)" ON)
# option(ENABLE_LIBRARYLINK "Enable use of the LibaryLink (default: yes)" ON)
option(ENABLE_MASS_ERROR_CHECK "Check mass eigenvalue precision (default: no)" OFF)
option(ENABLE_META "Create model class (default: yes)" ON)
option(ENABLE_ODEINT "Enable use of Boost's odeint function (default: yes)" ON)
option(ENABLE_SILENT "Suppress all command line output (default: no)" OFF)
cmake_dependent_option(ENABLE_VERBOSE "Verbose messages (default: no)" OFF
  "NOT ENABLE_SILENT" OFF)
option(ENABLE_TWO_SCALE_SOLVER "Enable use of the two-scale solver" ON)
option(ENABLE_SEMI_ANALYTIC_SOLVER "Enable use of the semi-analytic solver" OFF)
option(ENABLE_LATTICE_SOLVER "Enable use of the lattice solver" OFF)
option(ENABLE_TEST "Enable unit tests" OFF)
option(ENABLE_THREADS "Enable use of multi-threading" ON)
option(WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)

# enable/disable external packages
option(ENABLE_FFLITE "Enable use of the FFLite loop functions (default: no)" OFF)
option(ENABLE_HIMALAYA "Enable use of the Himalaya library (default: no)" OFF)
option(ENABLE_LAPACK "Enable use of the LAPACK library (default: yes)" ON)
option(ENABLE_LOOPTOOLS "Enable use of the LoopTools library (default: no)" OFF)
option(ENABLE_SQLITE "Enable use of the SQLite3 library (default: ON)" ON)
