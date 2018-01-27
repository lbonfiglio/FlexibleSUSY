# define additional build options
include(CMakeDependentOption)

option(ENABLE_COLORS "Colored output (default: no)" OFF)
option(ENABLE_COMPILE "Compile the source code (default: yes)" ON)
option(ENABLE_COMPILER_WARNINGS "Enable compiler warnings (default: no)" OFF)
option(ENABLE_MASS_ERROR_CHECK "Check mass eigenvalue precision (default: no)" OFF)
option(ENABLE_META "Create model class (default: yes)" ON)
option(ENABLE_SILENT "Suppress all command line output (default: no)" OFF)
cmake_dependent_option(ENABLE_VERBOSE "Verbose messages (default: no)" OFF
  "NOT ENABLE_SILENT" OFF)
option(ENABLE_TWO_SCALE_SOLVER "Enable use of the two-scale solver" ON)
option(ENABLE_SEMI_ANALYTIC_SOLVER "Enable use of the semi-analytic solver" OFF)
option(ENABLE_LATTICE_SOLVER "Enable use of the lattice solver" OFF)
option(WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)
