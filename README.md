FlexibleSUSY
============

FlexibleSUSY provides Mathematica and C++ code to create spectrum
generators for non-minimal supersymmetric models.  It is designed for
generating fast and modular C++ code, allowing for easy modification,
extension and reuse.

 * Homepage:                https://flexiblesusy.hepforge.org
 * Mailing list:            flexiblesusy@projects.hepforge.org
 * Source code repository:  https://github.com/FlexibleSUSY
 * Bug reports:             https://github.com/FlexibleSUSY/FlexibleSUSY/issues
 * References: [[CPC 190 (2015) 139-172 (arxiv:1406.2319)](https://arxiv.org/abs/1406.2319),
                [JHEP 1701 (2017) 079 (arXiv:1609.00371)](https://arxiv.org/abs/1609.00371)
                [arXiv:1710.03760](https://arxiv.org/abs/1710.03760)]


Requirements
============

 * C++ compiler (g++ >= 4.8.4 or clang++ >= 3.8.1 or icpc >= 15.0.0)
 * Fortran compiler (gfortran, ifort)
 * Mathematica (version 7.0 or higher)
 * SARAH (version 4.11.0 or higher)   [http://sarah.hepforge.org]
 * Boost (version 1.37.0 or higher)   [http://www.boost.org]
 * Eigen 3 (version 3.1 or higher)    [http://eigen.tuxfamily.org]
 * GNU scientific library             [http://www.gnu.org/software/gsl/]

Optional:

 * BLAS                               [http://www.netlib.org/blas/]
 * LAPACK                             [http://www.netlib.org/lapack/]
 * LoopTools (version 2.8 or higher)  [http://www.feynarts.de/looptools/]
 * Himalaya                           [https://github.com/Himalaya-Library/Himalaya]


Installation of SARAH
=====================

FlexibleSUSY requires SARAH to be installed and to be loadable with
the Needs["SARAH`"] command from inside Mathematica.  We recommend the
following setup:

     cd ~/.Mathematica/Applications/
     wget https://www.hepforge.org/archive/sarah/SARAH-4.12.3.tar.gz
     tar -xf SARAH-4.12.3.tar.gz
     ln -s $PWD/SARAH-4.12.3/ SARAH

     cd ~/.Mathematica/Kernel/
     echo "AppendTo[\$Path, \"$HOME/.Mathematica/Applications/SARAH/\"];" \
        >> init.m

All the above steps can be executed at once with the `install-sarah`
script.  Example:

     ./install-sarah

See `./install-sarah --help` for more options.


How to create a model
=====================

0. Before you setup a FlexibleSUSY model, you have to provide a SARAH
   model file.  To make it available in FlexibleSUSY, you can put it
   either into FlexibleSUSY's SARAH model directory

       FlexibleSUSY/sarah/<model>/

   or directly into SARAH's own model directly

       SARAH/Models/<model>/

   Here <model> is the name of your model (e.g. MSSM, NMSSM, etc.).
   Note, that there are already plenty of pre-installed model files in
   FlexibleSUSY's and SARAH's model directories that can be used.

1. Create a new / initialize an existing FlexibleSUSY model:

       ./createmodel --name=<model> # TODO

   See `./createmodel --help` for more details.  Afterwards there will
   be

   * a model directory models/<model>/
   * a makefile module models/<model>/module.mk
   * a Mathematica start script models/<model>/start.m
   * and a FlexibleSUSY model file models/<model>/FlexibleSUSY.m

   To modify the model details (input parameters, boundary conditions,
   etc.), edit the FlexibleSUSY model file
   models/<model>/FlexibleSUSY.m .

2. Create the Makefile and register your model(s):

       mkdir build
       cd build

       cmake -DWITH_MODELS=<model> ..

   Multiple models can be specified, separated by a comma.  See
   `./configure --help # TODO` for more options.

3. Compile FlexibleSUSY with your model:

       make

   Use `make -j<N>` to use <N> CPU cores.  When `make` is executed,
   Mathematica is called, which generates the C++ code for the
   specified models.  All C++ source files are written to the
   directory models/<model>/ .  When `make` has finished, the C++ code
   is compiled and the following spectrum generators are created for
   each specified model:

      models/<model>/run_<model>.x :  command line spectrum generator
      models/<model>/run_<model>.m :  Mathematica interface


Example
=======

    ./createmodel --name=HSSUSY # TODO

    mkdir build
    cd build
    cmake -DWITH_MODELS=HSSUSY ..
    make
    ./models/HSSUSY/run_HSSUSY.x \
       --slha-input-file=model_files/HSSUSY/LesHouches.in.HSSUSY


Creating the soucre code documentation
======================================

FlexibleSUSY's source code documentation (including the generated
source code files) can be generated with Doxygen in HTML or man
format.  To generate the HTML documentation please run:

    make doc-html

The generated HTML index file can then be found in doc/html/index.html
and can be viewed with any HTML browser, e.g.

    firefox doc/html/index.html

To generate the man documentation please run:

    make doc-man

The generated man pages can then be found in doc/man/man3/ and can be
viewed as

    man doc/man/man3/model_file_options.3


Creating only the source code files (no compilation)
====================================================

If you want to only create the C++ source files for your model, but do
not want to compile the code, you can set the ENABLE_COMPILE variable
to a false value:

    cmake -DENABLE_COMPILE=OFF ..
    make

Here, cmake will not check for installed compilers or libraries.  It
will only search for Mathematica and SARAH.  The execution of `make`
will stop as soon as all C++ source code files are generated.  See
below for how to export the generated source code.


Compile only (don't generate source code)
=========================================

If you want to only compile already created the C++ source files for
your model, you can use the ENABLE_META configure option:

    cmake -DENABLE_META=OFF ..
    make

Here, configure will only check for installed compilers or libraries.
It will not check for Mathematica and SARAH.

This option is useful if you want to generate the source code on one
computer and then transfer the generated code to another computer to
compile it.  This option can also be used with the pre-generated
FlexibleSUSY models, which are provided at the FlexibleSUSY home page.

Warning: Please make sure all C++ source files of your model are
available in the model directory models/<model>/ .  Otherwise the
compilation will fail.


Exporting the generated source code
===================================

The complete FlexibleSUSY source code, including the generated C++
code for the specified model(s) (but without the Mathematica meta
code), can be exported to a new directory.  The exported source code
is a complete standalone package, with it's own build system.  To
export the code, one has to set the target directory during
configuration via the `--with-install-dir=` option.  For example:

    ./configure --with-models=<models> \
       --with-install-dir=/path/to/export/directory # TODO

Afterwards

    make install-src

must be executed, which will copy the generated C++ source code for
all <models> to /path/to/export/directory , together with the
non-model specific source code from config/ , doc/ , slhaea/ and src/
.  Afterwards, the standalone package can be build like this:

    cd /path/to/export/directory
    ./configure
    make

It is also possible to create a "model package", which includes only
the generated source code for a given model, but does not contain the
whole FlexibleSUSY build system.  This is useful when the source code
for a model should be generated on one computer and later transferred
to another one to be compiled.  To create such a "model package" run

    make pack-<model>-src

where <model> is the name of the model whose generated source code
shall be packed.  After `make' has finished, the package file
<model>.tar.gz can be found in the working directory.


Dynamic libraries
=================

If you want to create dynamic model libraries (instead of static
libraries, which is the default) you need to set
`BUILD_SHARED_LIBS=ON` when calling cmake:

    cmake -DBUILD_SHARED_LIBS=ON ..


Statically linked executables
=============================

External libraries can be linked statically to the spectrum generator
executables by calling cmake as

    cmake -DBUILD_SHARED_LIBS=OFF \
       -DCMAKE_EXE_LINKER_FLAGS="-static" \
       -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ..

LoopTools
=========

It is possible to use LoopTools (http://www.feynarts.de/looptools/)
for calculating the loop functions, instead of using SOFTSUSY's loop
functions.  To enable LoopTools configure FlexibleSUSY via

    cmake -DENABLE_LOOPTOOLS=ON ..

To use the LoopTools library and header files from a specific
directory configure via

    cmake -DENABLE_LOOPTOOLS=ON -DLoopTools_BUILD_DIR=$HOME/packages/LoopTools-2.13/build ..

Note: LoopTools 2.8 or higher is required.


Plotting the mass spectrum and RG running
=========================================

The pole mass spectrum and the RG flow can be written to data files
for easy plotting.  In the MSSM for example these data files can be
generated via

    ./models/MSSM/run_MSSM.x \
       --slha-input-file=model_files/MSSM/LesHouches.in.MSSM \
       --rgflow-output-file=MSSM_rgflow.dat \
       --spectrum-output-file=MSSM_spectrum.dat

The generated files "MSSM_rgflow.dat" and "MSSM_spectrum.dat" can be
plotted with the gnuplot scripts in the model directory:

    gnuplot -persist -e "filename='MSSM_spectrum.dat'" \
       models/MSSM/MSSM_plot_spectrum.gnuplot

    gnuplot -persist -e "filename='MSSM_rgflow.dat'" \
       models/MSSM/MSSM_plot_rgflow.gnuplot

The gnuplot scripts are just for illustration and currently plot all
DR-bar parameters, regardless of mass dimension, so the resulting plot
is not particularly informative.  However, the user may adapt the
scripts to plot any chosen subset of the parameters.


Addons
======

A FlexibleSUSY addon is a program or library, which uses parts of the
FlexibleSUSY libraries or the generated models or is integrated into
FlexibleSUSY.  An example is GM2Calc [arXiv:1510.08071], which is
included in FlexibleSUSY in form of an addon.  An addon can be created
via

    ./createaddon --name=<addon> # TODO

where <addon> is the name of the addon.  The createaddon script
creates the directory addons/<addon>/ and the corresponding makefile
module addons/<addon>/module.mk .  If an addon has been created with
the above script, the user may edit the makefile module
(addons/<addon>/module.mk) to add source files in to the three
variables

    LIB@ADDON@_SRC   (list of source files to be included in library)
    EXE@ADDON@_SRC   (list of source files with a main())
    LIB@ADDON@_HDR   (list of header files)

Example:

    LIB@ADDON@_SRC := $(DIR)/file1.cpp
    EXE@ADDON@_SRC := $(DIR)/run.cpp
    LIB@ADDON@_HDR := $(DIR)/file1.hpp

To configure and compile the addon run

    cmake -DWITH_ADDONS=<addon> ..
    make

make compiles all source files and creates the addon library
addons/<addon>/lib<addon>.a (including the object file file1.o in the
above example) and an executable (addons/<addon>/run.x in the above
example).


Mathematica interface
=====================

FlexibleSUSY can be called from within Mathematica using Wolfram's
LibraryLink.  By default, FlexibleSUSY creates a LibraryLink library
for each spectrum genreator.  The generated library can be found in
models/<model>/<model>_librarylink.so, where <model> is the model
name.

Example for the CMSSM:

    Get["models/CMSSM/CMSSM_librarylink.m"];
 
    (* Create a handle to a model given the input parameters.
       See Options[FSCMSSMOpenHandle] for all default options. *)
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
          m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
    ];
 
    (* calculate pole mass spectrum *)
    FSCMSSMCalculateSpectrum[handle]
 
    (* calculate observables *)
    FSCMSSMCalculateObservables[handle]
 
    (* close the model handle *)
    FSCMSSMCloseHandle[handle];

For each model, FlexibleSUSY creates an example Mathematica script
which illustrates the use of the Mathematica interface.  The generated
example can be found in

    models/<model>/run_<model>.m

which can be run for example as

    math -run "<< \"models/<model>/run_<model>.m\""

Before running it, the model parameters in the script should be set to
reasonable values.  More advanced examples can be found in the
FlexibleSUSY documentation.

Note: In order to compile the library, Mathematica must be installed.
To disable the LibraryLink interface, run cmake with
`-DENABLE_LIBRARYLINK=OFF`.


Cleaning
========

There are several make targets to remove generated files, compiled
object files, libraries or executables:

    make clean      # deletes all .d .o .a .x files

    make distclean  # does `clean` and `clean-generated`
                    # and deletes in addition:
                    # Makefile flexiblesusy-config config.*
                    # config/list_sarah_model_files.sh

    make clean-dep  # deletes all .d files

    make clean-executables # deletes all .x files

    make clean-generated   # deletes generated files

    make clean-lib  # deletes all libraries

    make clean-obj  # deletes all .o files

For each model <model> or addon there are specific clean targets to
remove model-specific files:

    make clean-<model>     # deletes .d .o .a .x files

    make distclean-<model> # same as `make clean-<model> clean-<model>-src`

    make clean-<model>-dep # deletes .d files

    make clean-<model>-lib # deletes model library

    make clean-<model>-obj # deletes .o files

    make clean-<model>-src # deletes generated files


Package content
===============

In the following all sub-directories within the FlexibleSUSY package
are listed:

 * `addons/` contains addons for FlexibleSUSY, such as GM2Calc

 * `config/` contains helper scripts and makefile modules for the build
   system

 * `doc/` contains the FlexibleSUSY documentation

 * `examples/` contains examples how to build you own spectrum generator
   based on FlexibleSUSY

 * `fflite/` contains an alternative implementation of the
   Passarino-Veltman loop functions, based on FF

 * `meta/` contains the Mathematica meta code which generates the
   spectrum generators

 * `model_files/` contains default model files for some frequently used
   models (SM, SplitMSSM, MSSM, NMSSM, SMSSM, UMSSM, etc.)

 * `model_specific/` contains model-specific higher order corrections
   for the MSSM, NMSSM, SM and SplitMSSM from the literature

 * `models/` This is the output directory where the generated C++ code
   for the spectrum generators will be stored.

 * `Output/` contains SARAHs model-specific output files

 * `sarah/` contains SARAH model files shipped with FlexibleSUSY

 * `slhaea/` contains the SLHA reader library from
   [https://github.com/fthomas/slhaea]

 * `src/` contains non-model specific FlexibleSUSY C++ source code

 * `templates/` contains C++ template files for the spectrum generators

 * `test/` contains the FlexibleSUSY test suite

 * `utils/` contains some utility scripts to perform scans or extract
   data from SLHA files
