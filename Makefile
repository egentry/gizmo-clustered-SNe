#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=GIZMO
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype 
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype
#
###########
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   slighty by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly 
#   dealing with new files and filename conventions)
#
#############

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif
FINCL =


CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat GIZMO_config.h)

ifeq (FIRE_PHYSICS_DEFAULTS,$(findstring FIRE_PHYSICS_DEFAULTS,$(CONFIGVARS)))  # using 'fire default' instead of all the above
    CONFIGVARS += COOLING COOL_LOW_TEMPERATURES COOL_METAL_LINES_BY_SPECIES
    CONFIGVARS += GALSF METALS GALSF_SFR_MOLECULAR_CRITERION GALSF_SFR_VIRIAL_SF_CRITERION=0
    CONFIGVARS += GALSF_FB_GASRETURN GALSF_FB_HII_HEATING GALSF_FB_SNE_HEATING=1 GALSF_FB_RT_PHOTONMOMENTUM
    CONFIGVARS += GALSF_FB_LOCAL_UV_HEATING GALSF_FB_RPWIND_LOCAL GALSF_FB_RPROCESS_ENRICHMENT=6 GALSF_SFR_IMF_VARIATION
endif


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = #

GRACKLEINCL = -I./grackle/local/include
GRACKLELIBS = -L./grackle/local/lib -lgrackle -lhdf5

HDF5_HOME   = /usr
HDF5LIB     = -lhdf5
HDF5INCL    = -I/usr/include -DH5_USE_16_API

GSLLIB	    = -lgsl
GSLINCL	    = -I/usr/include


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
    FFTW_LIBNAMES =  #-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
    FFTW_LIBNAMES =  #-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
    FFTW_LIBNAMES =  #-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

# we only need fftw if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
# or if POWERSPEC_GRID is activated
ifeq (POWERSPEC_GRID, $(findstring POWERSPEC_GRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
  FFTW_LIBNAMES = #
endif

endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Curie")
CC       =  mpicc
CXX      =  mpicpc
FC	 =  $(CC)
OPTIMIZE = -O3 -funroll-loops -std=c99
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
HDF5_HOME   = /ccc/products/hdf5-1.18
GSL_HOME    = /ccc/products/gsl-1.14/
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += # -DUSE_MPI_IN_PLACE
endif

ifeq ($(SYSTYPE),"Horizon")
CC       =  mpicc
CXX      =  mpicpc
FC       =  $(CC)
OPTIMIZE = -O3 -funroll-loops -std=c99
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
HDF5_HOME   = /softs/hdf5/1.8.16
GSL_HOME    = /softs/gsl/2.1
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += # -DUSE_MPI_IN_PLACE
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Stampede")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O3 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
##MKL_LIBS = -L$(TACC_MKL_LIB) -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
#MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel mvapich2 gsl hdf5 fftw2
##  -- performance is very similar with impi (intel-mpi) instead of mpavich2, 
##   if preferred use that with MPICHLIB line uncommented
## newest version of code needed for compatibility with calls in MPI-2 libraries
##
endif



#----------------------------
ifeq ($(SYSTYPE),"MacBookPro")
CC       =  mpicc
CXX      =  mpiccxx
FC       =  mpifort
OPTIMIZE = -O1 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I/usr/local/include -I$(PORTINCLUDE)
GSL_LIBS = -L/usr/local/lib -L$(PORTLIB)
FFTW_INCL= -I/usr/local/include
FFTW_LIBS= -L/usr/local/lib
HDF5INCL = -I/usr/local/include -I$(PORTINCLUDE) -DH5_USE_16_API
HDF5LIB  = -L/usr/local/lib -L$(PORTLIB) -lhdf5 -lz
MPICHLIB = #
OPT     += #
##
## PFH: this is my own laptop installation (2013 MacBook Pro running Yosemite)
## --
## I have installed GSL and HDF5 through MacPorts (once you have it installed, just use:
## sudo port install gsl
## sudo port install hdf5
## then the shortcut PORTINCLUDE/PORTLIB are just my own links to the macports installation
##  directories. in my case they are the default:
## PORTLIB=/opt/local/lib
## PORTINCLUDE=/opt/local/include
## --
## Unfortunately, FFTW is more complicated, since macports, fink, and other repository systems
## do not support direct installation of the MPI version of FFTW2, which is what GIZMO needs
## if you want to run with PMGRID or POWERSPEC enabled (if not, it should just compile without
## FFTW just fine). Be sure to install FFTW 2.1.5: get it from http://www.fftw.org/
## then unpack it, go into the unpacked directory, and configure it with:
## ./configure --enable-mpi --enable-type-prefix --enable-float
## (this set of commands is important to install the correct version)
## then "make" and finally "sudo make install"
## that should install it to its default location, /usr/local/, which is where FFTW_INCL/FFW_LIBS
## are set to point (to the respective include and lib sub-directories). check to make sure you
## have the fftw libraries correctly installed.
## --
## With this done, and the code successfully compiled, you should be able to run it with
## mpirun -np X ./GIZMO 1>gizmo.out 2>gizmo.err &
## (here "X" is the number of processes you want to use, I'm assuming youre running from the
##  same directory with the code so ./GIZMO is just in the local directory, and GIZMO is the
##  compiled file, and the 1> and 2> commands route stdin and stderr to the desired files)
##--
## If you're having trouble, I recommend the excellent guides to installing GADGET-2 at:
## http://astrobites.org/2011/04/02/installing-and-running-gadget-2/
## and
## https://gauge.wordpress.com/2009/06/16/pitp-2009-installing-gadget2/
## (by Nathan Goldbaum and Javiera Guedes, respectively) -- the installation should be
## nearly identical here
##
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Zwicky")
CC       =  mpicc
CXX      =  mpicpc
FC       =  mpiifort -nofor_main
OPTIMIZE = -O3 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKL_HOME)/include
MKL_LIBS = -L$(MKL_HOME)/lib/em64t -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += # -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel/2011.4.191 impi/4.0.2.003 gsl/1.15-gcc HDF5 
##  -- the machine is quite picky, impi seems to be the only working mpi option right now
##  --  currently fftw2 isnt pre-installed, built library in my directory, with config flags:
##       ./configure --prefix=/home/phopkins/fftw --enable-mpi --enable-type-prefix --enable-float --with-gcc
##      linked via the above FFTW2_HOME=/home/phopkins/fftw (where the libraries are installed)
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"SciNet")
CC       =  mpicc     # sets the C-compiler
OPTIMIZE =  -O1 -xHost -funroll-loops -no-prec-div -fast-transcendentals -fp-model fast=2 -ipo  # speed
#OPTIMIZE += -openmp -parallel -par-num-threads=4  # for openmp mode
OPTIMIZE += -g -debug parallel -Wall  # warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE +=-openmp -parallel  # openmp required compiler flags
endif
FC       =  $(CC)
GSL_INCL =  -I${SCINET_GSL_INC}
GSL_LIBS =  -L${SCINET_GSL_LIB} #-limf
FFTW_INCL=  -I${SCINET_FFTW_INC}
FFTW_LIBS=  -L${SCINET_FFTW_LIB}
MPICHLIB =
HDF5INCL =  -I${SCINET_HDF5_INC} -DH5_USE_16_API
HDF5LIB  =  -L${SCINET_HDF5_LIB} -lhdf5 -lz
MPICHLIB =
##
## Notes:
## 
### benchmarking suggests these optimizations, 256 cores with omp=4 or 2, no DOUBLE, multidomain=16 or 32, altogether gives best performance in
###   simple galaxy merger experiment (6x speedup over 16-core run with old-but-highly-optimized code).
##
## module load intel use.experimental openmpi/intel/1.6.0 gsl fftw/2.1.5-intel-openmpi hdf5/intel-openmpi/1.8.9
## NOTYPEPREFIX_FFTW should not be set on this machine
## 
## flags: 
## OPT      += -DNOCALLSOFSYSTEM -DMPICH_IGNORE_CXX_SEEK -DNO_ISEND_IRECV_IN_DOMAIN
##   -- these used to be recommended, with new compiler settings they don't seem necessary, but may help
## If memory problems crash the code, recommend small-scale chunking: MPISENDRECV_SIZELIMIT=10-100 (but this costs speed!)
##
## old options with intelmpi (not as good as openmpi):
##    module load intel intelmpi gsl fftw/2.1.5-intel-intelmpi4 use.experimental hdf5/intelmpi/1.8.9
##    OPTIMIZE =  -O2 -m64 -mt_mpi -openmp -xhost -g -debug parallel -mcmodel=medium -funroll-loops -Wall
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Quest")
CC       =  mpiicc
CXX      =  mpiicpc
FC       =  $(CC)
OPTIMIZE = -O1 -funroll-loops
OPTIMIZE += -g -Wall -no-prec-div -ipo -heap-arrays
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I/software/gsl/1.16-intel/include
GSL_LIBS = -L/software/gsl/1.16-intel/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/software/FFTW/2.1.5-intel/include
FFTW_LIBS= -L/software/FFTW/2.1.5-intel/lib
HDF5INCL = -I/software/hdf5/1.8.12-serial/include -DH5_USE_16_API
HDF5LIB  = -L/software/hdf5/1.8.12-serial/lib -lhdf5 -lz
#MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## debugging:
#OPT     += -check_mpi -genv I_MPI_DEBUG 5
## modules to load:
##module load mpi/intel-mpi-4.1.0 gsl/1.16-intel hdf5/1.8.12-serial fftw/2.1.5-intel
endif




#----------------------------------------------------------------------------------------------
ifeq (Pleiades,$(findstring Pleiades,$(SYSTYPE)))
CC       =  icc -lmpi
CXX      =  icc -lmpi -lmpi++
FC       =  ifort -nofor_main -lmpi
ifeq ($(SYSTYPE),"Pleiades-Haswell")
OPTIMIZE = -O3 -ip -funroll-loops -no-prec-div -fp-model fast=2 -xCORE-AVX2 # Haswell cores
endif
ifeq ($(SYSTYPE),"Pleiades-SIBridge")
OPTIMIZE = -O3 -ip -funroll-loops -no-prec-div -fp-model fast=2 -xAVX # Sandy or Ivy-Bridge cores
endif
OPTIMIZE += -Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -qopenmp
endif
GMP_INCL =
GMP_LIBS =
GSL_INCL =
GSL_LIBS =
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5)/lib -lhdf5 -lz -L/nasa/szip/2.1/lib -lsz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
endif
##
## Notes:
##   1. modules to load:
##          module load comp-intel mpi-sgi/mpt hdf5/1.8.3/intel/mpt gsl python/2.7.9 szip
##   2. make sure you set the correct core-type: runs submitted to the wrong cores will not run
##   3. FFTW2: the pre-existing installation on Pleiades is incomplete and problematic.
##      you will need to install your own in your home directory. when building the library, use
##          ./configure --prefix=$HOME/fftw --enable-mpi --enable-type-prefix --enable-float
##      where "$HOME/fftw" can be renamed but is the install director (should be your home directory);
##      then you need to define the variable (here or in your bashrc file)
##          FFTW2_HOME=$HOME/fftw
##      (matching what you used for the installation) so that the code can find fftw2
##   4. in your job submission file, it is recommended for certain core types that additional settings
##      are used. for Sandy Bridge, they recommend:
##          setenv MPI_DSM_DISTRIBUTE 0
##          setenv KMP_AFFINITY disabled
##      before the actual lines submitting your job
##


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xO -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -parallel -openmp  # openmp
OPTIMIZE += -g -Wall -debug parallel # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =      # must be empty if using openmpi
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
##
## Notes:
## 
## include the following in your .bashrc file (there is no default fftw2 module):
## module load intel/10.1 openmpi/1.2.4 gmp gsl hdf5 #now have to add fftw2 manually
## export TACC_FFTW2_INC=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/include
## export TACC_FFTW2_LIB=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/lib
## export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
##
## Options
## OPT += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
##   are not necessary, but may improve stability in some cases
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc 
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -tp barcelona-64 -fast -Mipa=fast,inline -Munroll -Mvect -O4
OPTIMIZE += -mp -Mconcur  # openmp
OPTIMIZE += -Wall  # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPT     += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
## 
## Notes:
##
## include the following in your .bashrc file:
##   module load pgi mvapich gmp gsl fftw2 hdf5
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
## 
endif


ifeq ($(SYSTYPE),"Quest")
CC       =  mpiicc
CXX      =  mpiicpc
FC       =  $(CC)
OPTIMIZE = -O2 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2 
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I/software/gsl/1.16-intel/include
GSL_LIBS = -L/software/gsl/1.16-intel/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/software/FFTW/2.1.5-intel/include
FFTW_LIBS= -L/software/FFTW/2.1.5-intel/lib
HDF5INCL = -I/software/hdf5/1.8.12-serial/include -DH5_USE_16_API
HDF5LIB  = -L/software/hdf5/1.8.12-serial/lib -lhdf5 -lz
OPT     += -DUSE_MPI_IN_PLACE
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK 
FC       =  $(CC)
OPTIMIZE = -g -O2 -Wall -Wno-unused-but-set-variable
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
CXX      =  mpiCC  -m64
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
OPTIMIZE =   -g -Wall -fopenmp
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
OPT     +=  -DOLD_HDF5
endif

#----------------------------
ifeq ($(SYSTYPE),"denali")
CC       =  mpicc
CXX      =  mpicxx
FC       =  mpifort
OPTIMIZE = -O1 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lz
OPT     += #


tmp := $(shell cd grackle && csh ./configure)
tmp := $(shell echo "CONFIG_MACHINE = darwin" > grackle/src/clib/Make.config.machine)

tmp := $(shell git checkout grackle/src/clib/Make.mach.darwin)
tmp := $(shell echo "MACH_INSTALL_PREFIX = $(PWD)/grackle/local" >> grackle/src/clib/Make.mach.darwin)

endif

#----------------------------
ifeq ($(SYSTYPE),"hyades")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -g -Wall
GSL_INCL = -I$(HOME)/local/gsl-2.1/include
GSL_LIBS = -L$(HOME)/local/gsl-2.1/lib
FFTW_INCL= -I$(HOME)/local/fftw-2.1.5/include
FFTW_LIBS= -L$(HOME)/local/fftw-2.1.5/lib
HDF5INCL = -I$(HOME)/local/miniconda3/envs/hdf/include -DH5_USE_16_API
HDF5LIB  = -L$(HOME)/local/miniconda3/envs/hdf/lib -lhdf5 -lz
MPICHLIB = #

tmp := $(shell cd grackle && csh ./configure)

tmp := $(shell echo "CONFIG_MACHINE = linux-gnu " > grackle/src/clib/Make.config.machine)
tmp := $(shell git checkout grackle/src/clib/Make.mach.linux-gnu)
tmp := $(shell echo "MACH_INSTALL_PREFIX = $(PWD)/grackle/local" >> grackle/src/clib/Make.mach.linux-gnu)

endif


#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------



ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif


ifeq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
OPT += -DCONFIG_BFLOAT_8
else
GRACKLEINCL =
GRACKLELIBS =
endif



SYSTEM_OBJS =   system/system.o system/allocate.o system/mymalloc.o system/parallel_sort.o \
                system/peano.o system/parallel_sort_special.o system/mpi_util.o

GRAVITY_OBJS  = gravity/forcetree.o gravity/cosmology.o gravity/pm_periodic.o gravity/potential.o \
                gravity/gravtree.o gravity/forcetree_update.o gravity/pm_nonperiodic.o gravity/longrange.o \
                gravity/ags_hsml.o

HYDRO_OBJS = hydro/hydra_master.o hydro/density.o hydro/gradients.o eos/eos.o


L3_OBJS =


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = GIZMO

OBJS  =  main.o accel.o  timestep.o init.o restart.o io.o \
         predict.o global.o begrun.o run.o allvars.o read_ic.o \
         domain.o driftfac.o kicks.o ngb.o compile_time_info.o merge_split.o
FOBJS =

OBJS	+= $(GRAVITY_OBJS) $(HYDRO_OBJS) $(SYSTEM_OBJS)
OBJS	+= $(L3_OBJS)

INCL    += allvars.h proto.h gravity/forcetree.h domain.h system/myqsort.h kernel.h eos/eos.h Makefile \


ifeq (GALSF_SUBGRID_VARIABLEVELOCITY_DM_DISPERSION,$(findstring GALSF_SUBGRID_VARIABLEVELOCITY_DM_DISPERSION,$(CONFIGVARS)))
OBJS    += galaxy_sf/dm_dispersion_hsml.o
endif

ifeq (GRAIN_FLUID,$(findstring GRAIN_FLUID,$(CONFIGVARS)))
OBJS    += solids/grain_physics.o
endif

ifeq (GALSF,$(findstring GALSF,$(CONFIGVARS)))
OBJS    += galaxy_sf/sfr_eff.o
endif

ifeq (GALSF_FB_HII_HEATING,$(findstring GALSF_FB_HII_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/hII_heating.o
endif

ifeq (TWOPOINT_FUNCTION_COMPUTATION_ENABLED,$(findstring TWOPOINT_FUNCTION_COMPUTATION_ENABLED,$(CONFIGVARS)))
OBJS    += structure/twopoint.o
endif

ifeq (GALSF_FB_SNE_HEATING,$(findstring GALSF_FB_SNE_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/mechanical_fb.o
endif

ifeq (GALSF_FB_LUPI,$(findstring GALSF_FB_LUPI,$(CONFIGVARS)))
OBJS    += galaxy_sf/lupi_fb.o
endif

ifeq (CREASEY,$(findstring CREASEY,$(CONFIGVARS)))
OBJS    += creasey.o
endif

ifeq (GENTRY_FB,$(findstring GENTRY_FB,$(CONFIGVARS)))
OBJS    += gentry_fb.o
endif

ifeq (GALSF_FB_RPWIND_LOCAL,$(findstring GALSF_FB_RPWIND_LOCAL,$(CONFIGVARS)))
OBJS    += galaxy_sf/rp_localwinds.o
endif

ifeq (BLACK_HOLES,$(findstring BLACK_HOLES,$(CONFIGVARS)))
OBJS    += galaxy_sf/blackholes/blackhole.o
OBJS    += galaxy_sf/blackholes/blackhole_util.o
OBJS    += galaxy_sf/blackholes/blackhole_environment.o
OBJS    += galaxy_sf/blackholes/blackhole_feed.o
OBJS    += galaxy_sf/blackholes/blackhole_swallow_and_kick.o
INCL    += galaxy_sf/blackholes/blackhole.h
endif

ifeq (BH_LUPI,$(findstring BH_LUPI,$(CONFIGVARS)))
OBJS	+= galaxy_sf/blackholes/blackhole_lupi.o
OBJS	+= galaxy_sf/blackholes/blackhole_feed_lupi.o
endif

ifeq (SINGLE_STAR,$(findstring SINGLE_STAR,$(CONFIGVARS)))
OBJS	+= radiation/rt_utilities.o radiation/rt_CGmethod.o radiation/rt_source_injection.o radiation/rt_chem.o radiation/rt_cooling.o
OBJS    += galaxy_sf/sfr_eff.o galaxy_sf/hII_heating.o galaxy_sf/mechanical_fb.o galaxy_sf/rp_localwinds.o
OBJS    += galaxy_sf/blackholes/blackhole.o galaxy_sf/blackholes/blackhole_util.o galaxy_sf/blackholes/blackhole_environment.o galaxy_sf/blackholes/blackhole_feed.o galaxy_sf/blackholes/blackhole_swallow_and_kick.o
INCL    += galaxy_sf/blackholes/blackhole.h
endif



ifeq (SCFPOTENTIAL,$(findstring SCFPOTENTIAL,$(CONFIGVARS)))
OBJS    += modules/potentials/scf.o modules/potentials/scf_util.o
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += structure/fof.o
INCL	+= structure/fof.h
endif

ifeq (OUTPUTLINEOFSIGHT,$(findstring OUTPUTLINEOFSIGHT,$(CONFIGVARS)))
OBJS    += structure/lineofsight.o
endif

ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling/cooling.o
INCL	+= cooling/cooling.h
endif

ifeq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
OBJS    += cooling/grackle.o
endif

ifeq (BUBBLES,$(findstring BUBBLES,$(CONFIGVARS)))
OBJS    += modules/bubbles/bubbles.o
endif

ifeq (EOS_HELMHOLTZ,$(findstring EOS_HELMHOLTZ,$(CONFIGVARS)))
OBJS    += eos/eos_interface.o
INCL    += eos/helmholtz/helm_wrap.h
FOBJS   += eos/helmholtz/helm_impl.o eos/helmholtz/helm_wrap.o
FINCL   += eos/helmholtz/helm_const.dek eos/helmholtz/helm_implno.dek eos/helmholtz/helm_table_storage.dek eos/helmholtz/helm_vector_eos.dek
endif

ifeq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
OBJS	+= system/pinning.o
endif

ifeq (DISTORTIONTENSORPS,$(findstring DISTORTIONTENSORPS,$(CONFIGVARS)))
OBJS	+= modules/phasespace/phasespace.o modules/phasespace/phasespace_math.o
endif

ifeq (RT_,$(findstring RT_,$(CONFIGVARS)))
OBJS	+= radiation/rt_utilities.o radiation/rt_CGmethod.o radiation/rt_source_injection.o radiation/rt_chem.o radiation/rt_cooling.o
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o subfind/subfind_vars.o subfind/subfind_collective.o subfind/subfind_serial.o subfind/subfind_so.o subfind/subfind_cont.o \
	subfind/subfind_distribute.o subfind/subfind_findlinkngb.o subfind/subfind_nearesttwo.o subfind/subfind_loctree.o subfind/subfind_alternative_collective.o subfind/subfind_reshuffle.o \
	subfind/subfind_potential.o subfind/subfind_density.o
INCL	+= subfind/subfind.h
endif

ifeq (SIDM,$(findstring SIDM,$(CONFIGVARS)))
OBJS    +=  sidm/sidm_core.o sidm/sidm_allvars.o
INCL    +=  sidm/sidm_proto.h
endif

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS	+=  nuclear/nuclear_network_solver.o nuclear/nuclear_network.o
INCL	+=  nuclear/nuclear_network.h
endif

ifeq (TURB_DRIVING,$(findstring TURB_DRIVING,$(CONFIGVARS)))
OBJS	+= turb/turb_driving.o turb/turb_powerspectra.o
endif

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL) $(GRACKLEINCL)

ifeq (VIP,$(findstring VIP,$(CONFIGVARS)))
FFLAGS = $(FOPTIONS)
else
FFLAGS = $(OPTIONS)
endif


ifeq (ALTERNATIVE_PSORT,$(findstring ALTERNATIVE_PSORT,$(CONFIGVARS)))
OBJS  += fof_alt_psort.o modules/psort-1.0/error_handling.o
CXXFLAGS = $(CFLAGS)
FC    = $(CXX)
endif

FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 


LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW) $(GRACKLELIBS)

ifeq (OMP_NUM_THREADS,$(findstring OMP_NUM_THREADS,$(CONFIGVARS))) 
LIBS   +=  -lpthread
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c

$(FOBJS): %.o: %.f90
	$(FC) $(OPTIMIZE) -c $< -o $@

compile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~ compile_time_info.c GIZMO_config.h

all:
	make -C grackle/src/clib/ HDF5_HOME="$(HDF5_HOME)"
	mkdir -p grackle/local
	mkdir -p grackle/local/include
	mkdir -p grackle/local/lib
	make -C grackle/src/clib/ install
	make $(EXEC)

cleanall:
	make -C grackle/src/clib clean
	make clean
