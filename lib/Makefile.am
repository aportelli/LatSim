COM_CXXFLAGS = -Wall -Wconversion
if CXX_GNU
    COM_CXXFLAGS += -W -pedantic
else
if CXX_INTEL
    COM_CXXFLAGS +=
endif
endif

include eigen_files.mk
nobase_dist_pkginclude_HEADERS = $(eigen_files)

lib_LTLIBRARIES = libLatSim.la
	
libLatSim_la_SOURCES = \
    Global.cpp         \
    Layout.cpp         \
    Logger.cpp         \
    includes.hpp

libLatSim_ladir = $(pkgincludedir)
libLatSim_la_HEADERS =    \
    Expressions.hpp       \
    Global.hpp            \
    Hmc.hpp               \
    Integrator.hpp        \
    Lattice.hpp           \
    Layout.hpp            \
    LatOp.hpp             \
    Logger.hpp            \
    MpiTypes.hpp          \
    Phi4.hpp              \
    OmelyanIntegrator.hpp \
    VerletIntegrator.hpp

libLatSim_la_CXXFLAGS = $(COM_CXXFLAGS)

ACLOCAL_AMFLAGS = -I .buildutils/m4
