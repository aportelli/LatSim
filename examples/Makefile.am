COM_CXXFLAGS = -Wall -Wconversion
if CXX_GNU
    COM_CXXFLAGS += -W -pedantic
else
if CXX_INTEL
    COM_CXXFLAGS +=
endif
endif

noinst_PROGRAMS =           \
    exMatLibBench           \
    exMatOpBench            \
    exScalOpBench

exMatLibBench_SOURCES  = exMatLibBench.cpp
exMatLibBench_CXXFLAGS = $(COM_CXXFLAGS)
exMatLibBench_LDFLAGS  = -L../lib/.libs -lLatSim

exMatOpBench_SOURCES  = exMatOpBench.cpp
exMatOpBench_CXXFLAGS = $(COM_CXXFLAGS)
exMatOpBench_LDFLAGS  = -L../lib/.libs -lLatSim

exScalOpBench_SOURCES  = exScalOpBench.cpp
exScalOpBench_CXXFLAGS = $(COM_CXXFLAGS)
exScalOpBench_LDFLAGS  = -L../lib/.libs -lLatSim

ACLOCAL_AMFLAGS = -I .buildutils/m4
