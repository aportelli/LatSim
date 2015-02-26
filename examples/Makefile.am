COM_CXXFLAGS = -Wall -Wconversion
if CXX_GNU
    COM_CXXFLAGS += -W -pedantic
else
if CXX_INTEL
    COM_CXXFLAGS +=
endif
endif

noinst_PROGRAMS =           \
    exIndex                 \
    exMatCommBench          \
    exMatLibBench           \
    exMatOpBench            \
    exScalarSim             \
    exScalCommBench         \
    exScalOpBench

exIndex_SOURCES  = exIndex.cpp
exIndex_CXXFLAGS = $(COM_CXXFLAGS)

exMatCommBench_SOURCES  = exMatCommBench.cpp
exMatCommBench_CXXFLAGS = $(COM_CXXFLAGS)

exMatLibBench_SOURCES  = exMatLibBench.cpp
exMatLibBench_CXXFLAGS = $(COM_CXXFLAGS)

exMatOpBench_SOURCES  = exMatOpBench.cpp
exMatOpBench_CXXFLAGS = $(COM_CXXFLAGS)

exScalarSim_SOURCES  = exScalarSim.cpp
exScalarSim_CXXFLAGS = $(COM_CXXFLAGS)

exScalCommBench_SOURCES  = exScalCommBench.cpp
exScalCommBench_CXXFLAGS = $(COM_CXXFLAGS)

exScalOpBench_SOURCES  = exScalOpBench.cpp
exScalOpBench_CXXFLAGS = $(COM_CXXFLAGS)

ACLOCAL_AMFLAGS = -I .buildutils/m4
