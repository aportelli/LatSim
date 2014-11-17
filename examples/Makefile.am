if CC_GNU
    COM_CFLAGS = -Wall -W -pedantic
else
if CC_INTEL
    COM_CFLAGS = -Wall
endif
endif

if CXX_GNU
    COM_CXXFLAGS = -Wall -W -pedantic
else
if CXX_INTEL
    COM_CXXFLAGS = -Wall
endif
endif

noinst_PROGRAMS =           \
    exOp

exOp_SOURCES = exOp.cpp
exOp_CFLAGS  = -g -O2
exOp_LDFLAGS  = -L../lib/.libs -lLatSim

ACLOCAL_AMFLAGS = -I .buildutils/m4
