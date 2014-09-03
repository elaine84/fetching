O = 2
ifneq ($(wildcard /usr/include/python2.7),)
PY = 1
else ifneq ($(wildcard /software/linux/x86_64/python-2.7.1/include/python2.7),)
PY = 2
endif

CFLAGS := -g -W -Wall -O$(O) -I.
LDFLAGS :=
LIBS := -lgsl -lgslcblas -lboost_mpi -lboost_serialization

SAMPLER_OBJS := samplers/normal.o 

SAMPLER_OBJS += samplers/shell.o
LIBS += -lboost_iostreams

ifndef MPICXX
ifneq ($(shell which mpic++ 2>/dev/null),)
MPICXX := mpic++
else ifneq ($(shell which mpicxx 2>/dev/null),)
MPICXX := mpicxx
else ifneq ($(shell which openmpic++ 2>/dev/null),)
MPICXX := openmpic++
else ifneq ($(shell which mpicxx-openmpi-mp 2>/dev/null),)
MPICXX := mpicxx-openmpi-mp
else
MPICXX := openmpicxx
endif

ifneq ($(origin CXX),default)
MPICXX := CXX=${CXX} ${MPICXX}
else
ifneq ($(shell which gcc-mp-4.7 2>/dev/null),)
CC = gcc-mp-4.7
CXX = g++-mp-4.7
MPICXX := CXX=${CXX} ${MPICXX}
endif
endif
endif

# check whether to supply dependency flags
DEPSDIR := .deps
ifneq ($(shell ${MPICXX} --showme | grep -e '^[^ ]*\(clang\|gcc\|g[+][+]\|icpc\)'),)
DEPCFLAGS = -MD -MF $(DEPSDIR)/$(subst /,-,$*).d -MP
endif

ifneq ($(shell ${MPICXX} --version | grep ^icpc),)
# disable stupid Intel compiler warnings
# #383: value copied to temporary, reference to temporary used
# #981: operands are operated  in unspecified order
# #1418: external definition with no prior declaration
CFLAGS := $(subst -W ,,$(CFLAGS))
CFLAGS += -wd383,981,1418
endif

ifneq ($(wildcard /opt/local/*),)
CFLAGS += -I/opt/local/include -I/opt/local/include/openmpi
LDFLAGS += -L/opt/local/lib
endif

ifneq ($(wildcard ${HOME}/gsl/*),)
CFLAGS += -I${HOME}/gsl/include
LDFLAGS += -L${HOME}/gsl/lib
else ifneq ($(wildcard ${HOME}/include/gsl),)
CFLAGS += -I${HOME}/include
LDFLAGS += -L${HOME}/lib
endif

ifneq ($(wildcard ${HOME}/prefix/*),)
CFLAGS += -I${HOME}/prefix/include
LDFLAGS += -L${HOME}/prefix/lib
endif

ifeq ($(PY),1)
LIBS += -lboost_python -lpython2.7
SAMPLER_OBJS += samplers/python.o
CFLAGS += -I/usr/include/python2.7
else ifeq ($(PY),2)
LDFLAGS += -L/software/linux/x86_64/python-2.7.1/lib
LIBS += -lboost_python -lpython2.7
SAMPLER_OBJS += samplers/python.o
CFLAGS += -I/software/linux/x86_64/python-2.7.1/include/python2.7
endif

ifeq ($(shell uname),Darwin)
LIBS := $(subst -lboost_mpi,-lboost_mpi-mt,$(LIBS))
LIBS := $(subst -lboost_serialization,-lboost_serialization-mt,$(LIBS))
LIBS := $(subst -lboost_iostreams,-lboost_iostreams-mt,$(LIBS))
LIBS := $(subst -lboost_python,-lboost_python-mt,$(LIBS))
CFLAGS := $(subst -W ,,$(CFLAGS))
endif

ifdef BOOST_INC_DIR
CFLAGS += -I${BOOST_INC_DIR}
endif

ifdef BOOST_LIB_DIR
LDFLAGS += -L${BOOST_LIB_DIR}
endif

CLEAN = fetching *~ *.o *.so samplers/*.o pysamplers/*.pyc

all: fetching testmaster testmpi

fetching: master.o stype.o tree.o $(SAMPLER_OBJS)
	$(MPICXX) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

testmaster: testmaster.o master.o stype.o
	$(MPICXX) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

testmpi: testmpi.o
	$(MPICXX) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cc $(DEPSDIR)/stamp
	$(MPICXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

$(DEPSDIR)/stamp:
	mkdir -p $(dir $@)
	touch $@

clean:
	$(RM) $(wildcard $(CLEAN))
	$(RM) -r $(DEPSDIR)

DEPFILES := $(wildcard $(DEPSDIR)/*.d)
ifneq ($(DEPFILES),)
include $(DEPFILES)
endif

.PHONY: all clean
.SUFFIXES:
