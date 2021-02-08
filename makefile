CC = g++
RM = echo rm -rf
IDIR = -I include/

PLATFORM=$(shell uname -s)

CFLAGS = -g -O3 -std=c++11

ifeq ($(PLATFORM),Linux)
  LFLAGS = -lboost_program_options --static
else
  LFLAGS = -lboost_program_options 
endif

TOOLS = include/enum.h include/vector.h include/mxObject.h include/mxUtil.h include/rbtree.h include/set.h include/indexedHeap.h
BASICS = include/subindex.h include/VarSet.h include/Factor.h src/Factor.cpp
FULL   = include/FactorDense.h  src/FactorDense.cpp
SPARSE = include/FactorSparse.h src/FactorSparse.cpp

# Exact solvers: conditioned variable elimination (PR) or junction tree (MAR)
exactFull: ui/exact.cpp src/factorgraph.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/exact.cpp src/Factor.cpp src/factorgraph.cpp src/graphmodel.cpp $(LFLAGS) -o bin/exactFull

# Exact solver using "sparse" factor implementation
exactSparse: ui/exact.cpp src/factorgraph.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(SPARSE)
	$(CC) $(CFLAGS) $(IDIR) -DSPARSE_FACTORS ui/exact.cpp src/Factor.cpp src/factorgraph.cpp src/graphmodel.cpp $(LFLAGS) -o bin/exactSparse

# Basic Mini-bucket
mbe: ui/mbe.cpp include/mbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/mbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/mbe

# Basic weighted mini-bucket
wmb_basic: ui/wmb_basic.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/wmb_basic.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/wmb

# Weighted mini-bucket with importance sampling refinement
wmb_is: ui/wmb_is.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/wmb_is.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/wmb_is

# Weighted mini-bucket with search-based refinement (very simple; see Qi Lou's code for better versions)
#wmb_search: ui/wmb_search.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
#	$(CC) $(CFLAGS) $(IDIR) ui/wmb_search.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/wmb_search
#
# Qi Lou's solver (search + sampling)
qlou: ui/qlou.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp \
          src/mmap.cpp include/mmap.h src/mmapIS.cpp include/mmapIS.h src/randheur.cpp include/randheur.h \
          src/mcts.cpp include/mcts.h src/searchsample.cpp include/searchsample.h src/wmbsearch.cpp include/wmbsearch.h \
          $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/qlou.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp \
          src/mmap.cpp src/mmapIS.cpp src/mcts.cpp src/randheur.cpp src/searchsample.cpp src/wmbsearch.cpp \
          $(LFLAGS) -o bin/qlou

#roijers: ui/roijers.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
#	$(CC) $(CFLAGS) $(IDIR) ui/roijers.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/roijers

# Simple MPLP-like decomposition bound
mplp: ui/mplp.cpp src/Factor.cpp src/factorgraph.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/mplp.cpp src/Factor.cpp src/factorgraph.cpp src/graphmodel.cpp $(LFLAGS) -o bin/mplp


# UAI competition solvers
#
uai14: ui/uai14.cpp src/factorgraph.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(SPARSE)
	$(CC) $(CFLAGS) $(IDIR) ui/uai14.cpp src/Factor.cpp src/factorgraph.cpp src/graphmodel.cpp $(LFLAGS) -o bin/uai14

uai16_wmbis: ui/uai16_wmbis.cpp src/wmbe.cpp include/wmbe.h src/Factor.cpp src/graphmodel.cpp $(TOOLS) $(BASICS) $(FULL)
	$(CC) $(CFLAGS) $(IDIR) ui/uai16_wmbis.cpp src/wmbe.cpp src/Factor.cpp src/graphmodel.cpp $(LFLAGS) -o bin/uai16_wmbis



##########
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS)
	-@echo ' '

