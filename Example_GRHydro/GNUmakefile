# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

ebase := Main_PerfectFluid

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Matter  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/InitialConditions/ScalarFields \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/BlackHoles \
            ./GRHlib            

include $(CHOMBO_HOME)/mk/Make.test
