#!/bin/tcsh
#
# Makefile for IMPACTA
# This should work on Nyx with the 
# modules petsc/3.3.0-p2-intel-opt 
# and openmpi/1.6.0-intel and the
# default intel compilers
# or with GNU - tested on UBUNTU 10.04

COMPILER = mpicxx

# GNU compilers
#OPTFLAGS = -O4 
#DBGFLAGS = -v -g -O0 -Wall -pedantic

# INTEL compilers
OPTFLAGS = -std=c++14 -v
DBGFLAGS = -g -O0 -Wall -pedantic

SOURCES = impacta_main.cpp 

OBJECTS = ${addsuffix .o, ${basename ${SOURCES}}}

EXEC = impacta

SRCDIR = ../include
EXEDIR = ../bin
OBJDIR = ../tmp
INPDIR = ../input

BUILDSOURCES = ${addprefix ${SRCDIR}/,${SOURCES}}
BUILDOBJECTS = ${OBJECTS}
BUILDEXEC = ${addprefix ${EXEDIR}/,${EXEC}}


all : ${BUILDEXEC}
	mv ${OBJECTS} ${OBJDIR}/.

.PHONY: clean exec cleanerror

clean ::
	rm -f ${BUILDOBJECTS}
	rm -f ${BUILDEXEC}
	rm -f ${OBJDIR}/${BUILDOBJECTS}

exec:
	mpirun -np 2 ${EXEDIR}/${EXEC} -read ${INPDIR}/imstdin

cleanerror : 
	rm *.o*

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


${BUILDEXEC} : ${BUILDOBJECTS}
	-${COMPILER} ${OPTFLAGS} -o ${BUILDEXEC} ${PETSC_CC_INCLUDES} \
	${BUILDOBJECTS} ${PETSC_KSP_LIB} -L${BOOST_LIB} -I${BOOST_INCLUDE}

debug :
	${COMPILER} ${DBGFLAGS} -o ${BUILDEXEC} \
	${BUILDSOURCES} ${PETSC_CC_INCLUDES} ${PETSC_KSP_LIB} -L${BOOST_LIB} -I${BOOST_INCLUDE}


${BUILDOBJECTS} : ${BUILDSOURCES}
	cd ${OBJDIR}
	-${COMPILER} ${OPTFLAGS} ${PETSC_CC_INCLUDES} -c ${BUILDSOURCES}$ \
		${PETSC_KSP_LIB} -L${BOOST_LIB} -I${BOOST_INCLUDE}

