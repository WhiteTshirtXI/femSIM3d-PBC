## =================================================================== ##
## this is file bubble3d, created at 10-Jun-2009                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##

LIBDIR = ../lib
CXX = g++
CXXFLAGS = -O1 -g
LIBS += -lgsl -lm
LIBS += -L. -L$(HOME)/Programs/tetgen/1.4.3 -ltet
INCLUDES += -I. -I$(LIBDIR) 
INCLUDES += -I${PETSC_DIR}/include
INCLUDES += -I$(HOME)/Programs/tetgen/1.4.3

src += $(LIBDIR)/clVector.cpp
src += $(LIBDIR)/clMatrix.cpp
src += $(LIBDIR)/clDMatrix.cpp
src += $(LIBDIR)/PCGSolver.cpp
src += $(LIBDIR)/GMRes.cpp
src += $(LIBDIR)/PetscSolver.cpp
src += $(LIBDIR)/FEMLinElement3D.cpp
src += $(LIBDIR)/FEMMiniElement3D.cpp
src += $(wildcard ./*.cpp)

obj = $(src:%.cpp=%.o)

all: step bubble 2bubble diskNuC diskNuCte diskNuZ

diskNuC: ./scripts/mainDiskNuC.o libtest
	$(CXX) $< -L. -ltest -o $@

diskNuCte: ./scripts/mainDiskNuCte.o libtest
	$(CXX) $< -L. -ltest -o $@

diskNuZ: ./scripts/mainDiskNuZ.o libtest
	$(CXX) $< -L. -ltest -o $@

diskSurf: ./scripts/mainDiskSurf.o libtest
	$(CXX) $< -L. -ltest -o $@

2bubble: ./scripts/main2Bubble.o libtest
	$(CXX) $< -L. -ltest -o $@

bubble: ./scripts/mainBubble.o libtest
	$(CXX) $< -L. -ltest -o $@

step: ./scripts/mainStep.o libtest
	$(CXX) $< -L. -ltest -o $@

#--------------------------------------------------
# libtest: $(obj)
# 	$(CXX) -shared $(LIBS) ${PETSC_KSP_LIB} $(INCLUDES) $(obj) -o $@
# 	ln -s libtest libtest.so
#-------------------------------------------------- 

libtest: $(obj)
	$(CXX) -dynamiclib $(LIBS) ${PETSC_KSP_LIB} $(INCLUDES) $(obj) -o $@
	ln -s libtest libtest.dylib

%.o: %.cpp $(wildcard *.h)
	$(CXX) -fPIC $(INCLUDES) -c $< $(CXXFLAGS) -o $@
	
# Petsc new config
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

deepclean:
	@rm -f step bubble 2bubble diskNuC diskNuCte diskNuZ
	@rm -f libtest*
	@rm -f core
	@find $(LIBDIR) -name "*.o" -exec rm {} \;
	@find . -name "*.o" -exec rm {} \;
	@find . -name "*~" -exec rm {} \;
	@rm -f ./vtk/*.vtk
	@rm -f ./sim/vk*.dat
	@rm -f ./sim/sim*.dat
	@rm -f ./bin/*.bin
	@rm -f ./bin/*.dat
	@rm -f ./*.dat
	@rm -f ./relatorio.dat
	@rm -f ./info.dat

