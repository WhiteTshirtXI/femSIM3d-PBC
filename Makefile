## =================================================================== ##
## this is file bubble3d, created at 10-Jun-2009                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##

LIBDIR = ../lib
CXX = g++
CXXFLAGS = -O1 -g -fPIC
LIBS += -L/urs/lib -lgsl -lgslcblas -lm
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

diskNuC: ./scripts/mainDiskNuC.o $(obj)
	 -${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuCte: ./scripts/mainDiskNuCte.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuZ: ./scripts/mainDiskNuZ.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskSurf: ./scripts/mainDiskSurf.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

2bubbles: ./scripts/main2Bubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

bubble: ./scripts/mainBubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

#--------------------------------------------------
# step: ./scripts/mainStep.o libtest.so
# 	$(CXX) -L. -ltest -o $@ $<
#-------------------------------------------------- 

step: ./scripts/mainStep.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

libtest.so: $(obj)
	$(CXX) -shared $(LIBS) ${PETSC_KSP_LIB} $(INCLUDES) $(obj) -o $@

#--------------------------------------------------
# libtest: $(obj)
# 	$(CXX) -dynamiclib $(LIBS) ${PETSC_KSP_LIB} $(INCLUDES) $(obj) -o $@
# 	ln -s libtest libtest.dylib
#-------------------------------------------------- 

%.o: %.cpp $(wildcard *.h)
	$(CXX) $(INCLUDES) -c $< $(CXXFLAGS) -o $@
	
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

