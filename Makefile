## =================================================================== ##
## this is file bubble3d, created at 10-Jun-2009                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##

LIBDIR = ../lib
CXX = clang
CXXFLAGS = -O1 -g -fPIC
LIBS += -lgsl -lgslcblas -lm
LIBS += -L. -L${TETGEN_DIR} -ltet
INCLUDES += -I. -I$(LIBDIR) 
INCLUDES += -I${PETSC_DIR}/include
INCLUDES += -I${TETGEN_DIR}

src += $(LIBDIR)/clVector.cpp
src += $(LIBDIR)/clMatrix.cpp
src += $(LIBDIR)/clDMatrix.cpp
src += $(LIBDIR)/PCGSolver.cpp
src += $(LIBDIR)/GMRes.cpp
src += $(LIBDIR)/PetscSolver.cpp
src += $(LIBDIR)/FEMLinElement3D.cpp
src += $(LIBDIR)/FEMMiniElement3D.cpp
#src += $(LIBDIR)/FEMQuadElement3D.cpp
src += $(wildcard ./*.cpp)

obj = $(src:%.cpp=%.o)

all: step bubble 2bubbles diskNuC diskNuCte diskNuZ

diskNuC: ./script/mainDiskNuC.o $(obj)
	 -${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuCte: ./script/mainDiskNuCte.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuZ: ./script/mainDiskNuZ.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskSurf: ./script/mainDiskSurf.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

2bubbles: ./script/main2Bubbles.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

staticDroplet: ./script/mainStaticDroplet.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

oscillating: ./script/mainOscillating.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

bubble: ./script/mainBubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

#--------------------------------------------------
# step: ./script/mainStep.o libtest.so
# 	$(CXX) -L. -ltest -o $@ $<
#-------------------------------------------------- 

step: ./script/mainStep.o $(obj)
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

erase:
	@rm -f core
	@find . -name "*~" -exec rm {} \;
	@rm -f ./vtk/*.vtk
	@rm -f ./msh/*.msh
	@rm -f ./sim/vk*.dat
	@rm -f ./sim/sim*.dat
	@rm -f ./bin/*.bin
	@rm -f ./dat/*.dat

deepclean: 
	@rm -f staticDroplet step bubble 2bubble diskNuC diskNuCte diskNuZ
	@rm -f 2bubbles diskSurf
	@rm -f libtest*
	@rm -f core
	@find $(LIBDIR) -name "*.o" -exec rm {} \;
	@find . -name "*.o" -exec rm {} \;
	@find . -name "*~" -exec rm {} \;
	@rm -f ./vtk/*.vtk
	@rm -f ./msh/*.msh
	@rm -f ./sim/vk*.dat
	@rm -f ./sim/sim*.dat
	@rm -f ./bin/*.bin
	@rm -f ./dat/*.dat

