## =================================================================== ##
## this is file bubble3d, created at 10-Jun-2009                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##

CXX = clang
CXXFLAGS = -g -fPIC
LIBS += -lgsl -lgslcblas -lm
LIBS += -L. -L${TETGEN_DIR} -ltet
INCLUDES += -I. -I${FEMLIB_DIR}
INCLUDES += -I${PETSC_DIR}/include
INCLUDES += -I${TETGEN_DIR}

src += ${FEMLIB_DIR}/clVector.cpp
src += ${FEMLIB_DIR}/clMatrix.cpp
src += ${FEMLIB_DIR}/clDMatrix.cpp
src += ${FEMLIB_DIR}/PCGSolver.cpp
src += ${FEMLIB_DIR}/GMRes.cpp
src += ${FEMLIB_DIR}/PetscSolver.cpp
src += ${FEMLIB_DIR}/FEMLinElement3D.cpp
src += ${FEMLIB_DIR}/FEMMiniElement3D.cpp
#src += ${FEMLIB_DIR}/FEMQuadElement3D.cpp
src += $(wildcard ./*.cpp)

obj = $(src:%.cpp=%.o)

all: step bubble 2bubbles diskNuC diskNuCte diskNuZ \
     diskSurf curvature curvatureAndPressure \
	 staticDroplet staticTorus sessileDrop \
	 oscillating fallingDrop micro 

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

curvature: ./script/mainCurvature.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressure: ./script/mainCurvatureAndPressure.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

staticDroplet: ./script/mainStaticDroplet.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

staticTorus: ./script/mainStaticTorus.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

sessileDrop: ./script/mainSessileDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

oscillating: ./script/mainOscillating.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

fallingDrop: ./script/mainFallingDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

bubble: ./script/mainBubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

micro: ./script/mainMicro.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

step: ./script/mainStep.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

#--------------------------------------------------
# step: ./script/mainStep.o libtest.so
# 	$(CXX) -L. -ltest -o $@ $<
#-------------------------------------------------- 

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
	@rm -f ./vtk/*.vtk ./vtk/*.vtu
	@rm -f ./msh/*.msh
	@rm -f ./sim/vk*
	@rm -f ./sim/sim*.dat
	@rm -f ./bin/*.bin
	@rm -f ./dat/*.dat
	@rm -f ./dat/vk*

deepclean: 
	@rm -f staticDroplet step bubble 2bubble diskNuC 
	@rm -f 2bubbles diskSurf staticTorus sessileDrop 
	@rm -f curvature fallingDrop diskNuCte diskNuZ
	@rm -f curvatureAndPressure
	@rm -f oscillating micro
	@rm -f libtest*
	@rm -f core
	@find ${FEMLIB_DIR} -name "*.o" -exec rm {} \;
	@find . -name "*.o" -exec rm {} \;
	@find . -name "*~" -exec rm {} \;
	@rm -f ./vtk/*.vtk
	@rm -f ./msh/*.msh
	@rm -f ./sim/vk*.dat
	@rm -f ./sim/sim*.dat
	@rm -f ./bin/*.bin
	@rm -f ./*.dat
	@rm -f ./dat/*.dat
	@rm -f ./dat/vk*

