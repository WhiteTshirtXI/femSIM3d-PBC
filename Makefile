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
src += $(wildcard ${FEM3D_DIR}/*.cpp)

obj = $(src:%.cpp=%.o)

all: step bubble 2bubbles diskNuC diskNuCte diskNuZ \
     diskSurf curvature curvatureAndPressure \
	 staticDroplet staticTorus sessileDrop \
	 oscillating fallingDrop micro 

diskNuC: ${FEM3D_DIR}/script/mainDiskNuC.o $(obj)
	 -${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuCte: ${FEM3D_DIR}/script/mainDiskNuCte.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuZ: ${FEM3D_DIR}/script/mainDiskNuZ.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskSurf: ${FEM3D_DIR}/script/mainDiskSurf.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

2bubbles: ${FEM3D_DIR}/script/main2Bubbles.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvature: ${FEM3D_DIR}/script/mainCurvature.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressure: ${FEM3D_DIR}/script/mainCurvatureAndPressure.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

staticDroplet: ${FEM3D_DIR}/script/mainStaticDroplet.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

staticTorus: ${FEM3D_DIR}/script/mainStaticTorus.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

sessileDrop: ${FEM3D_DIR}/script/mainSessileDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

oscillating: ${FEM3D_DIR}/script/mainOscillating.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

fallingDrop: ${FEM3D_DIR}/script/mainFallingDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

bubble: ${FEM3D_DIR}/script/mainBubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

micro: ${FEM3D_DIR}/script/mainMicro.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

step: ${FEM3D_DIR}/script/mainStep.o $(obj)
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

