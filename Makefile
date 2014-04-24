## =================================================================== ##
## this is file Makefile, created at 10-Jun-2009                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##

# Compilers
#CXX = g++ 
CXX = clang++

# Flags
CXXFLAGS = -g -fPIC

# Libraries and includes
LIBS += -lm 
LIBS += -L${TETGEN_DIR} -ltet
INCLUDES += -I. -I${FEMLIB_DIR} 
INCLUDES += -I${PETSC_DIR}/include
INCLUDES += -I${TETGEN_DIR}
FEM3D_DIR = .

# Sorces
src += ${FEMLIB_DIR}/clVector.cpp
src += ${FEMLIB_DIR}/clMatrix.cpp
src += ${FEMLIB_DIR}/clDMatrix.cpp
src += ${FEMLIB_DIR}/PCGSolver.cpp
src += ${FEMLIB_DIR}/PetscSolver.cpp
src += ${FEMLIB_DIR}/FEMLinElement3D.cpp
src += ${FEMLIB_DIR}/FEMMiniElement3D.cpp
#src += ${FEMLIB_DIR}/FEMQuadElement3D.cpp
src += $(wildcard ${FEM3D_DIR}/*.cpp)

# Rules
obj = $(src:%.cpp=%.o)

all: single-phase two-phase two-phaseHT

single-phase: diskNuC diskNuZ diskNuCte diskSurf finiteDisk step stepALE \
              sphereNuCte 

single-phasePBC: channelPBC taylorVortexPBC 

two-phase: sphere cylinder torus curvatureSphere curvatureCylinder \
           curvatureHyperboloid curvatureTorus curvatureAndPressureSphere \
	       curvatureAndPressureCylinder curvatureAndPressureTorus\
		   hyperboloid curvatureAndPressureHyperboloid \
		   sessileDrop oscillatingDrop fallingDrop risingBubble \
		   2Bubbles micro zalesak vortex curvatureTest shear movingFrame \
		   
two-phasePBC: mainChannelPBC mainMovingFramePBC

two-phaseHT: risingBubbleHT sphereHMT


# --------------------<< backward step (single-phase) >>-------------------- #
#                                                                            #
step: ${FEM3D_DIR}/script/mainStep.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

stepALE: ${FEM3D_DIR}/script/mainStepALE.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

#                                                                            #
# -------------------------------------------------------------------------- #


# ------------------------------<< periodic >>------------------------------ #
#
channelPBC: ${FEM3D_DIR}/script/mainChannelPBC.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

taylorVortexPBC: ${FEM3D_DIR}/script/mainTaylorVortexPBC.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

mainMovingFramePBC: ${FEM3D_DIR}/script/mainMovingFramePBC.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
	#


# ------------------<< rotating disk (single-phase) >>---------------------- #
#                                                                            #
diskNuC: ${FEM3D_DIR}/script/mainDiskNuC.o $(obj)
	 -${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuCte: ${FEM3D_DIR}/script/mainDiskNuCte.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

sphereNuCte: ${FEM3D_DIR}/script/mainSphereNuCte.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

finiteDisk: ${FEM3D_DIR}/script/mainFiniteDisk.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

diskNuZ: ${FEM3D_DIR}/script/mainDiskNuZ.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #


# -----------------<< free surface disk (single-phase) >>------------------- #
#                                                                            #
# free surface (single-phase)
diskSurf: ${FEM3D_DIR}/script/mainDiskSurf.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #


# --------<< benchmarks for sphere,cylinder and torus (two-phase) >>-------- #
#                                                                            #
sphere: ${FEM3D_DIR}/script/mainSphere.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

zalesak: ${FEM3D_DIR}/script/mainZalesak.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

vortex: ${FEM3D_DIR}/script/mainVortex.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

shear: ${FEM3D_DIR}/script/mainShear.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

annular: ${FEM3D_DIR}/script/mainAnnular.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

cylinder: ${FEM3D_DIR}/script/mainCylinder.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

hyperboloid: ${FEM3D_DIR}/script/mainHyperboloid.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

torus: ${FEM3D_DIR}/script/mainTorus.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureTest: ${FEM3D_DIR}/script/mainCurvatureTest.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureSphere: ${FEM3D_DIR}/script/mainCurvatureSphere.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureCylinder: ${FEM3D_DIR}/script/mainCurvatureCylinder.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureHyperboloid: ${FEM3D_DIR}/script/mainCurvatureHyperboloid.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureTorus: ${FEM3D_DIR}/script/mainCurvatureTorus.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressureSphere: ${FEM3D_DIR}/script/mainCurvatureAndPressureSphere.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressureCylinder: ${FEM3D_DIR}/script/mainCurvatureAndPressureCylinder.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressureHyperboloid: ${FEM3D_DIR}/script/mainCurvatureAndPressureHyperboloid.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

curvatureAndPressureTorus: ${FEM3D_DIR}/script/mainCurvatureAndPressureTorus.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

sessileDrop: ${FEM3D_DIR}/script/mainSessileDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

oscillatingDrop: ${FEM3D_DIR}/script/mainOscillatingDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

fallingDrop: ${FEM3D_DIR}/script/mainFallingDrop.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

risingBubble: ${FEM3D_DIR}/script/mainRisingBubble.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

movingFrame: ${FEM3D_DIR}/script/mainMovingFrame.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

movingFramePBC: ${FEM3D_DIR}/script/mainMovingFramePBC.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

test: ${FEM3D_DIR}/script/mainTest.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

restart: ${FEM3D_DIR}/script/mainRestart.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #


# --------------------------<< misc (two-phase) >>-------------------------- #
#                                                                            #
2Bubbles: ${FEM3D_DIR}/script/main2Bubbles.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
2AxiBubbles: ${FEM3D_DIR}/script/main2AxiBubbles.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
helmholtz: ${FEM3D_DIR}/script/mainHelmholtz.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #


# -----------------------<< microchannels (two-phase) >>-------------------- #
#                                                                            #
micro: ${FEM3D_DIR}/script/mainMicro.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #


# ---------------------<< heat-transfer (two-phase) >>---------------------- #
#                                                                            #
risingBubbleHT: ${FEM3D_DIR}/script/mainRisingBubbleHT.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@

sphereHMT: ${FEM3D_DIR}/script/mainSphereHMT.o $(obj)
	-${CLINKER} $(obj) $(LIBS) ${PETSC_KSP_LIB} $< -o $@
#                                                                            #
# -------------------------------------------------------------------------- #

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
	@rm -f ./dat/edge.*

deepclean: 
	@find . -type f -perm -o+rx -exec rm {} \;
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
	@rm -f ./sim/vk*
	@rm -f ./dat/edge.*

