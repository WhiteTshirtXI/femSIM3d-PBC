## =================================================================== ##
## this is file Makefile, created at 10-Jun-2007                       ##
## maintained by Gustavo Rabello dos Anjos                             ##
## e-mail: gustavo.rabello@gmail.com                                   ##
## =================================================================== ##


TARGET = ns3d
DIR = .
LIBDIR = ../lib
CXX = g++
CXXFLAGS = -O3
LIBS += -L/opt/local/lib -lgsl -lgslcblas -lm 
LIBS += -L$(HOME)/Programs/tetgen/1.4.0 -ltet
INCLUDES += -I$(DIR) 
INCLUDES += -I$(LIBDIR) 
INCLUDES += -I$(HOME)/Programs/tetgen/1.4.0
INCLUDES += -I/opt/local/include

OBJECTS += $(LIBDIR)/clVector.o
OBJECTS += $(LIBDIR)/clMatrix.o
OBJECTS += $(LIBDIR)/clDMatrix.o
OBJECTS += $(LIBDIR)/CGSolver.o
OBJECTS += $(LIBDIR)/PCGSolver.o
OBJECTS += $(LIBDIR)/CGSSolver.o
OBJECTS += $(LIBDIR)/GMRes.o
OBJECTS += $(LIBDIR)/GSLSolver.o
OBJECTS += $(LIBDIR)/FEMLinElement3D.o
OBJECTS += $(LIBDIR)/FEMMiniElement3D.o
OBJECTS += Solver.o
OBJECTS += Model3D.o
OBJECTS += TElement.o
OBJECTS += Galerkin.o
OBJECTS += Interface3D.o
OBJECTS += MeshSmooth.o
OBJECTS += SemiLagrangean.o
OBJECTS += Simulator3D.o
OBJECTS += InOut.o
OBJECTS += main.o

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o $(TARGET) 

%.o : %.cpp $(wildcard *.h)
	$(CXX) $(INCLUDES) -c $< $(CXXFLAGS) -o $@

.PHONY: clean

clean:
	@rm -f core
	@rm -f $(TARGET)
	@rm -f ns3dDiskNuCte 
	@rm -f ns3dDiskNuZ
	@rm -f ns3dDiskNuC
	@rm -f ns3dSurf
	@rm -f ns3dBubble
	@find . -name "*~" -exec rm {} \;
	@rm -f ./oscillating*
	@rm -f ./mesh*
	@rm -f ./vtk/*.{vtk,vtu}
	@rm -f ./sim/vk?.*
	@rm -f ./bin/*.bin
	@rm -f ./dat/*.dat

deepclean:
	@rm -f core
	@find $(LIBDIR) -name "*.o" -exec rm {} \;
	@find . -name "*.o" -exec rm {} \;
	@find . -name "*~" -exec rm {} \;
	@rm -f $(TARGET)
	@rm -f ./vtk/*.vtk
	@rm -f ./sim/vk?.*
	@rm -f ./sim/sim*.dat
	@rm -f ./relatorio.dat
	@rm -f ./info.dat

# -- BUILD SPECIFIC SIMULATION -- #

#--------------------------------------------------
# nuCte:
# 	@rm -f main.cpp mainReserv.cpp mainDiskSurf.cpp mainDiskNuZ.cpp
# 	@rm -f reserv diskSurf diskNuZ
# 	@rm -f ./malhas/step* ./malhas/reserv*
# 	@mv diskNuCte Makefile
# 
# nuZ:
# 	@rm -f main.cpp mainReserv.cpp mainDiskSurf.cpp mainDiskNuCte.cpp
# 	@rm -f reserv diskSurf diskNuCte
# 	@rm -f ./malhas/step* ./malhas/reserv*
# 	@mv diskNuZ Makefile
# 
# reserv:
# 	@rm -f main.cpp mainDiskNuCte.cpp mainDiskNuZ.cpp mainDiskSurf.cpp
# 	@rm -f diskNuCte diskSurf diskNuZ
# 	@rm -f ./malhas/step* ./malhas/disk*
# 	@mv reserv Makefile
# 
# surf:
# 	@rm -f main.cpp mainReserv.cpp mainDiskNuCte.cpp mainDiskNuZ
# 	@rm -f reserv diskNuCte diskNuZ
# 	@rm -f ./malhas/step* ./malhas/reserv*
# 	@mv diskSurf Makefile
#-------------------------------------------------- 

# makefile help
# $@ is the name of the file to be made
# $? is the names of the changed dependents
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
