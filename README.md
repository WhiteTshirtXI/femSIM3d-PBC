# Single and Two-Phase flow Simulator
#### Module: femSIM3d-PBC
---

*DISCLAIMER:* `femSIM2D - femSIM3D - multiphaseSIM2D - multiphaseSIM3D` 
are an in-house CFD toolsuite intended to academic use. The code is under 
ongoing development in collaboration with the following universities/labs:
- [GESAR Lab](http://www.gesar.uerj.br) at *UERJ*: State University of Rio 
  de Janeiro/Group for Environmental Studies in Reservoirs;
- [LTCM Lab](http://www.epfl.ch) at *EPFL*: Ecole Polytechnique Fédérale de
  Lausanne/Laboratoire de Chaleur et de Masse;

### Metadata
---
- Authors and team:;
  - Gustavo R. ANJOS, Ph.D. (main developer): [mailto](mailto:gustavo.anjos@uerj.br);
  - Norberto MANGIAVACCHI, Ph.D. (GESAR Lab coordinator): [mailto](mailto:norberto.mangiavacchi@gmail.com)
  - José PONTES, D.Sc. (senior researcher): [mailto](mailto:jose.pontes@uerj.br)
  - John R. THOME, Ph.D. (LTCM Lab coordinator): [mailto](mailto:john.thome@epfl.ch)
  - Gustavo PEIXOTO DE OLIVEIRA, D.Sc. (developer):[mailto](mailto:gustavo.oliveira@uerj.br)
  - Erik GROS, Ph.D. candidate (developer): [mailto](mailto:erik.gros@epfl.ch)

- Method: classical finite element

- Elements: 
  - Taylor-Hood MINI element (4 nodes);
  - Quadratic QUAD element (6 nodes);

- Time discretization methods:
  - Semi-Lagrangian: convection term;
  - Galerkin: all the remnant terms;

- Linear system solution:  
  - In-house solvers: PCG, GMRES;
  - [PETSc](http://www.mcs.anl.gov/petsc) library solvers 

- Meshing tools:
  - [TRIANGLE](http://www.cs.cmu.edu/~quake/triangle.html), by J.R. Shewchuk
  - [GMSH](http://geuz.org/gmsh/), by Christophe Geuzaine and Jean-François Remacle

### Version notes
---
TAG 1.0

- Centralized version of single-phase flow simulator:
  - mainDiskNuCte
  - mainDiskNuZ
  - mainDiskNuC
  - mainStep
  
TAG 2.0

- Single-phase flow simulator tested and working. 
Two-phase flow simulator with 3D point re-meshing. 
Treatment of surface points still under development.

TAG 3.0

- Sigle and Two-Phase flow simulator tested and working with some 
improvements and suggestions as follow:
  - `mainDiskNuCte`   (minimum suggested mesh: 6-20-40)
  - `mainDiskNuZ`     (minimum suggested mesh: 6-20-40)
  - `mainDiskNuC`     (minimum suggested mesh: 6-20-60)
  - `mainBubble`      (optional: Laplace3D to distribute smooth points)
  - `mainOscillating` (add setSphereToEllipsoid to disturb surface)
  - `mainCurvature`   (add setBiggerSphere to perturb sphere radius)
  - `mainFallingDrop` (add test case: not tested yet)
  - `main2Bubbles`    (2 bubbles test case - Norberto-Portela)
  - `mainMicro`       (micro channel test case)


### Reference papers 
---

For an overview of the methodology applied in the code and peer-reviwed
accepted results, see:

- G. P. Oliveira; G. Anjos; J. Pontes; N. Mangiavacchi; J. R. Thome,
  ALE/finite element modeling of an unconfined bubble plume in periodic
  domain - bubble shape and oscillation analysis, Journal of the
  Brazilian Society of Mechanical Sciences and Engineering, 2015.
- Anjos, G.R.; Borhani, N.; Mangiavacchi, N.; Thome, J.R., 3D Moving
  Mesh Finite Element Method for Two-Phase Flows, Journal of
  Computational Physics, 2014.
- Anjos, G.R. e Mangiavacchi, N. e Pontes, J, Three-dimensional finite
  element method for rotating disk flows, Journal of the Brazilian
  Society of Mechanical Sciences and Engineering, 2014.

### Suggested compilation instructions
---

*Remark:* guidelines for Unix-based OSs (Mac/Linux).


### Proposed directory layout
---

``` bash
$HOME/projects/cpp/femSIM3d
$HOME/projects/cpp/lib
$HOME/Programs/petsc/petsc-version
$HOME/Programs/tetgen/tetgen-version
$HOME/Programs/boost/boost-version
$HOME/Programs/triangle/triangle-version
```


### Useful env. variables (Linux example)
---

```bash
export PETSC_ARCH=linux-gnu-cxx-debug
export PETSC_DIR=$HOME/Programs/petsc/petsc-dev
export DATA_DIR=$HOME/projects/db
export FEM3D_DIR=$HOME/projects/cpp/femSIM3d
export FEMLIB_DIR=$HOME/projects/cpp/lib
``` 


### Download and install [PETSc](http://www.mcs.anl.gov/petsc/download/index.html)
---

- Define environmental variables:

  - Mac OS:
``` bash
echo 'export PETSC_DIR=$HOME/Programs/petsc/petsc-version' > ~/.bashrc 
echo 'export PETSC_ARCH=darwin13.1.0-debug' > ~/.bashrc 
source ~/.bashrc
```
  - Linux:
``` bash
echo 'export PETSC_DIR=$HOME/Programs/petsc/petsc-version' > ~/.bashrc 
echo 'export PETSC_ARCH=linux-gnu-cxx-debug' > ~/.bashrc 
source ~/.bashrc
```
where `petsc-version` is the directory where you have installed PETSc.

*Remark*: On terminal, type `uname -a` to check your architecture so as to
set PETSC_ARCH. 

- Configure PETSc:

``` bash
cd $PETSC_DIR && sh petsc.sh
```

Content of file `petsc.sh` (compiler names might be changed depending on your
OS version):
``` bash
./config/configure.py --with-cc=gcc --with-fc=gfortran --with-cxx=g++ \
--download-f-blas-lapack=1 --download-mpich=1 --with-clanguage=cxx \
--download-parmetis=1 --download-blacs=1 --download-scalapack=1 \
--download-mumps=1; make all test
``` 

- If `make all test` fails, check your settings and try again. 

- Issues: 
 - Compiler on Mac OS Yosemite: `gcc` links to `clang` and not more
   available through XCode. Install `gcc` via e.g. Homebrew or MacPort.


### Download and install [TRIANGLE](http://www.cs.cmu.edu/~quake/triangle.html)
---
- Set env. variable:
``` bash
echo 'export TRIANGLE_DIR=$HOME/Programs/triangle' > ~/.bashrc
source ~/.bashrc
```

- Compile:
``` bash
cd $TRIANGLE_DIR
make trilibrary
```
> Change type `REAL` to `REALL` due to incompatibility with libTypes.h.


### Optionals 
---


#### download and install [BOOST](http://www.boost.org)
check if python-dev package is installed

``` bash
./bootstrap.sh --prefix=$HOME/Programs/boost/1.52.0/;./bjam install
export BOOST_DIR=$HOME/Programs/boost/1.52.0
export DYLD_LIBRARY_PATH=${BOOST_DIR}/lib:DYLD_LIBRARY_PATH (Mac OS X)
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:LD_LIBRARY_PATH (Linux)
```

#### download and install [PETSC4PY](http://code.google.com/p/petsc4py)
``` bash
python setup.py install --home=$HOME/Programs/python/
```

### Pre-processing meshing tools
---

#### Download and install [GMSH](http://geuz.org/gmsh/)
---
-- install depending on your O.S. --

#### Mesh directory examples
---
See: 
- [2D](http://github.com/gustavorabello/gmsh-2d); 
- [3D](http://github.com/gustavorabello/gmsh-3d).


### Post-processing and visualization tool
---

#### download and install [PARAVIEW](http://www.paraview.org/)
---
-- install depending on your O.S. --


#### How to compile the code to use MINI/QUAD elements
---

All the main files can be used with both MINI and QUAD elements, 
but for *two-phase flows* the QUAD element is *not* tested and probably 
won't work. To change the elements along the code, 2 steps are required:

- Makefile: comment/uncomment FEMMiniElement2D/FEMQuadElement2D as:
``` 
src += ${FEMLIB_DIR}/FEMMiniElement3D.cpp
#src += ${FEMLIB_DIR}/FEMQuadElement3D.cpp
```
- `TElement.h`: change the definition of NUMGLEU: 5 (MINI) or 10 (QUAD) 

``` cpp
#define NUMGLEU 4    ///< number of velocity nodes per element (MINI)
// #define NUMGLEU 6 ///< number of velocity nodes per element (QUAD)
```

### Strange behaviors
---

- Viscosity term: matrix *K*
``` cpp
Simulator3D::unCoupled()
```

#### BUBBLE MODE:
---
 - static test
	fint = GTilde*heaviside -> not working p=0
  	fint = G*heaviside -> pressure is wrong
 - oscillating
    fint = GTilde*heaviside -> does not work
    fint = G*heaviside -> does not work
 - rising ( Mrho*g )
    fint = GTilde*heaviside -> work with overshooting and bubble deform
    fint = G*heaviside -> work with no overshooting --BEST--
         ( Mrho-M )*g
    fint = GTilde*heaviside -> overshooting, bubble deform and osc velocity
    fint = G*heaviside -> no overshooting, bubble ok and osc velocity

CALL: `uv = uTilde + dt*invMLumped*fint + dt*invMrhoLumped*gravity;`


#### DROPLET MODE:
---
 - static test
    fint = GTilde*heaviside -> not working p=0
	fint = G*heaviside -> pressure correct
 - oscillating test
    fint = GTilde*heaviside -> don't work
 	fint = G*heaviside -> work (not so good at the end), vel osc --BEST--
 - rising test: gravity = ( Mrho*g )
    fint = GTilde*heaviside -> overshooting, bubble deform
    fint = G*heaviside -> vel osc,
 - rising test: gravity = ( Mrho-M )*g
    fint = GTilde*heaviside -> osc velo, bubble deform, overshooting
    fint = G*heaviside -> osc velo, no overshooting, bubble ok

CALL: `uv = uTilde + invA*fint + invA*gravity;`
