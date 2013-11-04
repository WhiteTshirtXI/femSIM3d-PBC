
#include "Galerkin.h"

Galerkin::Galerkin(Model3D &_m,clVector &_uSol,
                               clVector &_vSol,
							   clVector &_wSol,
							   clVector &_cSol,
							   clMatrix &_gx,
							   clMatrix &_gy,
							   clMatrix &_gz)
{
 outflow = _m.getOutflow();
 numVerts = _m.getNumVerts();
 numNodes = _m.getNumNodes();
 numElems = _m.getNumElems();
 uSol = _uSol;
 vSol = _vSol;
 wSol = _wSol;
 cSol = _cSol;
 gx = _gx;
 gy = _gy;
 gz = _gz;
}

void Galerkin::compute(double dt)
{
 clDMatrix diagU(uSol);
 clDMatrix diagV(vSol);
 clDMatrix diagW(wSol);
 clMatrix conv(numNodes,numVerts);

 // cria matrix de conveccao u*nabla
 conv = (diagU * gx) + (diagV * gy) + (diagW * gz); 

 // copia somente as velocidades dos vertices 
 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 clVector cVert(numVerts);
 uSol.CopyTo(0,uVert);
 vSol.CopyTo(0,vVert);
 wSol.CopyTo(0,wVert);
 cSol.CopyTo(0,cVert);

 convU = (-1)*conv * uVert; 
 convV = (-1)*conv * vVert;
 convW = (-1)*conv * wVert;
 convC = (-1)*conv * cVert;

 // condicao de contorno
 convU = convU.MultVec(*outflow);
 convV = convV.MultVec(*outflow);
 convW = convW.MultVec(*outflow);
 //convC = convC.MultVec(*outflow);

} // fecha metodo compute 

clVector* Galerkin::getConvU(){ return &convU; }
clVector* Galerkin::getConvV(){ return &convV; }
clVector* Galerkin::getConvW(){ return &convW; }
clVector* Galerkin::getConvC(){ return &convC; }

