// =================================================================== // 
// this is file MeshSmooth.cpp, created at 20-Ago-2009                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "MeshSmooth.h"

MeshSmooth::~MeshSmooth(){};

MeshSmooth::MeshSmooth(Model3D &_m,real _dt)
{
 m = &_m;
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 IEN = m->getIEN();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 neighbourVert = m->getNeighbourVert();
 neighbourPoint = m->getNeighbourPoint();
 surfMesh = m->getSurfMesh();
 inVert = m->getInVert();
 surface = m->getSurface();
 heaviside = m->getHeaviside();
 triEdge = m->getTriEdge();
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);
 uSmoothSurface.Dim(numVerts);
 vSmoothSurface.Dim(numVerts);
 wSmoothSurface.Dim(numVerts);
 dt = _dt;
}

clVector MeshSmooth::compute()
{
 //setSurface();

 stepSmooth();
 //stepSmoothSurface();
 setCentroid();

 clVector aux = uSmooth;
 aux.Append(vSmooth);
 aux.Append(wSmooth);
 return aux;
}

// calcula velocidade da malha em todos os vertices
void MeshSmooth::stepSmooth()
{
 real aux;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum;
 real size; // numero de elementos da lista
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 for (list<int>::iterator it=inVert->begin(); it!=inVert->end(); ++it)
 {
  plist = neighbourVert->at(*it);
  size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   xSum += X->Get(*vert)-X->Get(*it);
   ySum += Y->Get(*vert)-Y->Get(*it);
   zSum += Z->Get(*vert)-Z->Get(*it);
  }
  aux = (xSum/size)/dt; // velocidade USmooth
  uSmooth.Set(*it,aux);

  aux = (ySum/size)/dt; // velocidade VSmooth
  vSmooth.Set(*it,aux);

  aux = (zSum/size)/dt; // velocidade WSmooth
  wSmooth.Set(*it,aux);
  if( *it == 20 )
  {
   cout << xSum << " " << dt << " "  
	    << aux << endl;
  }
 }
} // fecha metodo stepSmooth

/*
 * Implicit Fairing of Irregular Meshes using Diffusion and Curvature
 * Flow
 * Mathieu Desbrun, Mark Meyer, Peter Schroder, Alan H. Barr
 * Caltech*
 *
 * Notes on Mesh Smoothing
 * Nicolas Bray
 *
 * OBS.: Theoretically uSmooth = (xSum/distSum)/dt, but I removed the dt
 * because the {u,v}Smooth velocities were too high and the
 * Simulator3D->dt based on the ALE velocity was decreasing too much for
 * values of order E-4. Now the parameter C3 on the scripts main may be
 * equal to 1 (in the past the max value was 0.01, that is in the same
 * scale of Simulator3D->dt, thus removing dt form the {u,v}Smooth we
 * can regularize the max ALE velocity.
 *
 * */
// calcula velocidade da malha em todos os vertices
void MeshSmooth::stepSmoothFujiwara()
{
 real aux;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum,distSum;
 //real size; // numero de elementos da lista
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 // loop over all the vertices except those belonging to the boundary
 for( int i=0;i<numVerts;i++ )
 {
  plist = neighbourVert->at(i);
  //size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  distSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   // distance between the vertex *it and all its neighbours
   real dist = distance( X->Get(*vert),Y->Get(*vert),Z->Get(*vert),
	                     X->Get(i),Y->Get(i),Z->Get(i) );
   // sum of distances
   distSum += dist;   

   // 
   xSum += ( X->Get(*vert)-X->Get(i) )*dist;
   ySum += ( Y->Get(*vert)-Y->Get(i) )*dist;
   zSum += ( Z->Get(*vert)-Z->Get(i) )*dist;
  }
  //aux = (1.0/distSum)*xSum; // USmooth
  aux = (1.0/distSum)*xSum/dt; // velocidade USmooth
  //aux = (2.0/distSum)*xSum/dt; // velocidade USmooth
  uSmooth.Set(i,aux);

  //aux = (1.0/distSum)*ySum; // VSmooth
  aux = (1.0/distSum)*ySum/dt; // velocidade VSmooth
  //aux = (2.0/distSum)*ySum/dt; // velocidade VSmooth
  vSmooth.Set(i,aux);

  //aux = (1.0/distSum)*zSum; //  WSmooth
  aux = (1.0/distSum)*zSum/dt; // velocidade WSmooth
  //aux = (2.0/distSum)*zSum/dt; // velocidade WSmooth
  wSmooth.Set(i,aux);
 }
} // fecha metodo stepSmooth

// calcula velocidade da malha em todos os vertices
void MeshSmooth::stepSurfaceSmoothFujiwara()
{
 real aux;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum,distSum;
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 for( int i=0;i<surface->Dim();i++ )
 {
  int surfaceNode = surface->Get(i);

  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  distSum = 0.0;

  int listSize = neighbourPoint->at(surfaceNode).size();
  list<int> plist = neighbourPoint->at(surfaceNode);
  list<int>::iterator vert=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   real P0x = X->Get(surfaceNode);
   real P0y = Y->Get(surfaceNode);
   real P0z = Z->Get(surfaceNode);

   int v1 = *vert;++vert;
   real P1x = X->Get(v1);
   real P1y = Y->Get(v1);
   real P1z = Z->Get(v1);

   real edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

   distSum += edgeLength;   
   xSum += ( P1x-P0x )*edgeLength;
   ySum += ( P1y-P0y )*edgeLength;
   zSum += ( P1z-P0z )*edgeLength;
  }
  aux = (1.0/distSum)*xSum/dt; // velocidade USmooth
  //aux = (2.0/distSum)*xSum/dt; // velocidade USmooth
  uSmoothSurface.Set(surfaceNode,aux);

  aux = (1.0/distSum)*ySum/dt; // velocidade USmooth
  //aux = (2.0/distSum)*ySum/dt; // velocidade USmooth
  vSmoothSurface.Set(surfaceNode,aux);

  aux = (1.0/distSum)*zSum/dt; // velocidade WSmooth
  //aux = (2.0/distSum)*zSum/dt; // velocidade WSmooth
  wSmoothSurface.Set(surfaceNode,aux);
 }
} // fecha metodo stepSmooth

void MeshSmooth::stepSmooth(clVector &_uVel,clVector &_vVel,clVector &_wVel)
{
 list<int> plist;
 list<int>::iterator vert;
 real uSum,vSum,wSum;
 int size; // numero de elementos da lista
 uSmoothSurface.Dim(numVerts);
 vSmoothSurface.Dim(numVerts);
 wSmoothSurface.Dim(numVerts);

 //for (list<int>::iterator it=inVert->begin(); it!=inVert->end(); ++it)
 for (int it=0;it<numVerts;it++ ) 
 {
  plist = neighbourVert->at(it);
  size = 0;
  uSum = 0.0;
  vSum = 0.0;
  wSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   if( heaviside->Get(*vert) == 0.5 )
   {
   uSum += _uVel.Get(*vert);
   vSum += _vVel.Get(*vert);
   wSum += _wVel.Get(*vert);
   size++;
   }
//--------------------------------------------------
//   if( it == 3361 )
//   {
//   cout << "-------- " << *vert << " --------" << endl;
//   cout << "u: " << _uVel.Get(*vert) << endl;
//   cout << "v: " << _vVel.Get(*vert) << endl;
//   cout << "w: " << _wVel.Get(*vert) << endl;
//   cout << "uSum: " << uSum << endl;
//   cout << "vSum: " << vSum << endl;
//   cout << "wSum: " << wSum << endl;
//   cout << "X: " << X->Get(it) << endl;
//   cout << "Y: " << Y->Get(it) << endl;
//   cout << "Z: " << Z->Get(it) << endl;
//   cout << "dt: " << dt << endl;
//   cout << "--------------------" << endl;
//   }
//-------------------------------------------------- 
  }
  if(size > 0 )
  {
  uSmooth.Set( it,uSum/size ); 
  vSmooth.Set( it,vSum/size ); 
  wSmooth.Set( it,wSum/size ); 
  }

 }
} // fecha metodo stepSmooth

void MeshSmooth::stepSmoothLonger(clVector &_uVel,clVector &_vVel,clVector &_wVel)
{
 list<int> plist;
 list<int>::iterator vert;
 real uSum,vSum,wSum;
 int size; // numero de elementos da lista
 uSmoothSurface.Dim(numVerts);
 vSmoothSurface.Dim(numVerts);
 wSmoothSurface.Dim(numVerts);

 //for (list<int>::iterator it=inVert->begin(); it!=inVert->end(); ++it)
 for (int it=0;it<numVerts;it++ ) 
 {
  plist = neighbourVert->at(it);
  size = plist.size();
  uSum = 0.0;
  vSum = 0.0;
  wSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   uSum += _uVel.Get(*vert);
   vSum += _vVel.Get(*vert);
   wSum += _wVel.Get(*vert);
  }
  uSmooth.Set( it,uSum/size ); 
  vSmooth.Set( it,vSum/size ); 
  wSmooth.Set( it,wSum/size ); 
 }
} // fecha metodo stepSmooth

void MeshSmooth::setCentroid()
{  
 // calculando os valores de uSmooth e vSmooth nos centroides
 // atraves da media dos valores nos vertices
 int v[5];
 real aux;

 for( int mele=0;mele<numElems;mele++ )
 {
  v[0] = (int) IEN->Get(mele,0);
  v[1] = (int) IEN->Get(mele,1);
  v[2] = (int) IEN->Get(mele,2);
  v[3] = (int) IEN->Get(mele,3);

  v[4] = (int) IEN->Get(mele,4);

  aux = ( uSmooth.Get(v[0])+uSmooth.Get(v[1])+
	      uSmooth.Get(v[2])+uSmooth.Get(v[3]) )/4.0;
  uSmooth.Set( v[4],aux );
  aux = ( vSmooth.Get(v[0])+vSmooth.Get(v[1])+
	      vSmooth.Get(v[2])+vSmooth.Get(v[3]) )/4.0;
  vSmooth.Set( v[4],aux );
  aux = ( wSmooth.Get(v[0])+wSmooth.Get(v[1])+
	      wSmooth.Get(v[2])+wSmooth.Get(v[3]) )/4.0;
  wSmooth.Set( v[4],aux );
 }
}

// calcula velocidade tangencial euleriana nos vertices da interface
// atraves da media dos vizinhos.
void MeshSmooth::stepSmoothSurface()
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 real aux;
 int surfaceNode;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum;
 real size; // numero de elementos da lista
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  plist = neighbourVert->at(surfaceNode);
  size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   xSum += X->Get(*vert);
   ySum += Y->Get(*vert);
   zSum += Z->Get(*vert);
  }
  real xAverage = xSum/size; // X medio
  aux = (xAverage - X->Get( surfaceNode ))/dt;
  uSmooth.Set( surfaceNode,aux );
//--------------------------------------------------
//   cout << surfaceNode << " " 
//        << xAverage << " "
//        << xAverage-X->Get(surfaceNode) << " "
//        << (xAverage-X->Get(surfaceNode)/dt) << " "
//        << X->Get(surfaceNode) << endl;
//-------------------------------------------------- 

  real yAverage = ySum/size; // Y medio
  aux = (yAverage - Y->Get( surfaceNode ))/dt;
  vSmooth.Set( surfaceNode,aux );

  real zAverage = zSum/size; // Z medio
  aux = (zAverage - Z->Get( surfaceNode))/dt;
  wSmooth.Set( surfaceNode,aux );
 }
} // fecha metodo stepSmoothSurface

// igual a stepSmoothSurface, porem pega apenas os vertices que estao
// localizados na superficie
void MeshSmooth::stepSmoothSurface2()
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 real aux;
 int surfaceNode;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum;
 real size; // numero de elementos da lista
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  plist = neighbourVert->at(surfaceNode);
  size = 0.0;
  xSum = 0.0;
  ySum = 0.0;
  zSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   if( cc->Get(*vert) == 0.5 )
   {
	xSum += X->Get(*vert);
	ySum += Y->Get(*vert);
	zSum += Z->Get(*vert);
	size++;
   }
  }
  real xAverage = xSum/size; // X medio
  aux = (xAverage - X->Get( surfaceNode ))/dt;
  uSmooth.Set( surfaceNode,aux );

  real yAverage = ySum/size; // Y medio
  aux = (yAverage - Y->Get( surfaceNode ))/dt;
  vSmooth.Set( surfaceNode,aux );

  real zAverage = zSum/size; // Z medio
  aux = (zAverage - Z->Get( surfaceNode ))/dt;
  wSmooth.Set( surfaceNode,aux );
 }
} // fecha metodo stepSmoothSurface2

// calcula velocidade tangencial nos vertices da interface
// atraves da media da posicao dos vizinhos.
void MeshSmooth::stepSmoothSurface(clVector &_uVel,
                                   clVector &_vVel,
								   clVector &_wVel)
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 int surfaceNode;
 list<int> plist;
 list<int>::iterator vert;
 real uSum,vSum,wSum;
 real size; // numero de elementos da lista
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  plist = neighbourVert->at(surfaceNode);
  size = plist.size();
  uSum = 0.0;
  vSum = 0.0;
  wSum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   uSum += _uVel.Get(*vert);
   vSum += _vVel.Get(*vert);
   wSum += _wVel.Get(*vert);
  }
  uSmooth.Set( surfaceNode,uSum/size ); 
  vSmooth.Set( surfaceNode,vSum/size ); 
  wSmooth.Set( surfaceNode,wSum/size ); 
 }
} // fecha metodo stepSmoothSurface

void MeshSmooth::setSurface()
{
 int aux;

 // xSurface representa os valores de X dos nos da interface
 // ySurface representa os valores de Y dos nos da interface
 xSurface.Dim(surface->Dim());
 ySurface.Dim(surface->Dim());
 zSurface.Dim(surface->Dim());
 for( int i=0;i<surface->Dim();i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get(aux));
  ySurface.Set(i,Y->Get(aux));
  zSurface.Set(i,Z->Get(aux));
 }
}

clVector* MeshSmooth::getUSmooth(){ return &uSmooth; }
clVector* MeshSmooth::getVSmooth(){ return &vSmooth; }
clVector* MeshSmooth::getWSmooth(){ return &wSmooth; }
clVector* MeshSmooth::getUSmoothSurface(){ return &uSmoothSurface; }
clVector* MeshSmooth::getVSmoothSurface(){ return &vSmoothSurface; }
clVector* MeshSmooth::getWSmoothSurface(){ return &wSmoothSurface; }

void MeshSmooth::stepSmoothFujiwaraByHeight()
{
 uSmooth.Dim(numVerts);
 vSmooth.Dim(numVerts);
 wSmooth.Dim(numVerts);
 for( int elem=0;elem<surfMesh->numElems;elem++ )
 {
  int v1 = surfMesh->IEN.Get(elem,0);
  int v2 = surfMesh->IEN.Get(elem,1);
  int v3 = surfMesh->IEN.Get(elem,2);

  // points
  real P1x = surfMesh->X.Get(v1);
  real P1y = surfMesh->Y.Get(v1);
  real P1z = surfMesh->Z.Get(v1);
  real P2x = surfMesh->X.Get(v2);
  real P2y = surfMesh->Y.Get(v2);
  real P2z = surfMesh->Z.Get(v2);
  real P3x = surfMesh->X.Get(v3);
  real P3y = surfMesh->Y.Get(v3);
  real P3z = surfMesh->Z.Get(v3);

  // centroid
  real xMid = (P1x+P2x+P3x)/3.0;
  real yMid = (P1y+P2y+P3y)/3.0;
  real zMid = (P1z+P2z+P3z)/3.0;

  // mid 12
  real xMid12 = (P1x+P2x)/2.0;
  real yMid12 = (P1y+P2y)/2.0;
  real zMid12 = (P1z+P2z)/2.0;

  // mid 13
  real xMid13 = (P1x+P3x)/2.0;
  real yMid13 = (P1y+P3y)/2.0;
  real zMid13 = (P1z+P3z)/2.0;

  // mid 23
  real xMid23 = (P2x+P3x)/2.0;
  real yMid23 = (P2y+P3y)/2.0;
  real zMid23 = (P2z+P3z)/2.0;

  // list of neighbouring points
  list<int> plist = neighbourVert->at(v1);
  for(list<int>::iterator vert=plist.begin(); vert != plist.end();++vert )
  {
   if( *vert > surfMesh->numVerts )
   {
	int vertID = surfMesh->vertIdRegion.Get(v1);

	real Pxvert = X->Get(*vert);
	real Pyvert = Y->Get(*vert);
	real Pzvert = Z->Get(*vert);

	real height1 = distance(P1x,P1y,P1z,Pxvert,Pyvert,Pzvert);
	real height2 = distance(P2x,P2y,P2z,Pxvert,Pyvert,Pzvert);
	real height3 = distance(P3x,P3y,P3z,Pxvert,Pyvert,Pzvert);
	// centroid
	real height4 = distance(xMid,yMid,zMid,Pxvert,Pyvert,Pzvert);
	// xMid12 
	real height5 = distance(xMid12,yMid12,zMid12,Pxvert,Pyvert,Pzvert);
	// xMid13 
	real height6 = distance(xMid13,yMid13,zMid13,Pxvert,Pyvert,Pzvert);
	// xMid23 
	real height7 = distance(xMid23,yMid23,zMid23,Pxvert,Pyvert,Pzvert);

	real minHeight = min(height1,height2);
	minHeight = min(minHeight,height3);
	minHeight = min(minHeight,height4);
	minHeight = min(minHeight,height5);
	minHeight = min(minHeight,height6);
	minHeight = min(minHeight,height7);

	if( minHeight < 0.4*triEdge[vertID] )
	{
	 list<int> plist2 = neighbourVert->at(*vert);
	 real xSum = 0.0;
	 real ySum = 0.0;
	 real zSum = 0.0;
	 real distSum = 0.0;
	 for(list<int>::iterator viz=plist2.begin(); viz != plist2.end(); ++viz )
	 {
	  // distance between the vertex *it and all its neighbours
	  real dist = distance( X->Get(*viz),Y->Get(*viz),Z->Get(*viz),
	                    	X->Get(*vert),Y->Get(*vert),Z->Get(*vert) );
	  // sum of distances
	  distSum += dist;   

	  // 
	  xSum += ( X->Get(*viz)-X->Get(*vert) )*dist;
	  ySum += ( Y->Get(*viz)-Y->Get(*vert) )*dist;
	  zSum += ( Z->Get(*viz)-Z->Get(*vert) )*dist;
	 }
	 //aux = (1.0/distSum)*xSum; // USmooth
	 real aux = (1.0/distSum)*xSum/dt; // velocidade USmooth
	 //aux = (2.0/distSum)*xSum/dt; // velocidade USmooth
	 uSmooth.Set(*vert,aux);

	 //aux = (1.0/distSum)*ySum; // VSmooth
	 aux = (1.0/distSum)*ySum/dt; // velocidade VSmooth
	 //aux = (2.0/distSum)*ySum/dt; // velocidade VSmooth
	 vSmooth.Set(*vert,aux);

	 //aux = (1.0/distSum)*zSum; //  WSmooth
	 aux = (1.0/distSum)*zSum/dt; // velocidade WSmooth
	 //aux = (2.0/distSum)*zSum/dt; // velocidade WSmooth
	 wSmooth.Set(*vert,aux);
	}
   }
  }
 }
}
