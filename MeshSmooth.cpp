// =================================================================== // 
// this is file MeshSmooth.cpp, created at 20-Ago-2009                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "MeshSmooth.h"

MeshSmooth::~MeshSmooth(){};

MeshSmooth::MeshSmooth(Model3D &_m)
{
 m = &_m;
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
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
 surfaceViz = m->getSurfaceViz();
 inVert = m->getInVert();
 surface = m->getSurface();
 //nonSurface = m->getNonSurface();
 uSmooth.Dim(numNodes);
 vSmooth.Dim(numNodes);
 wSmooth.Dim(numNodes);
 uSmoothSurface.Dim(numNodes);
 vSmoothSurface.Dim(numNodes);
 wSmoothSurface.Dim(numNodes);
}

clVector MeshSmooth::compute(real _dt)
{
 dt = _dt;
 //setSurface();

 stepSmooth();
 //stepSmoothTangent();
 setCentroid();

 clVector aux = uSmooth;
 aux.Append(vSmooth);
 aux.Append(wSmooth);
 return aux;
}

// calcula velocidade euleriana em todos os vertices
void MeshSmooth::stepSmooth()
{
 real aux;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum,zSum;
 real size; // numero de elementos da lista
 xAverage.Dim(numVerts);
 yAverage.Dim(numVerts);
 zAverage.Dim(numVerts);
 uSmooth.Dim(numNodes);
 vSmooth.Dim(numNodes);
 wSmooth.Dim(numNodes);

 for (list<int>::iterator it=inVert->begin(); it!=inVert->end(); ++it)
 {
  plist = neighbourVert->at(*it);
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
  xAverage.Set( *it,xSum/size ); // X medio
  aux = (xAverage.Get(*it) - X->Get(*it))/dt; // velocidade USmooth
  uSmooth.Set(*it,aux);

  yAverage.Set( *it,ySum/size ); // Y medio
  aux = (yAverage.Get(*it) - Y->Get(*it))/dt; // velocidade VSmooth
  vSmooth.Set(*it,aux);

  zAverage.Set( *it,zSum/size ); // Y medio
  aux = (zAverage.Get(*it) - Z->Get(*it))/dt; // velocidade WSmooth
  wSmooth.Set(*it,aux);
//--------------------------------------------------
//   if( *it == 0 )
//   {
//   cout << "-------- " << *it << " --------" << endl;
//   cout << "uSmooth: " << uSmooth.Get(*it) << endl;
//   cout << "vSmooth: " << vSmooth.Get(*it) << endl;
//   cout << "wSmooth: " << wSmooth.Get(*it) << endl;
//   cout << "xAverage: " << xAverage << endl;
//   cout << "yAverage: " << yAverage << endl;
//   cout << "zAverage: " << zAverage << endl;
//   cout << "X: " << X->Get(*it) << endl;
//   cout << "Y: " << Y->Get(*it) << endl;
//   cout << "Z: " << Z->Get(*it) << endl;
//   cout << "dt: " << dt << endl;
//   cout << "--------------------" << endl;
//   }
//-------------------------------------------------- 
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
	      uSmooth.Get(v[2])+uSmooth.Get(v[3]) )*0.25;
  uSmooth.Set( v[4],aux );
  aux = ( vSmooth.Get(v[0])+vSmooth.Get(v[1])+
	      vSmooth.Get(v[2])+vSmooth.Get(v[3]) )*0.25;
  vSmooth.Set( v[4],aux );
  aux = ( wSmooth.Get(v[0])+wSmooth.Get(v[1])+
	      wSmooth.Get(v[2])+wSmooth.Get(v[3]) )*0.25;
  wSmooth.Set( v[4],aux );
 }
}

// calcula velocidade tangencial euleriana nos vertices da interface
// atraves da media dos vizinhos.
void MeshSmooth::stepSmoothTangent()
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 real aux,surfaceNode;

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);

  // ******************************************************** //
  // eh necessario pensar e refazer a estrategia!!!
  // ******************************************************** //
  real closer = search( i,xAverage.Get(surfaceNode),
	                      yAverage.Get(surfaceNode),
						  zAverage.Get(surfaceNode) );

  // criando variaveis para nomear pontos e calcular vetor unitario
  real P0x = X->Get( surfaceNode );
  real P0y = Y->Get( surfaceNode );
  real P0z = Z->Get( surfaceNode );
  real P1x = X->Get( closer );
  real P1y = Y->Get( closer );
  real P1z = Z->Get( closer );

  // calculo do vetor unitario na direcao 0-1
  // x1Unit e y1Unit sao coordenadas do vetor unitario com origem em
  // 0,0. Para recoloca-lo em sua posicao original eh necessario 
  // somar P0x e P0y ao vetor unitario; 
  real dist01 = sqrt( (P1x-P0x)*(P1x-P0x)+
                      (P1y-P0y)*(P1y-P0y)+
					  (P1z-P0y)*(P1z-P0y) );
  real x1Unit = (P1x-P0x)/dist01;
  real y1Unit = (P1y-P0y)/dist01;
  real z1Unit = (P1z-P0z)/dist01;

  // velocidade suavizada ja calculada acima
  real velX = uSmooth.Get(surfaceNode); // velocidade na direcao x
  real velY = vSmooth.Get(surfaceNode); // velocidade na direcao y
  real velZ = wSmooth.Get(surfaceNode); // velocidade na direcao y

  // calculo da projecao da velocidade v na direcao da aresta 0-1
  // proj = | v.u^ | . u^
  aux = velX*x1Unit + velY*y1Unit + velZ*z1Unit; // v.u^ (produto escalar)

  //real xProj =  aux*x1Unit; 
  //real yProj =  aux*y1Unit;
  //real zProj =  aux*z1Unit;

  uSmooth.Set(surfaceNode,0.0);
  vSmooth.Set(surfaceNode,0.0);
  wSmooth.Set(surfaceNode,0.0);
  //uSmooth.Set(surfaceNode,xProj);
  //vSmooth.Set(surfaceNode,yProj);
  //wSmooth.Set(surfaceNode,zProj);
 }
} // fecha metodo stepSmoothTangent

// a velocidade eh calculada tomando como base a media das posicoes dos
// vertices vizinhos pertencentes a interface. A distribuicao de pontos
// na interface fica mais uniforme que na do metodo stepSmoothTangent.
void MeshSmooth::stepSmoothTangent2()
{
 real surfaceNode,aux;
 list<int> plist;
 list<int>::iterator vert;
 real xSum,ySum;
 real size; // numero de elementos da lista


 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);

  list<int> plist;
  list<int>::iterator vert;
  plist = surfaceViz->at(i);
  size = plist.size();
  xSum = 0.0;
  ySum = 0.0;
  for( vert=plist.begin(); vert != plist.end(); ++vert )
  {
   xSum = xSum + X->Get(*vert);
   ySum = ySum + Y->Get(*vert);
  }
  xAverage.Set( surfaceNode,xSum/size ); // X medio
  aux = (xAverage.Get(surfaceNode) - X->Get(surfaceNode))/dt; // velocidade USmooth
  uSmooth.Set(surfaceNode,aux);

  yAverage.Set( surfaceNode,ySum/size ); // Y medio
  aux = (yAverage.Get(surfaceNode) - Y->Get(surfaceNode))/dt; // velocidade VSmooth
  vSmooth.Set(surfaceNode,aux);

  real closer = search( i,xAverage.Get(surfaceNode),yAverage.Get(surfaceNode),zAverage.Get(surfaceNode) );

  // criando variaveis para nomear pontos e calcular vetor unitario
  real P0x = X->Get( surfaceNode );
  real P0y = Y->Get( surfaceNode );
  real P1x = X->Get( closer );
  real P1y = Y->Get( closer );

  // calculo do vetor unitario na direcao 0-1
  // x1Unit e y1Unit sao coordenadas do vetor unitario com origem em
  // 0,0. Para recoloca-lo em sua posicao original eh necessario 
  // somar P0x e P0y ao vetor unitario; 
  real dist01 = sqrt( (P1x-P0x)*(P1x-P0x)+(P1y-P0y)*(P1y-P0y) );
  real x1Unit = (P1x-P0x)/dist01;
  real y1Unit = (P1y-P0y)/dist01;

  // velocidade suavizada ja calculada acima
  real velX = uSmooth.Get(surfaceNode); // velocidade na direcao x
  real velY = vSmooth.Get(surfaceNode); // velocidade na direcao y

  // calculo da projecao da velocidade v na direcao da aresta 0-1
  // proj = | v.u^ | . u^
  aux = velX*x1Unit + velY*y1Unit; // v.u^ (produto escalar)
  real xProj =  aux*x1Unit; 
  real yProj =  aux*y1Unit;

  //uSmooth.Set(surfaceNode,0.0);
  //vSmooth.Set(surfaceNode,0.0);
  uSmooth.Set(surfaceNode,xProj);
  vSmooth.Set(surfaceNode,yProj);
 }
} // fecha metodo stepSmoothTangent2

// calcula velocidade tangencial euleriana nos vertices da interface
// atraves da media dos vizinhos. Utiliza vetor de normais dos pontos da
// superficie
void MeshSmooth::stepSmoothTangent3()
{
 // velocidade da interface eh igual a velocidade do fluido
 // para isso precisa-se encontrar os vertices da interface e impor a
 // velocidade do fluido (calculada pelo Semi-lagrangeano)
 real aux;

 // loop nos vertices da interface
 for( int i=0;i<surface->Dim();i++ )
 {
  real surfaceNode = surface->Get(i);

  real xNormal = xNormalSurface->Get(surfaceNode);
  real yNormal = yNormalSurface->Get(surfaceNode);
  real zNormal = zNormalSurface->Get(surfaceNode);

  // calculo do vetor unitario na direcao tangencial
  // x1Unit e y1Unit sao coordenadas do vetor unitario com origem em
  // 0,0. Para recoloca-lo em sua posicao original eh necessario 
  // somar P0x e P0y ao vetor unitario; 
  real len = sqrt( xNormal*xNormal+yNormal*yNormal+zNormal*zNormal );
  real xNormalUnit = xNormal/len;
  real yNormalUnit = yNormal/len;
  real zNormalUnit = zNormal/len;
	
  // velocidade suavizada ja calculada em stepSmooth
  real velX = uSmooth.Get(surfaceNode); // velocidade na direcao x
  real velY = vSmooth.Get(surfaceNode); // velocidade na direcao y
  real velZ = wSmooth.Get(surfaceNode); // velocidade na direcao y

  // calculo da projecao da velocidade v na direcao da aresta 0-1
  // proj = | v.u^ | . u^
  aux = velX*xNormalUnit + velY*yNormalUnit + velZ*zNormalUnit;
  real uProj =  aux*xNormalUnit; 
  real vProj =  aux*yNormalUnit;
  real wProj =  aux*zNormalUnit;

  real uTang = velX-uProj;
  real vTang = velY-vProj;
  real wTang = velZ-wProj;

//--------------------------------------------------
//   cout << surfaceNode << " " 
//        << velX << " " << velY << " " 
// 	   << uProj << " " << vProj << " "
// 	   << xNormalUnit << " " << yNormalUnit << " "
//        << uTang << " " << vTang << endl;
//-------------------------------------------------- 

  uSmooth.Set(surfaceNode,0.0);
  vSmooth.Set(surfaceNode,0.0);
  wSmooth.Set(surfaceNode,0.0);
  // subtraindo componente normal
  //uSmooth.Set(surfaceNode,uTang);
  //vSmooth.Set(surfaceNode,vTang);
  //wSmooth.Set(surfaceNode,wTang);
 }
} // fecha metodo stepSmoothTangent3


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

int MeshSmooth::search(int node,real _XI, real _YI, real _ZI)
{
 list<int> plist;
 list<int>::iterator vert;
 real dist, dmin;
 plist = surfaceViz->at(node);

 // o primeiro elemento da lista eh a coordenada de trabalho
 vert=plist.begin();
 ++vert; // indo para o 2o. elemento da lista

 // pega o 2o. elemento da lista para comparar
 dmin = (_XI-X->Get(*vert))*(_XI-X->Get(*vert))+
        (_YI-Y->Get(*vert))*(_YI-Y->Get(*vert))+
        (_ZI-Z->Get(*vert))*(_ZI-Z->Get(*vert));

 int vmin=*vert; // atribuindo o vmin do primeiro elemento
 //cout << "vmin1: " << vmin << endl;

 ++vert; // indo para o 3o. elemento da lista
 //cout << "vmin2: " << *vert << endl;

 for( vert=vert; vert != plist.end(); ++vert )
 {
 dist = (_XI-X->Get(*vert))*(_XI-X->Get(*vert))+
        (_YI-Y->Get(*vert))*(_YI-Y->Get(*vert))+
        (_ZI-Z->Get(*vert))*(_ZI-Z->Get(*vert));

  if (dist<dmin)
  {
   dmin=dist;
   vmin=*vert; 
  }
 }
 return vmin;
}


