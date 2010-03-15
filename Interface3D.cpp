// =================================================================== //
// this is file Interface3D.cpp, created at 26-Mar-2009                //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Interface3D.h"
#include "InOut.h"

Interface3D::Interface3D(Model3D &_m)
{
 cc = _m.getPointerCC();
 X = _m.getPointerX();
 Y = _m.getPointerY();
 Z = _m.getPointerZ();
 IEN = _m.getPointerIEN();
 numVerts = _m.getNumVerts();
 numNodes = _m.getNumNodes();
 numElems = _m.getNumElems();
 neighbourElem = _m.neighbourElem;
 neighbourVert = _m.neighbourVert; 
 neighbourFace = _m.neighbourFace; 
 elemSurface = _m.elemSurface; 
 neighbourFaceVert = _m.neighbourFaceVert;
 surfaceViz = _m.surfaceViz;
 xSurfaceViz = _m.xSurfaceViz;
 ySurfaceViz = _m.ySurfaceViz;
 zSurfaceViz = _m.zSurfaceViz;
 surface = _m.getPointerSurface();
 closer.Dim(numNodes);
 xCloser.Dim(numNodes);
 yCloser.Dim(numNodes);
 zCloser.Dim(numNodes);
 closerViz.Dim(numNodes);
 distance.Dim(numVerts);
 xCenter = _m.getXCenter();
 yCenter = _m.getYCenter();
 zCenter = _m.getZCenter();
 bubbleRadius = _m.getBubbleRadius();
 setSolverSmooth(new PCGSolver());

 //saveSurfaceVTK();
}

clVector Interface3D::curvature1()
{
 setCloser();
 real aux;

 // esta rotina retorna todos os valores em X,Y e Z dos nos da interface
 // vinculados pela menor distancia de todos os nos da malha.
 // calculo da distancia atraves da norma 2
 for( int i=0;i<numVerts;i++ )
 {
  aux = sqrt( (X->Get(i)-xCloser.Get(i))*(X->Get(i)-xCloser.Get(i))+
              (Y->Get(i)-yCloser.Get(i))*(Y->Get(i)-yCloser.Get(i))+
              (Z->Get(i)-zCloser.Get(i))*(Z->Get(i)-zCloser.Get(i)) );
  distance.Set(i,fabs(aux));
 }

 // sign distance function
 clVector ccAux = *cc;ccAux-0.5;ccAux=2*ccAux;

 // sign distance function (phi)
 distance = distance.MultVec(ccAux);

 return distance;
} // fecha metodo curvature1 

// imposicao da curvatura --> este metodo nao eh um calculo, so serve
// para testes simples do codigo
clVector Interface3D::curvature2()
{
 setCloser();
 real aux,aux1,aux2,aux3;
 clVector xAux(numNodes);
 clVector yAux(numNodes);
 clVector zAux(numNodes);

 //rotina para calculo da distancia dos nos da malha para seu
 // respectivo no da interface atraves do teorema de pitagoras
 for( int i=0;i<numVerts;i++ )
 {
  aux1 = (X->Get(i)-xCenter)*(X->Get(i)-xCenter);
  aux2 = (Y->Get(i)-yCenter)*(Y->Get(i)-yCenter);
  aux3 = (Z->Get(i)-zCenter)*(Z->Get(i)-zCenter);
  aux = sqrt( aux1+aux2+aux3 );
  aux = fabs(aux - bubbleRadius);
  distance.Set(i,aux);
 }

 // sign distance function
 clVector ccAux = *cc;ccAux-0.5;ccAux=2*ccAux;

 // sign distance function (phi)
 distance = distance.MultVec(ccAux);

 return distance;
} // fecha metodo curvature2

// calculo do Kappa geometricamente. Utiliza neighbourVert
clVector Interface3D::computeKappa()
{
 int surfaceNode;
 real area;
 listElem plist,plist2,plist3;
 list<int>::iterator face,vert;

 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);
  real fx = 0;
  real fy = 0;
  real fz = 0;
  real Px,Py,Pz;

  plist = neighbourVert.at(surfaceNode);
  for( vert=plist.begin();vert!=plist.end();++vert )
  {
   for( int j=0;j<surface->Dim();j++ )
   {
	if( *vert == surface->Get(j) )
	{
	 Px = X->Get(*vert);
	 Py = Y->Get(*vert);
	 Pz = Z->Get(*vert);
	 real dist = sqrt( (Px-P0x)*(Px-P0x)+(Py-P0y)*(Py-P0y)+(Pz-P0z)*(Pz-P0z) );
	 real xUnit = (Px-P0x)/dist;
	 real yUnit = (Py-P0y)/dist;
	 real zUnit = (Pz-P0z)/dist;
	 fx = fx + xUnit;
	 fy = fy + yUnit;
	 fz = fz + zUnit;
	 //cout << surfaceNode << " " << fx << " " << fy << " " << fz << endl;
	}
   }
  }
  //cout << "----------------------------"<<endl;

  area = 0;
  plist2 = elemSurface.at (surfaceNode); 
  for( face=plist2.begin();face!=plist2.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos da face
   plist3 = neighbourFaceVert.at (*face);
   vert=plist3.begin();
   real P1x = X->Get(*vert);
   real P1y = Y->Get(*vert);
   real P1z = Z->Get(*vert);
   ++vert;
   real P2x = X->Get(*vert);
   real P2y = Y->Get(*vert);
   real P2z = Z->Get(*vert);
   // comprimento das arestas a, b e c para aplicacao da formula de
   // Heron
   real a = sqrt( (P1x-P0x)*(P1x-P0x)+(P1y-P0y)*(P1y-P0y)+(P1z-P0z)*(P1z-P0z) );
   real b = sqrt( (P2x-P0x)*(P2x-P0x)+(P2y-P0y)*(P2y-P0y)+(P2z-P0z)*(P2z-P0z) );
   real c = sqrt( (P2x-P1x)*(P2x-P1x)+(P2y-P1y)*(P2y-P1y)+(P2z-P1z)*(P2z-P1z) );
   real s = (a+b+c)/2; //semiperimeter
   area = area + sqrt( s*(s-a)*(s-b)*(s-c) );
  }
  real force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );
  real L = sqrt(area);
  real pressure = force/L;
  //cout << surfaceNode << " " << area << " " << force << " " << pressure << endl;
  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa

// calculo do Kappa geometricamente. Utiliza neighbourVert
clVector Interface3D::computeKappa2()
{
 int surfaceNode;
 real area;
 listElem plist,plist2,plist3;
 list<int>::iterator face,vert;

 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);
  real fx = 0;
  real fy = 0;
  real fz = 0;
  real Px,Py,Pz;

  plist = neighbourVert.at(surfaceNode);
  for( vert=plist.begin();vert!=plist.end();++vert )
  {
   for( int j=0;j<surface->Dim();j++ )
   {
	if( *vert == surface->Get(j) )
	{
	 Px = X->Get(*vert);
	 Py = Y->Get(*vert);
	 Pz = Z->Get(*vert);
	 real dist = sqrt( (Px-P0x)*(Px-P0x)+(Py-P0y)*(Py-P0y)+(Pz-P0z)*(Pz-P0z) );
	 real xUnit = (Px-P0x)/dist;
	 real yUnit = (Py-P0y)/dist;
	 real zUnit = (Pz-P0z)/dist;
	 fx = fx + xUnit;
	 fy = fy + yUnit;
	 fz = fz + zUnit;
	 //cout << surfaceNode << " " << fx << " " << fy << " " << fz << endl;
	}
   }
  }
  //cout << "----------------------------"<<endl;

  area = 0;
  plist2 = elemSurface.at (surfaceNode); 
  for( face=plist2.begin();face!=plist2.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos da face
   plist3 = neighbourFaceVert.at (*face);
   vert=plist3.begin();
   real P1x = X->Get(*vert);
   real P1y = Y->Get(*vert);
   real P1z = Z->Get(*vert);
   ++vert;
   real P2x = X->Get(*vert);
   real P2y = Y->Get(*vert);
   real P2z = Z->Get(*vert);
   // comprimento das arestas a, b e c para aplicacao da formula de
   // Heron
   real a = sqrt( (P1x-P0x)*(P1x-P0x)+(P1y-P0y)*(P1y-P0y)+(P1z-P0z)*(P1z-P0z) );
   real b = sqrt( (P2x-P0x)*(P2x-P0x)+(P2y-P0y)*(P2y-P0y)+(P2z-P0z)*(P2z-P0z) );
   real c = sqrt( (P2x-P1x)*(P2x-P1x)+(P2y-P1y)*(P2y-P1y)+(P2z-P1z)*(P2z-P1z) );
   real s = (a+b+c)/2; //semiperimeter
   area = area + sqrt( s*(s-a)*(s-b)*(s-c) );
  }
  real force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );
  real L = sqrt(area);
  real pressure = force/L;
  cout << surfaceNode << " " << area << " " << force << " " << pressure << endl;
  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa2

// calculo do Kappa geometricamente. Utiliza neighbourVert
clVector Interface3D::computeKappa3()
{
 setCloser();
 int surfaceNode;
 real force,sumForce,sumArea,sumLength,fx,fy,fz;
 listElem plist,plist2;
 list<int>::iterator face,vert;

 // loop sobre todos os nos da superficie 
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);

  int c1 = 0;
  fx = 0;
  fy = 0;
  fz = 0;
  force = 0;
  sumForce = 0;
  sumArea = 0;
  sumLength = 0;

  // loop sobre os vizinhos do vertice i
  plist = elemSurface.at (surfaceNode); 
  for( face=plist.begin();face!=plist.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos do triangulo)
   plist2 = neighbourFaceVert.at (*face);
   vert=plist2.begin();
   int v1 = *vert;++vert;
   int v2 = *vert;
   vert=plist2.end();

   real P1x = X->Get(v1);
   real P1y = Y->Get(v1);
   real P1z = Z->Get(v1);
   real P2x = X->Get(v2);
   real P2y = Y->Get(v2);
   real P2z = Z->Get(v2);

   // ponto medio aresta 01
   real Pm01x = P0x + (P1x-P0x)/2.0; 
   real Pm01y = P0y + (P1y-P0y)/2.0; 
   real Pm01z = P0z + (P1z-P0z)/2.0; 

   // ponto medio aresta 02
   real Pm02x = P0x + (P2x-P0x)/2.0; 
   real Pm02y = P0y + (P2y-P0y)/2.0; 
   real Pm02z = P0z + (P2z-P0z)/2.0; 

   // distance do ponto 0 ate metade do segmento 01
   real a = sqrt( (P0x-Pm01x)*(P0x-Pm01x)+
	              (P0y-Pm01y)*(P0y-Pm01y)+
				  (P0z-Pm01z)*(P0z-Pm01z) );

   // distance do ponto 0 ate metade do segmento 02
   real b = sqrt( (P0x-Pm02x)*(P0x-Pm02x)+
	              (P0y-Pm02y)*(P0y-Pm02y)+
			      (P0z-Pm02z)*(P0z-Pm02z) );

   // distance da metade do segmento 01 ate metade do segmento 02
   real c = sqrt( (Pm02x-Pm01x)*(Pm02x-Pm01x)+
	              (Pm02y-Pm01y)*(Pm02y-Pm01y)+
			      (Pm02z-Pm01z)*(Pm02z-Pm01z) );

   // calculando area do trianglo 0-01medio-02medio
   real s = (a+b+c)/2; //semiperimeter
   real area = sqrt( s*(s-a)*(s-b)*(s-c) );

   // vetores unitarios
   real x1Unit = (Pm01x-P0x)/a;
   real y1Unit = (Pm01y-P0y)/a;
   real z1Unit = (Pm01z-P0z)/a;

   real x2Unit = (Pm02x-P0x)/b;
   real y2Unit = (Pm02y-P0y)/b;
   real z2Unit = (Pm02z-P0z)/b;

   real xRetaUnit = (Pm02x-Pm01x)/c;
   real yRetaUnit = (Pm02y-Pm01y)/c;
   real zRetaUnit = (Pm02z-Pm01z)/c;

   // soma dos vetores 1Unit + 2Unit = resultante
   real xUnit = x1Unit+x2Unit;
   real yUnit = y1Unit+y2Unit;
   real zUnit = z1Unit+z2Unit;

   // produto escalar --> projecao do vetor Unit no segmento de reta
   // | Unit.RetaUnit | . RetaUnit
   // resultado = vetor tangente a reta situado na superficie
   real prod = xUnit*xRetaUnit + yUnit*yRetaUnit + zUnit*zRetaUnit;
   real xTang = xRetaUnit*prod;
   real yTang = yRetaUnit*prod;
   real zTang = zRetaUnit*prod;

   // subtraindo vetor tangente do vetor unitario para encontrar as
   // coordenadas do vetor normal situada na superficie
   real xNormal = xUnit - xTang;
   real yNormal = yUnit - yTang;
   real zNormal = zUnit - zTang;

   real len = sqrt( (xNormal*xNormal)+(yNormal*yNormal)+(zNormal*zNormal) ); 

   // Unitario do vetor resultante situado no plano do triangulo
   // combinacao linear dos vetores unitarios das arestas do triangulo
   real xNormalUnit = xNormal/len;
   real yNormalUnit = yNormal/len;
   real zNormalUnit = zNormal/len;

   // normal integrada na distancia (MOD) dos 2 vertices medianos
   // force = resultante das componentes * tamanho da aresta que sera
   // usada como referencia no calculo da area do triangulo
   fx += xNormalUnit*c;
   fy += yNormalUnit*c;
   fz += zNormalUnit*c;

   sumArea += area;
   sumLength += c;

//--------------------------------------------------
//    if( surfaceNode == 2185 )
//    {
// 	cout << "Triangulo: ------------------------- " << c1 << endl;
// 	cout << "v1 = " << v1 << " " << "v2 = " << v2 << endl;
// 	cout << "Unit = " << xUnit << " " << yUnit << " " << zUnit << endl;
// 	cout << "RetaUnit = " << xRetaUnit << " " 
// 	                      << yRetaUnit << " " 
// 						  << zRetaUnit << endl;
// 	cout << "c = " << c << endl;
// 	cout << "Normal = " << xNormal << " " 
// 	                    << yNormal << " " 
// 						<< zNormal << endl;
// 	cout << "area = " << area << endl;
// 	cout << "force = " << force << endl;
// 	cout << "sumArea = " << sumArea << endl;
// 	cout << "fx = " << fx << endl;
// 	cout << "fy = " << fy << endl;
// 	cout << "fz = " << fz << endl;
// 	//cout << "sumForce = " << sumForce << endl;
// 	//cout << "pressure = " << sumForce/sumArea << endl;
// 	c1++;
//    }
//-------------------------------------------------- 

  }
  force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );
  real pressure = force/sumArea;
  //real pressure = force;

//--------------------------------------------------
//   if( surfaceNode == 1268 )
//   {
//   cout << surfaceNode << endl;
//   cout << "  " << "force = " << force << endl; 
//   cout << "  " << "fx = " << fx << endl; 
//   cout << "  " << "fy = " << fy << endl;
//   cout << "  " << "fz = " << fz << endl;
//   cout << "  " << "sumArea = " << sumArea << endl;
//   cout << "  " << "pressure = " << pressure << endl;
//   }
//-------------------------------------------------- 

  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa3

// calculo do Kappa geometricamente. Utiliza neighbourVert
clVector Interface3D::computeKappa4()
{
 setCloser();
 int surfaceNode;
 real force,sumForce,sumArea,sumLength,fx,fy,fz;
 listElem plist,plist2;
 list<int>::iterator face,vert;

 // loop sobre todos os nos da superficie 
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  real P0x = X->Get(surfaceNode);
  real P0y = Y->Get(surfaceNode);
  real P0z = Z->Get(surfaceNode);

  int c1 = 0;
  fx = 0;
  fy = 0;
  fz = 0;
  force = 0;
  sumForce = 0;
  sumArea = 0;
  sumLength = 0;

  // loop sobre os vizinhos do vertice i
  plist = elemSurface.at (surfaceNode); 
  for( face=plist.begin();face!=plist.end();++face )
  {
   // 3D: 2 pontos pois a face em 3D pertencente a superficie contem 
   // 3 pontos (P0 - surfaceNode, P1 e P2 que sao pontos do triangulo)
   plist2 = neighbourFaceVert.at (*face);
   vert=plist2.begin();
   int v1 = *vert;++vert;
   int v2 = *vert;
   vert=plist2.end();

   real P1x = X->Get(v1);
   real P1y = Y->Get(v1);
   real P1z = Z->Get(v1);
   real P2x = X->Get(v2);
   real P2y = Y->Get(v2);
   real P2z = Z->Get(v2);

   // ponto medio aresta 01
   real Pm01x = P0x + (P1x-P0x)/2.0; 
   real Pm01y = P0y + (P1y-P0y)/2.0; 
   real Pm01z = P0z + (P1z-P0z)/2.0; 

   // ponto medio aresta 02
   real Pm02x = P0x + (P2x-P0x)/2.0; 
   real Pm02y = P0y + (P2y-P0y)/2.0; 
   real Pm02z = P0z + (P2z-P0z)/2.0; 

   // distance do ponto 0 ate metade do segmento 01
   real a = sqrt( (P0x-P1x)*(P0x-P1x)+
	              (P0y-P1y)*(P0y-P1y)+
				  (P0z-P1z)*(P0z-P1z) );

   // distance do ponto 0 ate metade do segmento 02
   real b = sqrt( (P0x-P2x)*(P0x-P2x)+
	              (P0y-P2y)*(P0y-P2y)+
			      (P0z-P2z)*(P0z-P2z) );

   // distance da metade do segmento 01 ate metade do segmento 02
   real c = sqrt( (P2x-P1x)*(P2x-P1x)+
	              (P2y-P1y)*(P2y-P1y)+
			      (P2z-P1z)*(P2z-P1z) );

   // distance da metade do segmento 01 ate metade do segmento 02
   real d = sqrt( (Pm02x-Pm01x)*(Pm02x-Pm01x)+
	              (Pm02y-Pm01y)*(Pm02y-Pm01y)+
			      (Pm02z-Pm01z)*(Pm02z-Pm01z) );

   // calculando area do trianglo 0-01medio-02medio
   real s = (a+b+c)/2; //semiperimeter
   real area = sqrt( s*(s-a)*(s-b)*(s-c) );

   // vetores unitarios
   real x1Unit = (P1x-P0x)/a;
   real y1Unit = (P1y-P0y)/a;
   real z1Unit = (P1z-P0z)/a;

   real x2Unit = (P2x-P0x)/b;
   real y2Unit = (P2y-P0y)/b;
   real z2Unit = (P2z-P0z)/b;

   real xRetaUnit = (P2x-P1x)/c;
   real yRetaUnit = (P2y-P1y)/c;
   real zRetaUnit = (P2z-P1z)/c;

   // soma dos vetores 1Unit + 2Unit = resultante
   real xUnit = x1Unit+x2Unit;
   real yUnit = y1Unit+y2Unit;
   real zUnit = z1Unit+z2Unit;

   // produto escalar --> projecao do vetor Unit no segmento de reta
   // | Unit.RetaUnit | . RetaUnit
   // resultado = vetor tangente a reta situado na superficie
   real prod = xUnit*xRetaUnit + yUnit*yRetaUnit + zUnit*zRetaUnit;
   real xTang = xRetaUnit*prod;
   real yTang = yRetaUnit*prod;
   real zTang = zRetaUnit*prod;

   // subtraindo vetor tangente do vetor unitario para encontrar as
   // coordenadas do vetor normal situada na superficie
   real xNormal = xUnit - xTang;
   real yNormal = yUnit - yTang;
   real zNormal = zUnit - zTang;

   real len = sqrt( (xNormal*xNormal)+(yNormal*yNormal)+(zNormal*zNormal) ); 

   // Unitario do vetor resultante situado no plano do triangulo
   // combinacao linear dos vetores unitarios das arestas do triangulo
   real xNormalUnit = xNormal/len;
   real yNormalUnit = yNormal/len;
   real zNormalUnit = zNormal/len;

   // tamanho do comprimento da aresta do triangulo que eh metade do
   // elemento
   real e = sqrt( (c*c)/2.0 );

   // normal integrada na distancia (MOD) dos 2 vertices medianos
   // force = resultante das componentes * tamanho da aresta que sera
   // usada como referencia no calculo da area do triangulo
   fx += xNormalUnit*e;
   fy += yNormalUnit*e;
   fz += zNormalUnit*e;

   sumArea += area/2.0;
   sumLength += c;
   force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );

//--------------------------------------------------
//    if( surfaceNode == 173 )
//    {
// 	cout << "Triangulo: ------------------------- " << c1 << endl;
// 	cout << "v1 = " << v1 << " " << "v2 = " << v2 << endl;
// 	cout << "Unit = " << xUnit << " " << yUnit << " " << zUnit << endl;
// 	cout << "RetaUnit = " << xRetaUnit << " " 
// 	                      << yRetaUnit << " " 
// 						  << zRetaUnit << endl;
// 	cout << "c = " << c << endl;
// 	cout << "Normal = " << xNormal << " " 
// 	                    << yNormal << " " 
// 						<< zNormal << endl;
// 	cout << "area = " << area << endl;
// 	cout << "force = " << force << endl;
// 	cout << "sumArea = " << sumArea << endl;
// 	cout << "fx = " << fx << endl;
// 	cout << "fy = " << fy << endl;
// 	cout << "fz = " << fz << endl;
// 	//cout << "sumForce = " << sumForce << endl;
// 	//cout << "pressure = " << sumForce/sumArea << endl;
// 	c1++;
//    }
//-------------------------------------------------- 

  }
  real pressure = force/sumArea;
  //real pressure = force;

//--------------------------------------------------
//   if( surfaceNode == 1268 )
//   {
//   cout << surfaceNode << endl;
//   cout << "  " << "force = " << force << endl; 
//   cout << "  " << "fx = " << fx << endl; 
//   cout << "  " << "fy = " << fy << endl;
//   cout << "  " << "fz = " << fz << endl;
//   cout << "  " << "sumArea = " << sumArea << endl;
//   cout << "  " << "pressure = " << pressure << endl;
//   }
//-------------------------------------------------- 

  distance.Set( surfaceNode,pressure );
 }
 return distance;

} // fecha metodo computeKappa3


clVector Interface3D::smoothing(clMatrix &_AcTilde,clVector &_b1cTilde)
{
 clVector cTilde(_b1cTilde.Dim());
 cout << " ----> calculando escalar ----------- " << endl;
 solverC->solve(1E-08,_AcTilde,cTilde,_b1cTilde);
 cout << " ------------------------------------ " << endl;
 return cTilde;

}

void Interface3D::setCloser()
{
 int aux;

 // xSurface representa os valores de X dos nos da interface
 // ySurface representa os valores de Y dos nos da interface
 // zSurface representa os valores de Z dos nos da interface
 xSurface.Dim( surface->Dim() );
 ySurface.Dim( surface->Dim() );
 zSurface.Dim( surface->Dim() );
 for( int i=0;i<surface->Dim();i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get( aux ));
  ySurface.Set(i,Y->Get( aux ));
  zSurface.Set(i,Z->Get( aux ));
 }

 // closer=surface(dsearchn(X(surface),Y(surface),X,Y));
 // esta funcao retorna o noh da interface (surface) mais 
 // proximo de cada noh da malha (vertices)
 closer = dsearchn(xSurface,ySurface,zSurface,*X,*Y,*Z);
 for( int i=0;i<closer.Dim();i++ )
 {
  aux = closer.Get(i);
  closer.Set(i,surface->Get(aux)); // alterando os valores de closer(i)
  xCloser.Set(i,X->Get( closer.Get(i) ));
  yCloser.Set(i,Y->Get( closer.Get(i) ));
  zCloser.Set(i,Z->Get( closer.Get(i) ));
 }
}

clDMatrix Interface3D::setKappaSurface(clVector &_kappaAux)
{
 int aux;
 clDMatrix kappa(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  aux = closer.Get(i);
  kappa.Set(i,_kappaAux.Get(aux));
  kappa.Set(i+numNodes,_kappaAux.Get(aux));
  kappa.Set(i+2*numNodes,_kappaAux.Get(aux));
 }
 return kappa;
}

clVector Interface3D::dsearchn(clVector _X,clVector _Y,clVector _Z,
                               clVector &_XI,clVector &_YI,clVector &_ZI)
{
 clVector vmin(_XI.Dim());
 vmin.SetAll(0);
 real dist, dmin;

 // loop sobre todos os pontos a serem procurados
 for( int i=0;i<_XI.Dim();i++ )
 {

  // distancia de um ponto qualquer x(0),y(o) para o ponto informado xi,yi
  dmin = (_XI.Get(i)-_X.Get(0))*(_XI.Get(i)-_X.Get(0))+
         (_YI.Get(i)-_Y.Get(0))*(_YI.Get(i)-_Y.Get(0))+
         (_ZI.Get(i)-_Z.Get(0))*(_ZI.Get(i)-_Z.Get(0));

  // loop sobre todos os pontos da malha (menos o ponto inicial) para
  // saber qual o ponto apresenta a menor distancia com relacao ao ponto
  // informado xi,yi
  for (int j=1;j<_X.Dim();j++)
  {
   //calculo da distancia de outro ponto da malha x(i),y(i)
   dist = (_XI.Get(i)-_X.Get(j))*(_XI.Get(i)-_X.Get(j))+
	      (_YI.Get(i)-_Y.Get(j))*(_YI.Get(i)-_Y.Get(j))+
	      (_ZI.Get(i)-_Z.Get(j))*(_ZI.Get(i)-_Z.Get(j));
   if (dist<dmin)
   {
	dmin=dist;
	vmin.Set(i,j); 
   }
  }
 }

 // retorna o indice do ponto da malha
 return vmin;
}

void Interface3D::setSolverSmooth(Solver *s){ solverC = s; }

// rotina nao funciona pois cada vertice da superfice tem o mapeamento
// das faces vizinhas e com isso as faces se repetem ao longo da
// impressao das celulas, impossibilitando a construcao da malha da
// bolha
void Interface3D::saveSurfaceVTK()
{
 int v0,v1,v2,surfaceNode;
 listElem plist,plist2;
 list<int>::iterator face,vert;

 const char* filename = "vtk/surface.vtk";
 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Bubble Surface 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "POINTS " << surface->Dim() << " float" << endl;

 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  vtkFile << X->Get(surfaceNode) << " " 
          << Y->Get(surfaceNode) << " " 
		  << Z->Get(surfaceNode) << endl;
 }
 vtkFile << endl;
 
 vtkFile << "CELLS " << neighbourFaceVert.size() << " " 
         << 4*neighbourFaceVert.size() << endl;

 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  plist = elemSurface.at (surfaceNode); 
  for( face=plist.begin();face!=plist.end();++face )
  {
   v0 = surfaceNode;
   plist2 = neighbourFaceVert.at (*face);
   vert=plist2.begin();
   v1 = *vert;
   ++vert;
   v2 = *vert;
   vtkFile << "3 " << v0 << " " << v1 << " " << v2 << endl;
  }
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << neighbourFaceVert.size() << endl;
 for( int i=0;i<neighbourFaceVert.size();i++ )
  vtkFile << "5 ";

 vtkFile << endl;

 vtkFile <<  "SURFACE NODES " << endl;
 for( int i=0;i<surface->Dim();i++ )
 {
  surfaceNode = surface->Get(i);
  vtkFile << surfaceNode << endl;
 }
 vtkFile << endl;
 vtkFile.close();

 cout << "malha da superfie da bolha gravada em VTK" << endl;

} // fecha metodo saveSurfaceVTK

clVector Interface3D::getCloser(){ return closer; }
