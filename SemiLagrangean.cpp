// =================================================================== // 
// this is file SemiLagrangean.cpp, created at 20-Ago-2007             //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "SemiLagrangean.h"

SemiLagrangean::SemiLagrangean( Model3D &_m )
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
 idbcc = m->getIdbcc();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 neighbourElem = m->getNeighbourElem(); 
 oFace = m->getOFace();

 convLin.Dim(numVerts,numVerts);
 convQuad.Dim(numNodes,numNodes);
}

SemiLagrangean::SemiLagrangean(Model3D &_m,clVector &_uSol,
                                           clVector &_vSol,
										   clVector &_wSol,
										   clVector &_velU,
										   clVector &_velV,
										   clVector &_velW,
										   clVector &_cSol)
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
 idbcc = m->getIdbcc();
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 uSol = _uSol;
 vSol = _vSol;
 wSol = _wSol;
 velU = _velU;
 velV = _velV;
 velW = _velW;
 cSol = _cSol;
 neighbourElem = m->getNeighbourElem(); 
 oFace = m->getOFace();

 convLin.Dim(numVerts,numVerts);
 convQuad.Dim(numNodes,numNodes);
}

void SemiLagrangean::compute(real dt)
{
 if( NUMGLEU == 5 )
 {
  getDepartElem(dt); // procura de elemento

  clVector uVert(numVerts);
  clVector vVert(numVerts);
  clVector wVert(numVerts);
  uSol.CopyTo(0,uVert);
  vSol.CopyTo(0,vVert);
  wSol.CopyTo(0,wVert);

  uParticle = convLin*uVert; 
  vParticle = convLin*vVert; 
  wParticle = convLin*wVert;
  cParticle = convLin*cSol;

  clVector zerosConv(numNodes-numVerts);
  uParticle.Append(zerosConv);
  vParticle.Append(zerosConv);
  wParticle.Append(zerosConv);
  setCentroid();
 }
 else 
 {
  getDepartElemQuad(dt); // procura de elemento

  uParticle = convQuad*uSol; 
  vParticle = convQuad*vSol; 
  wParticle = convQuad*wSol;
  cParticle = convLin*cSol;
 }

 setBC();
} // fim do metodo compute

void SemiLagrangean::computeFreeSurface(real dt)
{
 getDepartElem2(dt); // procura de elemento

 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 uSol.CopyTo(0,uVert);
 vSol.CopyTo(0,vVert);
 wSol.CopyTo(0,wVert);

 uParticle = convLin*uVert; 
 vParticle = convLin*vVert; 
 wParticle = convLin*wVert;
 cParticle = convLin*cSol;

 clVector zerosConv(numNodes-numVerts);
 uParticle.Append(zerosConv);
 vParticle.Append(zerosConv);
 wParticle.Append(zerosConv);

 setBC();
 setCentroid();

} // fim do metodo computeFreeSurface

void SemiLagrangean::setCentroid()
{
 int v[NUMGLEU];
 real aux;

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  aux = ( uParticle.Get(v[0])+
	      uParticle.Get(v[1])+
		  uParticle.Get(v[2])+
		  uParticle.Get(v[3]) )/4.0;
  uParticle.Set(v[4],aux);

  aux = ( vParticle.Get(v[0])+
          vParticle.Get(v[1])+
          vParticle.Get(v[2])+
	 	  vParticle.Get(v[3]) )/4.0;
  vParticle.Set(v[4],aux);

  aux = ( wParticle.Get(v[0])+
	      wParticle.Get(v[1])+
          wParticle.Get(v[2])+
		  wParticle.Get(v[3]) )/4.0;
  wParticle.Set(v[4],aux);
 }

}// fim do metodo compute -> setCentroid

//--------------------------------------------------
// void SemiLagrangean::setQuad()
// {
//  int v[NUMGLEU];
//  real aux;
// 
//  for( int mele=0;mele<numElems;mele++ )
//  {
//   for( int n=0;n<NUMGLEU;n++ )
//    v[n] = (int) IEN->Get(mele,n);
// 
//   // v[4] (v[0]-v[1])
//   aux = ( uParticle.Get(v[0])+
// 		  uParticle.Get(v[1]) )*0.5;
//   uParticle.Set(v[4],aux);
// 
//   aux = ( vParticle.Get(v[0])+
// 		  vParticle.Get(v[1]) )*0.5;
//   vParticle.Set(v[4],aux);
// 
//   aux = ( wParticle.Get(v[0])+
// 		  wParticle.Get(v[1]) )*0.5;
//   wParticle.Set(v[4],aux);
// 
//   // v[5] (v[1]-v[2])
//   aux = ( uParticle.Get(v[1])+
// 		  uParticle.Get(v[2]) )*0.5;
//   uParticle.Set(v[5],aux);
// 
//   aux = ( vParticle.Get(v[1])+
// 		  vParticle.Get(v[2]) )*0.5;
//   vParticle.Set(v[5],aux);
// 
//   aux = ( wParticle.Get(v[1])+
// 		  wParticle.Get(v[2]) )*0.5;
//   wParticle.Set(v[5],aux);
// 
//   // v[6] (v[0]-v[2])
//   aux = ( uParticle.Get(v[0])+
// 		  uParticle.Get(v[2]) )*0.5;
//   uParticle.Set(v[6],aux);
// 
//   aux = ( vParticle.Get(v[0])+
// 		  vParticle.Get(v[2]) )*0.5;
//   vParticle.Set(v[6],aux);
// 
//   aux = ( wParticle.Get(v[0])+
// 		  wParticle.Get(v[2]) )*0.5;
//   wParticle.Set(v[6],aux);
// 
//   // v[7] (v[0]-v[3])
//   aux = ( uParticle.Get(v[0])+
// 		  uParticle.Get(v[3]) )*0.5;
//   uParticle.Set(v[7],aux);
// 
//   aux = ( vParticle.Get(v[0])+
// 		  vParticle.Get(v[3]) )*0.5;
//   vParticle.Set(v[7],aux);
// 
//   aux = ( wParticle.Get(v[0])+
// 		  wParticle.Get(v[3]) )*0.5;
//   wParticle.Set(v[7],aux);
// 
//   // v[8] (v[1]-v[3])
//   aux = ( uParticle.Get(v[1])+
// 		  uParticle.Get(v[3]) )*0.5;
//   uParticle.Set(v[8],aux);
// 
//   aux = ( vParticle.Get(v[1])+
// 		  vParticle.Get(v[3]) )*0.5;
//   vParticle.Set(v[8],aux);
// 
//   aux = ( wParticle.Get(v[1])+
// 		  wParticle.Get(v[3]) )*0.5;
//   wParticle.Set(v[8],aux);
// 
//   // v[0]0 (v[2]-v[3])
//   aux = ( uParticle.Get(v[2])+
// 		  uParticle.Get(v[3]) )*0.5;
//   uParticle.Set(v[9],aux);
// 
//   aux = ( vParticle.Get(v[2])+
// 		  vParticle.Get(v[3]) )*0.5;
//   vParticle.Set(v[9],aux);
// 
//   aux = ( wParticle.Get(v[2])+
// 		  wParticle.Get(v[3]) )*0.5;
//   wParticle.Set(v[9],aux);
//  }
// 
// }// fim do metodo compute -> setQuad
//-------------------------------------------------- 

void SemiLagrangean::getDepartElem(real dt)
{
 real xP,yP,zP;

 clVector xParticle = *X - velU*dt; 
 clVector yParticle = *Y - velV*dt;
 clVector zParticle = *Z - velW*dt;

 list<int> plist;
 list<int>::iterator mele;

 for (int ii = 0; ii < numVerts; ii++)
 {
  plist = neighbourElem->at(ii);
  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = xParticle.Get(ii);
  yP = yParticle.Get(ii);
  zP = zParticle.Get(ii);

  jumpToElem(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElem

void SemiLagrangean::getDepartElemQuad(real dt)
{
 real xP,yP,zP;

 clVector xParticle = *X - velU*dt; 
 clVector yParticle = *Y - velV*dt;
 clVector zParticle = *Z - velW*dt;

 list<int> plist;
 list<int>::iterator mele;

 for (int ii = 0; ii < numNodes; ii++)
 {
  plist = neighbourElem->at(ii);
  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = xParticle.Get(ii);
  yP = yParticle.Get(ii);
  zP = zParticle.Get(ii);

  jumpToElemQuad(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElemQuad

// rotina usada em simulacao FREE SURFACE
void SemiLagrangean::getDepartElem2(real dt)
{
 real xP,yP,zP;

 //clVector xParticle = *X; 
 //clVector yParticle = *Y - vSol*dt;
 //clVector zParticle = *Z - wSol*dt;

 //clVector xParticle = *X - uSol*dt; 
 //clVector yParticle = *Y;
 //clVector zParticle = *Z - wSol*dt;

 clVector xParticle = *X - uSol*dt; 
 clVector yParticle = *Y - vSol*dt;
 clVector zParticle = *Z;

 list<int> plist;
 list<int>::iterator mele;

 for (int ii = 0; ii < numVerts; ii++)
 {
  plist = neighbourElem->at(ii);
  mele=plist.begin(); // pega o primeiro elemento da lista, seja la qual for

  xP = xParticle.Get(ii);
  yP = yParticle.Get(ii);
  zP = zParticle.Get(ii);

  //cout << "vertice de origem = " << ii << endl;
  jumpToElem(*mele,ii,xP,yP,zP);
 } 
} // fim do metodo compute -> getDepartElem2

// verifica se esta dentro do elemento, se estiver calcula os
// coeficientes de interpolacao (convLin). Se nao encontrar e tiver
// vizinho, pula para o novo vizinho. Se nao, interpola na face
void SemiLagrangean::jumpToElem(int destElem,int iiVert,real R2X,real R2Y,real R2Z)
{
 real l1,l2,l3,l4;
 real Bl1,Bl2,Bl3;
 int v[4],v1,v2,v3,v4,vjump;
 int ib1=0;
 int ib2=0;
 int ib3=0;
//--------------------------------------------------
//  cout << iiVert << " " << R2X << " " 
//                        << R2Y << " " 
// 					   << R2Z << endl;
//-------------------------------------------------- 
 if( testElement(destElem,iiVert,R2X,R2Y,R2Z,&l1,&l2,&l3,&l4) )
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);
  convLin.Set(iiVert,v[0],l1);
  convLin.Set(iiVert,v[1],l2);
  convLin.Set(iiVert,v[2],l3);
  convLin.Set(iiVert,v[3],l4);

  //cout << " Ponto localizado " << iiVert << endl;
 } 
 else 
 {
  v[0] = v1 = (int) IEN->Get(destElem,0);
  v[1] = v2 = (int) IEN->Get(destElem,1);
  v[2] = v3 = (int) IEN->Get(destElem,2);
  v[3] = v4 = (int) IEN->Get(destElem,3);

  if((l1<=l2) && (l1<=l3) && (l1<=l4)){ vjump=0; ib1=v2; ib2=v3; ib3=v4; };
  if((l2<=l1) && (l2<=l3) && (l2<=l4)){ vjump=1; ib1=v1; ib2=v3; ib3=v4; };
  if((l3<=l1) && (l3<=l2) && (l3<=l4)){ vjump=2; ib1=v1; ib2=v2; ib3=v4; };
  if((l4<=l1) && (l4<=l2) && (l4<=l3)){ vjump=3; ib1=v1; ib2=v2; ib3=v3; };

  if ( oFace->Get(destElem,vjump)!=-1)
  {
   jumpToElem( (int) oFace->Get(destElem,vjump), iiVert,R2X,R2Y,R2Z );
  }
  else
  {
   computeIntercept( iiVert,R2X,R2Y,R2Z,ib1,ib2,ib3,&Bl1,&Bl2,&Bl3 );
   convLin.Set( iiVert,ib1,Bl1 );
   convLin.Set( iiVert,ib2,Bl2 );  
   convLin.Set( iiVert,ib3,Bl3 );  
  }
 }
}

// verifica se esta dentro do elemento, se estiver calcula os
// coeficientes de interpolacao (convLin). Se nao encontrar e tiver
// vizinho, pula para o novo vizinho. Se nao, interpola na face
void SemiLagrangean::jumpToElemQuad(int destElem,int iiVert,real R2X,real R2Y,real R2Z)
{
 real l1,l2,l3,l4;
 real N1,N2,N3,N4,N5,N6,N7,N8,N9,N10;
 real Bl1,Bl2,Bl3;
 real BN1,BN2,BN3,BN4,BN5,BN6;
 int v[NUMGLEU],vjump;
 int ib1=0;
 int ib2=0;
 int ib3=0;
 int ib4=0;
 int ib5=0;
 int ib6=0;

 for( int n=0;n<NUMGLEU;n++ )
  v[n] = (int) IEN->Get(destElem,n);

 if( testElement(destElem,iiVert,R2X,R2Y,R2Z,&l1,&l2,&l3,&l4) )
 {
  N1 = (2*l1-1.0)*l1;
  N2  = (2*l2-1.0)*l2;
  N3  = (2*l3-1.0)*l3;
  N4  = (2*l4-1.0)*l4; 
  N5  = 4*l1*l2;
  N6  = 4*l2*l3;
  N7  = 4*l1*l3;
  N8  = 4*l1*l4;
  N9  = 4*l2*l4;
  N10 = 4*l3*l4;

  if( iiVert < numVerts )
  {
   convLin.Set(iiVert,v[0],l1);
   convLin.Set(iiVert,v[1],l2);
   convLin.Set(iiVert,v[2],l3);
   convLin.Set(iiVert,v[3],l4);
  }

  convQuad.Set(iiVert,v[0],N1);
  convQuad.Set(iiVert,v[1],N2);
  convQuad.Set(iiVert,v[2],N3);
  convQuad.Set(iiVert,v[3],N4);
  convQuad.Set(iiVert,v[4],N5);
  convQuad.Set(iiVert,v[5],N6);
  convQuad.Set(iiVert,v[6],N7);
  convQuad.Set(iiVert,v[7],N8);
  convQuad.Set(iiVert,v[8],N9);
  convQuad.Set(iiVert,v[9],N10);

  //cout << " Ponto localizado " << iiVert << endl;
 } 
 else 
 {
  if((l1<=l2) && (l1<=l3) && (l1<=l4))
  { 
   vjump=0; ib1=v[1]; ib2=v[2]; ib3=v[3]; ib4=v[5]; ib5=v[9]; ib6=v[8]; 
  };
  if((l2<=l1) && (l2<=l3) && (l2<=l4))
  { 
   vjump=1; ib1=v[0]; ib2=v[2]; ib3=v[3]; ib4=v[6]; ib5=v[9]; ib6=v[7]; 
  };
  if((l3<=l1) && (l3<=l2) && (l3<=l4))
  { 
   vjump=2; ib1=v[0]; ib2=v[1]; ib3=v[3]; ib4=v[4]; ib5=v[8]; ib6=v[7]; 
  };
  if((l4<=l1) && (l4<=l2) && (l4<=l3))
  { 
   vjump=3; ib1=v[0]; ib2=v[1]; ib3=v[2]; ib4=v[4]; ib5=v[5]; ib6=v[6]; 
  };

  if ( oFace->Get(destElem,vjump)!=-1)
  {
   jumpToElemQuad( (int) oFace->Get(destElem,vjump), iiVert,R2X,R2Y,R2Z );
  }
  else
  {
   computeIntercept( iiVert,R2X,R2Y,R2Z,ib1,ib2,ib3,&Bl1,&Bl2,&Bl3 );

   BN1 = (2*Bl1-1.0)*Bl1;
   BN2 = (2*Bl2-1.0)*Bl2;
   BN3 = (2*Bl3-1.0)*Bl3;
   BN4 = 4*Bl1*Bl2;
   BN5 = 4*Bl2*Bl3;
   BN6 = 4*Bl1*Bl3;

   if( iiVert < numVerts )
   {
	convLin.Set( iiVert,ib1,Bl1 );
	convLin.Set( iiVert,ib2,Bl2 );  
	convLin.Set( iiVert,ib3,Bl3 );  
   }

   convQuad.Set( iiVert,ib1,BN1 );
   convQuad.Set( iiVert,ib2,BN2 );  
   convQuad.Set( iiVert,ib3,BN3 );  
   convQuad.Set( iiVert,ib4,BN4 );  
   convQuad.Set( iiVert,ib5,BN5 );  
   convQuad.Set( iiVert,ib6,BN6 );  
  }
 }
}

// calcula ponto de intersecao resolvendo um sistema linear 3x3 atraves
// da regra de Cramer
//
// | a1 b1 c1 |   | x1 | = | d1 |
// | a2 b2 c2 | . | x2 | = | d2 |
// | a3 b3 c3 |   | x3 | = | d3 |
//
//       detx          dety
// x1 = ------   x2 = ------   x3 = 1.0-x1-x2
//       det            det
void SemiLagrangean::computeIntercept(int ii,real R2X,real R2Y,real R2Z,
  int ib1,int ib2,int ib3,real *Bl1,real *Bl2,real *Bl3)
{
 real R1X = X->Get(ii); real R1Y = Y->Get(ii); real R1Z = Z->Get(ii);

 real B1X = X->Get(ib1); real B1Y = Y->Get(ib1); real B1Z = Z->Get(ib1);
 real B2X = X->Get(ib2); real B2Y = Y->Get(ib2); real B2Z = Z->Get(ib2);
 real B3X = X->Get(ib3); real B3Y = Y->Get(ib3); real B3Z = Z->Get(ib3);

 real a1 = B1X-B3X; real b1 = B2X-B3X; real c1 = R1X-R2X; real d1 = R1X-B3X;
 real a2 = B1Y-B3Y; real b2 = B2Y-B3Y; real c2 = R1Y-R2Y; real d2 = R1Y-B3Y;
 real a3 = B1Z-B3Z; real b3 = B2Z-B3Z; real c3 = R1Z-R2Z; real d3 = R1Z-B3Z;

 real det  = (a1*b2*c3)+(a3*b1*c2)+(a2*b3*c1)-
             (a3*b2*c1)-(a1*b3*c2)-(a2*b1*c3);
 real detx = (d1*b2*c3)+(d3*b1*c2)+(d2*b3*c1)-
             (d3*b2*c1)-(d1*b3*c2)-(d2*b1*c3);
 real dety = (a1*d2*c3)+(a3*d1*c2)+(a2*d3*c1)-
             (a3*d2*c1)-(a1*d3*c2)-(a2*d1*c3);

 real x1 = detx/det;
 real x2 = dety/det;

 *Bl1 = x1;
 *Bl2 = x2;
 *Bl3 = 1.0-*Bl1-*Bl2;
}

bool SemiLagrangean::testElement(int mele,int ii,real xP,real yP,real zP, real *l1,real *l2,real *l3,real *l4)
{
 int v[NUMGLE];
 real V,V1,V2,V3;
 real EPSlocal = 10e-6;

  for( int n=0;n<NUMGLE;n++ )
   v[n] = (int) IEN->Get(mele,n);

 V = (1.0/6.0) * (+1*( (X->Get(v[1])*Y->Get(v[2])*Z->Get(v[3])) 
	                    +(Y->Get(v[1])*Z->Get(v[2])*X->Get(v[3])) 
	                    +(Z->Get(v[1])*X->Get(v[2])*Y->Get(v[3])) 
	                    -(Y->Get(v[1])*X->Get(v[2])*Z->Get(v[3])) 
	                    -(X->Get(v[1])*Z->Get(v[2])*Y->Get(v[3])) 
	                    -(Z->Get(v[1])*Y->Get(v[2])*X->Get(v[3])) )
	        -X->Get(v[0])*( +Y->Get(v[2])*Z->Get(v[3])
		                 +Y->Get(v[1])*Z->Get(v[2]) 
		                 +Z->Get(v[1])*Y->Get(v[3])
		                 -Y->Get(v[1])*Z->Get(v[3])
		 	  		     -Z->Get(v[2])*Y->Get(v[3]) 
			 		     -Z->Get(v[1])*Y->Get(v[2]) )
	        +Y->Get(v[0])*( +X->Get(v[2])*Z->Get(v[3])
	 	 			     +X->Get(v[1])*Z->Get(v[2])
	 				     +Z->Get(v[1])*X->Get(v[3])
		                 -X->Get(v[1])*Z->Get(v[3])
					     -Z->Get(v[2])*X->Get(v[3]) 
					     -Z->Get(v[1])*X->Get(v[2]) )
		    -Z->Get(v[0])*( +X->Get(v[2])*Y->Get(v[3])
			             +X->Get(v[1])*Y->Get(v[2]) 
					     +Y->Get(v[1])*X->Get(v[3])
		                 -X->Get(v[1])*Y->Get(v[3])
				         -Y->Get(v[2])*X->Get(v[3]) 
					     -Y->Get(v[1])*X->Get(v[2]) ) );

 V1 = (1.0/6.0) * (+1*( (X->Get(v[1])*Y->Get(v[2])*Z->Get(v[3])) 
	                     +(Y->Get(v[1])*Z->Get(v[2])*X->Get(v[3])) 
	                     +(Z->Get(v[1])*X->Get(v[2])*Y->Get(v[3])) 
	                     -(Y->Get(v[1])*X->Get(v[2])*Z->Get(v[3])) 
	                     -(X->Get(v[1])*Z->Get(v[2])*Y->Get(v[3])) 
	                     -(Z->Get(v[1])*Y->Get(v[2])*X->Get(v[3])) )
               -xP*( +Y->Get(v[2])*Z->Get(v[3])
		                  +Y->Get(v[1])*Z->Get(v[2]) 
			              +Z->Get(v[1])*Y->Get(v[3])
		                  -Y->Get(v[1])*Z->Get(v[3])
				   	      -Z->Get(v[2])*Y->Get(v[3]) 
					      -Z->Get(v[1])*Y->Get(v[2]) )
	           +yP*( +X->Get(v[2])*Z->Get(v[3])
			 		      +X->Get(v[1])*Z->Get(v[2])
					      +Z->Get(v[1])*X->Get(v[3])
		                  -X->Get(v[1])*Z->Get(v[3])
					      -Z->Get(v[2])*X->Get(v[3]) 
					      -Z->Get(v[1])*X->Get(v[2]) )
		       -zP*( +X->Get(v[2])*Y->Get(v[3])
			              +X->Get(v[1])*Y->Get(v[2]) 
					      +Y->Get(v[1])*X->Get(v[3])
		                  -X->Get(v[1])*Y->Get(v[3])
				          -Y->Get(v[2])*X->Get(v[3]) 
					      -Y->Get(v[1])*X->Get(v[2]) ) );

 V2 = (1.0/6.0) * (+1*( (xP*Y->Get(v[2])*Z->Get(v[3])) 
	                    +(yP*Z->Get(v[2])*X->Get(v[3])) 
	                    +(zP*X->Get(v[2])*Y->Get(v[3])) 
	                    -(yP*X->Get(v[2])*Z->Get(v[3])) 
	                    -(xP*Z->Get(v[2])*Y->Get(v[3])) 
	                    -(zP*Y->Get(v[2])*X->Get(v[3])) )
	                -X->Get(v[0])*( +Y->Get(v[2])*Z->Get(v[3])
		                 +yP*Z->Get(v[2]) 
			             +zP*Y->Get(v[3])
		                 -yP*Z->Get(v[3])
					             -Z->Get(v[2])*Y->Get(v[3]) 
					     -zP*Y->Get(v[2]) )
	                +Y->Get(v[0])*( +X->Get(v[2])*Z->Get(v[3])
					     +xP*Z->Get(v[2])
					     +zP*X->Get(v[3])
		                 -xP*Z->Get(v[3])
					             -Z->Get(v[2])*X->Get(v[3]) 
					     -zP*X->Get(v[2]) )
		            -Z->Get(v[0])*( +X->Get(v[2])*Y->Get(v[3])
			             +xP*Y->Get(v[2]) 
					     +yP*X->Get(v[3])
		                 -xP*Y->Get(v[3])
				                 -Y->Get(v[2])*X->Get(v[3]) 
					     -yP*X->Get(v[2]) ) );

 V3 = (1.0/6.0) * (+1*( (X->Get(v[1])*yP*Z->Get(v[3])) 
	                    +(Y->Get(v[1])*zP*X->Get(v[3])) 
	                    +(Z->Get(v[1])*xP*Y->Get(v[3])) 
	                    -(Y->Get(v[1])*xP*Z->Get(v[3])) 
	                    -(X->Get(v[1])*zP*Y->Get(v[3])) 
	                    -(Z->Get(v[1])*yP*X->Get(v[3])) )
	        -X->Get(v[0])*( +yP*Z->Get(v[3])
	                     +Y->Get(v[1])*zP 
	                     +Z->Get(v[1])*Y->Get(v[3])
	                     -Y->Get(v[1])*Z->Get(v[3])
						 -zP*Y->Get(v[3]) 
						 -Z->Get(v[1])*yP )
			+Y->Get(v[0])*( +xP*Z->Get(v[3])
			             +X->Get(v[1])*zP
						 +Z->Get(v[1])*X->Get(v[3])
						 -X->Get(v[1])*Z->Get(v[3])
						 -zP*X->Get(v[3]) 
						 -Z->Get(v[1])*xP )
			-Z->Get(v[0])*( +xP*Y->Get(v[3])
			             +X->Get(v[1])*yP 
						 +Y->Get(v[1])*X->Get(v[3])
						 -X->Get(v[1])*Y->Get(v[3])
						 -yP*X->Get(v[3]) 
						 -Y->Get(v[1])*xP ) );


 *l1 = V1/V;
 *l2 = V2/V;
 *l3 = V3/V;
 *l4 = 1.0 - *l3 - *l2 - *l1;

 // retorna verdadeiro quando a particula cai dentro do elemento
 return ( (*l1>=0.0-EPSlocal) && (*l1<=1.0+EPSlocal) && 
          (*l2>=0.0-EPSlocal) && (*l2<=1.0+EPSlocal) && 
		  (*l3>=0.0-EPSlocal) && (*l3<=1.0+EPSlocal) && 
		  (*l4>=0.0-EPSlocal) && (*l4<=1.0+EPSlocal) ); 
}

void SemiLagrangean::setBC()
{
 int aux;
 for( int i=0;i<idbcu->Dim();i++ )
 {
  aux = (int) idbcu->Get(i);
  uParticle.Set( aux,uc->Get(aux) ); 
 }

 for( int i=0;i<idbcv->Dim();i++ )
 {
  aux = (int) idbcv->Get(i);
  vParticle.Set( aux,vc->Get(aux) ); 
 }

 for( int i=0;i<idbcw->Dim();i++ )
 {
  aux = (int) idbcw->Get(i);
  wParticle.Set( aux,wc->Get(aux) ); 
 }

 /* APENAS UTILIZADO NAS SIMULACOES DISK (NUCTE, NUC, NUZ) */
 /* ------------- caso com c.c. livre em w --------------- */
//--------------------------------------------------
//  for( int i=0;i<idbcv->Dim();i++ )
//  {
//   aux = (int) idbcv->Get(i);
//   if( Z->Get(aux) > 1.0 )
//    wParticle.Set( aux,wSol.Get(aux) ); 
//  }
//-------------------------------------------------- 
 /* ------------------------------------------------------ */

 for( int i=0;i<idbcc->Dim();i++ )
 {
  aux = (int) idbcc->Get(i);
  cParticle.Set( aux,cc->Get(aux) ); 
 }
}

clVector* SemiLagrangean::getUSL(){ return &uParticle; }
clVector* SemiLagrangean::getVSL(){ return &vParticle; }
clVector* SemiLagrangean::getWSL(){ return &wParticle; }
clVector* SemiLagrangean::getCSL(){ return &cParticle; }

