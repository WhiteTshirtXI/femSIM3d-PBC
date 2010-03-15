/**=========================================================================
 *\file Id:TElement3D.h, GESAR
 *\author Gustavo Rabello dos Anjos
 *\date   20-jun-2007
 *\comment Classe TElement 3D para o modulo HYDRO
 *\references
 *\version 1.0
 *\updates
 	version    date         author             comment
	1.0        08/07/2007   Gustavo            finalizacao da classe
 =========================================================================**/

#ifndef TElement_H
#define TElement_H

#include "clMatrix.h"
#include "Model3D.h"

class TElement
{
 public:
  
  TElement( Model3D &_m );
  virtual ~TElement();
  int numGLEU;                  ///< numero de graus de liberdade de vel.
  int numGLEP;                  ///< numero de graus de liberdade de pressao
  int numGLEC;                  ///< numero de graus de liberdade de pressao
  Model3D *m;

  double massele[5][5];         ///< matriz de massa do elemento           
  double kxx[5][5];             ///< matriz do laplaciano do bloco 11
  double kxy[5][5];             ///< matriz do laplaciano do bloco 12
  double kxz[5][5];             ///< matriz do laplaciano do bloco 13
  double kyx[5][5];             ///< matriz do laplaciano do bloco 21
  double kyy[5][5];             ///< matriz do laplaciano do bloco 22
  double kyz[5][5];             ///< matriz do laplaciano do bloco 23
  double kzx[5][5];             ///< matriz do laplaciano do bloco 31
  double kzy[5][5];             ///< matriz do laplaciano do bloco 32
  double kzz[5][5];             ///< matriz do laplaciano do bloco 33
  double gxele[5][4];           ///< matriz do oper. grad do bloco 1
  double gyele[5][4];           ///< matriz do oper. grad do bloco 2
  double gzele[5][4];           ///< matriz do oper. grad do bloco 3
  double dxele[4][5];           ///< matriz do oper. div do bloco 1
  double dyele[4][5];           ///< matriz do oper. div do bloco 2
  double dzele[4][5];           ///< matriz do oper. div do bloco 3
  double masselec[4][4];         ///< matriz de massa do elemento           
  double masselec2[4][4];         ///< matriz de massa do elemento           
  double kelec[4][4];             ///< matriz do laplaciano do bloco 11
  double kxxc[4][4];             ///< matriz do laplaciano do bloco 11
  double kyyc[4][4];             ///< matriz do laplaciano do bloco 11
  double kzzc[4][4];             ///< matriz do laplaciano do bloco 11
  double gxelec[4][4];           ///< matriz do oper. grad do bloco 1
  double gyelec[4][4];           ///< matriz do oper. grad do bloco 2
  double gzelec[4][4];           ///< matriz do oper. grad do bloco 3
  double dxelec[4][4];           ///< matriz do oper. div do bloco 1
  double dyelec[4][4];           ///< matriz do oper. div do bloco 2
  double dzelec[4][4];           ///< matriz do oper. div do bloco 3

};

#endif

