// =================================================================== // 
// this is file TElement.cpp, created at 21-Jun-2007                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "TElement.h"

TElement::~TElement(){}

TElement::TElement( Model3D &_m )
{
 m = &_m;
 numGLEU = _m.getNumGLEU(); 
 numGLEP = _m.getNumGLEP(); 
 numGLEC = _m.getNumGLEC(); 
}




