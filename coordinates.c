/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"
#include <emmintrin.h>


#define sub(a,b) (_mm_sub_ps(a,b))
#define mul(a,b) (_mm_mul_ps(a,b))
#define add(a,b) (_mm_add_ps(a,b))
#define sq(a)	 (_mm_sqrt_ps(a))
#define set(a)	 (_mm_set1_ps(a))
/************/

int gord( float position , float grid_span , int grid_size ) {

  int ordinate ;

  float one_span = grid_span / (float)grid_size ;

  ordinate = (int)( position / one_span ) + ( grid_size / 2 ) ;

  if( position < 0 ) ordinate -= 1 ;

  return ordinate ;

}

/************/

float pythagoras( float x1 , float y1 , float z1 , float x2 , float y2 , float z2 ) {
  return sqrt( ( ( x1 - x2 ) * ( x1 - x2 ) ) + ( ( y1 - y2 ) * ( y1 - y2 ) ) + ( ( z1 - z2 ) * ( z1 - z2 ) ) ) ;

}

float pythagoras2( float x1 , float y1 , float z1 , float x2 , float y2 , float z2 ) {

  return sqrt( ( ( x1 - x2 ) * ( x1 - x2 ) ) + ( ( y1 - y2 ) * ( y1 - y2 ) ) + ( ( z1 - z2 ) * ( z1 - z2 ) ) ) ;

}




void pythagorasVectCore2Duo2( float *x1 , float *y1 , float *z1 , float x2 , float y2 , float z2, float *d) {
	__m128 *vx1, *vy1, *vz1;
	__m128 r1, r2, r3;
	vx1 = (__m128*)x1;
	vy1 = (__m128*)y1;
	vz1 = (__m128*)z1;
	r1 = sub(*vx1,set(x2));
	r1 = mul(r1,r1);
	r2 = sub(*vy1,set(y2));
	r2 = mul(r2,r2);
	r3 = sub(*vz1,set(z2));
	r3 = mul(r3,r3);
	*((__m128*)d) = sq(add(add(r1,r2),r3));
}

