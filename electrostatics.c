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
#include <stdio.h>

#define sub(a,b) (_mm_sub_ps(a,b))
#define mul(a,b) (_mm_mul_ps(a,b))
#define add(a,b) (_mm_add_ps(a,b))
#define sq(a)	 (_mm_sqrt_ps(a))
#define set(a)	 (_mm_set1_ps(a))

#define pythagorasVectCore2Duo(x1, y1, z1, x2, y2, z2, d) {\
	__m128 *vx1, *vy1, *vz1;			\
	__m128 r1, r2, r3;				\
	vx1 = (__m128*)x1;				\
	vy1 = (__m128*)y1;				\
	vz1 = (__m128*)z1;				\
	r1 = sub(*vx1,set(x2));				\
	r1 = mul(r1,r1);				\
	r2 = sub(*vy1,set(y2));				\
	r2 = mul(r2,r2);				\
	r3 = sub(*vz1,set(z2));				\
	r3 = mul(r3,r3);				\
	*((__m128*)d) = sq(add(add(r1,r2),r3));		\
}

void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int	residue , atom ;

/************/

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      This_Structure.Residue[residue].Atom[atom].charge = 0.0 ;

      /* peptide backbone */

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " N  " ) == 0 ) {
        if( strcmp( This_Structure.Residue[residue].res_name , "PRO" ) == 0 ) {
          This_Structure.Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          This_Structure.Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 ) This_Structure.Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " O  " ) == 0 ) {
        This_Structure.Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == This_Structure.length  ) This_Structure.Residue[residue].Atom[atom].charge = -1.00 ;
      }

      /* charged residues */

      if( ( strcmp( This_Structure.Residue[residue].res_name , "ARG" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "ASP" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "GLU" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "LYS" ) == 0 ) && ( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NZ " ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  1.00 ;

    }
  }

/************/

}


inline void allocatedata_electric_field(float **charge, float **coord1, float **coord2, float **coord3, float **epsilon, float **distance, int maxTotalElements) {
  if ((posix_memalign((void**)charge, 16, maxTotalElements*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
	if ((posix_memalign((void**)coord1, 16, maxTotalElements*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
	if ((posix_memalign((void**)coord2, 16, maxTotalElements*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
	if ((posix_memalign((void**)coord3, 16, maxTotalElements*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
  if ((posix_memalign((void**)distance, 16, 4*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
  if ((posix_memalign((void**)epsilon, 16, 4*sizeof(float))!=0)) {
    printf("No memory.\n");
    exit(-1);
  }
}


/************************/
#define calcEpsilon(k) {																		\
	if (distance[k] < 2.0 ) distance[k] = 2.0 ;               \
	if (distance[k] >= 8.0 ) {                                \
		epsilon[k] = 80 ;                                       \
	} else {                                                  \
		if (distance[k] <= 6.0 ) {                              \
			epsilon[k] = 4;                                       \
		} else {                                                \
			epsilon[k] = ( 38 * distance[k] ) - 224;              \
		}                                                       \
	}                                                         \
}

//Asumeix que el block es multiple de 4
#define computeBlock(start, stop) {			  											\
	for (i = start; i < stop; i+=4) {                             \
		pythagorasVectCore2Duo(&coord1[i], &coord2[i], &coord3[i],  \
											x_centre, y_centre, z_centre, distance);  \
		                                                            \
		for (k = 0; k < 4; k++) calcEpsilon(k);					            \
                                                                \
		_phiSet=_mm_div_ps(*((__m128 *)&charge[i]),                 \
			_mm_mul_ps(*((__m128*)epsilon),*((__m128*)distance)));    \
		phi += phiSet[0] + phiSet[1] + phiSet[2] + phiSet[3];       \
	}                                                             \
}

//No assumeix que es block es multiple de 4
#define computeEndBlock(start, stop) {													\
	computeBlock(start, (stop-3));                                \
                                                                \
	for (;i < stop; i++) {                                        \
	  distance[0] = pythagoras2( coord1[i], coord2[i],            \
					coord3[i], x_centre, y_centre, z_centre);             \
		calcEpsilon(0);                                             \
		phi += ( charge[i] / ( epsilon[0] * distance[0] ) ) ;       \
	}                                                             \
}

void electric_field( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {
  /* Counters */
  int	residue , atom ;
  /* Co-ordinates */
  int	x, y, z, i, k;
  float		x_centre , y_centre , z_centre ;
  /* Variables */
  float		phi;
	float *distance, *epsilon;
  /* Vectorization stuff */
  __m128 _phiSet;
	float *phiSet = (float *) &_phiSet;
	/* Blocking stuff */
	int block_size = 512;
	/* Threads stuff */
	int num_threads = 2; 
 
  int maxTotalElements = 0;
  int totalElements = 0;
	float *charge;
	float *coord1, *coord2, *coord3;

	for (residue = 1; residue <= This_Structure.length; residue++)
		maxTotalElements += This_Structure.Residue[residue].size;

	allocatedata_electric_field(&charge, &coord1, &coord2, &coord3, &epsilon, &distance, maxTotalElements);

	for (residue = 1; residue <= This_Structure.length; residue++)
		for (atom = 1; atom <= This_Structure.Residue[residue].size; atom++)
			if (This_Structure.Residue[residue].Atom[atom].charge != 0){
				charge[totalElements] = This_Structure.Residue[residue].Atom[atom].charge;
				coord1[totalElements] = This_Structure.Residue[residue].Atom[atom].coord[1]; 
				coord2[totalElements] = This_Structure.Residue[residue].Atom[atom].coord[2]; 
				coord3[totalElements] = This_Structure.Residue[residue].Atom[atom].coord[3]; 
				totalElements++;
			}


  for (x = 0; x < grid_size; x++)
    for (y = 0; y < grid_size; y++)
      for (z = 0; z < grid_size; z++)
        grid[gaddress(x,y,z,grid_size)] = (fftw_real)0;

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;
  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;


	int block_ini, block_fin;

	for (block_ini = 0; block_ini < totalElements-(block_size-1); block_ini += block_size) { 
		block_fin = block_ini + block_size;
	  for( x = 0 ; x < grid_size ; x ++ ) {
	    printf( "." ) ;
	    x_centre  = gcentre( x , grid_span , grid_size ) ;
	    for( y = 0 ; y < grid_size ; y ++ ) {
	      y_centre  = gcentre( y , grid_span , grid_size ) ;
	      for( z = 0 ; z < grid_size ; z ++ ) {
	        z_centre  = gcentre( z , grid_span , grid_size ) ;
	        phi = 0 ;
	
					computeBlock(block_ini, block_fin);
	    	  grid[gaddress(x,y,z,grid_size)] += (fftw_real)phi ;
	      }
	    }
	  }
	}

	if (block_ini != totalElements) { //ultima passada
	  for( x = 0 ; x < grid_size ; x ++ ) {
	    printf( "." ) ;
	    x_centre  = gcentre( x , grid_span , grid_size ) ;
	    for( y = 0 ; y < grid_size ; y ++ ) {
	      y_centre  = gcentre( y , grid_span , grid_size ) ;
	      for( z = 0 ; z < grid_size ; z ++ ) {
	        z_centre  = gcentre( z , grid_span , grid_size ) ;
	        phi = 0 ;
	
					computeEndBlock(block_ini, totalElements);
	    	  grid[gaddress(x,y,z,grid_size)] += (fftw_real)phi ;
	      }
	    }
	  }
	}
  printf( "\n" ) ;


  return ;

}



/************************/



void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	x_low , x_high , y_low , y_high , z_low , z_high ;

  float		a , b , c ;
  float		x_corner , y_corner , z_corner ;
  float		w ;

  /* Variables */

  float		one_span ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/

  one_span = grid_span / (float)grid_size ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( This_Structure.Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( This_Structure.Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( This_Structure.Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        x_high = x_low + 1 ;
        y_high = y_low + 1 ;
        z_high = z_low + 1 ;

        a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        for( x = x_low ; x <= x_high  ; x ++ ) {
 
          x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;

          for( y = y_low ; y <= y_high  ; y ++ ) {

            y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;

            for( z = z_low ; z <= z_high  ; z ++ ) {

              z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;

              w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;

              grid[gaddress(x,y,z,grid_size)] += (fftw_real)( w * This_Structure.Residue[residue].Atom[atom].charge ) ;

            }
          }
        }

      }

    }
  }

/************/

  return ;

}



/************************/



void electric_field_zero_core( int grid_size , fftw_real *elec_grid , fftw_real *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x , y , z ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        if( surface_grid[gaddress(x,y,z,grid_size)] == (fftw_real)internal_value ) elec_grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/

  return ;

}
