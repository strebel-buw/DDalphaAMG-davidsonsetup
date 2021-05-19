#include "includes.h"

extern int nx, ny, nz, nt;
extern int sites_on_node;
extern su3_matrix *tlink;
double plaq=0;

// if there are endian problems ...
void byteswap( void *in_p )
{
  //      char *in=(char *) in_p;
  //      char tmp[4];
  //      tmp[0] = in[3];
  //      tmp[1] = in[2];
  //      tmp[2] = in[1];
  //      tmp[3] = in[0];
  //  
  //      in[0] = tmp[0];
  //      in[1] = tmp[1];
  //      in[2] = tmp[2];
  //      in[3] = tmp[3];
}

/* Copy a su3 matrix:  b <- a   */
void su3mat_copy( su3_matrix * a, su3_matrix * b )
{
  register int i, j;
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ )
    {
      b->ROWCOL( i, j ).real = a->ROWCOL( i, j ).real;
      b->ROWCOL( i, j ).imag = a->ROWCOL( i, j ).imag;
    }
}

void dump_mat( su3_matrix * m )
{
  int i, j;
  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 3; j++ )
      printf( "(%.5e,%.2e)\t", m->ROWCOL( i, j ).real, m->ROWCOL( i, j ).imag );
    printf( "\n" );
  }
  printf( "\n" );
}

void special_su3mat( double number, su3_matrix *dest )
{
  register int i,j;
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      dest->e[i][j].real = dest->e[i][j].imag = 0.0;
      if (i==j) 
        dest->e[i][j].real=number;	
    }
}

/* b <- a^T */
void su3mat_transpose( su3_matrix * a, su3_matrix * b )
{
  register int i,j;
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      b->e[i][j].real = a->e[j][i].real;
      b->e[i][j].imag = a->e[j][i].imag;
    }
}



int memsize=0;
#define MEMALIGN( variable, type, size ) do{ (variable) = ( type * ) memalign( 16, (size) * sizeof( type ) ); if( (variable) == NULL ) { printf( "ERROR %s:%d in %s(): memory allocation failed\n",__FILE__,__LINE__,__FUNCTION__);  exit(1);}; memsize+=(size)*sizeof(type); }while(0)
#define FREE( variable, type, size ) do{ free(variable); variable=NULL; memsize-=(size)*sizeof(type); }while(0)

int node_index( int x, int y, int z, int t )
{
  register int i, xr, yr, zr;
  int coords[4];
  coords[XUP] = x;
  coords[YUP] = y;
  coords[ZUP] = z;
  coords[TUP] = t;
  xr = coords[XUP] % nx;
  yr = coords[YUP] % ny;
  zr = coords[ZUP] % nz;
  i = xr + nx * (yr + ny * (zr + nz * coords[TUP]));
  
  if( ( x + y + z + t ) % 2 == 0 )
  { /* even site */
    return ( i / 2 );
  }
  else
  {
    return ( ( i + sites_on_node ) / 2 );
  }
}



void dump_conf( void ) {
  
  int x, y, z, t, ind, dir;
  
  printf("\n\nconfiguration dump:\n\n");
  
  for ( t = 0; t < nt; t++ )
    for ( z = 0; z < nz; z++ )
      for ( y = 0; y < ny; y++ )
        for ( x = 0; x < nx; x++ )
        {
          ind = node_index( x, y, z, t );
          for ( dir = TUP; dir >= XUP; dir-- )
          {
            printf("(t=%i,z=%i,y=%i,x=%i) ", t, z, y, x );
            if ( dir == XUP )
              printf("dir=+x\n");
            else if ( dir == YUP )
              printf("dir=+y\n");
            else if ( dir == ZUP )
              printf("dir=+z\n");
            else if ( dir == TUP )
              printf("dir=+t\n");
            
            dump_mat( &(tlink[4 * ind + dir ]) );
          }
        }
}



/* save lattice */
void save_lattice( int flag, char *filename )
{
  FILE *f=NULL;
  int x, y, z, t, ind;
  su3_matrix mat_transpose;
  su3_matrix_comp mat_comp;
  su3_matrix_comp *alllinks_comp, *tlink_comp;
  int size, nd, dir, i;
  int myx, myy, myz;
  
  
  switch (flag)
  {
    case SAVE_SERIAL:
      printf( "Writing lattice %s ...\n", filename );
      if( ( f = fopen( filename, "w" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        exit( 1 );
      }
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      fwrite( &nx, sizeof( int ), 1, f );
      fwrite( &ny, sizeof( int ), 1, f );
      fwrite( &nz, sizeof( int ), 1, f );
      fwrite( &nt, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      
      for ( i = 0; i < 4 * sites_on_node * 2; i++ )
        byteswap( ( char * ) ( ( float * ) tlink + i ) );
      
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              ind = node_index( x, y, z, t );
              for ( dir = XUP; dir <= TUP; dir++ )
              {
                fwrite( &(tlink[4 * ind + dir ]), sizeof( su3_matrix ), 1, f );
              }
            }
      }
      fclose( f );
      
      break;
    case SAVE_SERIAL_COMP:
      printf( "Writing compressed lattice %s...\n", filename );
      fflush(0);
      if( ( f = fopen( filename, "w" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        exit( 1 );
      }
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      fwrite( &nx, sizeof( int ), 1, f );
      fwrite( &ny, sizeof( int ), 1, f );
      fwrite( &nz, sizeof( int ), 1, f );
      fwrite( &nt, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      fflush(0);
      
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              ind = node_index( x, y, z, t );
              for ( dir = XUP; dir <= TUP; dir++ )
              {
                //dump_mat ( &(tlink[4 * ind + dir]));
                su3_to_comp( &(tlink[4 * ind + dir]), &mat_comp );
                for ( i = 0; i < 8; i++ ){
                  //printf ("%e->", mat_comp.a[i]); 
                  byteswap( (char * ) &mat_comp.a[i] );	
                  //printf ("%e.\n", mat_comp.a[i]); 
                }	
                fwrite( &mat_comp, sizeof( su3_matrix_comp ), 1, f );
              }
            }
      }
      fclose( f );
      
      break;
    case SAVE_MG4QCD:
      // 4 integers: lattice size (T,Z,Y,X)
      // 1 double: plaquette
      
      // t slowest running index
      // x fastest running index
      // all positive directions
      // ordering: +T,+Z,+Y,+X
      // SU3 matrices stored in row major format
      
      printf( "Writing lattice %s in MG4QCD format...\n", filename );
      fflush(0);
      if( ( f = fopen( filename, "wb" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        exit( 1 );
      }
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      fwrite( &nt, sizeof( int ), 1, f );
      fwrite( &nz, sizeof( int ), 1, f );
      fwrite( &ny, sizeof( int ), 1, f );
      fwrite( &nx, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      fflush(0);
      plaq = 0;
      fwrite( &plaq, sizeof( double ), 1, f );
      
      for ( t = 0; t < nt; t++ ) {
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ ) {
              
              ind = node_index( x, y, z, t );
              
              for ( dir = TUP; dir >= XUP; dir-- ) {
                su3mat_transpose( &(tlink[4 * ind + dir ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
              }
              
            }
      }
      fclose( f );
      
      break;
      
    case SAVE_DDHMC:
      // 4 integers: lattice size (T,X,Y,Z)
      // 1 double: plaquette
      
      // t slowest running index
      // x fastest running index
      // all odd points: all positive directions and negative directions
      // ordering: +T,-T,+X,-X,+Y,-Y,+Z,-Z
      //     ^
      //     |
      // --> o -->
      //     ^
      //     |
      // SU3 matrices stored in row major format
      
      printf( "Writing lattice %s in DDHMC format...\n", filename );
      fflush(0);
      if( ( f = fopen( filename, "wb" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        exit( 1 );
      }
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      fwrite( &nt, sizeof( int ), 1, f );
      fwrite( &nz, sizeof( int ), 1, f );
      fwrite( &ny, sizeof( int ), 1, f );
      fwrite( &nx, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      fflush(0);
      
      fwrite( &plaq, sizeof( double ), 1, f );
      
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              if ( ( x+y+z+t ) % 2 != 0 ) // for all odd points
              {
                
                // plus T direction
                ind = node_index( x, y, z, t );
                su3mat_transpose( &(tlink[4 * ind + TUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                // minus T direction (not daggered)
                ind = node_index( x, y, z, (t-1+nt)%nt );
                su3mat_transpose( &(tlink[4 * ind + TUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                
                // plus Z direction
                ind = node_index( x, y, z, t );
                su3mat_transpose( &(tlink[4 * ind + ZUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                // minus Z direction (not daggered)
                ind = node_index( x, y, (z-1+nz)%nz, t );
                su3mat_transpose( &(tlink[4 * ind + ZUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                
                // plus Y direction
                ind = node_index( x, y, z, t );
                su3mat_transpose( &(tlink[4 * ind + YUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                // minus Y direction (not daggered)
                ind = node_index( x, (y-1+ny)%ny, z, t );
                su3mat_transpose( &(tlink[4 * ind + YUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );				
                
                // plus X direction
                ind = node_index( x, y, z, t );
                su3mat_transpose( &(tlink[4 * ind + XUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                // minus X direction (not daggered)
                ind = node_index( (x-1+nx)%nx, y, z, t );
                su3mat_transpose( &(tlink[4 * ind + XUP ]), &mat_transpose );
                fwrite( &mat_transpose, sizeof( su3_matrix ), 1, f );
                
              }
            }
      }
      fclose( f );
      
      break;
      
      
  }
  
  return;
}

/* reads lattice */
void reload_lattice( int flag, char *filename )
{
  FILE *f;
  int i, j, x, y, z, t, ind, dir;
  double max_deviation;
  su3_matrix mat_transpose;
  su3_matrix_comp *alllinks_comp;
  su3_matrix_comp mat_comp;
  su3_matrix mat;
  
  switch ( flag )
  {
    case RELOAD_SERIAL: /* read binary lattice serially */
      printf( "Restoring lattice %s ...\n", filename );
      fflush( 0 );
      if( ( f = fopen( filename, "r" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        fflush( 0 );
        exit ( 1 );
      }
      fread( &nx, sizeof( int ), 1, f );
      fread( &ny, sizeof( int ), 1, f );
      fread( &nz, sizeof( int ), 1, f );
      fread( &nt, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      sites_on_node=nx*ny*nz*nt;
      MEMALIGN(tlink, su3_matrix, sites_on_node*4);
      
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t);
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              ind = node_index( x, y, z, t );
              for ( dir = XUP; dir <= TUP; dir++ )
              {
                fread( &mat, sizeof( su3_matrix ), 1, f );
                su3mat_copy( &mat, &( tlink[4 * ind + dir] ) );
                
              }
            }
      }
      for ( i = 0; i < 4 * sites_on_node * 2; i++ )
        byteswap( ( char * ) ( ( float * ) tlink + i ) );
      fclose( f );
      printf( "done.\n" );
      break;
      
    case RELOAD_SERIAL_COMP:
      printf( "Restoring lattice %s...\n", filename );
      fflush(0);
      if( ( f = fopen( filename, "r" ) ) == 0 )
      {
        printf( "ERROR reload_lattice: Unable to open file %s\n", filename );
        fflush(0);
        exit( 1 );
      }
      
      fread( &nx, sizeof( int ), 1, f );
      fread( &ny, sizeof( int ), 1, f );
      fread( &nz, sizeof( int ), 1, f );
      fread( &nt, sizeof( int ), 1, f );
      
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      sites_on_node=nx*ny*nz*nt;
      MEMALIGN(tlink, su3_matrix, sites_on_node*4);
      
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      fflush(0);
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              //printf ("."); 
              ind = node_index( x, y, z, t );
              for ( dir = XUP; dir <= TUP; dir++ )
              {
                fread( &mat_comp, sizeof( su3_matrix_comp ), 1, f );
                
                for ( i = 0; i < 8; i++ ){
                  //printf ("%e->", mat_comp.a[i]); 
                  byteswap( (char * ) &mat_comp.a[i] );	
                  //printf ("%e.\n", mat_comp.a[i]); 
                }	
                //printf("%i %i\n", ind, dir); fflush(0);
                comp_to_su3( &mat_comp, &(tlink[4 * ind + dir]) );
              }
            }
      }
      fclose( f );
      printf( "done.\n" );
      fflush(0);
      break;
      
      
      case RELOAD_MG4QCD:
        // 4 integers: lattice size (T,Z,Y,X)
        // 1 double: plaquette
        
        // t slowest running index
        // x fastest running index
        // all positive directions
        // ordering: +T,+Z,+Y,+X
        // SU3 matrices stored in row major format
        
        printf( "reading lattice %s in MG4QCD format...\n", filename );
        fflush(0);
        if( ( f = fopen( filename, "rb" ) ) == 0 )
        {
          printf( "Unable to open file %s\n", filename );
          exit( 1 );
        }
        byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
        fread( &nt, sizeof( int ), 1, f );
        fread( &nz, sizeof( int ), 1, f );
        fread( &ny, sizeof( int ), 1, f );
        fread( &nx, sizeof( int ), 1, f );
        byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
        printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
        fflush(0);
        plaq = 0;
        fread( &plaq, sizeof( double ), 1, f );
        
        sites_on_node=nx*ny*nz*nt;
        MEMALIGN(tlink, su3_matrix, sites_on_node*4);
        
        for ( t = 0; t < nt; t++ ) {
          printf("t=%i\n", t); 
          for ( z = 0; z < nz; z++ )
            for ( y = 0; y < ny; y++ )
              for ( x = 0; x < nx; x++ ) {
                ind = node_index( x, y, z, t );
                
                for ( dir = TUP; dir >= XUP; dir-- ) {
                  fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                  su3mat_transpose( &mat_transpose, &(tlink[4 * ind + dir ]) );
                }
                
              }
        }
        fclose( f );
        
        break;
        
        
case RELOAD_DDHMC:
      // 4 integers: lattice size (T,X,Y,Z)
      // 1 double: plaquette
      
      // t slowest running index
      // x fastest running index
      // all odd points: all positive directions and negative directions
      // ordering: +T,-T,+X,-X,+Y,-Y,+Z,-Z
      //     ^
      //     |
      // --> o -->
      //     ^
      //     |
      // SU3 matrices stored in row major format
      
      printf( "Writing lattice %s in DDHMC format...\n", filename );
      fflush(0);
      if( ( f = fopen( filename, "r" ) ) == 0 )
      {
        printf( "Unable to open file %s\n", filename );
        exit( 1 );
      }
      fread( &nx, sizeof( int ), 1, f );
      fread( &ny, sizeof( int ), 1, f );
      fread( &nz, sizeof( int ), 1, f );
      fread( &nt, sizeof( int ), 1, f );
      byteswap(&nx); byteswap(&ny); byteswap(&nz); byteswap(&nt);
      printf( "nx=%i, ny=%i, nz=%i, nt=%i.\n", nx, ny, nz, nt );
      fflush(0);
      
      sites_on_node=nx*ny*nz*nt;
      MEMALIGN(tlink, su3_matrix, sites_on_node*4);
      
      fread( &plaq, sizeof( double ), 1, f );
      
      for ( t = 0; t < nt; t++ ){
        printf("t=%i\n", t); 
        for ( z = 0; z < nz; z++ )
          for ( y = 0; y < ny; y++ )
            for ( x = 0; x < nx; x++ )
            {
              if ( ( x+y+z+t ) % 2 != 0 ) // for all odd points
              {
                
                // plus T direction
                ind = node_index( x, y, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );                
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + TUP ]) );
                // minus T direction (not daggered)
                ind = node_index( x, y, z, (t-1+nt)%nt );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + TUP ]) );
                
                // plus Z direction
                ind = node_index( x, y, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + ZUP ]) );
                // minus Z direction (not daggered)
                ind = node_index( x, y, (z-1+nz)%nz, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + ZUP ]) );
                
                // plus Y direction
                ind = node_index( x, y, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + YUP ]) );
                // minus Y direction (not daggered)
                ind = node_index( x, (y-1+ny)%ny, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );   
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + YUP ]) );  
                
                // plus X direction
                ind = node_index( x, y, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + XUP ]) );
                // minus X direction (not daggered)
                ind = node_index( (x-1+nx)%nx, y, z, t );
                fread( &mat_transpose, sizeof( su3_matrix ), 1, f );
                su3mat_transpose( &mat_transpose, &(tlink[4 * ind + XUP ]) );
              }
            }
      }
      fclose( f );
      
      break;
      
      
  }
  
  return;
}
