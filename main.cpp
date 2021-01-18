#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "intfit.h"


int make_data( long* size_, double** x_, double** y_, double xl, double xh )
{
   long size = 100;
   double *x,*y;

   x = (double *) malloc( ((size_t) size) * sizeof(double) );
   y = (double *) malloc( ((size_t) size) * sizeof(double) );

   for(long n=0;n<size;++n) {
      x[n] = ((double) n)/((double) size) * (xh - xl) + xl;
      y[n] = 100.0 + x[n]/100.0 + 10.0*sin( x[n]/1000.0 );
   }

   FILE *fp = fopen( "data.dat", "w" );
   for(long n=0;n<size;++n) {
      fprintf( fp, " %16.9e %16.9e \n", x[n], y[n] );
   }
   fclose( fp );

   *size_ = size;
   *x_ = x;
   *y_ = y;

   return 0;
}


int main( int argc, char *argv[] )
{
   double xlow = 1.0;
   double xhigh = 20000.0;

   long size;
   double *x, *y;
   (void) make_data( &size, &x, &y, xlow, xhigh );

//while(1) {

   inTFit_MultiFit mfit;

   mfit.setBounds( xlow, xhigh );

   inTFit_Fit fit1;
   fit1.setBounds( xlow, xhigh/2.0 );
// fit1.setTerms( "X X^2  X^3 X^4 LN(X) X^5 X^-1" );
   fit1.setTerms( "1 X X^2 X^3 X^4 X^5 X^6 " );
// fit1.setTerms( "1 X X^2 X^3 X^4 X^5 X^6 X^7 X^8 X^9 X^10 " );
   fit1.finalize();
   printf("Number of terms: %d \n", fit1.getNumTerms() );
   mfit.addFit( fit1 );
   fit1.clear();

/**/
   fit1.setBounds( xhigh/2.0, xhigh );
// fit1.setTerms( "X   X^2 LOG(X)  " );
   fit1.setTerms( "X X^2 " );
   fit1.finalize();
   printf("Number of terms: %d \n", fit1.getNumTerms() );
   mfit.addFit( fit1 );
   fit1.clear();
/**/

   mfit.finalize();

   mfit.storeData( size, x, y );

   mfit.compute();

  printf("========================= LOOP ==========================\n");
//}

   return(0);
}


