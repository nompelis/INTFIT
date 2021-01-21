#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "intfit.h"


int read_data( long* size_, double** x_, double** y_ )
{
   long size;
   double *x=NULL,*y=NULL;
   size_t mem_size=0;

   FILE *fp = fopen( "data.dat", "r" );
   if( fp == NULL ) {
      fprintf( stdout, "Error opening data file \n");
      return 1;
   }

   size = 0;
   while( 1 ) {
      char string[100];

      if( size == ((long) mem_size) ) {
         double *p;
         p = (double *) realloc( x, (mem_size+100)*sizeof(double) );
         if( p == NULL ) {
            fprintf( stdout, "Error (re)allocating memory \n");
            free( x );
            free( y );
            return -1;
         }
         x = p;

         p = (double *) realloc( y, (mem_size+100)*sizeof(double) );
         if( p == NULL ) {
            fprintf( stdout, "Error (re)allocating memory \n");
            free( x );
            free( y );
            return -2;
         }
         y = p;

         mem_size += 100;
         fprintf( stdout, "Reallocation made (size: %ld) \n",size);
      }

      fgets( string, 100, fp );
      if( feof( fp ) ) break;
   // fprintf( stdout, "STRING: \"%s\"\n", string );
      sscanf( string, "%lf %lf", &(x[size]), &(y[size]) );
      ++size;
   }
   fclose( fp );




   *size_ = size;
   *x_ = x;
   *y_ = y;

   return 0;
}


int main( int argc, char *argv[] )
{
   long size;
   double *x, *y;
   (void) read_data( &size, &x, &y );

   inTFit_MultiFit mfit;

   double xlow, xhigh;
   xlow = x[0], xhigh = x[size-1];    // pick these numbers from the data
   mfit.setBounds( xlow, xhigh );

   inTFit_Fit fit1;
   fit1.setBounds( xlow, xhigh );
   fit1.setTerms( "1 X X^2 X^3 X^4 X^5 " );
   fit1.finalize();
   printf("Number of terms: %d \n", fit1.getNumTerms() );
   mfit.addFit( fit1 );
   fit1.clear();

   // form the derivative on the left numerically end from the data
   double dval = (y[1] - y[0])/(x[1] - x[0]);
   mfit.addConstraint( CONSTRAINT_DERIVATIVE, 0, dval, xlow );
   // form the derivative on the right numerically end from the data
          dval = (y[size-1] - y[size-2])/(x[size-1] - x[size-2]);
   mfit.addConstraint( CONSTRAINT_DERIVATIVE, 0, dval, xhigh );
   // add constrain on the left side of the fit to match the data
   dval = y[0];
   mfit.addConstraint( CONSTRAINT_VALUE, 0, dval, xlow );

   mfit.finalize();

   mfit.storeData( size, x, y );

   mfit.compute();
   mfit.dumpFitData( (char*) "data", 1 );
   mfit.clear();

   return(0);
}

