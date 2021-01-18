/******************************************************************************
 Copyright (c) 2021, Ioannis Nompelis
 All rights reserved.

 Redistribution and use in source and binary forms, with or without any
 modification, are permitted provided that the following conditions are met:
 1. Redistribution of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistribution in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. All advertising materials mentioning features or use of this software
    must display the following acknowledgement:
    "This product includes software developed by Ioannis Nompelis."
 4. Neither the name of Ioannis Nompelis and his partners/affiliates nor the
    names of other contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.
 5. Redistribution or use of source code and binary forms for profit must
    have written permission of the copyright holder.
 
 THIS SOFTWARE IS PROVIDED BY IOANNIS NOMPELIS ''AS IS'' AND ANY
 EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL IOANNIS NOMPELIS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


#include "intfit.h"

#ifdef __cplusplus
extern "C" {
#endif


//
// A function to perform Gaussian elimination with pivoting
// (This is taken from IN's "inUtils" collrection. This is _not_ Saad's
// form, but the textbook version.)
//

int inUtils_GaussianElimination( int ne, double* a, double* c, double* d )
{
   double rm,ra,ba,rl;
   int i,j,k,il;

// loop over columns and skip previous rows

   for(j=0;j<ne;++j) {                       // loop over columns

      rm = 0.0;                              // init to minimum
      il = j;                                // initialize to running top
      for(i=j;i<ne;++i) {                    // loop over remaining rows
         ra = fabs(a[i*ne + j]);             // get magnitude of term
         if(ra > rm) {                       // new value is larger
            il = i;                          // store row number
            rm = ra;                         // update maximum magnitude
         }
      }

      for(k=0;k<ne;++k) {
         d[k] = a[il*ne + k];                // store elements of large row
      }
      ba      = c[il];                       // store RHS of large row
      for(k=0;k<ne;++k) {
         a[il*ne + k] = a[j*ne + k];         // replace elem. of large w/ top
      }
      c[il]      = c[j];                     // replace RHS of large row w/ t
      for(k=0;k<ne;++k) {
         a[j*ne + k] = d[k];                 // put elem. of large row to top
      }
      c[j]      = ba;                        // put RHS of large row to top

      ra = 1.0/d[j];                         // get inverse of large row's nz
      for(i=j+1;i<ne;++i) {                  // work on remaining rows
         rl = a[i*ne + j]*ra;                // form equation multiplier
         a[i*ne + j] = 0.0;                  // zero identically
         for(k=j+1;k<ne;++k) {               // loop over lower rows elements
            a[i*ne + k] = a[i*ne + k] - rl*d[k];  // subtract large row
         }
         c[i] = c[i] - rl*ba;                // subtract large RHS from RHS
      }
   }

// back-substitution
   c[ne-1] = c[ne-1]/a[(ne-1)*ne + ne-1];    // invert for last element
   for(i=ne-2;i>=0;--i) {                    // eliminate backwards
      ra = 0.0;                              // init the accumulated sum
      for(k=ne-1;k>=i+1;--k) {               // backward sum of solved elems
         ra = ra + a[i*ne + k]*c[k];         // accumulate sum of solved
      }
      c[i] = (c[i] - ra)/a[i*ne + i];        // solve for element
   }

   return 0;
}


//
// Term class
//

inTFit_Term::inTFit_Term( inTerm_type type_, int order_ )
{
   type = type_;
   order = order_;
}


inTFit_Term::~inTFit_Term()
{

}


inTerm_type inTFit_Term::getType( void ) const
{
   return type;
}


int inTFit_Term::getOrder( void ) const
{
   return order;
}


double inTFit_Term::eval( double t ) const
{
   double result = 0.0;
   switch( type ) {
    case( TERM_NULL ):
       result = 0.0;
    break;
    case( TERM_CONST ):
       result = 1.0;
    break;
    case( TERM_MONOMIAL ):
       result = pow( t, (double) order );
    break;
    case( TERM_EXP ):
       result = exp( t );
    break;
    case( TERM_LOG ):
       result = log( t );
    break;
   }

   return result;
}


double inTFit_Term::evalDer( double t ) const
{
   double result = 0;

   return result;
}


double inTFit_Term::evalInt( double t ) const
{
   double result = 0;

   return result;
}


//
// Fit class
//

inTFit_Fit::inTFit_Fit()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Constructor \n",CLASS);
#endif

   num_terms = 0;
   xs = -9.9e99;
   xe = -9.9e99;
   ibound = 0;
}


inTFit_Fit::~inTFit_Fit()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Deconstructing \n",CLASS);
#endif

}


int inTFit_Fit::setBounds( double xstart, double xend )
{
   if( ibound != 0 ) {
      fprintf( stdout, " Error: bounds have already been set \n");
      fprintf( stdout, "        xs = %16.9e, xe = %16.9e \n", xs, xe );
      return 1;
   }

   xs = xstart;
   xe = xend;
   ibound = 1;
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Bounds: %16.9e to %16.9e \n",CLASS, xs,xe);
#endif

   return 0;
}


int inTFit_Fit::getBounds( double & xstart, double & xend )
{
   if( ibound != 1 ) {
      fprintf( stdout, " Error: bounds of fit have NOT been set \n");
      return 1;
   }

   xstart = xs;
   xend = xe;

   return 0;
}


int inTFit_Fit::setTerms( const char* term_string )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Setting terms via string \n",CLASS);
   fprintf( stdout, "    Terms: \"%s\" \n", term_string );
#endif

   int iret = parseTermsString( (const char*) term_string );
   if( iret != 0 ) {
      fprintf( stdout, " Error: could not parse terms! \n");
      return 1;
   }

   return 0;
}


int inTFit_Fit::getNumTerms( void ) const
{
   return num_terms;
}


int inTFit_Fit::finalize()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Method \"finalize()\" invoked \n",CLASS);
#endif

   if( ibound != 1 ) {
      fprintf( stdout, " Error: Fit structure has no bounds set!\n");
      return 1;
   }

   if( num_terms <= 0 ) {
      fprintf( stdout, " Error: Fit structure has no terms! \n");
      return 2;
   }
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Fit has %d terms \n",CLASS, num_terms );
#endif

   // check for whether monomials repeat and whether other terms are unique
   int iret=0;
   int ihave_const=0, ihave_log=0, ihave_exp=0;
   for(int n=0;n<num_terms && iret == 0;++n) {
      if( terms[n].getType() == TERM_MONOMIAL ) {
         int order1 = terms[n].getOrder();
         if( order1 == 0 ) ++ihave_const;
         for(int m=0;m<n && iret == 0;++m) {
            if( terms[m].getType() == TERM_MONOMIAL ) {
               int order2 = terms[m].getOrder();
               if( order1 == order2 ) iret = 1;
            }
         }
      }
      if( terms[n].getType() == TERM_CONST ) ++ihave_const;
      if( terms[n].getType() == TERM_LOG ) ++ihave_log;
      if( terms[n].getType() == TERM_EXP ) ++ihave_exp;
   }

   if( iret != 0 ) {
      fprintf( stdout, " Error: fit contains duplicate monomials \n");
      return 100;
#ifdef _DEBUG_
   } else {
      fprintf( stdout, " [%s]  Monomials are unique \n",CLASS);
#endif
   }

   if( ihave_const > 1 ) {
      fprintf( stdout, " Error: fit contains duplicate cosnt terms \n");
      return 200;
#ifdef _DEBUG_
   } else {
      fprintf( stdout, " [%s]  Constant term is unique \n",CLASS);
#endif
   }

   if( ihave_log > 1 ) {
      fprintf( stdout, " Error: fit contains duplicate log() terms \n");
      return 300;
#ifdef _DEBUG_
   } else {
      fprintf( stdout, " [%s]  Logarithm term is unique \n",CLASS);
#endif
   }

   if( ihave_exp > 1 ) {
      fprintf( stdout, " Error: fit contains duplicate exp() terms \n");
      return 400;
#ifdef _DEBUG_
   } else {
      fprintf( stdout, " [%s]  Exponential term is unique \n",CLASS);
#endif
   }

   return 0;
}


void inTFit_Fit::clear( void )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Clearing...\n",CLASS);
#endif
   terms.clear();
   num_terms = 0;
// values.clear();

   xs = -9.9e99;
   xe = -9.9e99;
   ibound = 0;
}


int inTFit_Fit::parseTermsString( const char* string )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Parsing terms: \"%s\" \n",CLASS, string );
#endif

   if( string == NULL ) return 1;

   int ilen = strlen( string );
   char *data = (char *) malloc( (size_t) ilen+1 );
   if( data == NULL ) return -1;
   strncpy( data, string, (size_t) ilen+1 );

   // parse terms in the string
   int iloop = 0, iret;
   char *s1,*s2,*save;
   for( s1 = data, s2 = data; s2 != NULL && iloop == 0; s1 = NULL ) {
#ifdef _DEBUG2_
      printf("data= \"%s\"  \n",data);//HACK
#endif
      s2 = strtok_r( s1, " ", &save );
#ifdef _DEBUG2_
      printf("S1= \"%s\"   S2= \"%s\" \n", s1,s2 );   // HACK
usleep(200000);//HACK
#endif

      if( s2 != NULL ) {
         inTerm_type type = TERM_NULL;
         int order = -99999;
         iret = parseTerm( s2, &type, &order );
         if( iret < 0 ) {
            iloop = 2;        // to potentially trap errors this way
         } else {
            iret = registerTerm( type, order );
            if( iret != 0 ) {
               iloop = 3;     // top potentially trap errors this way
            }
         }
      } else {
         iloop = 1;        // done parsing terms
      }
   }

   free( data );

   return 0;
}


int inTFit_Fit::parseTerm( const char* string,
                           inTerm_type* type_, int* order_ )
{
#ifdef _DEBUG2_
   fprintf( stdout, " [%s]  Parsing term: \"%s\" \n",CLASS,string);
#endif
   if( string == NULL ) return -100;

   int ilen = strlen( string );

   inTerm_type type = TERM_NULL;
   int order = -9999;

///////////////////////////////////////////////////////////////////////////////
///////   DANGER   ////////////////////////////////////////////////////////////
/////// Avoiding the nesting of code blocks by throwing GOTOs.
/////// This is so that I can move on for now and I will come fix this later.
/////// Violating my own rules of not ever using GOTOs, ever!
///////////////////////////////////////////////////////////////////////////////
   // level one query is for the type of term
goto level_one;   // avoid compiler warnings
level_one:
   if( string[0] == '1' ) {
      type = TERM_CONST;
#ifdef _DEBUG2_
   fprintf( stdout, "   Identified as \"CONST\" \n");
#endif
      goto level_two_const;
   } else if( string[0] == 'X' ) {
      type = TERM_MONOMIAL;
#ifdef _DEBUG2_
   fprintf( stdout, "   Identified as \"MONOMIAL\" \n");
#endif
      goto level_two_monomial;
   } else if( string[0] == 'L' ) {
      type = TERM_LOG;
#ifdef _DEBUG2_
   fprintf( stdout, "   Identified as \"NATURAL LOGARITHM\" \n");
#endif
      goto level_two_log;
   } else if( string[0] == 'E' ) {
      type = TERM_EXP;
#ifdef _DEBUG2_
   fprintf( stdout, "   Identified as \"EXPONENTIAL\" \n");
#endif
      goto level_two_exp;
   } else {
      return -200;       // did not find something recognizeable
   }

   // level two queries are for consistency of terms
level_two_const:
   if( ilen > 1 ) {   // specified something longer than "1"
      return -301;
   }
   goto level_three;
   return -999;   // guard from any GOTO fuckup by throwing bad return...

level_two_monomial:
   // I cannot of any consistency checks that I can use here...
   goto level_three;
   return -999;   // guard from any GOTO fuckup by throwing bad return...

level_two_log:
   if( ilen < 5 ) {   // acceptable is a minimum of "LN(X)", which is 5 chars
      return -501;
   }
   if( string[1] == 'N' || string[1] == 'O' ) {   // only accept "LN" and "LOG"
      // for now assume that it all went really well...
      goto level_three;
   } else {
      return -502;
   }
   return -999;   // guard from any GOTO fuckup by throwing bad return...

level_two_exp:
   if( ilen < 5 ) {   // acceptable is a minimum of "EXP(X)", which is 6 chars
      return -601;
   }
   if( strncmp( "EXP(X)", string, 6 ) == 0 ) {
      // for now assume that it all went really well...
      goto level_three;
   } else {
      return -602;
   }
   return -999;   // guard from any GOTO fuckup by throwing bad return...

   // level three queries are to pick up metadata for select terms
level_three:      // we got here because a known term type was identified
   if( type == TERM_MONOMIAL ) {
      if( ilen == 1 ) {
         order = 1;     // it must be "X"
      } else if( ilen < 3 ) {    // expect at least "X^N"
         return -701;
      } else {
         if( string[1] == '^' ) {
            // without guarding for crap following, do this mad thing:
            order = atoi( (const char*) &( string[2] ) );
#ifdef _DEBUG2_
            fprintf( stdout, "   Order \"%d\" \n", order );
#endif
         } else {
            return -702;
         }
      }
      goto level_four;
   } else {
      goto level_four;
   }
   return -999;   // guard from any GOTO fuckup by throwing bad return...

level_four:       // final stage (jsut in case)
   *type_ = type;
   *order_ = order;

   return 0;
}


int inTFit_Fit::registerTerm( inTerm_type type_, int order_ )
{
   inTFit_Term t( type_, order_ );

   terms.push_back( t );
   ++num_terms;

   return 0;
}


__inline double inTFit_Fit::evalProd( int i, int j, double t ) const
{
   double result = terms[i].eval( t ) * terms[j].eval( t );

   return result;
}


__inline double inTFit_Fit::evalTerm( int i, double t ) const
{
   return terms[i].eval( t );
}


//
// Multi-fit class
//

inTFit_MultiFit::inTFit_MultiFit()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Constructor \n",CLASS);
#endif

   xs = -9.9e99;
   xe = -9.9e99;
   ibound = 0;

   num_ranges = 0;

   num_data = 0;

   ne = 0;
   amat = NULL;
   rhs = NULL;
   tp = NULL;
}


inTFit_MultiFit::~inTFit_MultiFit()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Deconstructing \n",CLASS);
#endif
   if( amat == NULL ) free( amat );
   if( rhs == NULL ) free( rhs );
   if( tp == NULL ) free( tp );
}


int inTFit_MultiFit::setBounds( double xstart, double xend )
{
   if( ibound != 0 ) {
      fprintf( stdout, " Error: bounds have already been set \n");
      fprintf( stdout, "        xs = %16.9e, xe = %16.9e \n", xs, xe );
      return 1;
   }

   xs = xstart;
   xe = xend;
   ibound = 1;
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Bounds: %16.9e to %16.9e \n",CLASS, xs,xe);
#endif

   return 0;
}


int inTFit_MultiFit::addFit( inTFit_Fit & fit_ )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Method \"addFit()\" invoked \n",CLASS);
#endif
   if( ibound != 1 ) {
      fprintf( stdout, " Error: bounds of multi-fit have NOT been set \n");
      return 1;
   }

   double xstart,xend;
   fit_.getBounds( xstart, xend );

   if( xstart < xs || xe < xend ) {
      fprintf( stdout, " Error: incompatible bounds! \n");
      fprintf( stdout, "        Multi-fit bounds: %16.9e : %16.9e \n", xs, xe );
      fprintf( stdout, "              Fit bounds: %16.9e : %16.9e \n",
               xstart, xend );
      return 2;
#ifdef _DEBUG_
   } else {
      fprintf( stdout, " [%s]  Bounds are compatible! \n",CLASS);
#endif
   }

   fits.push_back( fit_ );
   lows.push_back( xstart );
   highs.push_back( xend );
   ++num_ranges;

   return 0;
}


int inTFit_MultiFit::finalize()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Method \"finalize()\" invoked \n",CLASS);
#endif
   ne=0;
   for(int n=0;n<num_ranges;++n) {
      ne += fits[n].getNumTerms();
   }
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Using %d terms \n",CLASS, ne );
#endif

   // allocate memory for the linear system
   amat = (double *) malloc( ((size_t) (ne*ne)) * sizeof(double) );
   rhs = (double *) malloc( ((size_t) ne) * sizeof(double) );
   tp = (double *) malloc( ((size_t) ne) * sizeof(double) );
   if( amat == NULL || rhs == NULL || tp == NULL ) {
      if( amat != NULL ) free( amat );
      amat = NULL;
      if( rhs != NULL ) free( rhs );
      rhs = NULL;
      if( tp != NULL ) free( tp );
      tp = NULL;
      return -1;
   }

   return 0;
}


int inTFit_MultiFit::storeData( long size,
                                const double* x_, const double* y_ )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Storing data \n",CLASS);
#endif
   x.clear();
   y.clear();

   if( size <= 0 ) return 1;

   x.resize( size );
   y.resize( size );
   for(long n=0;n<size;++n) {
      x[n] = x_[n];
      y[n] = y_[n];
   }
   num_data = size;

   (void) findIndices();

   return 0;
}


int inTFit_MultiFit::compute( void )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Computing...\n",CLASS);
#endif
   if( ne <= 0 || amat == NULL || rhs == NULL || tp == NULL ) return 1;

   (void) form();
   (void) inUtils_GaussianElimination( ne, amat, rhs, tp );
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Solution \n",CLASS);
   for(int n=0;n<ne;++n) fprintf( stdout, " n=%d   %16.9e \n", n, rhs[n] );
#endif


   return 0;
}


void inTFit_MultiFit::clear( void )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Clearing...\n",CLASS);
#endif
   fits.clear();
   lows.clear();
   highs.clear();
   idx_lo.clear();
   idx_hi.clear();
   num_ranges = 0;

   xs = -9.9e99;
   xe = -9.9e99;
   ibound = 0;

   x.clear();
   y.clear();
   num_data = 0;
}


int inTFit_MultiFit::findIndices()
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Finding indices for ranges \n",CLASS);
#endif

   if( num_data < 2 ) return 1;

   idx_lo.resize( num_ranges );
   idx_hi.resize( num_ranges );

   for(int n=0;n<num_ranges;++n) {
      double xlo = lows[n], xhi = highs[n];
#ifdef _DEBUG2_
fprintf( stdout, " Iteration for fit %d (low) \n",n);
#endif
      long il = 0, ih = num_data - 1;
      long itmp = il + (ih - il)/2;
      while( il != ih-1 ) {
#ifdef _DEBUG2_
fprintf( stdout, " XL: %.6ld   ITMP: %.6ld   XH: %.6ld \n",il,itmp,ih);
usleep(300000);
#endif
         if( xlo < x[ itmp ] ) {
            ih = itmp;
         } else {
            il = itmp;
         }
         itmp = il + (ih - il)/2;
      }
#ifdef _DEBUG2_
fprintf( stdout, " XL: %.6ld   ITMP: %.6ld   XH: %.6ld \n",il,itmp,ih);
usleep(300000);
#endif
      if( x[ il ] >= lows[n] ) {
         idx_lo[n] = il;
      } else {
         idx_lo[n] = ih;
      }

#ifdef _DEBUG2_
fprintf( stdout, " Iteration for fit %d (high) \n",n);
#endif
      il = 0, ih = num_data - 1;
      itmp = il + (ih - il)/2;
      while( il != ih-1 ) {
#ifdef _DEBUG2_
fprintf( stdout, " XL: %.6ld   ITMP: %.6ld   XH: %.6ld \n",il,itmp,ih);
usleep(300000);
#endif
         if( x[ itmp ] < xhi ) {
            il = itmp;
         } else {
            ih = itmp;
         }
         itmp = il + (ih - il)/2;
      }
#ifdef _DEBUG2_
fprintf( stdout, " XL: %.6ld   ITMP: %.6ld   XH: %.6ld \n",il,itmp,ih);
usleep(300000);
#endif
      if( x[ ih ] <= highs[n] ) {
         idx_hi[n] = ih;
      } else {
         idx_hi[n] = il;
      }

   }
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Summary of indices \n",CLASS);
#endif
   for(int n=0;n<num_ranges;++n) {
      fprintf( stdout, "  %d:  [%6ld <-> %6ld ] \n", n, idx_lo[n], idx_hi[n] );
   }

   return 0;
}


int inTFit_MultiFit::form( void )
{
#ifdef _DEBUG_
   fprintf( stdout, " [%s]  Forming...\n",CLASS);
#endif
   // clean-up data
   for(int n=0;n<ne*ne;++n) amat[n] = 0.0;
   for(int n=0;n<ne;++n) rhs[n] = 0.0;

   // sweep over ranges (fit segments)
   int ntt=0;    // term offset
   for(int n=0;n<num_ranges;++n) {
      int nt = fits[n].getNumTerms();

      // sweep over the points of this segment
      inTFit_Fit* fit = &( fits[n] );
      long is = idx_lo[n], ie = idx_hi[n];
#ifdef _DEBUG3_
      fprintf( stdout, "    Segment %d, index bounds %ld, %ld \n", n, is, ie );
#endif
      for(long i=is;i<ie;++i) {
         double t = x[i];

         for(int m=0;m<nt;++m) {
#ifdef _DEBUG2_
            if( i == is )
               for(int k=0;k<ntt;++k)
                  fprintf( stdout, "- ");
#endif
            for(int k=0;k<nt;++k) {
#ifdef _DEBUG2_
               if( i == is )
                  fprintf( stdout, "X ");
#endif
               // assign matrix element
               double ee = fit->evalProd( m, k, t );
               amat[ ntt*ne + ne*m + ntt + k ] += ee;
            }
            // assign RHS vector element
            double rr = y[i] * fit->evalTerm( m, t );
            rhs[ ntt + m ] += rr;
#ifdef _DEBUG2_
            if( i == is )
               for(int k=ntt+nt;k<ne;++k)
                  fprintf( stdout, "- ");
            if( i == is )
               fprintf( stdout, "\n");
#endif
         }

      }

      // advance term offset
      ntt += nt;
   }
#ifdef _DEBUG3_
   fprintf( stdout, " [%s]  Linear system \n",CLASS);
   for(int n=0;n<ne;++n) {
      for(int m=0;m<ne;++m) {
         fprintf( stdout, " %16.9e", amat[ n*ne + m ] );
      }
      fprintf( stdout, "    %16.9e\n", rhs[ n ] );
   }
#endif

   return 0;
}


#ifdef __cplusplus
}
#endif

