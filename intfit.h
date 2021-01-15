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

#ifndef _INTFIT_H_
#define _INTFIT_H_

#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <vector>
#include <map>
#include <string>


#ifdef __cplusplus
extern "C" {
#endif

enum inTerm_type {
   TERM_NULL,
   TERM_CONST,
   TERM_MONOMIAL,
   TERM_LOG,
   TERM_EXP,
};


//
// object representing a term in the approximation
//

class inTFit_Term {
 public:
   inTFit_Term( inTerm_type type_, int order_ );
   ~inTFit_Term();

   int getOrder( void ) const;
   inTerm_type getType( void ) const;
   double eval( double t ) const;
   double evalDer( double t ) const;
   double evalInt( double t ) const;

 protected:

 private:
   inTerm_type type;
   int order;

};


//
// object representing a fit to data that comprises of the degrees of
// freedom and the terms in the approximation
//

class inTFit_Fit {
 public:
   inTFit_Fit();
   ~inTFit_Fit();

   int setBounds( double xstart, double xend );
   int getBounds( double & xstart, double & xend );
   int setTerms( const char* term_string );
   int getNumTerms( void ) const;
   int finalize();
   void clear( void );
   double evalProd( int i, int j, double t ) const;

 protected:

 private:
   double xs,xe;
   int ibound;

   std::vector< inTFit_Term > terms;
   int num_terms;
// std::vector< double > values;

   // methods
   int parseTermsString( const char* );
   int parseTerm( const char*, inTerm_type* type_, int* order_ );
   int registerTerm( inTerm_type type_, int order_ );

#ifdef _DEBUG_
   char* CLASS = (char*) "inTFit_Fit";
#endif
};


//
// object representing a multi-segment fit
//

class inTFit_MultiFit {
 public:
   inTFit_MultiFit();
   ~inTFit_MultiFit();

   int setBounds( double xstart, double xend );
   int addFit( inTFit_Fit & fit_ );
   int finalize();
   void clear( void );
   int storeData( long size, const double* x_, const double* y_ );

   int form( void );

 protected:

 private:
   double xs,xe;
   int ibound;

   std::vector< inTFit_Fit > fits;
   std::vector< double > lows;
   std::vector< double > highs;
   std::vector< long > idx_lo, idx_hi;
   int num_ranges;

   long num_data;
   std::vector< double > x,y;

   int findIndices();

   int ne;
   double *amat, *rhs;

#ifdef _DEBUG_
   char* CLASS = (char*) "inTFit_MultiFit";
#endif
};




#ifdef __cplusplus
}
#endif

#endif

