# INTFIT

This is a code to generate some kind of fit to some kind of data; nothing
to see here. Data is read from a file to produce the fit, and the data must
be formatted in two columns per row (like a scatter plot). How the data is
fit in terms of the functional form is done programmatically (this is the
power of this library code). The data can also be fitted with multiple fit
segments, and the segments can further be constrained (among eachother or
independently) to produce a very specific global fit over the domain.

- Motivation

I am tired of modifying old code to do new things when it comes to fitting
funcitonal forms to thermodynamic data that adopt the format of the NASA
CEA database. This is an abstraction that allows to programmatically
synthesize a functional form that comprises of segments, each wtih specific
terms, to use for fitting one-dimensional data. It is meant to produce a
"best fit" in a least-squares sense. Specific constraints can also be added
to the system to achieve specific behavriour, for example to make the global
fit continuous and to have a continuous derivative everywhere.

- What this software does and how to use it

You should not use it, because you have no idea what this is for, and what
the underlying assumptions are. (You should contact me if you really need
to use it.) But should you insist, look at the example provided in the
driver. You must adjust the filename for providing the data in columnar
form, you must adjust the range in the domain over which the data is valid,
you need to programmatically specify the segmented fit object by adding
multiple segments as needed, you must specify constraints if you have any,
and you can finally look at the approximation that is created.

You are supposed to try and understand what the software does at the high
level and then extract the coefficients that represent the fit. (In the
future I _may_ augment the code to spit out this information in text form.)

- How it works

The software requires that one or more local fits are specified and then they
are all made part of a multi-segment fit. Each of the segments is specified
by establishing a range (on the x axis) where the curve will fit the data (on
the y axis), and a number of functions that will become the terms of the
approximation for the segment. The terms are provided via a string with
(well-structured) terms that are recognizeable by the parser. Those can
presently be: a constant (what multiplies "1"), monomials of any degree
--including inverses-- like ("X^N" type), the natural logarithm "log(x)".
and the exponential "exp(x)".

Constraints for making the global curve and its derivative continuous
yet been implemented. Constraints can presently be specified for forming
an equation that is meant to constrain either the value of the fit at a
given point or the derivative. For example, say you have data that is to
be fitted in a least-squares sense over two ranges (splitting your domain
in two). Further suppose that you want the values of each fit to evaluate
to exacty the same value at a specific point (implying continuity of the
fit at that point). You can do so by specifying a single "value" constraint
that involves fit "0" and fit "1". In this case you will also need to specify
the "sign" to be "-1" and the "right-hand side" to be "0.0. This is because
the constraints are imposed via an equation that we should read as "the value
of fit 0 at point x plus the sign-multiplied value of the fit 1 at point x
is quual to right-hand side". That is: "fit0_at_x - sign * fit1_at_x = rhs".
This allows for jumps in the value or the derivative via the inhomogeneous
term, as well as other tricks.


TO BE CONTINUED...

IN 2021/01/19

