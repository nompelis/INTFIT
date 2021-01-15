# INTFIT

This is a code to generate some kind of fit to some kind of data; nothing
to see here. (THIS IS NOT A WORKING VERSION EITHER!)

- Motivation

I am tired of modifying old code to do new things when it comes to fitting
funcitonal forms to thermodynamic data that adopt the format of the NASA
CEA database.

- What this software does and how to use it

You should not use it, because you have no idea what this is for, and what
are the underlying assumptions.

- How it works

The software requires that one or more local fits are specified and then they
are all made part of a multi-segment fit. Each of the segmetns is specified
by establishing a range (on the x axis) where the fit will fit the data (on
the y axis). Then a number of functions that will become the terms of the
approximation are to be specified for a given segment fit; this is done by
providing a string with (well-structured) terms. Those can presently be:
a constant (what multiplies "1"), monomials of any degree, including inverses,
like ("X^N" type), the natural logarithm "log(x)". and the exponential "exp(x)".




TO BE CONTINUED...

IN 2021/01/13

