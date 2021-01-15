# INTFIT

This is a code to generate some kind of fit to some kind of data; nothing
to see here. (THIS IS NOT A WORKING VERSION EITHER!)

- Motivation

I am tired of modifying old code to do new things when it comes to fitting
funcitonal forms to thermodynamic data that adopt the format of the NASA
CEA database. This is an abstraction that allows to programmatically
synthesize a functional form that comprises of segments, each wtih specific
terms, to use for fitting one-dimensional data. It is meant to produce a
"best fit" in a least-squares sense. Specific constraints can also be added
to the system to achieve specific behavriour.

- What this software does and how to use it

You should not use it, because you have no idea what this is for, and what
the underlying assumptions are. (You should contact me if you really need
to use it.)

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

The constraints for making the global curve and its derivative continuous
have not yet been implemented. At present, the software can build a multi-fit
with many fits, but those are independent and not constrained. It also does
not provide any meangful output, yet.



TO BE CONTINUED...

IN 2021/01/15

