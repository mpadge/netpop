# netpop
Simluations of populations reproducing and moving throughout different kinds of
networks of 2, 3, 4, and 25 nodes.

Two nodes are theoretical values calculated with an R script; all others are C++
simulations. Three-node networks are triangular; there are four different
four-node networks; and there are two different kinds of 25-node networks:
dendritic and lattice, with their own separate routines.

### NOTE:

These are the final versions used to produce the figures for initial submission
(24th Feb, 2016) - commit#5ec4c1c

### build:
1. cd ./build  
2. cmake ..  
3. make


[![Build
Status](https://travis-ci.org/mpadge/netpop.svg?branch=master)](https://travis-ci.org/mpadge/netpop)
