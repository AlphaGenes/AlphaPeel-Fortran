# README #

### What is this repository for? ###

This is a multilocus peeler that (probably) implements Meuwissen and Goddard (2010). It extends geneprob by estimating segregation information at each locus and passing that segregation information to the next locus.

### Known Issues

1. Inconsistency between macOS and Linux and with using different number of threads.  

   * Linux and macOS; Fortran compiler 19 update 5; PURE functions fixed: different numbers of threads give **different** results
   * Fortran compiler 19 update 5; PURE functions fixed; `-O3`: macOS and Linux versions give **different** results
   * Fortran compiler 19 update 5; PURE functions fixed; `-O0`: macOS and Linux versions give **identical** results for same number of threads