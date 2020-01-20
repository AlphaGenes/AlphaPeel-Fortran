# README #

### What is this repository for? ###

This is a multilocus peeler that (probably) implements Meuwissen and Goddard (2010). It extends geneprob by estimating segregation information at each locus and passing that segregation information to the next locus.

### Known Issues

1. Inconsistency between macOS and Linux and with using different number of threads.  

   * Linux and macOS; Fortran compiler 19 update 5; PURE functions fixed: different numbers of threads give **different** results
   * Fortran compiler 19 update 5; PURE functions fixed; `-O3`: macOS and Linux versions give **different** results
   * Fortran compiler 19 update 5; PURE functions fixed; `-O0`: macOS and Linux versions give **identical** results for same number of threads
   
   With default compilation options (`-O3`, `-qopenmp` etc). Fixed PURE functions for 2019 compiler only. Values from `compare_files.r`:
   * Linux; Fortran compiler 17 update 4; 1 thread vs 4 threads: max diff 0.9886; av diff -7.16544850498341e-06
   * Linux; Fortran compiler 19 update 5; 1 thread vs 4 threads: max diff 0.2273; av diff -1.61046511627906e-06
   * Mac; Fortran compiler 19 update 5; 1 thread vs 4 threads: max diff 0.3794; av diff 1.04086378737541e-06