
<!-- README.md is generated from README.Rmd. Please edit that file -->
MODclustsig
===========

This package is a redistribution of the original [clustsig package](https://cran.r-project.org/web/packages/clustsig/index.html) written by [Douglas Whitaker](http://www.douglaswhitaker.com) and Mary Christman.

The changes I made were the following:

1.  I included the Hellinger Distance ("hellinger") as a method to calculate the distance matrix. In the original version of the package you can pass any function to the method.distance, but since we will need to work with Hellinger Distance quite often I thought that it was a good idea to have it included by default.
2.  The Hellinger Distance matrix is calculated in C++ to speed up the process a little bit. This means that the package now imports form Rcpp evalCpp.
3.  I made a class simprof to handle the result of the simprof function and implemented the simprof.plot as a method for this class (instead of just a function as before).

Everything else remains the same. If you intend to use this redistribution for anything else but passing the Hellinger Distance to hclust, I think it would be better to use the package original versiona. Anyways, you are free to use it if you want. I'll be available to answer any questions that you may have.

Thank you Douglas Whitaker and Mary Christman, you save my day.

Installation
------------

Install the current development from github via:

``` r
remotes::install_github("fcorra/mod_clustsig")
```
