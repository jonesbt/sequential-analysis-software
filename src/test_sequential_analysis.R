#!/usr/bin/env Rscript

## This file contains test code for the R wrapper and also serves as an
## incomplete reference to the R bindings. To run the tests, execute the below
## commands from within R.
##
## library(testthat)
## test_file('test_sequential_analysis.R')

library(testthat)

dyn.load(paste("sequential_analysis_wrap", .Platform$dynlib.ext, sep=""))
source("sequential_analysis_wrap.R")
cacheMetaData(1)

## Create the origins.
origins = StringVector()
StringVector_push_back(origins, "1")
StringVector_push_back(origins, "2")
StringVector_push_back(origins, "3")
## Create the destinations.
destinations = StringVector()
StringVector_push_back(destinations, "Mid-Maine")
StringVector_push_back(destinations, "Three States")
StringVector_push_back(destinations, "Mass Bay")
StringVector_push_back(destinations, "Nantucket")
StringVector_push_back(destinations, "Other")
## Create the matrix of counts.
counts = vector(mode='list', len=StringVector_size(origins))
counts[[1]] = c(5, 3, 15, 2, 142)
counts[[2]] = c(0, 1, 30, 1, 135)
counts[[3]] = c(0, 0, 57, 2, 107)
## Create the connectivity matrix.
cm = ConnectivityMatrix(origins, destinations)
## Check that the initial counts are correct.
new_counts = ConnectivityMatrix_get_counts(cm)
expect_that(length(new_counts), equals(3))
for(i in seq(length(new_counts)))
    expect_that(length(new_counts[[i]]), equals(5))
for(i in seq(length(new_counts)))
    for(j in seq(length(new_counts[[i]])))
        expect_that(new_counts[[i]][j], equals(0))
ConnectivityMatrix_set_obj_fn_cv_args(cm, 0.005, 0.05)
## Allocate the particles and check that the allocation is correct.
alloc = ConnectivityMatrix_allocate(cm, 500)
expect_that(length(alloc), equals(3))
## Update the matrix and check that the update was correct.
ConnectivityMatrix_update(cm, counts)
new_counts = ConnectivityMatrix_get_counts(cm)
for(i in seq(length(counts)))
    expect_that(new_counts[[i]], is_equivalent_to(counts[[i]]))
## Reallocate particles and test that the allocation is correct.
alloc = ConnectivityMatrix_allocate(cm, 500)
expect_that(alloc, is_equivalent_to(c(41, 229, 230)))



