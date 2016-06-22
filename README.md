# MSBWT-IS
## Introduction
MSBWT-IS is a tool for building BWTs of string collections using induced sorting.  In general, the time and memory usage of MSBWT-IS scales linearly with the length of the input string collection.  As a result, characteristics of the string collection like string length have low impact on the performance of the algorithm.  This implementation excels at building BWTs for genomic, long-read sequencing datasets like those produced by PacBio or nanopore sequencing technologies.  Additionally, the source code can be adapted to suit any type of string collection (aka, any alphabet) with relative ease and without impacting performance.

## Installation and Setup
First, download the latest version of MSBWT-IS and unzip it.  Then simply make the program and run it with no parameters to verify it installed.

    cd msbwt-is
    make
    ./msbwtis -h

## Building a BWT with MSBWT-IS
MSBWT-IS currently requires an output BWT directory and a list of FASTA formatted files as input.

    ./msbwtis <outputDirectory> <input1.fa> ...

## Using the BWT output
MSBWT-IS outputs the BWT in an uncompressed numpy format that is recognized by the [*msbwt*](https://github.com/holtjma/msbwt) python package.  Once built, we recommend visiting the [*msbwt* wiki pages](https://github.com/holtjma/msbwt/wiki) for information on compressing or querying the BWT.