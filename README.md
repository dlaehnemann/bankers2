# bankers2
A set of tools to create plotting files from subset counts in banker's sequenc order. I.e. helper scripts for the visualisation of count data across all possible subsets of a sample set.

## Install

1. clone the repository
2. move to the directory you cloned to
3. type `make`
4. make the created binaries of the tools mentioned below available in your `$PATH`

## bankers2circos
This tool automatically generates all the files necessary for visualising the counts of all subsets in a file in banker's sequence order. It automatically assigns color using a color profile from the [Brewer palettes](http://colorbrewer2.org/), optimizes the layout of subset connection ribbons to minimize overlap and generates tick marks and tick mark labels in useful intervals. Darker shades of green signify an increasing number of samples sharing the respective counts' instances, with the lightest green with no ribbons giving instances unique to that sample. Optionally, counts of missing instances (e.g. a missing genotype when using input generated by 

### Usage / Help message
```
About:   Create files for a circos plot representing all subsets/intersections
         contained in a banker's sequence order file.

Usage:   bankers2circos [options] <banker's-seq-subsets-file>

Options:
    -p, --prefix <prefix>    prefix for circos file names [default: prefix of input without path]
    -m, --missing            if set, include count of missing genotypes per sample in
                             output, including the circos files, if they are requested
    -h, --help               this help message

Input:
    1) Header line specifying samples: '@SMPS SMP1,SMP2,...' (required)
         This header line needs to start with @SMPS, followed by a tab or a whitespace. 
         Then comes a list of sample names in the order they appear, separated by comma, tab or whitespace.
    2) Count lines just contain one number per line. (required)
         They should appear in banker's sequence order with regard to the @SMPS sample order.
         If one missing count per sample is included, these values are in the first #samples lines.
    3) Comment lines starting with '#'. (optional)
         These lines are meant to document details of how the counts were generated and are ignored.
```
### Example
With [circos](http://circos.ca/) [installed](http://circos.ca/tutorials/lessons/configuration/distribution_and_installation/), using the following commands and the file `test/view.gtisec.out`, will create the figure:

    # create all the necessary circos files with prefix view.gtisec.
    bankers2circos view.gtisec.out
    # plot using circos
    circos -conf view.gtisec.circos.conf -file view.gtisec

![alt tag](/test/view.gtisec.png)
