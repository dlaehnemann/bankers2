# bankers2
A set of tools to create plotting files from subset counts in banker's sequenc order. I.e. helper scripts for the visualisation of count data across all possible subsets of a sample set.

## bankers2circos
About:   Create files for a circos plot representing all subsets/intersections contained in a banker's sequence order file.

Usage:   bankers2circos [options] <banker's-seq-subsets-file>

Options:
 -p, --prefix <prefix>    prefix for circos file names [default: prefix of input without path]
 -m, --missing            if set, include count of missing genotypes per sample in output, including the circos files, if they are requested
 -h, --help               this help message

Input:
 1) Header line specifying samples: '@SMPS SMP1,SMP2,...' (required)
     This header line needs to start with @SMPS, followed by a tab or a whitespace. Then comes a list of sample names in the order they appear, separated by comma, tab or whitespace.
 2) Count lines just contain one number per line. (required)
     They should appear in banker's sequence order with regard to the @SMPS sample order. If one missing count per sample is included, these values are in the first #samples lines.
 3) Comment lines starting with '#'. (optional)
     These lines are meant to document details of how the counts were generated and are ignored.
