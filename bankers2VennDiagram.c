
#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include <inttypes.h>
#include <error.h>
#include <errno.h>

typedef struct
{
    char *prefix; /*! prefix of VennDiagram plotting file */
    int nsmp; /*! number of samples, determined from number of single-sample lines at beginning of file */
    int nsmpp2; /*! 2^(nsmp) (is needed multiple times) */
    char **smpns; /*! list of samples as specified in banker's sequence file of intersection counts */
    char *palette; /*! palette prefix for VennDiagram samples */
/*    uint32_t *bankers; /*! array to store banker's sequence for all possible sample subsets for
                                programmatic indexing into smp_is for output printing, e.g. for three
                                samples A, B and C this would be the following order:
                                [   C,   B,   A,  CB,  CA,  BA, CBA ]
                                [ 100, 010, 001, 110, 101, 011, 111 ]
                                */
    uint64_t *quick; /*! array to store lookup table of n choose k values from choose() */
    uint8_t missing; /*! flag, whether missing values should be extracted from banker's sequence file */
    uint64_t *smp_is; /*! array to track all possible intersections between
                 samples, with each bit in the index integer belonging to one
                 sample. E.g. for three samples A, B and C, count would be in
                 the following order:
                 [   A,   B,  AB,   C,  AC,  BC, ABC ]
                 [ 001, 010, 011, 100, 101, 110, 111 ]
                 */
    FILE *vdout; /*! VennDiagram plotting file handle */
}
args_t;


static void init_data(args_t *args)
{
    args->nsmpp2 = pow( 2, args->nsmp );
//    args->bankers = (uint32_t*) calloc( args->nsmpp2, sizeof(uint32_t) );
    args->quick = (uint64_t*) calloc((args->nsmp * (args->nsmp + 1)) / 4, sizeof(unsigned long));
    if ( args->nsmp > 5 )
    {
        fprintf(stderr, "Warning: VennDiagram can plot a maximum of five samples.\n");
        exit(EXIT_FAILURE);
    } else {
        args->palette = malloc(24);
        sprintf(args->palette, "brewer.pal(%i, \"Set1\" )[", args->nsmp);
    }
    char fn[256];
    sprintf(fn, "%s.R", args->prefix);
    args->vdout = fopen(fn, "w");
    args->smp_is = (uint64_t*) calloc( args->nsmpp2-1, sizeof(uint64_t));

    /* Compute banker's sequence for following systematic input reading. */
/*    int i;
    for ( i = 0; i < args->nsmpp2; i++ )
    {
        args->bankers[i] = compute_bankers(i, args);
    }
*/
}

static void destroy_data(args_t *args)
{
    free(args->smp_is);
//    free(args->bankers);
    free(args->quick);
    int i;
    for ( i = 0; i <= args->nsmp; i++ )
    {
        free(args->smpns[i]);
    }
    free(args->smpns);
    free(args->prefix);
    free(args->palette);
    fclose(args->vdout);
    free(args);
}


/* ADAPTED CODE FROM CORIN LAWSON (START)
 * https://github.com/au-phiware/bankers/blob/master/c/bankers.c
 * who implemented ideas of Eric Burnett:
 * http://www.thelowlyprogrammer.com/2010/04/indexing-and-enumerating-subsets-of.html
 */

/*
 * Compute the binomial coefficient of `n choose k'.
 * Use the fact that binom(n, k) = binom(n, n - k).
 * Use a lookup table (triangle, actually) for speed.
 * Otherwise it's dumb (heart) recursion.
 * Added relative to Corin Lawson:
 * * Passing in of sample number through pointer to args struct
 * * quick lookup table is now maintained externally in args struct
 */
/*
uint64_t choose(unsigned int n, unsigned int k, args_t *args) {
    if (n == 0)
        return 0;
    if (n == k || k == 0)
        return 1;
    if (k > n / 2)
        k = n - k;

    unsigned int i = (n * (n - 1)) / 4 + k - 1;
    if (args->quick[i] == 0)
        args->quick[i] = choose(n - 1, k - 1, args) + choose(n - 1, k, args);

    return args->quick[i];
}
*/
/*
 * Returns the Banker's number at the specified position, a.
 * Derived from the recursive bit flip method.
 * Added relative to Corin Lawson:
 * * Uses same lookup table solution as choose function, just
 *   maintained externally to persist across separate function calls.
 * * Uses bitwise symmetry of banker's sequence to use bitwise inversion
 *   instead of recursive bit flip for second half of sequence.
 */
/*
uint32_t compute_bankers(unsigned long a, args_t *args)
{
    if (a == 0)
        return 0;

    if ( args->bankers[a] == 0 )
    {
        if ( a >= (args->nsmpp2 / 2) )
        {
            return args->bankers[a] = ( compute_bankers(args->nsmpp2 - (a+1), args) ^ (args->nsmpp2 - 1) ); // use bitwise symmetry of bankers sequence
        }
        unsigned int c = 0;
        uint32_t n = args->nsmp;
        uint64_t e = a, binom;
        binom = choose(n, c, args);
        do {
            e -= binom;
        } while ((binom = choose(n, ++c, args)) <= e);

        do {
            if (e == 0 || (binom = choose(n - 1, c - 1, args)) > e)
                c--, args->bankers[a] |= 1;
            else
                e -= binom;
        } while (--n && c && ((args->bankers[a] <<= 1) || 1));
        args->bankers[a] <<= n;
    }

    return args->bankers[a];
}
*/
// CODE BY CORIN LAWSON END

/* DEBUG PRINTING OF TWO-DIMENSIONAL ARRAY is_cmp_cnt
static inline print_isc(args_t *args, uint64_t *is_smp_cnt)
{
    int s1, s2;
    fprintf(stderr, "Samples");
    for ( s2 = 0; s2 < args->nsmp; s2++)
    {
        fprintf(stderr, "\t%s", args->smpns[s2]);
    }
    fprintf(stderr, "\n");
    for ( s1 = 0; s1 < args->nsmp; s1++ )
    {
        fprintf(stderr, "%s", args->smpns[s1]);
        for ( s2 = 0; s2 < args->nsmp; s2++ )
        {
            fprintf(stderr, "\t%"PRIu64"", *(is_smp_cnt + s1*args->nsmp + s2));
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}
*/

int bankersMapping(int n, args_t *args)
{
    switch(args->nsmp)
    {
    case 3:
        switch(n)
        {
        case 0: return 0; break;
        case 1: return 3; break;
        case 2: return 1; break;
        case 3: return 4; break;
        case 4: return 6; break;
        case 5: return 5; break;
        case 6: return 2; break;
        } break;
    case 4:
        switch (n)
        {
        case  0: return  2; break;
        case  1: return  9; break;
        case  2: return  3; break;
        case  3: return  5; break;
        case  4: return 12; break;
        case  5: return 14; break;
        case  6: return 13; break;
        case  7: return  8; break;
        case  8: return  0; break;
        case  9: return  6; break;
        case 10: return 11; break;
        case 11: return 10; break;
        case 12: return  7; break;
        case 13: return  1; break;
        case 14: return  4; break;
        } break;
    case 5:
        switch (n)
        {
        case  0: return  0; break;
        case  1: return  1; break;
        case  2: return  2; break;
        case  3: return  3; break;
        case  4: return  4; break;
        case  5: return 13; break;
        case  6: return  8; break;
        case  7: return  7; break;
        case  8: return  5; break;
        case  9: return 11; break;
        case 10: return  9; break;
        case 11: return  6; break;
        case 12: return 12; break;
        case 13: return 10; break;
        case 14: return 14; break;
        case 15: return 24; break;
        case 16: return 19; break;
        case 17: return 20; break;
        case 18: return 16; break;
        case 19: return 17; break;
        case 20: return 22; break;
        case 21: return 15; break;
        case 22: return 18; break;
        case 23: return 21; break;
        case 24: return 23; break;
        case 25: return 29; break;
        case 26: return 28; break;
        case 27: return 27; break;
        case 28: return 26; break;
        case 29: return 25; break;
        case 30: return 30; break;
        } break;
    }
    return -1;
}

void print_VennDiagram_file(args_t *args)
{
    int n;
    char *vn;

    fprintf(args->vdout, "library(\"RColorBrewer\")\n");
    fprintf(args->vdout, "library(\"VennDiagram\")\n");
    switch(args->nsmp)
    {
    case 1: vn = "single"; break;
    case 2: vn = "pairwise"; break;
    case 3: vn = "triple"; break;
    case 4: vn = "quad"; break;
    case 5: vn = "quintuple"; break;
    }
    fprintf(args->vdout, "venn.plot <- draw.%s.venn( ",vn);
    fprintf(args->vdout, "direct.area = TRUE, area.vector = c( %"PRIu64"", args->smp_is[ bankersMapping(0, args)] );
    for ( n=1; n < args->nsmpp2-1; n++ )
    {
        fprintf(args->vdout, ", %"PRIu64"", args->smp_is[ bankersMapping(n, args) ] );
    }
    fprintf(args->vdout, " ), \n    category = c( \"%s\"", args->smpns[0]);
    for ( n=1; n < args->nsmp; n++ )
    {
        fprintf(args->vdout, ", \"%s\"", args->smpns[n] );
    }
    fprintf(args->vdout, " ), \n    lwd = rep(2, %i", args->nsmp);
    fprintf(args->vdout, " ), \n    lty = rep(\"solid\", %i", args->nsmp);
    fprintf(args->vdout, " ), \n    col = c( %s1]", args->palette );
    for ( n=1; n < args->nsmp; n++ )
    {
        fprintf(args->vdout, ", %s%i]", args->palette, n+1 );
    }
    fprintf(args->vdout, " ), \n    fill = c( %s1]", args->palette );
    for ( n=1; n < args->nsmp; n++ )
    {
        fprintf(args->vdout, ", %s%i]", args->palette, n+1 );
    }
    fprintf(args->vdout, " ), \n    cat.col = c( %s1]", args->palette );
    for ( n=1; n < args->nsmp; n++ )
    {
        fprintf(args->vdout, ", %s%i]", args->palette, n+1 );
    }
    fprintf(args->vdout, " ), \n    alpha = rep(0.5, %i", args->nsmp);
    fprintf(args->vdout, " ), \n    label.col = rep(\"black\", %i", args->nsmpp2-1);
    fprintf(args->vdout, " ), \n    cex = rep(.9, %i", args->nsmpp2-1);
    fprintf(args->vdout, " ), \n    fontface = rep(\"plain\", %i", args->nsmpp2-1);
    fprintf(args->vdout, " ), \n    fontfamily = rep(\"sans\", %i", args->nsmpp2-1);
    switch(args->nsmp)
    {
    case 1: vn = "0 "; break;
    case 2: vn = "-50, 50 "; break;
    case 3: vn = "-40, 40, 180 "; break;
    case 4: vn = "-15, 15, 0, 0 "; break;
    case 5: vn = "0, 287.5, 215, 145, 70 "; break;
    }
    fprintf(args->vdout, " ), \n    cat.pos = c( %s", vn );
    switch(args->nsmp)
    {
    case 1: vn = "c( 0.025"; break;
    case 2: vn = "c( 0.025, 2"; break;
    case 3: vn = "c( 0.05, 0.05, 0.025"; break;
    case 4: vn = "c( 0.22, 0.22, 0.11, 0.11"; break;
    case 5: vn = "rep( 0.1, 5"; break;
    }
    fprintf(args->vdout, " ), \n    cat.dist = %s", vn );
    fprintf(args->vdout, " ), \n    cat.cex = rep(1.5, %i", args->nsmp);
    fprintf(args->vdout, " ), \n    cat.just = rep(list(c(0.5, 0.5)), %i", args->nsmp);
    fprintf(args->vdout, " ), \n    cat.fontface = rep(\"plain\", %i", args->nsmp);
    fprintf(args->vdout, " ), \n    cat.fontfamily = rep(\"sans\", %i", args->nsmp);
    fprintf(args->vdout, " ), \n    rotation.degree = 0");
    fprintf(args->vdout, ", \n    rotation.centre = c(0.5, 0.5 ");
    fprintf(args->vdout, " ), \n    print.mode = c(\"raw\", \"percent\" ");
    fprintf(args->vdout, " ), \n    sigdigs = 2");
    fprintf(args->vdout, ", \n    ind = FALSE");
    fprintf(args->vdout, ", \n    cex.prop = NULL \n");
    fprintf(args->vdout, ")\n");
    fprintf(args->vdout, "\npdf(\"%s.pdf\", width = 12, height = 12)", args->prefix);
    fprintf(args->vdout, "\ngrid.draw(venn.plot)");
    fprintf(args->vdout, "\ndev.off()\n");
}

static int usage(int status)
{
    fprintf(stderr, "\n"
                    "About:   Create file for a plot representing all subsets/intersections contained in\n"
                    "         a banker's sequence order file, using the R package VennDiagram.\n"
                    "\n"
                    "Usage:   bankers2VennDiagram [options] <banker's-seq-subsets-file>\n"
                    "\n"
                    "Options:\n"
                    "    -p, --prefix <prefix>    prefix for circos file names [default: prefix of input without path]\n"
                    "    -h, --help               this help message\n"
                    "\n"
                    "Input:\n"
                    "    1) Header line specifying samples: '@SMPS SMP1,SMP2,...' (required)\n"
                    "         This header line needs to start with @SMPS, followed by a tab or a whitespace. \n"
                    "         Then comes a list of sample names in the order they appear, separated by comma, tab or whitespace.\n"
                    "    2) Count lines just contain one number per line. (required)\n"
                    "         They should appear in banker's sequence order with regard to the @SMPS sample order.\n"
                    "         If one missing count per sample is included, these values are in the first #samples lines.\n"
                    "    3) Comment lines starting with '#'. (optional)\n"
                    "         These lines are meant to document details of how the counts were generated and are ignored.\n"
                    "\n");
    exit(status);
}

int main(int argc, char *argv[])
{
    int c;

    if ( argc < 2 ) usage(EXIT_FAILURE);

    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->prefix = NULL;
    args->smpns = calloc(32, sizeof(char*));

    static struct option loptions[] =
    {
        {"help",        no_argument,      0,'h'},
        {"prefix",      required_argument,0,'p'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "p:h",loptions,NULL)) >= 0) {
        switch (c) {
            case 'p': args->prefix = optarg; break;
            case 'h':
            case '?': usage(EXIT_SUCCESS);
            default: error(EXIT_FAILURE, EINVAL, "Unknown argument: %s\n", optarg);
        }
    }

    char *fname = argv[optind];
    if ( !args->prefix )
    {
        char *bn = basename(fname);
        char *o = strrchr(bn, '.');
        args->prefix = calloc(o - bn + 1, sizeof(char));
        strncpy(args->prefix, bn, o - bn);
    }
    if ( argc>optind+1 )  usage(EXIT_FAILURE);  // too many files given
        FILE *in = fopen(fname, "r");
        if ( in == NULL ) error(EXIT_FAILURE, ENOENT, "Could not open file %s\n", fname);

        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        uint64_t j = 0;
//        int m = 0;
        char *end;
        while ( ( read = getline(&line, &len, in)) != -1 )
        {
            if ( line[0] == '#' )
            {
                continue; // skip comment lines
            } else if ( line[0] == '@' && line[1] == 'S' && line[2] == 'M' && line[3] == 'P' && line[4] == 'S' ) // sample names line
            {
                line[strcspn(line, "\r\n")] = 0;
                char *p;
                int c = 0;
                p = strtok(line, " \t"); // chop off leading "@SMPS[ \t]"
                p = strtok(NULL, ", \t"); // split on following commas, spaces or tabs
                while ( p != NULL )
                {
                    if ( c > 5 ) error(EXIT_FAILURE, EIO, "Too many samples. A maximum of 5 is supported.\n");
                    args->smpns[c] = calloc(strlen(p)+1, sizeof(char));
                    strcpy(args->smpns[c], p);
                    p = strtok(NULL, ", \t"); // split on following commas, spaces or tabs
                    c++;
                }
                if ( c < 3 ) error( EXIT_FAILURE, EIO, "Not enough samples. A minimum of 3 is required.\n");
                args->nsmp = c;
                init_data(args); // we need to know the number of samples for args initialization
            } else {
                args->smp_is[j] = strtoul(line, &end, 10);
/*
                uint64_t curr = strtoul(line, &end, 10);
                int k;
                uint64_t s = 0, e = 1;
                for ( k = 1; k <= args->nsmp; k++ ) // iterate over subset sizes
                {

                    uint64_t p; // position in banker's sequence
                    s = e;
                    uint64_t subs = choose(args->nsmp, k, args); // get current subset size
                    e += subs; // get end position of current subset size
                    for ( p = s; p < e; p++ ) // iterate over all subsets of current size k
                    {
                        if ( args->smp_is[p] > 0 ) // Nothing to see here, otherwise
                        {
                            uint8_t m, n;
                            uint64_t c = 0;
                            n = 0; // current banker's subset inspected
                            c = args->bankers[p]; // get current sample int from banker's sequence
                            for ( n = 1; n < s; n++ ) // go through all subsets smaller than current cardinality
                            {
                                if ( (n|c) == c ) args->smp_is[n] += curr;
                            }
                            args->smp_is[j] += curr;
                        }
                    }
                } */
                j++;
                if ( j > args->nsmpp2-1 )
                {
                    fprintf(stderr, "\nError: Too many entries in input file.\n\n");
                    usage(EXIT_FAILURE);
                }
            }
        }
        if ( j < args->nsmpp2-1 )
        {
            fprintf(stderr, "\nError: Not enough entries in input file.\n\n");
            usage(EXIT_FAILURE);
        }
        free(line);
        fclose(in);

    print_VennDiagram_file(args);

    destroy_data(args);

    return 0;
}
