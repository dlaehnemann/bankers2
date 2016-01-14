
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
    char *prefix; /*! prefix of circos plotting files, iff they should be output */
    int nsmp; /*! number of samples, determined from number of single-sample lines at beginning of file */
    int nsmpp2; /*! 2^(nsmp) (is needed multiple times) */
    char **smpns; /*! list of samples as specified in banker's sequence file of intersection counts */
    uint64_t total; /*! total number of counts in all samples, needed to determine tick marks */
    int label_multipl; /*! label multiplier, any multiple of 3 for the decimal powers */
    char *label_suff; /*! suffix for tick mark label suffix, corresponding to label_multipl, i.e. k, M, G or T */
    double tick_spacing[3]; /*! values for the labelled tick marks and the two smaller distance tick mark categories */
    char *palette[2]; /*! palette prefix for plotting circos files, iff they should be output
                               First palette: Banding pattern encoding sample numbers.
                               Second palette: Differentiate sets within sample numbers. */
    uint32_t *bankers; /*! array to store banker's sequence for all possible sample subsets for
                                programmatic indexing into smp_is for output printing, e.g. for three
                                samples A, B and C this would be the following order:
                                [   C,   B,   A,  CB,  CA,  BA, CBA ]
                                [ 100, 010, 001, 110, 101, 011, 111 ]
                                */
    uint64_t *quick; /*! array to store lookup table of n choose k values from choose() */
    uint8_t missing; /*! flag, whether missing values should be extracted from banker's sequence file */
    uint64_t *missing_gts; /*! array to count missing genotypes of each sample */
    uint64_t *smp_is; /*! array to track all possible intersections between
                 samples, with each bit in the index integer belonging to one
                 sample. E.g. for three samples A, B and C, count would be in
                 the following order:
                 [   A,   B,  AB,   C,  AC,  BC, ABC ]
                 [ 001, 010, 011, 100, 101, 110, 111 ]
                 */
    FILE **couts; /*! array of circos files to output with following entries:
                       0: <prefix>.circos.conf
                       1: <prefix>.karyotype.txt
                       2: <prefix>.links.2smps.txt
                       3: <prefix>.links.3smps.txt
                       ...
                       nsmp: <prefix>.links.nsmps.txt */
}
args_t;


static void init_data(args_t *args)
{
    args->nsmpp2 = pow( 2, args->nsmp );
    args->bankers = (uint32_t*) calloc( args->nsmpp2, sizeof(uint32_t) );
    args->quick = (uint64_t*) calloc((args->nsmp * (args->nsmp + 1)) / 4, sizeof(unsigned long));
    if ( args->missing ) args->missing_gts = (uint64_t*) calloc( args->nsmp, sizeof(uint64_t));
    if ( args->nsmp > 9 )
    {
        fprintf(stderr, "Warning: For more than 9 samples, no good standard Color Brewer palette is \n"
                        "         available. You will have to adjust 'color's in the LINK FILE section of\n"
                         "         %s.circos.conf, in the %s.links.nsmps.tab files\n"
                        "         and band colors in %s.karyotype.txt manually (search\n"
                        "         and replace IDEOGRAM-COLORS-<smp_n> and SUBSET-RIBBON-COLORS-<smp_n>).\n"
                        "         However, your circos image will probably not be very informative with\n"
                        "         so many samples...\n"
                        "\n",
                        args->prefix,
                        args->prefix,
                        args->prefix);
        args->palette[0] = "IDEOGRAM-COLORS-";
        args->palette[1] = "SUBSET-RIBBON-COLORS-";
    } else {
        args->palette[0] = "ylgn-9-seq-";
        args->palette[1] = "set1-9-qual-";
    }
    args->couts = (FILE**) calloc( args->nsmp+3, sizeof(FILE*) );
    char fn[256];
    sprintf(fn, "%s.circos.conf", args->prefix);
    args->couts[0] = fopen(fn, "w");
    sprintf(fn, "%s.karyotype.txt", args->prefix);
    args->couts[1] = fopen(fn, "w");
    int i;
    for ( i = 2; i <= args->nsmp; i++ )
    {
        sprintf(fn, "%s.links.%ismps.tab", args->prefix, i);
        args->couts[i] = fopen(fn, "w");
    }
    args->smp_is = (uint64_t*) calloc( args->nsmpp2, sizeof(uint64_t));
}

static void destroy_data(args_t *args)
{
    if (args->missing)
    {
        free(args->missing_gts);
    }
    free(args->smp_is);
    free(args->bankers);
    free(args->quick);
    int i;
    for ( i = 0; i <= args->nsmp; i++ )
    {
        free(args->smpns[i]);
        fclose(args->couts[i]);
    }
    free(args->smpns);
    free(args->prefix);
    free(args->couts);
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

/*
 * Returns the Banker's number at the specified position, a.
 * Derived from the recursive bit flip method.
 * Added relative to Corin Lawson:
 * * Uses same lookup table solution as choose function, just
 *   maintained externally to persist across separate function calls.
 * * Uses bitwise symmetry of banker's sequence to use bitwise inversion
 *   instead of recursive bit flip for second half of sequence.
 */
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

// CODE BY CORIN LAWSON END

void determine_ticks(args_t *args)
{
    uint64_t ld = args->total / args->nsmp;
    int p10;
    for ( p10 = 0; ld >= 10; )
    {
        ld /= 10;
        p10++;
    }

    int multiplier = p10 / 3;
    args->label_multipl = multiplier;
    switch (multiplier)
    {
    case 0: args->label_suff = ""; break;
    case 1: args->label_suff = "k"; break;
    case 2: args->label_suff = "M"; break;
    case 3: args->label_suff = "G"; break;
    case 4: args->label_suff = "T"; break;
    case 5: args->label_suff = "P"; break;
    case 6: args->label_suff = "E"; break;
    case 7: args->label_suff = "Z"; break;
    default:    args->label_suff = "Unknown_Unit";
                break;
    }
    p10 %= 3;
    if ( ld <= 1 )                                          // working example with label_power = 0
    {
        args->tick_spacing[0] = pow(10, p10) / 2;           // 0.5
        args->tick_spacing[1] = args->tick_spacing[0] / 5;  // 0.1
        args->tick_spacing[2] = args->tick_spacing[0] / 10; // 0.05
    }
    else if ( ld > 1 && ld < 5 )
    {
        args->tick_spacing[0] = pow(10, p10);               // 1.0
        args->tick_spacing[1] = args->tick_spacing[0] / 2;  // 0.5
        args->tick_spacing[2] = args->tick_spacing[0] / 10; // 0.1
    }
    else if ( ld > 4 && ld < 10 )
    {
        args->tick_spacing[0] = pow(10, p10) * 2;           // 2.0
        args->tick_spacing[1] = args->tick_spacing[0] / 2;  // 1.0
        args->tick_spacing[2] = args->tick_spacing[0] / 4;  // 0.5
    }
}

void print_circos_conf(args_t *args)
{
    fprintf(args->couts[0],
            "# INCLUDE PREDEFINED COLORS, FONTS AND PATTERNS\n"
            "# most importantly, the Brewer color palettes\n"
            "<<include colors_fonts_patterns.conf>>\n"
            "\n"
            "\n"
            "# IDEOGRAM CONFIGURATION\n"
            "\n"
            "<ideogram>\n"
            "\n"
            "<spacing>\n"
            "default = 0.01r\n"
            "break   = 0.5r\n"
            "</spacing>\n"
            "\n"
            "# Ideogram positioning configuration\n"
            "radius           = 0.90r\n"
            "thickness        = 100p\n"
            "fill             = no\n"
            "# fill_color       = black\n"
            "stroke_thickness = 1\n"
            "stroke_color     = lgrey\n"
            "\n"
            "# Ideogram label configuration\n"
            "show_label       = yes\n"
            "label_font       = condensed\n"
            "label_radius     = dims(ideogram,radius) + 0.075r\n"
            "label_with_tag   = yes\n"
            "label_size       = 36\n"
            "label_parallel   = yes\n"
            "#label_case       = lower\n"
            "label_format     = eval(sprintf(\"%%s\",var(label)))\n"
            "\n"
            "# Ideogram bands configuration\n"
            "show_bands            = yes\n"
            "fill_bands            = yes\n"
            "band_stroke_thickness = 0\n"
            "#band_stroke_color     = lgrey\n"
            "band_transparency     = 0\n"
            "\n"
            "</ideogram>\n"
            "\n"
            "\n"
            "# TICKS CONFIGURATION\n"
            "\n"
            "show_ticks          = yes\n"
            "show_tick_labels    = yes\n"
            "\n"
            "<ticks>\n"
            "\n"
            "radius           = dims(ideogram,radius_outer)\n"
            "orientation      = out\n"
            "label_multiplier = %i%s\n"
            "color            = black\n"
            "force_display    = yes\n"
            "show_label       = yes\n"
            "size             = 20p\n"
            "thickness        = 3p\n"
            "label_offset     = 5p\n"
            "#skip_last_label  = yes # uncomment this line, if the end label of a sample overlaps the last regular label\n"
            "\n"
            "<tick>\n"
            "spacing        = %fu\n"
            "show_label     = no\n"
            "size           = 0.05r\n"
            "</tick>\n"
            "\n"
            "<tick>\n"
            "spacing        = %fu\n"
            "show_label     = no\n"
            "size           = 0.1r\n"
            "</tick>\n"
            "\n"
            "<tick>\n"
            "spacing        = %fu\n"
            "size           = 0.2r\n"
            "show_label     = yes\n"
            "label_size     = 24p\n"
            "format         = %%d %s\n"
            "</tick>\n"
            "\n"
            "<tick>\n"
            "show           = yes\n"
            "position       = end\n"
            "size           = 0.2r\n"
            "label_size     = 24p\n"
            "format         = %%d %s\n"
            "</tick>\n"
            "\n"
            "</ticks>\n",
            args->label_multipl,
            args->label_multipl == 0 ? "" : "e-3",
            args->tick_spacing[2],
            args->tick_spacing[1],
            args->tick_spacing[0],
            args->label_suff,
            args->label_suff);
    fprintf(args->couts[0],
            "\n"
            "\n"
            "<image>\n"
            "<<include etc/image.conf>>\n"
            "</image>\n"
            "\n"
            "\n"
            "# INCLUDE KARYOTYPE DEFINITION DERIVED FROM SAMPLE AND INTERSECTION COUNTS\n"
            "\n"
            "karyotype   = %s.karyotype.txt\n",
            args->prefix);
    fprintf(args->couts[0],
            "\n"
            "chromosomes_units = %.0f\n",
            pow(10, args->label_multipl * 3) );
    fprintf(args->couts[0],
            "## If you want to restrict plotting to certain samples, uncomment the\n"
            "## two following lines and adjust the sample list to contain those wanted.\n"
            "#chromosomes       = %s",
            args->smpns[0] );
    int i;
    for ( i = 1; i < args->nsmp; i++) fprintf(args->couts[0], ";%s", args->smpns[i] );
    fprintf(args->couts[0],
            "\n"
            "#chromosomes_display_default = no\n"
            "\n"
            "# If you adjust the radius of the ideograms, links incident\n"
            "# on these ideograms will inherit the new radius.\n"
            "#chromosomes_radius = hs2:0.9r;hs3:0.8r\n"
            "\n"
            "\n"
            "# LINK FILE INCLUDES (SAMPLE INTERSECTIONS) AND LINK CONFIGURATIONS\n"
            "\n"
            "# Links (bezier curves or straight lines) are defined in <links> blocks.\n"
            "# Each link data set is defined within a <link> block. \n"
            "#\n"
            "# As with highlights, parameters defined\n"
            "# in the root of <links> affect all data sets and are considered\n"
            "# global settings. Individual parameters value can be refined by\n"
            "# values defined within <link> blocks, or additionally on each\n"
            "# data line within the input file.\n"
            "\n"
            "<links>\n"
            "\n"
            "radius = 0.99r\n"
            "crest  = 1\n"
            "ribbon           = yes\n"
            "flat             = yes\n"
            "stroke_color     = dgrey\n"
            "stroke_thickness = 1\n"
            "color            = grey_a2 # this is a fallback value, if it is used,\n"
            "                           # something went wrong in the color\n"
            "                           # assignment of the circos file generation\n"
            "bezier_radius        = 0r\n"
            "bezier_radius_purity = 0.5\n"
            "\n");
    for ( i = 2; i <= args->nsmp; i++)
    {
        fprintf(args->couts[0],
            "<link>\n"
            "file          = %s.links.%dsmps.tab\n"
            "# radius        = 0.95r\n"
            "# color         = %s%i_a2\n"
            "z             = %i\n"
            "bezier_radius_purity = %f\n"
            "\n"
            "# Curves look best when this value is small (e.g. 0.1r or 0r)\n"
            "# bezier_radius = 0.1r\n"
            "# thickness     = 1\n"
            "\n"
            "# These parameters have default values. To unset them\n"
            "# use 'undef'\n"
            "#crest                = undef\n"
            "#bezier_radius_purity = undef\n"
            "\n"
            "# Limit how many links to read from file and draw\n"
            "#record_limit  = 2000\n"
            "\n"
            "</link>\n"
            "\n",
            args->prefix,
            i,
            args->palette[0],
            i,
            ( i == args->nsmp ) ? 1 : i,
            (1.0/args->nsmp)*(i-1)
            );
    }
    fprintf(args->couts[0],
            "</links>\n"
            "\n"
            "\n"
            "# STANDARD INCLUDE FROM CIRCOS' etc DIR\n"
            "<<include etc/housekeeping.conf>>\n"
            "data_out_of_range* = trim\n"
            "\n"
            "# If you want to turn off all track default values, \n"
            "# uncomment the line below. This overrides the\n"
            "# value of the parameter imported from etc/housekeeping.conf\n"
            "\n"
            "#track_defaults* = undef\n"
            "# The defaults for links are\n"
            "#\n"
            "ribbon           = yes\n"
            "# color            = black\n"
            "# thickness        = 1\n"
            "radius           = 0.95r\n"
            "bezier_radius    = 0.1r\n"
            "# crest                = 0.5\n"
            "# bezier_radius_purity = 0.75\n"
            "#\n"
            "# See etc/tracks/link.conf in Circos distribution\n");
}

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

/*!
 * Print karyotype file and link files for circos plotting.
 *
 * The karyotype file contains one "chromosome" per sample, with one band for
 * each subset size k (from 1, which are counts unique to that sample, to nsmp,
 * which are the counts that all samples share). Subset sizes are colored green,
 * with darker green signifying a larger subset size. An extra light red band
 * is created for each sample's missing value count, if requested via -m.
 * For each subset size of at least 2 (i.e. where counts are shared with at
 * least one other sample), one link file is created. For each count, ribbons
 * of the same color connect all samples that share it. To make different
 * counts as distinguishable as possible, Ribbon colors are as diverse as
 * possible and within each subset size, the outer starting positions of the
 * counts differ, covering cases where similar or same color ribbons and next
 * to each other. To decrease ribbon overlap, smaller subset sizes have their
 * ribbons pushed to the edge, while larger subset sizes have their ribbons
 * pushed to the center.
 */
void print_circos_karyotype_links(args_t *args)
{
    /* counts of how many shared genotypes a sample (first index) has with how
       many other samples (second index) */
    uint64_t is_smp_cnt[args->nsmp][args->nsmp];
    memset(is_smp_cnt, 0, sizeof(is_smp_cnt) );

    // LINK FILE PRINTING AND KARYOTYPE DATA COLLECTION
    uint8_t k;
    uint64_t s = 0, e = 1; /* start and end of current subset size in banker's sequence */
    uint64_t smp_s[args->nsmp]; /* current start position tracking for all samples
                                   during the link file printing */
    memset(smp_s, 0, sizeof(smp_s) );
    for ( k = 1; k <= args->nsmp; k++ ) // iterate over subset sizes
    {
        uint64_t p; /* position in banker's sequence */
        s = e;
        uint64_t subs = choose(args->nsmp, k, args); // get current subset size
        uint64_t subs_i = 0; // subset index within current subset size
        e += subs; // get end position of current subset size
        for ( p = s; p < e; p++, subs_i++ ) // iterate over all subsets of current size k
        {
            if ( args->smp_is[p] > 0 ) // Nothing to see here, otherwise
            {
                uint8_t m, n;
                uint64_t c = 0;
                n = 0; /* current sample inspected */
                m = 0; /* samples seen */
                uint8_t smps[k]; /* keep track of encountered samples for link file printing */
                c = args->bankers[p]; // get current sample int from banker's sequence
                while ( m < k && n < args->nsmp )
                {
                    if ( c & 1 ) // collect stats from each sample's perspective
                    {
                        is_smp_cnt[args->nsmp-1 - n][k-1] += args->smp_is[p];
                        smps[m] = args->nsmp-1 - n;
                        m++;
                    }
                    c >>= 1;
                    n++;
                }

                // link file printing
                if ( k >= 2 )
                {
                    uint64_t i1, i2;
                    i1 = 0;
                    while ( i1 < k-1 )
                    {
                        i2 = i1 + 1;
                        while ( i2 < k )
                        {
                            fprintf(args->couts[k],
                                    "%s_id\t%"PRIu64"\t%"PRIu64"\t%s_id\t%"PRIu64"\t%"PRIu64"\tcrest=%f,color=%s%"PRIu64"_a%i,radius=%fr\n",
                                    args->smpns[ smps[i1] ],
                                    smp_s[ smps[i1] ],
                                    smp_s[ smps[i1] ] + args->smp_is[p],
                                    args->smpns[ smps[i2] ],
                                    smp_s[ smps[i2] ],
                                    smp_s[ smps[i2] ] + args->smp_is[p],
                                    (1.2/args->nsmp)*(k-1) + ((3.0/args->nsmp)/subs)*subs_i,
                                    ( k == args->nsmp ) ? args->palette[0] : args->palette[1],
                                    ( k == args->nsmp ) ? 9 : (p % 7) + 1,
                                    ( k == args->nsmp ) ? 5 : 2,
                                    0.99 - 0.01*subs_i
                                    );
                            i2++;
                        }
                        i1++;
                    }
                }
                // update samples' start positions for link file printing
                uint8_t k1;
                for ( k1 = 0; k1 < k; k1++ ) smp_s[ smps[k1] ] += args->smp_is[p];
            }
        }
    }

    // KARYOTYPE PRINTING (SAMPLES AND BANDS)
    uint8_t n;
    for ( n = 0; n < args->nsmp; n++ )
    {
        uint64_t sum = 0; /* running sharing sum for current sample */
        const char *sn; /* name of the current sample */
        sn = args->smpns[n]; /* name of current sample */
        for ( k = 1; k <= args->nsmp; k++) // print sharing bands for current sample
        {
            if ( is_smp_cnt[n][k-1] > 0 ) // nothing to see here, otherwise
            {
                uint64_t e = sum + is_smp_cnt[n][k-1];
                fprintf(args->couts[1],
                        "band\t%s_id\t%s_%"PRIu8"\t%"PRIu8"_smps\t%"PRIu64"\t%"PRIu64"\t%s%"PRIu8"\n",
                        sn,
                        sn,
                        k,
                        k,
                        sum,
                        e,
                        args->palette[0],
                        9 - (args->nsmp - k)
                        );
                sum = e;
            }
        }
        if ( args->missing )
        {
            uint64_t e = sum + args->missing_gts[n];
            fprintf(args->couts[1], // print missing genotypes band of sample
                            "band\t%s_id\t%s_MV\tmissing_values\t%"PRIu64"\t%"PRIu64"\t%s\n",
                            sn,
                            sn,
                            sum,
                            e,
                            "rdylgn-11-div-1_a4"
                            );
            sum = e;
        }
        fprintf(args->couts[1], // print current sample
                "chr\t-\t%s_id\t%s\t0\t%"PRIu64"\twhite\n",
                sn,
                sn,
                sum
                );
    }
}


/*! array of circos files to output with following entries:
                       0: <prefix>.circos.conf
                       1: <prefix>.karyotype.txt
                       2: <prefix>.links.2smps.txt
                       3: <prefix>.links.3smps.txt
                       ...
                       nsmp-1+1: <prefix>.links.nsmps.txt */
/*
 * Function for printing out circos plot files.
 */
void print_circos_files(args_t *args)
{
    determine_ticks(args);
    print_circos_conf(args);
    print_circos_karyotype_links(args);
}


static int usage(int status)
{
    fprintf(stderr, "\n"
                    "About:   Create files for a circos plot representing all subsets/intersections\n"
                    "         contained in a banker's sequence order file.\n"
                    "\n"
                    "Usage:   bankers2circos [options] <banker's-seq-subsets-file>\n"
                    "\n"
                    "Options:\n"
                    "    -p, --prefix <prefix>    prefix for circos file names [default: prefix of input without path]\n"
                    "    -m, --missing            if set, include count of missing genotypes per sample in\n"
                    "                             output, including the circos files, if they are requested\n"
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
    args->missing = 0;
    args->total = 0;
    args->smpns = calloc(32, sizeof(char*));

    static struct option loptions[] =
    {
        {"help",        no_argument,      0,'h'},
        {"prefix",      required_argument,0,'p'},
        {"missing",     no_argument,      0,'m'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "p:mh",loptions,NULL)) >= 0) {
        switch (c) {
            case 'p': args->prefix = optarg; break;
            case 'm': args->missing = 1; break;
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
        uint64_t j = 1;
        int m = 0;
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
                    if ( c >= 32 ) error(EXIT_FAILURE, EIO, "Too many samples. A maximum of 32 is supported.\n");
                    args->smpns[c] = calloc(strlen(p)+1, sizeof(char));
                    strcpy(args->smpns[c], p);
                    p = strtok(NULL, ", \t"); // split on following commas, spaces or tabs
                    c++;
                }
                args->nsmp = c;
                init_data(args); // we need to know the number of samples for args initialization
            } else {
                if ( args->missing && m < args->nsmp )
                {
                    args->missing_gts[m] = strtoul(line, &end, 10);
                    args->total += args->missing_gts[m];
                    m++;
                    continue;
                }
                args->smp_is[j] = strtoul(line, &end, 10);
                args->total += args->smp_is[j];
                j++;
                if ( j > args->nsmpp2 )
                {
                    fprintf(stderr, "\nError: Too many entries in input file. Did you forget to set the missing values flag (-m)?\n\n");
                    usage(EXIT_FAILURE);
                }
            }
        }
        if ( j < args->nsmpp2 )
        {
            fprintf(stderr, "\nError: Not enough entries in input file. Did you accidentally set the missing values flag (-m)?\n\n");
            usage(EXIT_FAILURE);
        }
        free(line);
        fclose(in);

    uint32_t i;
    /* Compute banker's sequence for following circos file printing.
     */
    for ( i = 0; i < args->nsmpp2; i++ )
    {
        args->bankers[i] = compute_bankers(i, args);
    }

    print_circos_files(args);

    destroy_data(args);

    return 0;
}
