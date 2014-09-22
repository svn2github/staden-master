/*
 * Other ideas:
 *
 * 1. Identify all readings on inconsistent templates and strike them out of
 *    the contig. If holes remain, then verify assembly at that point.
 *
 * 2. Look for regions of high depth. If high depth and avg length of insert
 *    spanning the gap is shorter than normal, then probable collapsed repeat.
 *
 * 3. Use the diploid scanner and haplotype splitting code too as part of
 *    the auto-break mechanism.
 *
 * 4. Spot repeats and divergence => end of repeats. Check all divergence
 *    points to ensure consistency. (Need a decent local alignment algorithm)
 *
 * 5. Emulate neighbourhoods by using a word that is (eg) 6mer, 1mer gap, 6mer.
 *    Hence with all overlapping coordinates this allows 1 base to differ
 *    while still using the memory for a 12mer and not 13mer.
 *
 * 7. Ressurrect stop-detector. Clip after stop - it's probably wrong data.
 */

/*

  Given word size N, the number of combinations of R of N being GC
  is 2^N * N! / ((N-R)! R!)
  
  Ie for N=8:
  
  0/8 256
  1/8 2048
  2/8 7168
  3/8 14336
  4/8 17920
  5/8 14336
  6/8 7168
  7/8 2048
  8/8 256

  Given a GC content of x (0 <= x <= 1) we therefore expect
  (x/2)^R * ((1-x)/2)^(N-R) * 2^N * N! / ((N-R)!R!) matches.
*/

#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <io_lib/hash_table.h>

#include "tg_gio.h"
#include "misc.h"
#include "xalloc.h"
#include "text_output.h"
#include "filter_words.h"
#include "auto_break.h"
#include "gap_globals.h"
#include "qual.h"
#include "qualIO.h"
#include "dstring.h"
#include "consensus.h"

#define MIN3(a,b,c) (MIN(MIN((a),(b)),(c)))
#define MAX3(a,b,c) (MAX(MAX((a),(b)),(c)))

#ifndef WS
#    define WS 12
#endif
#define WS2 ((int)(WS/2))

#define ALLB(ws) ((1<<(2*(ws)))-1)

#define MIN_OVERLAP 5
#define MIN_CLIP 3
/* #define NORMALISE_FOR_GC 1 */
#define NORMALISE_STR_SCORES

#define CONTIG_END_IGNORE 200

/* Size to consider neighbouring problems over */
#define MAX_NEIGHBOUR 10000

/* #define GAPPED_WORDS */

#ifdef GAPPED_WORDS
#  define WORD2GAPPED(w) ( ((w) & ALLB(WS2)) | (((w) >> 2) & ~ALLB(WS2)) )
#else
#  define WORD2GAPPED(w) (w)
#endif

typedef unsigned short count_t;
#define MAX_COUNTS ((1<<16)-1)

//#define DEBUG
//#define DEBUG_SEQ

typedef struct contig_region {
    int start;
    int end;
    int deleted;
    int valid;

    // Array of records covering gap
    rangec_t *r;
    int nr;

    // scores, for adding to the tag
    int score;
    int num_unique_good;
    int num_good;
    int num_unique_bad;
    int num_bad;
    int num_unique_large;
    int num_large;
    int num_spanning;
    int num_unique_spanning;
    int num_single;
    int num_unique_single;
} contig_region_t;


typedef struct clip_pos {
    int left;
    int right;
} clip_pos_t;

void dump_gaps(Array gaps) {
    int i;
    puts("\n");
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	printf("Gap %d\t%d %d %d\n", i, 
	       gap->start, gap->end, gap->deleted);
    }
}

static count_t counts[1<<(2*WS)]; /* 4^WS */
static int lookup[256];
static int clookup[256];
static double probs[4];

/*
 * Turns an array of coverage data (0 => no coverage) to an Array of
 * contig_region_t structs.
 */
static Array coverage2contig_regions(char *valid, int *clips, int *depth,
				     int cstart, int clen, int end_skip_len) {
    int i, j, ngaps = 0;
    Array gaps;

    gaps = ArrayCreate(sizeof(contig_region_t), 0);
    if (NULL == gaps)
	return NULL;

    for (i = 0; i < clen; i++) {
	if (valid[i])
	    break;
    }
    if (i != 0) {
	contig_region_t *r;
	printf("START gap %d..%d\n", cstart, i+cstart);
	/* Skip start gap */
	r = (contig_region_t *)ArrayRef(gaps, ngaps);
	r->start = 1;
	r->end = i;
	r->deleted = 0;
	r->valid = 1;
	ngaps++;
    }

    for (; i < clen; i++) {
	if (clips[i] >= 3 && 100*clips[i]/depth[i] > 33) {
	    contig_region_t *r;

	    printf("SOFT_CLIP %d\n", i+cstart);

	    r = (contig_region_t *)ArrayRef(gaps, ngaps);
	    r->start = i+cstart;
	    r->end =   i+cstart;
	    r->deleted = 0;
	    r->valid = i < end_skip_len || clen-i < end_skip_len
		? 1 : 0;
	    ngaps++;
	}

	if (!valid[i]) {
	    contig_region_t *r;
	    for (j = i+1; j < clen; j++) {
		if (valid[j])
		    break;
	    }
	    if (j >= clen) {
		/* Skip end gaps */
		printf("END gap %d..%d\n", i+cstart, j-1+cstart);
	    } else {
		printf("INTERNAL gap %d..%d\n", i+cstart, j-1+cstart);
	    }
	    r = (contig_region_t *)ArrayRef(gaps, ngaps);
	    r->start = i+cstart;
	    r->end = j >= clen ? clen : j-1+cstart;
	    r->deleted = 0;
	    r->valid = i < end_skip_len || clen-i < end_skip_len
		? 1 : (j >= clen);
	    ngaps++;

	    /* Here simply to generate SOFT_CLIP messages */
	    for (i++; i <= j-1 && i < clen; i++) {
		if (clips[i] >= 3 && 100*clips[i]/depth[i] > 33) {
		    printf("SOFT_CLIP %d\n", i+cstart);
		    
		    r = (contig_region_t *)ArrayRef(gaps, ngaps);
		    r->start = i+cstart;
		    r->end =   i+cstart;
		    r->deleted = 0;
		    r->valid = i < end_skip_len || clen-i < end_skip_len
			? 1 : 0;
		    ngaps++;
		}
	    }

	    i = j-1;
	}
    }
    
    return gaps;
}

/*
 * Counts all WS-mers in the gap4 database to detect overly frequent words.
 * Only uses clipped data (to avoid vector etc).
 */
static void init_tables(void) {
    int i;

    /* Initialise lookup tables */
    for (i = 0; i < 256; i++) {
	lookup[i] = -1;
	clookup[i] = -1;
    }
    lookup['A']  = lookup['a'] = 0;
    lookup['C']  = lookup['c'] = 1;
    lookup['G']  = lookup['g'] = 2;
    lookup['T']  = lookup['t'] = 3;
#ifdef GAPPED_WORDS
    clookup['A'] = clookup['a'] = 3 << (2*WS);
    clookup['C'] = clookup['c'] = 2 << (2*WS);
    clookup['G'] = clookup['g'] = 1 << (2*WS);
    clookup['T'] = clookup['t'] = 0 << (2*WS);
#else
    clookup['A'] = clookup['a'] = 3 << (2*WS-2);
    clookup['C'] = clookup['c'] = 2 << (2*WS-2);
    clookup['G'] = clookup['g'] = 1 << (2*WS-2);
    clookup['T'] = clookup['t'] = 0 << (2*WS-2);
#endif

    memset(counts, 0, (1<<(2*WS)) * sizeof(*counts));
}

#ifdef NORMALISE_FOR_GC
static void init_gc_table(double gc) {
    probs[1] = probs[2] = gc/2;
    probs[0] = probs[3] = (1-gc)/2;
}
#endif


static char *word2str(int word) {
    static char str[WS+2];
    signed int i, j;
    
    for (j = WS, i = 0; i < WS; i++) {
#ifdef GAPPED_WORDS
	if (i == (int)(WS/2))
	    str[j--] = '.';
#endif
	str[j--] = "ACGT"[word & 3];
	word >>= 2;
    }
    str[WS+1] = 0;
    
#ifdef GAPPED_WORDS
    return str;
#else
    return str+1;
#endif
}

#if 0
static char *word2str2(int word) {
    static char str[WS+2];
    int i, j;

    for (j = WS, i = 0; i < WS+1; i++) {
	str[j--] = "ACGT"[word & 3];
	word >>= 2;
    }
    str[WS+1] = 0;

    return str;
}
#endif

double compute_prob(int word) {
    int i;
    double prob = 1;
    
    for (i = 0; i < WS; i++) {
	int base = (word >> (2*(WS-1)-2*i)) & 3;
	prob *= probs[base];
    }

    return prob;
}


/*
 * Computes the maximum theoretical redundancy of this word. Ie if we have
 * a GT repeat then word GTGTGTGTGTGT occurs 6 times more frequently within
 * a run than we'd expect.
 *
 * It does this by sliding word w along looking at the overlap. Ie:
 * ABCDEABC
 *  ABCDEABC          no match in overlap as BCDEABC != ABCDEAB
 *
 * ABCDEABC
 *      ABCDEABC      matches as ABC == ABC
 *
 * Note that on average even a random word will have 1/4 chance of the last
 * base matching the first base so claim to have a redundancy of WS/(WS-1).
 * Similarly 1/16 times it'll have redundancy of WS/(WS-2). We could
 * compensate for this by factoring this expectation into the equation, maybe
 * even utilising the GC content in there too, but for now I think it skews
 * things sufficiently little to not be too concerned.
 *
 * Returns the total number of new counts.
 */
#ifdef NORMALISE_STR_SCORES
int64_t normalise_str_scores(void) {
    int w;
    int64_t tc = 0;
    
    for (w = 0; w < 1<<(2*WS); w++) {
	int i, m = (1 << (2*(WS-1))) - 1;

	if (!counts[w])
	    continue;

	for (i = 1; i <= WS; i++, m >>= 2) {
	    if ( (w >> (2*i)) == (w & m) )
		break;
	}
	counts[w] /= (double)WS/i;
	if (counts[w] == 0)
	    counts[w] = 1; /* Min count of 1 */

	tc += counts[w];
    }

    return tc;
}
#endif

/*
 * Fills out the counts[] array based on all words of length WS in the
 * consensus sequences. 
 *
 * Returns the total number of words indexed.
 */
int64_t word_count_cons(GapIO *io, int argc, contig_list_t *argv) {
    int j, cnum;
    int64_t tw = 0, gc = 0, at = 0;
    char *cons = NULL;

    init_tables();

    for (cnum = 0; cnum < argc; cnum++) {
	contig_t *c;
	int clen;
	unsigned char *s;
	unsigned int word, cword;

	c = cache_search(io, GT_Contig, argv[cnum].contig);
	clen = c->end - c->start + 1;
	cons = xrealloc(cons, clen);

	calc_consensus(c->rec, c->start, c->end, CON_SUM, cons,
		       NULL, NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);

	/* printf("CONS =%d/%d: %.*s\n", cnum, NumContigs(io), clen, cons); */

	if (clen <= 2*CONTIG_END_IGNORE)
	    continue;

	s = (unsigned char *) cons + CONTIG_END_IGNORE;
	cons[clen-1-CONTIG_END_IGNORE] = 0;

	cword = word = 0;
	for (j = 0; *s; s++) {
	    if (*s == '*')
		continue;

	    switch (lookup[*s]) {
	    case -1:
		j = 0;
		break;
	    case 0:
	    case 3:
		at++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    case 1:
	    case 2:
		gc++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    }
	    
	    /*
	    printf("%c %s ", j >= WS ? *s : '.', word2str(word));
	    printf("%s %08x %08x\n", word2str(cword), cword, clookup[*s]);
	    */

	    if (j >= WS) {
		/* Remove middle nucleotide from word before storing */
		unsigned int w2, cw2;

		w2  = WORD2GAPPED(word);
		cw2 = WORD2GAPPED(cword);
		
		/* Count both original and complementary word */
		if (counts[ w2 & ALLB(WS)] < MAX_COUNTS)
		    counts[ w2 & ALLB(WS)]++;
		if (counts[cw2 & ALLB(WS)] < MAX_COUNTS)
		    counts[cw2 & ALLB(WS)]++;

		/*
		printf("%c %s (%d)",
		       *s, word2str(w2), counts[w2 & ALLB(WS)]);
		printf("\t%s (%d)\n",
		       word2str(cw2), counts[cw2 & ALLB(WS)]);
		*/
		tw += 2;
	    }
	}
    }

    xfree(cons);

    printf("Total words = %"PRId64", GC = %5.2f%%\n", tw, (100.0*gc)/(gc+at));

#ifdef NORMALISE_FOR_GC
    /* normalise counts by GC content */
    init_gc_table((double)gc/(gc+at));

    for (j = 0; j < (1<<(2*WS)); j++) {
	int c = counts[j];
	c /= (1<<(2*WS)) * compute_prob(j);
	if (c > MAX_COUNTS)
	    c = MAX_COUNTS;
	counts[j] = c;
    }
#endif

#ifdef NORMALISE_STR_SCORES
    tw = normalise_str_scores();
#endif

    return tw;
}

void print_counts(double min) {
    int i;
    for (i = 0; i < 1<<(2*WS); i++) {
	if (counts[i] >= min)
	    printf("%s %d\n", word2str(i), counts[i]);
    }
}

void print_bins(void) {
    int i, j;
    int bins[10000];
    
    memset(bins, 0, sizeof(*bins)*10000);
    for (i = 0; i < (1<<(2*WS)); i++) {
	if (counts[i] < 10000)
	    bins[counts[i]]++;
    }

    for (i = 0; i < 10000; i++)
	if (bins[i])
	    break;
    for (j = 9999; j >= 0; j--)
	if (bins[j])
	    break;
    for (; i <= j; i++) {
	printf("%d %d\n", i, bins[i]);
    }
}

int filter_common_words(char *seq, char *filt, size_t len, int tw,
			double depth, double score, char filter_char,
			int debug) {
    size_t i, j;
    unsigned int word = 0;
    int pads = 0;
    double rescale = 1;

    memcpy(filt, seq, len);

    /* Start with an entire word */
    for (i = j = 0; i < WS && i < len; i++) {
	if (seq[i] == '*') {
	    pads++;
	    continue;
	}

	word = (word << 2) | lookup[(unsigned char) seq[i]];
	j++;
    }
    
    /*
    printf("compare %f * %f (=%f) vs %f * %d / %d (=%f)\n",
	   score, depth, score * depth,
	   score, tw, 1<<(2*WS),
	   score * (double)tw / (1<<(2*WS)));
    */

    /*
     * A quick and easy hack, but not particularly robust:
     * If we've observed more words than there are buckets then we
     * expect the number of hits to be higher than the depth alone would
     * suggest.
     * Therefore we downscale a bit to compensate.
     */
    if (tw /  (1<<(2*WS)) > 1) {
	rescale = ((double)tw / (1<<(2*WS))) / depth;
    }

    /* Scan through looking for matches of a word */
    for (; i < len; i++) {
	if (seq[i] == '*'){
	    pads++;
	    continue;
	}

	word = (word << 2) | lookup[(unsigned char) seq[i]];
	
	if (debug)
	    printf("Seq pos %ld %.*s: => %d",
		   (long)i, WS, seq+i, counts[word & ((1<<(2*WS))-1)]);
	/* if (counts[word & ((1<<(2*WS))-1)] >= score * (double)tw / (1<<(2*WS))) { */
	if ((counts[WORD2GAPPED(word) & ALLB(WS)]) / rescale
	    >= score * depth) {
	    /* FIXME: ignores pads for now */
#ifdef GAPPED_WORDS
	    memset(&filt[i-WS], filter_char, WS+1);
#else
	    memset(&filt[i-(WS-1)], filter_char, WS);
#endif
	    if (debug) putchar('*');
	}

	if (debug) putchar('\n');
    }

    /* Merge filtered blocks together if they have <= 4 bases gap */
    /* FIXME: merge if lots of valid read-pairs span gaps.
     * Eg this means that ....gap...gap.... becomes one gap if we find
     * linkage info claiming that the sequence gap...gap looks as if it
     * forms a valid contig itself.
     */
    for (i = 1; i < len; i++) {
	if (filt[i-1] == filter_char && filt[i] != filter_char) {
	    int g = i;
	    while (i < len && filt[i] != filter_char)
		i++;
	    if (i-g <= 4) {
		while (g != i && g < len)
		    filt[g++] = filter_char;
	    }
	}
    }

    return 0;
}

/*
 * Identifies regions where the input data mismatches the consensus and
 * mask out N characters either side of this location.
 *
 * If using automatic repeat filtering, you probably want to supply
 * this parameter as WS-1.
 *
 * The logic is that these are base calling errors and as such will bias
 * the word-counting (as errors in repeats will make it appear to be a
 * unique word).
 */
void filter_consen_diffs(char *in, char *out, int len, char *cons, int N) {
    int i, j;
    //printf("%.*s\n", MIN(79,len), in);
    //printf("%.*s\n", MIN(79,len), cons);
    for (i = 0; i < len; i++) {
	if (in[i] == cons[i])
	    continue;
	if (toupper(in[i]) == cons[i])
	    continue;
	if (in[i] == '-' && cons[i] == 'N')
	    continue;
	for (j = i-N >= 0 ? i-N : 0;
	     j <= i+N && j < len;
	     j++) {
	    out[j] = '%';
	}
    }
}

/*
 * Loops through all readings in a contig looking to see that specific low
 * complexity regions have readings that span either end. If not then we
 * label this region as a suspicious join.
 *
 * Returns: an Array of type contig_region_t on success
 *          NULL on failure
 */
static Array suspect_joins(GapIO *io, tg_rec contig, int64_t tw,
			   double filter_score, int filter_consensus,
			   double avg_depth, int min_mq,
			   HashTable *clip_hash, int end_skip_len) {
    int i, clen, cstart, cend;
    char *valid = NULL, *cp;
    int *clip_depth = NULL, *total_depth = NULL;
    char legal_chars[256];
    Array gaps;
    char *cons = NULL;
    contig_iterator *ci;
    rangec_t *r;
    int skip_non_paired = 1;

    consensus_valid_range(io, contig, &cstart, &cend);
    clen = cend - cstart + 1;

    if (NULL == (valid = (char *)xcalloc(clen+1, 1)))
	goto error;

    if (NULL == (clip_depth = (int *)xcalloc(clen+1, sizeof(int))))
	goto error;
    if (NULL == (total_depth = (int *)xcalloc(clen+1, sizeof(int))))
	goto error;

    memset(legal_chars, 0, 256);
    for (cp = "ACGTacgt"; *cp; cp++)
	legal_chars[(unsigned char) *cp] = 1;

    /* Compute consensus */
    if (NULL == (cons = (char *)xmalloc(clen+1)))
	goto error;
    calc_consensus(contig, cstart, cend, CON_SUM, cons, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    
    /*
     * Loop through readings computing new clip points after masking.
     * Mark the remainder as 'valid' to identify the invalid bits.
     */
    ci = contig_iter_new(io, contig, 1, CITER_FIRST |CITER_ISTART |CITER_PAIR,
			 CITER_CSTART, CITER_CEND);
    while (NULL != (r = contig_iter_next(io, ci))) {
	seq_t *s = cache_search(io, GT_Seq, r->rec), *sorig = s;
	char *seq, *fseq;
	size_t len;
	int unique, first, last, firstq, lastq, p;
	int p_start, p_end;

	// if (skip_non_paired) {
	//     if (!r->library_rec || r->pair_rec == 0)
	// 	continue;
	// }

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	seq = s->seq;
	len = ABS(s->len);
	if (NULL == (fseq = malloc(len))) {
	    if (sorig != s)
		free(s);
	    continue;
	}
	memcpy(fseq, seq, len);

	/* Aggregate cutoff consensus */

	//if (s->left > 1) printf("L #%"PRIrec" %.*s\n", r->rec, s->left-1, seq);
	//if (s->right < len) printf("R #%"PRIrec" %.*s\n", r->rec, len - s->right, &seq[s->right]);

	if (s->left > MIN_CLIP && (p = r->start + s->left - cstart - 2) >= 0)
	    clip_depth[p]++;

	if (len - s->right-2 >= MIN_CLIP-1 &&
	    (p = r->start + s->right - cstart - 2) <= clen) {
	    if (p < 0)    p = 0;
	    if (p > clen) p = clen;
	    clip_depth[p]++;
	}

	p_start = r->start + s->left - cstart - 2;
	if (p_start < 0) p_start = 0;
	p_end = r->start + s->right - cstart - 1;
	if (p_end > clen) p_end = clen;

	for (p = p_start; p < p_end; p++) {
	    total_depth[p]++;
	}

	/* Filter over represented words in consensus */
	if (filter_consensus)
	    filter_common_words(seq, fseq, len, tw, avg_depth,
				filter_score, '~', 0);

	/* Filter low mapping quality */
	if (r->mqual < min_mq) {
	    memset(fseq, 'q', len);
	} else {
	    
	    /* Filter where the seq is low qual and disagrees with consen */
	    filter_consen_diffs(&seq[s->left-1], &fseq[s->left-1],
				s->right - s->left +1,
				&cons[r->start + s->left - cstart -1],
				WS2);

	    /* Filter specific low-complexity regions */
	    filter_words_local1(seq, fseq, len, "A",   12, 10, '#');
	    filter_words_local1(seq, fseq, len, "C",   12, 10, '#');
	    filter_words_local1(seq, fseq, len, "G",   12, 10, '#');
	    filter_words_local1(seq, fseq, len, "T",   12, 10, '#');
	    filter_words_local2(seq, fseq, len, "AC",  12, 10, '#');
	    filter_words_local2(seq, fseq, len, "AG",  12, 10, '#');
	    filter_words_local2(seq, fseq, len, "AT",  12, 10, '#');
	    filter_words_local2(seq, fseq, len, "CG",  12, 10, '#');
	    filter_words_local2(seq, fseq, len, "CT",  12, 10, '#');
	    filter_words_local2(seq, fseq, len, "GT",  12, 10, '#');
	}

#ifdef DEBUG_SEQ
	printf("SEQ #%"PRIrec" %d+%d %.*s\n",
	       s->rec, r->start + s->left-1, s->right - s->left -1,
	       s->right - s->left -1, seq+s->left+1);
	printf("OUT #%"PRIrec" %d+%d %.*s\n",
	       s->rec, r->start + s->left-1, s->right - s->left -1,
	       s->right - s->left -1, fseq+s->left+1);
#endif


	/* Check first good base */
	unique = 0;
	for (i = s->left+1; i < s->right && unique < MIN_OVERLAP; i++) {
#ifdef DEBUG
	    printf("<Pos %d: %c\n", i, fseq[i]);
#endif
	    if (legal_chars[(unsigned)fseq[i]])
		unique++;
	}
	first = i-1 + r->start;
	firstq = i-1;

	/* Check for last good base */
	unique = 0;
	for (i = s->right-1; i > s->left && unique < MIN_OVERLAP; i--) {
#ifdef DEBUG
	    printf(">Pos %d: %c\n", i, fseq[i]);
#endif
	    if (legal_chars[(unsigned)fseq[i]])
		unique++;
	}
	last = i + r->start;
	lastq = i;

#ifdef DEBUG
	printf("Seq %"PRIrec", first=%d, last=%d\n", s->rec, first, last);
#endif

	/*
	 * Adjust first/last to take into account mapping quality too. We want
	 * low quality regions to be flagged as gaps to validate, but
	 * we don't want to consider them in the clips[] array as
	 * it'll end up being the entire reading and is too blunt an
	 * instrument.
	 */
	while (firstq >= s->left+1 && fseq[firstq] == 'q')
	    firstq--;
	while (lastq < s->right && fseq[lastq] == 'q')
	    lastq++;
	if ((firstq += r->start) < r->start + MIN_OVERLAP+1)
	    firstq = r->start + MIN_OVERLAP+1;
	if ((lastq += r->start) > r->end - MIN_OVERLAP)
	    lastq = r->end - MIN_OVERLAP;

	if (clip_hash && (firstq != r->start + MIN_OVERLAP+1 ||
			  lastq  != r->end - MIN_OVERLAP)) {
	    clip_pos_t *c = malloc(sizeof(*c));
	    HashData hd;

	    if (!c)
		return NULL;

	    c->left  = firstq;
	    c->right = lastq;
	    hd.p = c;
	    HashTableAdd(clip_hash, (char *)&r->rec, sizeof(r->rec), hd, 0);
	}

	if (last >= first) {
	    if (!skip_non_paired || (r->library_rec && r->pair_rec)) {
#ifdef DEBUG_SEQ
		printf("Valid from %d to %d (cstart=%d)\n",
		       first-1, last, cstart);
#endif	    
		memset(&valid[first-1 -cstart], 1, last-first+1);
	    }
	}

	if (sorig != s)
	    free(s);

	xfree(fseq);
    }

    contig_iter_del(ci);

    gaps = coverage2contig_regions(valid, clip_depth, total_depth,
				   cstart, clen, end_skip_len);
    xfree(valid);
    xfree(cons);
    xfree(clip_depth);
    xfree(total_depth);

    return gaps;

 error:
    if (valid)
	xfree(valid);

    if (cons)
	xfree(cons);

    if (clip_depth)
	xfree(clip_depth);

    if (total_depth)
	xfree(total_depth);

    return NULL;
}

/* Cache the lib_type for this library rec so we can
 * identify whether it is expected to be pointing
 * inwards or outwards.
 */
static int compute_lib_type(GapIO *io, tg_rec library_rec, HashTable *lt_h, 
			    int *isize_min, int *isize_max, int *valid) {
    if (library_rec) {
	HashItem *hi;
	library_t *lib;

	if ((hi = HashTableSearch(lt_h, (char *)&library_rec,
				  sizeof(tg_rec)))) {
	    lib = hi->data.p;
	} else {
	    HashData hd;

	    update_library_stats(io, library_rec, 100,
				 NULL, NULL, NULL);
	    lib = cache_search(io, GT_Library, library_rec);
	    if (NULL == lib) goto give_up;
	    cache_incr(io, lib);
	    hd.p = lib;
	    HashTableAdd(lt_h, (char *)&library_rec,
			 sizeof(tg_rec),
			 hd, NULL);
	}

	/*
	 * FIXME, consider isize_min/max to be not gaussian, but just scan
	 * through the table to compute the actual size? Eg 99% of reads
	 * are between < and >, regardless of distribution shape?
	 *
	 * If so, this belongs in update_library_stats().
	 */
	*isize_max = *isize_min = lib->insert_size[lib->lib_type];
	*isize_min -= 3*lib->sd[lib->lib_type];
	*isize_max += 3*lib->sd[lib->lib_type];
	if (valid) *valid = (lib->flags & 2) ? 0 : 1;
	return lib->lib_type;
    } else {
    give_up:
	*isize_min = 20;
	*isize_max = 2000;
	if (valid) *valid = 0;
	return LIB_T_INWARD;
    }
}

/*
 * Determines whether a pair WITHIN the same contig should be
 * considered as consistent or not.
 *
 * Returns severity.
 *         1 for OK
 *         0 for warning/error, with *severity filled out if not NULL
 *
 * Severity ratings:
 *         0 for OK
 *         1 for OK, but slightly too large/small
 *         2 for contig spanning (shouldn't happen in this code)
 *         3 for OK
 */
static int consistent_pair(GapIO *io, rangec_t *r, HashTable *lt_h,
			   int *severity) {
    int isize_min, isize_max, isize, lib_type;

    lib_type = compute_lib_type(io, r->library_rec, lt_h,
				&isize_min, &isize_max, NULL);

    /* Check COMP1 vs COMP2 */
    switch (lib_type) {
    case LIB_T_INWARD:
	if (r->start < r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) == 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) != 0))
	    return severity?*severity=3,0:0;

	if (r->start > r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) != 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) == 0))
	    return severity?*severity=3,0:0;
	break;

    case LIB_T_OUTWARD:
	if (r->start > r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) == 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) != 0))
	    return severity?*severity=3,0:0;

	if (r->start < r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) != 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) == 0))
	    return severity?*severity=3,0:0;
	break;

    case LIB_T_SAME:
	if (((r->flags & GRANGE_FLAG_COMP1) == 0) !=
	    ((r->flags & GRANGE_FLAG_COMP2) == 0))
	    return severity?*severity=3,0:0;
	break;
    }


    /* Check insert size */
    isize = MAX(MAX(r->pair_start, r->pair_end), MAX(r->start, r->end))
	  - MIN(MIN(r->pair_start, r->pair_end), MIN(r->start, r->end));
    isize = ABS(isize);
    if (!(isize >= isize_min/3 && isize <= isize_max*2))
	return severity?*severity=3,0:0; // Far too big/small
    if (!(isize >= isize_min && isize <= isize_max))
	return severity?*severity=1,0:0; // Slightly too big/small

    return severity?*severity=0,1:1;
}

/*
 * Our range query brought back all data within any library insert size of
 * the problem region start..end.
 *
 * However many libraries may have smaller expected insert sizes, so discount
 * reads from these libraries if they are too far away from the problem
 * region to be considered as valid confirmation / denial.
 */
static int pair_in_range(GapIO *io, rangec_t *r, HashTable *lt_h,
			 int start, int end) {
    int isize_min, isize_max, isize, lib_type;

    lib_type = compute_lib_type(io, r->library_rec, lt_h,
				&isize_min, &isize_max, NULL);

    if (r->start      < start && start - r->start      < isize_max) return 1;
    if (r->end        < start && start - r->end        < isize_max) return 1;
    if (r->pair_start < start && start - r->pair_start < isize_max) return 1;
    if (r->pair_end   < start && start - r->pair_end   < isize_max) return 1;

    if (r->start      > end && r->start      - end < isize_max) return 1;
    if (r->end        > end && r->end        - end < isize_max) return 1;
    if (r->pair_start > end && r->pair_start - end < isize_max) return 1;
    if (r->pair_end   > end && r->pair_end   - end < isize_max) return 1;
    
    return 0;
}

struct range_loc {
    RB_ENTRY(range_loc) link;
    int start;
    int end;
};

int rl_cmp(struct range_loc *l1, struct range_loc *l2) {
    int d = l1->end - l2->end;
    // exact matches still need to return +/- otherwise RB_INSERT will fail
    return d ? d : (l1 < l2 ?1 :-1);
}


RB_HEAD(rlTREE, range_loc);
RB_PROTOTYPE(rlTREE, range_loc, link, rl_cmp);
RB_GENERATE(rlTREE, range_loc, link, rl_cmp);

/*
 * FIXME: implement per library (if sufficiently well sampled).
 */
static void dump_template_dist(GapIO *io, tg_rec contig) {
    HashTable *lt_h = NULL;
    int cstart, cend;
    int *tdist = NULL;
    contig_iterator *ci;
    rangec_t *r;
    int i, n, last_i;
    int64_t sum, sum_sq;
    struct range_loc *node, *next;

    struct rlTREE rltree = RB_INITIALIZER(&rltree);

    lt_h = HashTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (!lt_h)
	goto cleanup;

    if (consensus_valid_range(io, contig, &cstart, &cend) == -1)
	goto cleanup;

    if (!(tdist = calloc(cend-cstart+1, sizeof(*tdist))))
	goto cleanup;

    ci = contig_iter_new(io, contig, 1, CITER_FIRST |CITER_ISTART |CITER_PAIR,
			 CITER_CSTART, CITER_CEND);

    n = 0;
    sum = 0;
    sum_sq = 0;
    last_i = cstart;
    while (NULL != (r = contig_iter_next(io, ci))) {
	int st, en, severity;

	sequence_get_range_pair_position(io, r, contig, 0);

	if (contig != r->pair_contig)
	    continue;

	if (!consistent_pair(io, r, lt_h, &severity)) {
	    // Reject unless this read is overlapping the boundary of
	    // the template. Hence extreme long range templates don't count
	    // except for when we're at the actual ends.

	    //printf("Rec %"PRIrec" inconsistent pair\n", r->rec);
	    if (severity > 1)
		continue;
	}

	// Only do left end per pair
	if (r->pair_start < r->start ||
	    (r->pair_start == r->start && r->pair_rec < r->rec))
	    continue;
	
	en = MAX(MAX(r->pair_start, r->pair_end), MAX(r->start, r->end));
	st = MIN(MIN(r->pair_start, r->pair_end), MIN(r->start, r->end));

	//printf("Rec %"PRIrec" consistent pair %d..%d\n", r->rec, st, en);
	if (st < cstart) st = cstart;
	if (en > cend)   en = cend;
	for (i = st; i <= en; i++)
	    tdist[i-cstart]++;

	//printf("#%"PRIrec" %d..%d %d..%d\n", r->rec, r->start, r->end, st, en);

	/* Check tree for templates ending before start 'st' */
	//RB_FOREACH(node, rlTREE, &rltree) {
	for (node = RB_MIN(rlTREE, &rltree); node != NULL; node = next) {
	    next = RB_NEXT(rlTREE, &rltree, node);
	    if (node->end >= st)
		break;
	    //printf("    Node %d..%d\n", node->start, node->end);
	    if (node->end >= last_i)
		printf("  depth %d %d\t%d\t%f\t%f\n", last_i, node->end, n,
		       (double)sum/n,
		       31/sqrt(n));
		    //sqrt((double)sum_sq/n-((double)sum/n)*((double)sum/n)));
	    last_i = node->end+1;
	    RB_REMOVE(rlTREE, &rltree, node);
	    n--;
	    sum -= node->end - node->start;
	    sum_sq -= (node->end - node->start) * (node->end - node->start);
	    free(node);
	}

	/*
	 * Given N samples from a distribution with mean M and SD S,
	 * we can calculate the observed mean with an SD of S/sqrt(N).
	 * 
	 * This means we have a quick way of estimating whether an insert size
	 * adjustment is significant.
	 */
	
//	if (last_i == 6183) {
//	    RB_FOREACH(node, rlTREE, &rltree) {
//		printf("@6183 size %d\n", node->end - node->start + 1);
//	    }
//	}

	//if (last_i == 6183) {
	//if (last_i == 62321) {
	//if (last_i == 103605) {
	//if (st-1 >= last_i) {
	if (0) {
	    double chi2 = 0;
	    double count = 0;
	    int stats[500], N = 0, x;

	    memset(stats, 0, 500*sizeof(*stats));

	    library_t *lib = cache_search(io, GT_Library, 10); // FIXME
	    RB_FOREACH(node, rlTREE, &rltree) {
		int isz  = node->end - node->start + 1;
		int ibin = isize2ibin(isz);
		int iwid = ibin_width(isz);
		if (isz > 0 && isz < 500) {
		    stats[isz]++;
		    N++;
		}
		count += (double)lib->size_hist[0][ibin] / iwid;
	    }
	    //printf("@count at %d = %f, N=%d\n", last_i, count, N);
	    for (x = 0; x < 500; x++) {
		if (!stats[x])
		    continue;

		int isz  = x;
		int ibin = isize2ibin(x);
		int iwid = ibin_width(x);
		double e = ((double)lib->size_hist[0][ibin]/iwid) / count * N;
		//printf("@%d %d %f\n", x, stats[x], e);

		chi2 += (stats[x]-e)*(stats[x]-e)/e;
	    }
	    //printf("@chi2 = %f, N=%d\n", chi2, N);
	    
	    {
		int dof = N-1, i;
		double lgamma = 0, p;
		for (i = 1; i <= (dof/2)-1; i++)
		    lgamma += log(i);
		
		p = (dof/2.0 - 1)*log(chi2) - chi2/2 - (dof/2.0)*log(2)
		    - lgamma;
		//printf("@ln(p) = %f, p=%g\n", p, exp(p));

		printf("%f", p);
	    }
	}
	
	if (st-1 >= last_i)
	    printf("  Depth %d %d\t%d\t%f\t%f\n", last_i, st-1, n,
		   (double)sum/n,
		   31/sqrt(n));
		   //sqrt((double)sum_sq/n - ((double)sum/n)*((double)sum/n)));
	last_i = st;

	/* Insert new node */
	node = malloc(sizeof(*node));
	node->start = st;
	node->end   = en;
	RB_INSERT(rlTREE, &rltree, node);
	n++;
	sum += en-st;
	sum_sq += (en-st)*(en-st);
    }
    contig_iter_del(ci);

    printf("### n1=%d\n", n);
    n = 0;
    RB_FOREACH(node, rlTREE, &rltree)
	printf("Remaining %d..%d\n", node->start, node->end),n++;
    printf("### n2=%d\n", n);

    printf("\n\nDist for contig =%"PRIrec"\n", contig);
    for (i = cstart; i <= cend; i++) {
	printf("%d\t%d\n", i, tdist[i-cstart]);
    }
    printf("\n");

 cleanup:
    if (lt_h)
	HashTableDestroy(lt_h, 0);

    if (tdist)
	free(tdist);

    RB_FOREACH(node, rlTREE, &rltree) {
	RB_REMOVE(rlTREE, &rltree, node);
	free(node);
    }
}

static double gc(GapIO *io, tg_rec contig, int start, int end) {
    static char *cons = NULL;
    static tg_rec ctg = 0;
    static int cstart;
    int clen;
    contig_t *c;
    int i, gc = 0, n = 0;

    if (ctg != contig) {
	ctg  = contig;

	c = cache_search(io, GT_Contig, ctg);
	clen = c->end - (cstart = c->start) + 1;
	if (cons) free(cons);
	cons = malloc(clen);

	calc_consensus(ctg, c->start, c->end, CON_SUM, cons,
		       NULL, NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);
    }

    for (i = start; i <= end; i++) {
	switch(cons[i-cstart]) {
	case 'A':
	case 'T':
	    n++;
	    break;
	    
	case 'G':
	case 'C':
	    n++;
	    gc++;
	}
    }

    return n ? (double)gc/n : 0;
}


/*
 * Uses read-pair information to confirm whether a gap appears to be
 * valid (or invalid).
 *
 * 1. Find all read-pairs in the gap +/- the expected template sizes.
 * 2. Mark these with good, bad and unknown (spanning).
 * 3. Score... (also consider spanning as bad?)
 */
static void confirm_gaps(GapIO *io, tg_rec contig, Array gaps,
			 int unique_mqual,
			 int good_score, int good_unique_score,
			 int bad_score, int bad_unique_score,
			 int large_score, int large_unique_score,
			 int spanning_score, int spanning_unique_score,
			 int singleton_score, int singleton_unique_score,
			 int min_score) {
    int i, severity;
    HashTable *lt_h;
    HashIter *iter;
    HashItem *hi;
    int cstart, cend;
    int last_gap, next_gap;
    int first_pass = 1;
    
    lt_h = HashTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (!lt_h)
	return;

    iter = HashTableIterCreate();
    if (!iter)
	goto cleanup;

    if (consensus_valid_range(io, contig, &cstart, &cend) == -1)
	goto cleanup;

 second_pass:
    last_gap = INT_MIN;

    /* Now process gaps validating by read-pair */
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	contig_t *c;
	int score;
	int j, nr;
	rangec_t *r;
	HashTable *h;
	int num_unique_good = 0;
	int num_good = 0;
	int num_unique_bad = 0;
	int num_bad = 0;
	int num_unique_large = 0;
	int num_large = 0;
	int num_unique_spanning = 0;
	int num_spanning = 0;
	int num_unique_single = 0;
	int num_single = 0;

	if (gap->deleted) {
	    abort();
	    continue;
	}

	if (gap->valid)
	    continue;

	for (j = i+1; j < ArrayMax(gaps); j++) {
	    if (!(arrp(contig_region_t, gaps, j))->valid)
		break;
	}
	next_gap = j < ArrayMax(gaps) 
	    ? (arrp(contig_region_t, gaps, j))->start
	    : INT_MAX;

	printf("Gap %d..%d ", gap->start, gap->end);

	if (!(c = cache_search(io, GT_Contig, contig)))
	    continue;

//FIXME: use template_max_size(io) instead.
#define INS_SIZE 10000

	if (first_pass) {
	    r = contig_seqs_in_range(io, &c,
				     gap->start-INS_SIZE, gap->end+INS_SIZE,
				     CSIR_PAIR, &nr);
	    if (!r)
		continue;
	} else {
	    r = gap->r;
	    nr = gap->nr;
	}

	h = HashTableCreate(1024, HASH_DYNAMIC_SIZE |
			          HASH_POOL_ITEMS |
			          HASH_NONVOLATILE_KEYS);
	if (!h)
	    continue;

	/* Pair up all seqs */
	for (j = 0; j < nr; j++) {
	    HashData hd;
	    int new;

	    /*
	     * Reads without read-pairs that come from a paired end library
	     * are inherently suspicious. Why don't we have the other end?
	     * Presumably it didn't align, which typically means we have
	     * a contig missing.
	     */
	    if (!r[j].pair_rec) {
		int x, valid = 0, lib_type;
		int isize_min, isize_max, isize_mid;

		if (!r[j].library_rec)
		    continue;

		lib_type = compute_lib_type(io, r[j].library_rec, lt_h,
					    &isize_min, &isize_max,
					    &valid);
		if (!valid)
		    continue;

		/* 3 s.d. is a bit harsh for penalising things, so limit
		 * it to 1 s.d.
		 */
		isize_mid = (2*isize_min + isize_max) / 3;

		// Other end didn't align?
		switch (lib_type) {
		case LIB_T_INWARD:
		    if ((r[j].flags & GRANGE_FLAG_COMP1) == 0) {
			if (gap->start > r[j].start && 
			    gap->start - r[j].start < isize_mid) {
			    if (r[j].mqual >= unique_mqual)
				num_unique_single++;
			    else
				num_single++;
			}
		    } else {
			if (r[j].end > gap->end &&
			    r[j].end - gap->end < isize_mid) {
			    if (r[j].mqual >= unique_mqual)
				num_unique_single++;
			    else
				num_single++;
			}
		    }
		    break;

		case LIB_T_OUTWARD:
		    if ((r[j].flags & GRANGE_FLAG_COMP1) != 0) {
			if (gap->start > r[j].start &&
			    gap->start - r[j].start < isize_mid) {
			    if (r[j].mqual >= unique_mqual)
				num_unique_single++;
			    else
				num_single++;
			}
		    } else {
			if (r[j].end > gap->end &&
			    r[j].end - gap->end < isize_mid) {
			    if (r[j].mqual >= unique_mqual)
				num_unique_single++;
			    else
				num_single++;
			}
		    }
		    break;

		case LIB_T_SAME:
		    if (r[j].mqual >= unique_mqual)
			num_unique_single++;
		    else
			num_single++;
		    break;
		}
		continue;
	    }

	    /* Paired in this contig, so see if also spans the problem. */
	    if ((hi = HashTableSearch(h, (char *)&r[j].pair_rec,
				      sizeof(r[j].pair_rec)))) {
		rangec_t *p = (rangec_t *)hi->data.p;
		sequence_get_range_pair_position(io, p, contig, 0);

		if (MIN(p->start, p->pair_start) <= gap->start-MIN_OVERLAP &&
		    MAX(p->end, p->pair_end) >= gap->end+MIN_OVERLAP &&
		    /*(MIN(p->start, p->pair_start) >= last_gap ||
		      MAX(p->end, p->pair_end) <= next_gap) && */
		    pair_in_range(io, p, lt_h, gap->start, gap->end)) {
		    
		    /* Spans the gap */
		    if (consistent_pair(io, p, lt_h, &severity)) {
			int unique = (p->mqual >= unique_mqual &&
				      p->pair_mqual >= unique_mqual);
			if (unique)
			    num_unique_good++;
			else
			    num_good++;
		    } else {
			int unique = (p->mqual >= unique_mqual &&
				      p->pair_mqual >= unique_mqual);
			if (severity == 1) {
			    if (unique)
				num_unique_large++;
			    else
				num_large++;
			} else {
			    if (unique)
				num_unique_bad++;
			    else
				num_bad++;
			}
		    }
		}
		HashTableDel(h, hi, 0);
		continue;
	    }

	    hd.p = &r[j];
	    hi = HashTableAdd(h, (char *)&r[j].rec, sizeof(r[j].rec),
			      hd, &new);
	    if (!new) {
		fprintf(stderr, "Error: seq already in hash\n");
		HashTableDel(h, hi, 0);
	    }
	}

	/*
	 * Count unmatched pairs remaining in the hash table.
	 *
	 * Some of these will be due to the edges of our search
	 * window; they're in this contig, but pointing outwards. We
	 * discard those.
	 *
	 * The remainder are either singletons or contig spanning
	 * read-pairs. Singletons are counted as "unknowns" if it's a
	 * paired end library. Spanning read-pairs are counted as bad
	 * if and only if they're pointing towards the problem and
	 * aren't within the likely insert-size of the end of the
	 * contig.
	 */
	HashTableIterReset(iter);

	while (((hi = HashTableIterNext(h, iter)))) {
	    rangec_t *r = hi->data.p;
	    sequence_get_range_pair_position(io, r, contig, 0);
	    int isize_min, isize_max, isize_mid, lib_type, unique;

	    /*
	     * Paired end libraries with singletons should be
	     * considered as an unknown states.
	     *
	     * Single ended libraries obviously are considered as OK,
	     * although they probably should have been filtered out
	     * before reachign this point.
	     */
	    assert(r->pair_contig);
	    if (!r->pair_contig) {
		/* Shouldn't get here as we only added pairs */
		int x, valid;
		compute_lib_type(io, r->library_rec, lt_h, &x, &x, &valid);
		if (valid) {
		    if (r->mqual >= unique_mqual)
			num_unique_single++;
		    else
			num_single++;
		}
		continue;
	    }
	    /* Check consistency */
	    if (contig == r->pair_contig) {
		/* Make sure it spans the gap */
		if (!(MIN(r->start, r->pair_start) <= gap->start-MIN_OVERLAP &&
		      MAX(r->end, r->pair_end) >= gap->end+MIN_OVERLAP))
		    continue;

		/*
		if (!(MIN(r->start, r->pair_start) >= last_gap ||
		      MAX(r->end, r->pair_end) <= next_gap))
		    continue;
		*/

		if (!pair_in_range(io, r, lt_h, gap->start, gap->end))
		    continue;

		if (consistent_pair(io, r, lt_h, &severity)) {
		    int unique = (r->mqual >= unique_mqual &&
				  r->pair_mqual >= unique_mqual);
		    if (unique)
			num_unique_good++;
		    else
			num_good++;
		} else {
		    int unique = (r->mqual >= unique_mqual &&
				  r->pair_mqual >= unique_mqual);
		    if (severity == 1) {
			if (unique)
			    num_unique_large++;
			else
			    num_large++;
		    } else  {
			if (unique)
			    num_unique_bad++;
			else
			    num_bad++;
		    }
		}
		continue;
	    }

	    /* Spanning => bad, but only if it points towards this gap */
	    lib_type = compute_lib_type(io, r->library_rec, lt_h,
					&isize_min, &isize_max, NULL);
	    isize_mid = (2*isize_min + isize_max) / 3;

	    unique = (r->mqual >= unique_mqual &&
		      r->pair_mqual >= unique_mqual);

	    switch (lib_type) {
	    case LIB_T_INWARD:
		if ((r->flags & GRANGE_FLAG_COMP1) == 0) {
		    if (cend - r->start < isize_max)
			// pointing off end of contig, maybe real join
			break;

		    if (gap->start > r->start &&
			gap->start - r->start < isize_mid) {
			// just left of gap
			if (unique)
			    num_unique_spanning++;
			else
			    num_spanning++;
		    }
		} else {
		    if (r->end - cstart < isize_max)
			break;

		    if (r->end > gap->end && r->end - gap->end < isize_mid) {
			if (unique)
			    num_unique_spanning++;
			else
			    num_spanning++;
		    }
		}
		break;

	    case LIB_T_OUTWARD:
		if ((r->flags & GRANGE_FLAG_COMP1) != 0) {
		    if (cend - r->start < isize_max)
			break;

		    if (gap->start > r->start &&
			gap->start - r->start < isize_mid) {
			if (unique)
			    num_unique_spanning++;
			else
			    num_spanning++;
		    }
		} else {
		    if (r->end - cstart < isize_max)
			break;

		    if (r->end > gap->end && r->end - gap->end < isize_mid) {
			if (unique)
			    num_unique_spanning++;
			else
			    num_spanning++;
		    }
		}
		break;

	    case LIB_T_SAME:
		/* We can't tell which is 1st and 2nd read easily? */
		if (gap->start > r->start && gap->start - r->start < isize_mid){
		    if (unique)
			num_unique_spanning++;
		    else
			num_spanning++;
		} else if (r->end > gap->end && r->end - gap->end < isize_mid) {
		    if (unique)
			num_unique_spanning++;
		    else
			num_spanning++;
		}
		break;

	    default:
		if (unique)
		    num_unique_bad++;
		else
		    num_bad++;
	    }
	}

        score =           num_good * good_score
              +    num_unique_good * good_unique_score
              +            num_bad * bad_score
              +     num_unique_bad * bad_unique_score
              +          num_large * large_score
              +   num_unique_large * large_unique_score
              +       num_spanning * spanning_score
              +num_unique_spanning * spanning_unique_score
              +      num_single    * singleton_score
              +num_unique_single   * singleton_unique_score;

	//score -= (gap->end - gap->start) / 5;
 
	printf("%d/%d good, %d/%d bad, %d/%d small/large, %d/%d spanning, "
	       "%d/%d single => score %d\n",
	       num_unique_good,     num_good,
	       num_unique_bad,      num_bad,
	       num_unique_large,    num_large,
	       num_unique_spanning, num_spanning,
	       num_unique_single,   num_single,
	       score);
	gap->score = score;
	gap->num_unique_good = num_unique_good;
	gap->num_good = num_good;
	gap->num_unique_bad = num_unique_bad;
	gap->num_bad = num_bad;
	gap->num_unique_large = num_unique_large;
	gap->num_large = num_large;
	gap->num_unique_spanning = num_unique_spanning;
	gap->num_spanning = num_spanning;
	gap->num_unique_single = num_unique_single;
	gap->num_single = num_single;

	gap->valid = score >= min_score ? 1 : 0;

	HashTableDestroy(h, 0);

	if (gap->valid) {
	    free(r);
	} else {
	    gap->r  = r;
	    gap->nr = nr;

	    last_gap = gap->end;
	}

	fflush(stdout);
    }

//    if (first_pass) {
//	first_pass = 0;
//	goto second_pass; /* Sorry! It's a "temporary" trial */
//    }

 cleanup:
    /* Clean up cached library information */
    if (NULL != iter) {
	HashTableIterReset(iter);
	while (NULL != (hi = HashTableIterNext(lt_h, iter))) {
	    library_t *lib = hi->data.p;
	    if (lib) cache_decr(io, lib);
	}
	HashTableIterDestroy(iter);
    }
    HashTableDestroy(lt_h, 0);

    //uninit_template_checks(io, tarr);
}


/*
 * Merges gaps if they're close together, within min_distance apart.
 */
static void merge_gaps(Array gaps, int min_distance) {
    int i, j;
    contig_region_t *last_gap;

    if (ArrayMax(gaps) == 0)
	return;

    last_gap = arrp(contig_region_t, gaps, 0);
    
    for (i = 1; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);

	if (gap->start - last_gap->end < min_distance) {
	    printf("Merging gap %d..%d with %d..%d\n",
		   last_gap->start, last_gap->end,
		   gap->start, gap->end);

	    last_gap->end = MAX(last_gap->end, gap->end);
	    last_gap->valid |= gap->valid;
	    gap->deleted = 1;
	} else {
	    last_gap = gap;
	}
    }

    for (i = j = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	if (gap->deleted)
	    continue;

	if (i == j)
	    continue;

	arr(contig_region_t, gaps, j) = *gap;
	j++;
    }
    ArrayMax(gaps) = j;
}


/*
 * This code doesn't actually make the break itself, rather it just finds
 * where and how to break. It also produces the tags if required.
 *
 * Analyses gaps to work out where to break. We'll identify a left contig,
 * a right contig, and possibly multiple single-read contigs for readings
 * that are contained entirely within the gap section.
 *
 * Tag_params contains the parameters used for auto-break and are added to
 * annotation covering the break site. (It should end up in both contigs.)
 * It can be NULL, in which case no annotations are created.
 *
 * Fills out 'ds' with a list of breaks to make. Each item is in itself
 * another list with the first two items being a contig number and position
 * within the contig to run Break Contig at with the remainder of the list
 * being an array of reading numbers to run disassembly on (these are
 * reads that fall within the gap and should be taken out to form their
 * own contig).
 */
static void break_gaps(GapIO *io, tg_rec contig, Array gaps, char *tag_params,
		       HashTable *clip_hash, dstring_t *ds) {
    int i, j, new_start;
    tg_rec right_start = 0;

    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	int worst_score, worst_idx;

	if (gap->deleted)
	    continue;

	if (gap->valid) {
	    printf("Skipping as marked as valid\n");
	    continue;
	}

	/* Only break a problem if it is the worst problem within a
	 * stretch of invalid regions.
	 *
	 * This avoids the issue where a single problem can cause all
	 * scores within a region to reduce due to invalid large inserts
	 * spanning multiple potential problem sites.
	 */
	worst_score = INT_MAX;
	worst_idx = i;
	for (j = i; j < ArrayMax(gaps); j++) {
	    contig_region_t *gap2 = arrp(contig_region_t, gaps, j);
	    if (gap2->valid || gap2->start - gap->start > MAX_NEIGHBOUR)
		break;
	    if (worst_score > gap2->score) {
		worst_score = gap2->score;
		worst_idx = j;
	    }
	}

	/* For consistent reporting output only */
	for (j = i; j < ArrayMax(gaps); j++) {
	    contig_region_t *gap2 = arrp(contig_region_t, gaps, j);
	    if (gap2->valid || gap2->start - gap->start > MAX_NEIGHBOUR)
		break;
	    if (gap2->score > worst_score) {
		printf("Gap from %d to %d\n", gap2->start, gap2->end);
		printf("Skipping as poorer neighbouring problem; %d vs %d\n",
		       gap2->score, worst_score);
	    }
	}

	i = --j;
	gap = arrp(contig_region_t, gaps, worst_idx);

	printf("Gap from %d to %d\n", gap->start, gap->end);

	/* 
	 * Identify readings that after clipping are entirely contained
	 * within the region. These are candidates for disassembling
	 * completely
	 */
	dstring_appendf(ds, " {");
	
	printf("  New starting point for right contig = %d\n", gap->start+1);
	dstring_appendf(ds, "%"PRIrec" %d", contig, gap->start+1);

	new_start = INT_MAX;
	for (j = 0; j < gap->nr; j++) {
	    int left, right;
	    HashItem *hi;

	    hi = HashTableSearch(clip_hash, (char *)&gap->r[j].rec,
				 sizeof(tg_rec));
	    if (hi) {
		clip_pos_t *c = (clip_pos_t *)hi->data.p;
		left = c->left;
		right = c->right;
	    } else {
		left = gap->r[j].start + MIN_OVERLAP+1;
		right = gap->r[j].end - MIN_OVERLAP;
	    }

	    //printf("Read #%"PRIrec" %d..%d clipped %d..%d\n",
	    //	   gap->r[j].rec,
	    //	   gap->r[j].start, gap->r[j].end,
	    //	   left, right);

	    if (left >= gap->start && right <= gap->end) {
		printf("  Read #%"PRIrec" to self-contig\n", gap->r[j].rec);
		dstring_appendf(ds, " #%"PRIrec, gap->r[j].rec);
	    } else if (gap->r[j].end > gap->start &&
		       new_start > gap->r[j].start) {
		seq_t *s = cache_search(io, GT_Seq, gap->r[j].rec);
		seq_t *sorig = s;
		if ((s->len < 0) ^ gap->r[j].comp) {
		    s = dup_seq(s);
		    complement_seq_t(s);
		}
		if (gap->r[j].start + s->left-1 > gap->start) {
		    if (new_start > gap->r[j].start + s->left-1)
			new_start = gap->r[j].start + s->left-1;
		}

		if (sorig != s)
		    free(s);
	    }
	}

	dstring_appendf(ds, "}");

	if (tag_params && !io->read_only) {
	    dstring_t *ds = dstring_create(NULL);
	    if (new_start == INT_MAX)
		new_start = gap->end;
	    new_start = MAX(new_start, gap->end);
	    dstring_appendf(ds, "Score:           %d\n", gap->score);
	    dstring_appendf(ds, "#good pairs:     %d/%d\n",
			    gap->num_unique_good, gap->num_good);
	    dstring_appendf(ds, "#bad pairs:      %d/%d\n",
			    gap->num_unique_bad, gap->num_bad);
	    dstring_appendf(ds, "#large pairs:    %d/%d\n",
			    gap->num_unique_large, gap->num_large);
	    dstring_appendf(ds, "#spanning pairs: %d\n", gap->num_spanning);
	    dstring_appendf(ds, "#singletons:     %d\n", gap->num_single);
	    dstring_appendf(ds, "\nParams:\n%s", tag_params);
	    anno_ele_add(io, GT_Contig, contig, 0, str2type("ABRK"),
			 dstring_str(ds), gap->start, new_start, ANNO_DIR_NUL);
	    dstring_destroy(ds);
	}

	free (gap->r);
    }
}

void auto_break_single_contig(GapIO *io, tg_rec contig, int start, int end,
			      int end_skip_len,
			      int64_t tw, double filter_score,
			      int filter_consensus, double depth,
			      int min_mq, int min_score, int unique_mqual,
			      int good_score, int good_unique_score,
			      int bad_score, int bad_unique_score,
			      int large_score, int large_unique_score,
			      int spanning_score, int spanning_unique_score,
			      int singleton_score, int singleton_unique_score,
			      dstring_t *ds) {
    Array gaps;
    HashTable *clip_hash;
    dstring_t *p_ds = dstring_create(NULL);

    dstring_appendf(p_ds, "Minimum mapping quality:  %d\n", min_mq);
    dstring_appendf(p_ds, "Minimum break score:      %d\n", min_score);
    dstring_appendf(p_ds, "Weight for good pair:     %d/%d\n",
		    good_unique_score, good_score);
    dstring_appendf(p_ds, "Weight for bad  pair:     %d/%d\n",
		    bad_unique_score, bad_score);
    dstring_appendf(p_ds, "Weight for large pair:    %d/%d\n",
		    large_unique_score, large_score);
    dstring_appendf(p_ds, "Weight for spanning pair: %d/%d\n",
		    spanning_unique_score, spanning_score);
    dstring_appendf(p_ds, "Weight for singleton:     %d/%d\n",
		    singleton_unique_score, singleton_score);

    printf("\n=== Checking contig %"PRIrec" ===\n", contig);

    clip_hash = HashTableCreate(1024,
				HASH_POOL_ITEMS);

    printf("  = Identifying suspect joins\n");
    fflush(stdout);
    gaps = suspect_joins(io, contig, tw, filter_score, filter_consensus,
			 depth, min_mq, clip_hash, end_skip_len);
    fflush(stdout);

    printf("  = Merging gaps\n");
    merge_gaps(gaps, MIN_OVERLAP*2);
    fflush(stdout);
    //dump_gaps(gaps);

    printf("  = Confirming gaps\n");
    confirm_gaps(io, contig, gaps, unique_mqual,
		 good_score, good_unique_score,
		 bad_score, bad_unique_score,
		 large_score, large_unique_score,
		 spanning_score, spanning_unique_score,
		 singleton_score, singleton_unique_score,
		 min_score);

    printf("  = Finding break points\n"); // Doesn't actually do the break
    break_gaps(io, contig, gaps, dstring_str(p_ds), clip_hash, ds);

    HashTableDestroy(clip_hash, 1);
    ArrayDestroy(gaps);
    dstring_destroy(p_ds);
    //xfree(clips);
}

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      int end_skip_len,
			      double filter_score, int filter_consensus,
			      int min_mq, int min_score, int unique_mqual,
			      int good_score, int good_unique_score,
			      int bad_score, int bad_unique_score,
			      int large_score, int large_unique_score,
			      int spanning_score, int spanning_unique_score,
			      int singleton_score, int singleton_unique_score)
{
    int i;
    int64_t tw = 0;
    double gc;
    int depth = 1; // filter by consensus

    dstring_t *ds = dstring_create(NULL);

    if (filter_consensus) {
	tw = word_count_cons(io, argc, argv);
	//print_counts(10);
    }
    
    for (i = 0; i < argc; i++) {
	//dump_template_dist(io, argv[i].contig);

	auto_break_single_contig(io, argv[i].contig,
				 argv[i].start, argv[i].end, end_skip_len, 
				 tw, filter_score, filter_consensus, depth,
				 min_mq, min_score, unique_mqual,
				 good_score, good_unique_score,
				 bad_score, bad_unique_score,
				 large_score, large_unique_score,
				 spanning_score, spanning_unique_score,
				 singleton_score, singleton_unique_score,
				 ds);
    }

    return ds;
}
