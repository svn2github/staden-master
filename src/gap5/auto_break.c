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

#define MIN3(a,b,c) (MIN(MIN((a),(b)),(c)))
#define MAX3(a,b,c) (MAX(MAX((a),(b)),(c)))

#ifndef WS
#    define WS 12
#endif
#define WS2 ((int)(WS/2))

#define ALLB(ws) ((1<<(2*(ws)))-1)

#define MIN_OVERLAP 5
#define MAXTSIZE 10000
#define MIN_HARD_SCORE -9
#define MIN_SOFT_SCORE 0
#define MIN_SOFT_VALID 4
/* #define NORMALISE_FOR_GC 1 */

#define CONTIG_END_IGNORE 200

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
 * This checks if an inconsistent template appears to be capable of
 * spanning 'gap'.
 *
 * As there are so many reasons why a template can be inconsistent this
 * iterates through all readings on the template in 'contig' to see if the
 * range of possibilities for the "other end" could be on the opposite side
 * to the gap.
 *
 * eg where F is invalid and may span if repositioned.
 * <---F---.........---R--->   |gap|
 *  ---F--->........---R--->   |gap|
 *
 * eg where no span is possible unless both readings are misassembled.
 * <---F---........<---R---    |gap|
 *
 * Returns 1 for potentially spanning
 *         0 for not spanning
 */
#if 0
static int bad_template_span(GapIO *io, int contig, template_c *t,
			     contig_region_t *gap) {
    GTemplates te;
    item_t *item;

    template_read(io, t->num, te);

    for (item = head(t->gel_cont); item; item = item->next) {
	GReadings r;
	gel_cont_t *gc = (gel_cont_t *)(item->data);

	if (gc->contig != contig)
	    continue;

	gel_read(io, gc->read, r);

	/* ---> |gap| */
	if (r.position < gap->start &&
	    r.sense == GAP_SENSE_ORIGINAL &&
	    r.position + te.insert_length_max > gap->end) {
	    return 1;
	}

	/* |gap| <--- */
	if (r.position + r.sequence_length - 1 > gap->end &&
	    r.sense != GAP_SENSE_ORIGINAL &&
	    r.position + r.sequence_length - 1 - te.insert_length_max <
	        gap->start) {
	    return 1;
	}
    }

    return 0;
}
#endif

/*
 * Checks the region from 'start' to 'end' to see whether it can be
 * validated by nearby read-pairs.
 *
 * It does this by looking up to MAXTSIZE either size of start/end
 * identifying readings that point inwards to the hole. For these reads it
 * then verifies if they validate, have no impact, or reject the region.
 *
 * Returns a score, +ve for confirmed, -ve for denied and 0 for unknown.
 */
static int check_read_pairs(GapIO *io, tg_rec contig,
			    contig_region_t *gap) {
    int i;
    int valid = 0, invalid = 0;

#ifdef DEBUG    
    printf("\n***** Contig %"PRIrec", gap %d..%d *****\n",
	   contig, gap->start, gap->end);
#endif

    if (valid - invalid < MIN_HARD_SCORE ||
	(valid - invalid < MIN_SOFT_SCORE &&
	 valid < MIN_SOFT_VALID))
	gap->valid = 0;
    else
	gap->valid = 1;

    return valid - invalid;
}

/*
 * Turns an array of coverage data (0 => no coverage) to an Array of
 * contig_region_t structs.
 */
static Array coverage2contig_regions(char *valid, int cstart, int clen) {
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
	/*
	r = (contig_region_t *)ArrayRef(gaps, ngaps);
	r->start = 1;
	r->end = i;
	r->deleted = 0;
	r->valid = 1;
	ngaps++;
	*/
    }

    for (; i < clen; i++) {
	if (!valid[i]) {
	    contig_region_t *r;
	    for (j = i+1; j < clen; j++) {
		if (valid[j])
		    break;
	    }
	    if (j >= clen) {
		/* Skip end gaps */
		printf("END gap %d..%d\n", i+cstart, j-1+cstart);
		break;
	    } else {
		printf("INTERNAL gap %d..%d\n", i+cstart, j-1+cstart);
	    }
	    r = (contig_region_t *)ArrayRef(gaps, ngaps);
	    r->start = i;
	    r->end = j-1;
	    r->deleted = 0;
	    r->valid = 1;
	    ngaps++;
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

	word = (word << 2) | lookup[seq[i]];
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

	word = (word << 2) | lookup[seq[i]];
	
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
static Array suspect_joins(GapIO *io, tg_rec contig, int tw,
			   double filter_score,
			   double depth, HashTable *clip_hash) {
    int i, clen, cstart, cend;
    char *valid, *cp;
    char legal_chars[256];
    Array gaps;
    char *cons;
    contig_iterator *ci;
    rangec_t *r;
    int min_mq = 1;

    consensus_valid_range(io, contig, &cstart, &cend);
    clen = cend - cstart + 1;

    if (NULL == (valid = (char *)xcalloc(clen+1, 1)))
	return NULL;

    memset(legal_chars, 0, 256);
    for (cp = "ACGTacgt"; *cp; cp++)
	legal_chars[*cp] = 1;

    /* Compute consensus */
    if (NULL == (cons = (char *)xmalloc(clen+1))) {
	xfree(valid);
	return NULL;
    }
    calc_consensus(contig, cstart, cend, CON_SUM, cons, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    
    /*
     * Loop through readings computing new clip points after masking.
     * Mark the remainder as 'valid' to identify the invalid bits.
     */
    ci = contig_iter_new(io, contig, 1, CITER_FIRST | CITER_ISTART |
			 CITER_SMALL_BS, CITER_CSTART, CITER_CEND);
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = cache_search(io, GT_Seq, r->rec), *sorig = s;
	char *seq, *fseq;
	size_t len;
	int unique, first, last;

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

	/* Filter over represented words */
	//filter_common_words(seq, fseq, len, tw, depth, filter_score, '#', 0);

	/* Filter low mapping quality */
	if (r->mqual < min_mq) {
	    memset(fseq, '/', len);
	} else {
	    
	    /* Filter where the seq is low qual and disagrees with consen */
	    filter_consen_diffs(&seq[s->left-1], &fseq[s->left-1],
				s->right - s->left +1,
				&cons[r->start + s->left - cstart -1], 2);
	
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

#ifdef DEBUG
	printf("Seq %"PRIrec", first=%d, last=%d\n", s->rec, first, last);
#endif

	if (clip_hash && (first != r->start + MIN_OVERLAP+1 ||
			  last  != r->end - MIN_OVERLAP)) {
	    clip_pos_t *c = malloc(sizeof(*c));
	    HashData hd;

	    if (!c)
		return NULL;

	    c->left  = first;
	    c->right = last;
	    hd.p = c;
	    HashTableAdd(clip_hash, (char *)&r->rec, sizeof(r->rec), hd, 0);
	}

	if (last >= first) {
#ifdef DEBUG_SEQ
	    printf("Valid from %d to %d (cstart=%d)\n", first-1, last, cstart);
#endif	    
	    memset(&valid[first-1 -cstart], 1, last-first+1);
	}

	if (sorig != s)
	    free(s);

	xfree(fseq);
    }

    gaps = coverage2contig_regions(valid, cstart, clen);
    xfree(valid);
    xfree(cons);

    return gaps;
}

#define GOOD_SCORE +1
#define BAD_SCORE  -2
#define UNKNOWN_SCORE -1

/*
 * Determines whether a pair within the same contig should be
 * considered as consistent or not.
 * Returns 1 for yes;
 *         0 for no.
 */
static int consistent_pair(GapIO *io, rangec_t *r, HashTable *lt_h) {
    int isize_min, isize_max, isize, lib_type;

    /* Cache the lib_type for this library rec so we can
     * identify whether it is expected to be pointing
     * inwards or outwards.
     */
    if (r->library_rec) {
	HashItem *hi;
	library_t *lib;

	if ((hi = HashTableSearch(lt_h, (char *)&r->library_rec,
				  sizeof(tg_rec)))) {
	    lib = hi->data.p;
	} else {
	    HashData hd;

	    lib = cache_search(io, GT_Library, r->library_rec);
	    update_library_stats(io, lib->rec, 100,
				 NULL, NULL, NULL);
	    hd.p = lib;
	    HashTableAdd(lt_h, (char *)&r->library_rec,
			 sizeof(tg_rec),
			 hd, NULL);
	}

	lib_type = lib->lib_type;
	isize_max = isize_min = lib->insert_size[lib->lib_type];
	isize_min -= 2*lib->sd[lib->lib_type];
	isize_max += 2*lib->sd[lib->lib_type];
    } else {
	lib_type = LIB_T_INWARD;
	isize_min = 20;
	isize_max = 2000;
    }


    /* Check COMP1 vs COMP2 */
    switch (lib_type) {
    case LIB_T_INWARD:
	if (r->start < r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) == 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) != 0))
	    return 0;

	if (r->start > r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) != 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) == 0))
	    return 0;
	break;

    case LIB_T_OUTWARD:
	if (r->start > r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) == 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) != 0))
	    return 0;

	if (r->start < r->pair_start &&
	    !((r->flags & GRANGE_FLAG_COMP1) != 0 &&
	      (r->flags & GRANGE_FLAG_COMP2) == 0))
	    return 0;
	break;

    case LIB_T_SAME:
	if (((r->flags & GRANGE_FLAG_COMP1) == 0) !=
	    ((r->flags & GRANGE_FLAG_COMP2) == 0))
	    return 0;
	break;
    }


    /* Check insert size */
    isize = MAX(r->pair_start, r->pair_end) - MIN(r->start,r->end);
    isize = ABS(isize);
    if (!(isize >= isize_min && isize <= isize_max))
	return 0;

    return 1;
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
			 int good_score, int bad_score, int unknown_score) {
    int i;
    HashTable *lt_h;
    
    lt_h = HashTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (!lt_h)
	return;

    /* Now process gaps validating by read-pair */
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	contig_t *c;
	int score;
	int j, nr;
	rangec_t *r;
	HashTable *h;
	HashIter *iter;
	HashItem *hi;
	int num_good = 0;
	int num_bad = 0;
	int num_unknown = 0;

	if (gap->deleted)
	    continue;

	score = check_read_pairs(io, contig, gap);
	printf("Gap %d..%d,  validity %d\n",
	       gap->start, gap->end, score);

	if (!(c = cache_search(io, GT_Contig, contig)))
	    continue;

//FIXME: query library info
#define INS_SIZE 1000

	r = contig_seqs_in_range(io, &c,
				 gap->start-INS_SIZE, gap->end+INS_SIZE,
				 CSIR_PAIR, &nr);
	if (!r)
	    continue;

	h = HashTableCreate(1024, HASH_DYNAMIC_SIZE |
			          HASH_POOL_ITEMS |
			          HASH_NONVOLATILE_KEYS);
	if (!h)
	    continue;

	/* Pair up all seqs */
	for (j = 0; j < nr; j++) {
	    HashData hd;
	    int new;

	    if (!r[j].pair_rec)
		continue;

	    if ((hi = HashTableSearch(h, (char *)&r[j].pair_rec,
				      sizeof(r[j].pair_rec)))) {
		rangec_t *p = (rangec_t *)hi->data.p;
		sequence_get_range_pair_position(io, p);

		if (MIN(p->start, p->pair_start) <= gap->start &&
		    MAX(p->end, p->pair_end) >= gap->end) {
		    /* Spans the gap */
		    if (consistent_pair(io, p, lt_h))
			num_good++;
		    else
			num_bad++;
		}
		HashTableDel(h, hi, 0);
		continue;
	    }

	    hd.p = &r[j];
	    hi = HashTableAdd(h, (char *)&r[j].rec, sizeof(r[j].rec),
			      hd, &new);
	    if (!new) {
		fprintf(stderr, "Error: seq already in hash\n");
		num_unknown++;
		HashTableDel(h, hi, 0);
	    }
	}

	/* Count unmatched pairs remaining in the hash table */
	if (!(iter = HashTableIterCreate())) {
	    HashTableDestroy(h, 0);
	    free(r);
	    continue;
	}

	while (((hi = HashTableIterNext(h, iter)))) {
	    rangec_t *r = hi->data.p;
	    sequence_get_range_pair_position(io, r);

	    /* FIXME:
	     * Paired end libraries with singletons should be
	     * considered as an unknown states.
	     *
	     * Single ended libraries obviously are considered as OK,
	     * although they probably should have been filtered out
	     * before reachign this point.
	     */

	    if (!r->pair_contig) {
		num_unknown++;
		continue;
	    }

	    /* Make sure it spans the gap */
	    if (!(MIN(r->start, r->pair_start) <= gap->start &&
		  MAX(r->end, r->pair_end) >= gap->end))
		continue;

	    /* Check consistency */
	    if (contig == r->pair_contig) {
		if (consistent_pair(io, r, lt_h)) {
		    num_good++;
		    continue;
		}
	    }

	    /* Spanning => bad */
	    num_bad++;
	}

	HashTableIterDestroy(iter);

	score =  num_good * good_score
	        + num_bad * bad_score
	    + num_unknown * unknown_score;

	printf("GAP %d good, %d bad, %d unknown => score %d\n",
	       num_good, num_bad, num_unknown, score);

	/* FIXME: parameterise this via scoring function.
	 * Eg +1 for good, -2 for bad?
	 * Or just require a percentage?
	 */
	gap->valid = score >= 0 ? 1 : 0;

	HashTableDestroy(h, 0);

	if (gap->valid) {
	    free(r);
	} else {
	    gap->r  = r;
	    gap->nr = nr;
	}
    }

    HashTableDestroy(lt_h, 0);

    //uninit_template_checks(io, tarr);
}


/*
 * Merges gaps if they're close together, within min_distance apart.
 */
static void merge_gaps(Array gaps, int min_distance) {
    int i;
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

	    last_gap->end = gap->end;
	    gap->deleted = 1;
	} else {
	    last_gap = gap;
	}
    }
}


/*
 * Analyses gaps to work out where to break. We'll produce a left contig,
 * a right contig, and possibly multiple single-read contigs for readings
 * that are contained entirely within the gap section.
 *
 * Fills out 'ds' with a list of breaks to make. Each item is in itself
 * another list with the first two items being a contig number and position
 * within the contig to run Break Contig at with the remainder of the list
 * being an array of reading numbers to run disassembly on (these are
 * reads that fall within the gap and should be taken out to form their
 * own contig).
 */
static void break_gaps(GapIO *io, tg_rec contig, Array gaps,
		       HashTable *clip_hash, dstring_t *ds) {
    int i, j;
    tg_rec right_start = 0;

    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);

	if (gap->deleted)
	    continue;

	printf("Gap from %d to %d\n", gap->start, gap->end);

	if (gap->valid) {
	    printf("Skipping as marked as valid\n");
	    continue;
	}

	/* 
	 * Identify readings that after clipping are entirely contained
	 * within the region. These are candidates for disassembling
	 * completely
	 */
	dstring_appendf(ds, " {");
	
	printf("  New starting point for right contig = %d\n", gap->start+1);
	dstring_appendf(ds, "%"PRIrec" %d", contig, gap->start+1);

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
	    }
	}

	dstring_appendf(ds, "}");

	free (gap->r);
    }
}

void auto_break_single_contig(GapIO *io, tg_rec contig, int start, int end,
			      int tw, double filter_score, double depth,
			      int good_score, int bad_score,
			      int unknown_score, dstring_t *ds) {
    Array gaps;
    HashTable *clip_hash;

    printf("=== Checking contig %d ===\n", contig);

    clip_hash = HashTableCreate(1024,
				HASH_POOL_ITEMS);

    printf("  = Identifying suspect joins\n");
    gaps = suspect_joins(io, contig, tw, filter_score, depth, clip_hash);

    printf("  = Merging gaps\n");
    merge_gaps(gaps, MIN_OVERLAP*2);
    dump_gaps(gaps);

    printf("  = Confirming gaps\n");
    confirm_gaps(io, contig, gaps, good_score, bad_score, unknown_score);

    printf("  = Finding break points\n");
    break_gaps(io, contig, gaps, clip_hash, ds);

    HashTableDestroy(clip_hash, 1);
    ArrayDestroy(gaps);
    //xfree(clips);
}

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      double filter_score, int by_consensus) {
    int tw, i;
    double gc;
    int depth;

    dstring_t *ds = dstring_create(NULL);

    for (i = 0; i < argc; i++) {
	auto_break_single_contig(io, argv[i].contig, argv[i].start,
				 argv[i].end, tw, filter_score, depth,
				 GOOD_SCORE, BAD_SCORE, UNKNOWN_SCORE,
				 ds);
    }

    return ds;
}
