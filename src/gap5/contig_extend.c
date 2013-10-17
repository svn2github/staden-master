/*
 * Assemblies produced by NGS tools often collapse data to an internal
 * representation during assembly, and then map reads back to the 
 * consensus sequences produced during assembly.
 *
 * At that stage we see oddities though, such as a consensus ending early
 * and dozens of reads with cutoff that agrees perfectly. For example,
 * with lowercase being cutoff:
 *
 * Consensus: AGCTATGAGGCGATACGATCCGTAA
 * Read1:     AGCTATGAGGCGATACGATCCGTAAgata
 * Read2:     AGCTATGAGGCGATACGATCCGTAAgatat
 * Read3:     AGCTATGTGGCGATACGATCCGTaagatct
 * Read4:     AGCTATGAGGCGATACGATCCGTAAgtatg
 * Read5:     AGCTATGAGGCGATACG-TCCGTAAgatatgcgtt
 * Read6:     AGCTATGAGGCGATACGATCCGTAAgatatgcgttaa
 * Read7:     AGCTATGAGGCGATACGATCCGTAAgatatgcgttaagatcg
 *
 * In this case we'd be happy to extend to gatatgcgtt say (3 fold depth),
 * but leaving the trailing end clipped still as it's not got deep enough
 * coverage.
 *
 * Also we'd extend reads 1,2,3,5,6,7 but not read4 due to an indel.
 * Note read3 has a different clip point than the others, but that is fine
 * as it aligns still.
 *
 * Conversely, sometimes due to a false join by an assembler followed by
 * mapping a lot of sequences get clipped in the middle of contigs. We
 * can detect and break these, and extend the clipped data. However we also
 * need an automatic trimmer first to cut off any outlier reads which extend
 * beyond the new contig ends, so that the auto-extended can then extend the
 * true sequence instead.
 */

#include <string.h>
#include <stdio.h>
#include <tg_gio.h>
#include <ctype.h>

#include "dna_utils.h"
#include "contig_extend.h"
#include "consensus.h"
#include "text_output.h"
#include "gap4_compat.h"  /* complement_contig */

#define CSZ 1024 /* consensus size - how far back in contig we'll look */
#define ESZ 1024 /* extension size - how far we can extend */

/*
 * Extends the right hand end of a single contig.
 *
 * Min_depth is the minimum depth for extension. If lower then even if the
 * data matches we'll not extend further.
 *
 * Match_score (+ve) and mismatch_score (-ve) are accumulated during
 * extension to ensure that we don't extend into junk mismatching DNA.
 */
static int contig_extend_single(GapIO *io, tg_rec crec, int dir, int min_depth,
				int match_score, int mismatch_score) {
    int end;
    rangec_t *r;
    int nr, i;
    contig_t *c;
    char cons[CSZ], new_cons[ESZ];
    int freqs[ESZ][5], depth[ESZ];
    double score, best_score;
    int best_pos, nseq;

    vmessage("Processing contig #%"PRIrec", %s end\n",
	     crec, dir ? "left" : "right");

    for (i = 0; i < ESZ; i++) {
	freqs[i][0] = freqs[i][1] = freqs[i][2] = freqs[i][3] =freqs[i][4] = 0;
	depth[i] = 0;
    }

    c = cache_search(io, GT_Contig, crec);
    if (NULL == c) return -1;
    cache_incr(io, c);

    if (consensus_valid_range(io, crec, NULL, &end) != 0) {
	cache_decr(io, c);
	return -1;
    }

    calculate_consensus_simple(io, crec, end-(CSZ-1), end, cons, NULL);

    /* Start */
    /* Not implemented for now: rev complement and go again! */

    /* End */
    r = contig_seqs_in_range(io, &c, end, end, 0, &nr);
    if (!r) {
	cache_decr(io, c);
	return -1;
    }

    for (i = 0; i < nr; i++) {
	seq_t *s = cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int cstart, cend;
	int j, k, slen;

	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	cstart = r[i].start + s->left-1;
	cend   = r[i].start + s->right-1;

	/* Does cutoff extend to contig end, if so does it match cons? */
	if (cend < end) {
	    int mis = 0, len = 0;
	    if (end - cend >= CSZ) {
		/*
		fprintf(stderr,"Skipping #%"PRIrec" due to length of cutoff\n",
			r[i].rec);
		*/
		if (sorig != s)
		    free(s);
		r[i].rec = 0; /* Mark for removal */
		continue;
	    }

	    for (k = s->right, j = cend+1; j <= end; j++, k++) {
		//printf("%d: %c %c\n", j, s->seq[k], cons[j-(end-(CSZ-1))]);
		if (s->seq[k] != cons[j-(end-(CSZ-1))])
		    mis++;
	    }
	    len = end - cend;
	    if (100*mis/len > 5) {
		/*
		fprintf(stderr, "Skipping #%"PRIrec" due to high disagreement "
			"with consensus.\n", r[i].rec);
		*/
		if (sorig != s)
		    free(s);
		r[i].rec = 0;
		continue;
	    }
	}

	/* So we got here, let's accumulate extension stats */
	slen = ABS(s->len);
	for (k = 0, j = end+1 - r[i].start; j < slen && k < ESZ; j++, k++) {
	    //printf("%d: %c\n", j + r[i].start, s->seq[j]);
	    if(s->seq[j] == 'N')
		continue;

	    freqs[k][dna_lookup[(uint8_t) s->seq[j]]]++;
	    depth[k]++;
	}

	if (sorig != s)
	    free(s);
    }

    score = best_score = 0;
    best_pos = 0;
    
    for (i = 0; i < ESZ; i++) {
	int call = 4, best = 0, j;
	double dd;

	if (depth[i] < min_depth)
	    break;

	for (j = 0; j < 5; j++) {
	    if (best < freqs[i][j]) {
		best = freqs[i][j];
		call = j;
	    }
	}
	new_cons[i] = "ACGT*"[call];

	dd = (double)depth[i];
	switch (call) {
	case 0:
	    score +=  freqs[i][0] / dd;
	    score -= (freqs[i][1]+freqs[i][2]+freqs[i][3]+freqs[i][4]) / dd;
	    break;
	case 1:
	    score +=  freqs[i][1] / dd;
	    score -= (freqs[i][0]+freqs[i][2]+freqs[i][3]+freqs[i][4]) / dd;
	    break;
	case 2:
	    score +=  freqs[i][2] / dd;
	    score -= (freqs[i][0]+freqs[i][1]+freqs[i][3]+freqs[i][4]) / dd;
	    break;
	case 3:
	    score +=  freqs[i][3] / dd;
	    score -= (freqs[i][0]+freqs[i][1]+freqs[i][2]+freqs[i][4]) / dd;
	    break;
	case 4:
	    score +=  freqs[i][4] / dd;
	    score -= (freqs[i][0]+freqs[i][1]+freqs[i][2]+freqs[i][3]) / dd;
	    break;
	}

	if (best_score <= score) {
	    best_score = score;
	    best_pos = i+1;
	}
	/*
	printf("%3d %3d\t%c\t%3d %3d %3d %3d %7.1f\n",
	       i, depth[i], "ACGT"[call],
	       freqs[i][0], freqs[i][1], freqs[i][2], freqs[i][3],
	       score);
	*/
    }
    //printf("Best score is %f at %d\n", best_score, best_pos);

    /* Extend */
    nseq = 0;
    if (best_pos > 0) {
	int furthest_left = end;

	for (i = 0; i < nr; i++) {
	    seq_t *s;
	    int r_pos;
	    int score;

	    if (r[i].rec == 0)
		continue;

	    s = cache_search(io, GT_Seq, r[i].rec);
	    s = cache_rw(io, s);

	    if (furthest_left > r[i].start)
		furthest_left = r[i].start;

	    /*
	     * end + best_pos is the furthest right we can go, but this
	     * specific read may not be justified in reaching that far
	     * if it has too many disagreements.
	     */
	    if ((s->len > 0) ^ r[i].comp) {
		int best_r = 0, j, k;
		int len = ABS(s->len);

		//printf(">%s\t", s->name);

		r_pos = s->right;
		score = 0;
		//for (k = s->right, j = 0; j < best_pos && k < len; j++, k++) {
		for (k = end - r[i].start + 1, j = 0; j < best_pos && k < len; j++, k++) {
		    if (new_cons[j] == toupper(s->seq[k])) {
			score += match_score;
			if (best_r <= score) {
			    best_r  = score;
			    r_pos = k+1;
			}
		    } else {
			score += mismatch_score;
		    }

		    //putchar(new_cons[j] == toupper(s->seq[k])
		    //	    ? toupper(s->seq[k])
		    //        : tolower(s->seq[k]));
		}
		//putchar('\n');

		if (s->right != r_pos) {
		    s->right  = r_pos;
		    nseq++;
		}
	    } else {
		int best_r = 0, j, k;

		//printf("<%s\t", s->name);

		r_pos = s->left-2;
		score = 0;
		//for (k = s->left-2, j = 0; j < best_pos && k >= 0; j++, k--) {
		for (k = r[i].end - end - 1, j = 0; j < best_pos && k >= 0; j++, k--) {
		    char b = complement_base(s->seq[k]);
		    if (new_cons[j] == b) {
			score += match_score;
			if (best_r <= score) {
			    best_r  = score;
			    r_pos = k-1;
			}
		    } else {
			score += mismatch_score;
		    }

		    //putchar(new_cons[j] == toupper(b)
		    //	    ? toupper(b)
		    //	    : tolower(b));
		}
		//putchar('\n');

		if (s->left != r_pos+2) {
		    s->left  = r_pos+2;
		    nseq++;
		}
	    }
	}

	vmessage("    Extended by %d, adjusting %d sequence clip%s\n",
		 best_pos, nseq, nseq == 1 ? "" : "s");

	bin_invalidate_consensus(io, crec, furthest_left, end + best_pos);
    } else {
	vmessage("    Unable to extend contig\n");
    }
    free(r);

    cache_decr(io, c);
    cache_flush(io);
    return 0;
}

/*
 * The main contig_extend interface.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_extend(GapIO *io, tg_rec *contig, int ncontigs, int min_depth,
		  int match_score, int mismatch_score) {
    int i, err = 0;

    for (i = 0; i < ncontigs; i++) {
	/* Left end */
	UpdateTextOutput();
	complement_contig(io, contig[i]);
	err |= contig_extend_single(io, contig[i], 1, min_depth,
				    match_score, mismatch_score);

	/* Right end */
	UpdateTextOutput();
	complement_contig(io, contig[i]);
	err |= contig_extend_single(io, contig[i], 0, min_depth,
				    match_score, mismatch_score);
    }

    return err ? -1 : 0;
}

/*-----------------------------------------------------------------------------
 * Contig trimmer.
 *
 * This is the opposite to Contig Extend, but works in harmony with it.
 *
 * If we have faked scaffold consensus assembled into a project then breaking
 * a contig in two can leave a long overhanging fake read that makes FIJ and
 * auto_join fail.  Similarly if an incorrect join was made by an assembler
 * off the back of a single chimeric read.
 *
 * In both cases the junction of the incorrect join, given it was produced
 * by an assembler consensus and then remapped back to, will probably
 * have a lot of clipped sequences for the same reason that existing
 * contig ends do. Hence: break, trim, extend.
 */


/* Callback from consensus_pileup */
static int trim_func(GapIO *io, tg_rec contig, int pos,
		     consensus_t *cons,
		     pileup_base_t *p, int depth, void *data) {
    int used_depth = 0;
    int *depth_pos = (int *)data;
    pileup_base_t *porig;

    /* Compute depth */
    for (porig = p; p; p = p->next) {
	if (p->base_index >= p->s->right)
	    break;
	else if (p->base_index >= p->s->left-1)
	    used_depth++;

	if (p->base_index == ABS(p->s->len)-1)
	    break;
    }

    depth_pos[1] = pos;

    /* Trim if sufficient depth or end of sequence */
    if (p || used_depth >= depth_pos[0]) {
	for (p = porig; p; p = p->next) {
	    if (p->comp) {
		seq_t *s = cache_search(io, GT_Seq, p->r.rec);
		cache_rw(io, s);
		if (ABS(s->len) - s->right < p->base_index) {
		    s->right = ABS(s->len) - p->base_index;
		    if (s->right < s->left)
			s->right = s->left;
		}
		sequence_range_length(io, &s);
	    } else {
		cache_rw(io, p->s);
		if (p->s->left < p->base_index+1) {
		    p->s->left = p->base_index+1;
		    if (p->s->left > p->s->right)
			p->s->left = p->s->right;
		}
		sequence_range_length(io, &p->s);
	    }
	}
	
	/* Aborts pileup calling function */
	return 1;
    }

    return 0;
}

static int contig_trim_single(GapIO *io, tg_rec crec, int dir, int min_depth) {
    contig_t *c;
    int depth_pos[2] = {min_depth, 0};

    if (!(c = cache_search(io, GT_Contig, crec)))
	return -1;
    if (c->nseqs < min_depth)
	return 0;

    consensus_pileup(io, crec, CITER_CSTART, CITER_CEND,
		     0/*CONS_ALL*/, trim_func, (void *)depth_pos);
    vmessage("  Trimmed %s end to pos %d (from end)\n",
	     dir ? "right" : "left", depth_pos[1]);
    
    return 0;
}

/*
 * The main contig_trim interface.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_trim(GapIO *io, tg_rec *contig, int ncontigs, int min_depth) {
    int i, err = 0, skip_clip = 0;

    // Hack to avoid clipping consensus annotations now (we'll do later)
    if (ncontigs < 0) {
	ncontigs *= -1;
	skip_clip = 1;
    }

    for (i = 0; i < ncontigs; i++) {
	/* Left end */
	vmessage("Contig =%"PRIrec" (%d/%d)\n", contig[i], i+1, ncontigs);
	err |= contig_trim_single(io, contig[i], 0, min_depth);
	UpdateTextOutput();
	complement_contig(io, contig[i]);

	/* Right end */
	err |= contig_trim_single(io, contig[i], 1, min_depth);
	UpdateTextOutput();
	complement_contig(io, contig[i]);

	if (!skip_clip) {
	    contig_visible_start(io, contig[i], CITER_CSTART);
	    contig_visible_end(io, contig[i], CITER_CEND);
	}
    }

    return err ? -1 : 0;
}

/*-----------------------------------------------------------------------------
 * Joint trim + extend interface
 */
int contig_trim_and_extend(GapIO *io, tg_rec *contig, int ncontigs,
			   int do_trim, int do_extend,
			   int trim_depth, int ext_depth,
			   int ext_match_score, int ext_mismatch_score) {
    int i, err = 0;

    for (i = 0; i < ncontigs; i++) {
	vmessage("\n");

	if (do_trim)
	    err |= contig_trim(io, &contig[i], -1, trim_depth);

	if (do_extend)
	    err |= contig_extend(io, &contig[i], 1, ext_depth,
				 ext_match_score, ext_mismatch_score);

	if (do_trim) {
	    contig_visible_start(io, contig[i], CITER_CSTART);
	    contig_visible_end(io, contig[i], CITER_CEND);
	}
	vmessage("\n");
    }

    return err ? -1 : 0;
}
