/* ----------------------------------------------------------------------
 * General purpose fasta and fastq reading code.
 *
 * We provide next_fasta and next_fastq iterators on an opened FILE*.
 * Deallocate memory by passing in NULL to next_fast[aq] function.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#include "text_output.h"
#include "tg_gio.h"
#include "tg_index_common.h"
#include "zfio.h"
#include "fasta.h"

#define HDR_FASTA '>'
#define HDR_FASTQ '@'
#define DELIM_FASTA '>'
#define DELIM_FASTQ '+'

typedef struct {
    char *fn;
    unsigned long line;
    char *name;
    char *seq;
    char *qual;
    size_t  max_name_len;
    size_t  max_seq_len;
    size_t  max_qual_len;
    size_t  seq_len;
    char header, seq_delimiter;
} fastq_entry_t;


static inline int grow_char_string(char **str, size_t *max_len,
				   size_t desired_len) {
    char *new_str;
    size_t new_max = *max_len;
    
    while (new_max < desired_len) {
	new_max = new_max ? new_max * 2 : 128;
    }
    
    new_str = realloc(*str, new_max);
    if (!new_str) {
	verror(ERR_WARN, "grow_char_string", "Out of memory");
	return -1;
    }

    *str = new_str;
    *max_len = new_max;
    return 0;
}

/*
 * Reads a new fasta or fastq entry and returns it in *e.
 * For fasta, the quality string will be NULL
 *
 * Returns 0 on success
 *         1 on eof
 *        -1 on failure
 */

#define BLK_SIZE 8192
int fastaq_next(zfp *fp, fastq_entry_t *e) {
    char line[BLK_SIZE];
    size_t l, pos = 0;
    int discard = 0, more;

    /* Skip blank lines */
    line[BLK_SIZE - 1] = '*';
    do {
	if (NULL == zfgets(line, BLK_SIZE, fp)) return 1;
	e->line++;
    } while (*line == '\n');

    /* Read name */
    for (;;) {
	if (!discard) {
	    char *start = line, *cp;
	    if (0 == pos) { /* Check first char */
		if (e->header) {
		    if (*line != e->header) {
			verror(ERR_WARN, "fastaq_next",
			       "Error: sequence name does not start with '%c' "
			       "at %s line %lu.",
			       e->header, e->fn, e->line);
			return -1;
		    }
		} else {
		    if (*line != HDR_FASTA && *line != HDR_FASTQ) {
			verror(ERR_WARN, "fastaq_next",
			       "Error: sequence name does not start with "
			       "'%c' or '%c' at %s line %lu.",
			       HDR_FASTA, HDR_FASTQ, e->fn, e->line);
			return -1;
		    } else {
			e->header = *line;
			e->seq_delimiter = (*line == HDR_FASTA
					    ? DELIM_FASTA : DELIM_FASTQ);
		    }
		}
		start++;
	    }
	    cp = start;
	    while (*cp && !isspace(*cp)) cp++;
	    if (e->max_name_len < pos + cp - start + 1) {
		if (0 != grow_char_string(&e->name, &e->max_name_len,
					  pos + cp - start + 1)) {
		    return -1;
		}
	    }
	    memcpy(e->name + pos, start, cp - start);
	    pos += cp - start;
	    e->name[pos] = '\0';
	    if (isspace(*cp)) discard = 1;
	}
	if (line[BLK_SIZE - 1] == '*' || line[BLK_SIZE - 2] == '\n') break;
	line[BLK_SIZE - 1] = '*';
	if (NULL == zfgets(line, BLK_SIZE, fp)) {
	    verror(ERR_WARN, "fastaq_next",
		   "Error: Unexpected end-of-file while reading sequence name "
		   "at %s line %lu", e->fn, e->line);
	    return -1;
	}
    }
    if (*e->name == '\0') {
	verror(ERR_WARN, "fastaq_next",
	       "Error: Sequence entry with no name at %s line %lu",
	       e->fn, e->line);
	return -1;
    }

    /* Sequence */
    e->seq_len = 0;
    more = 0;
    while (!zfeof(fp) && (zfpeek(fp) != e->seq_delimiter || more)) {
	char *src = line, *dest;
	if (NULL == zfgets(line, BLK_SIZE, fp)) break;
	if (!more) e->line++;
	l = strlen(line);
	if (0 == l) {
	    verror(ERR_WARN, "fastaq_next",
		   "Error: Unexpected NUL byte at %s line %lu\n",
		   e->fn, e->line);
	    return -1;
	}
	if (e->max_seq_len < e->seq_len + l + 1) {
	    if (0 != grow_char_string(&e->seq, &e->max_seq_len,
				      e->seq_len + l + 1)) {
		return -1;
	    }
	}
	for (dest = e->seq + e->seq_len; *src; src++) {
	    if (!isspace(*src)) *dest++ = *src;
	}
	e->seq_len = dest - e->seq;
	more = line[l - 1] != '\n';
    }
    if (e->seq) e->seq[e->seq_len] = '\0';

    /* Done if fasta */
    if (e->header != HDR_FASTQ) return 0;

    /* + line: skip */
    line[BLK_SIZE - 1] = '*'; /* Sentinal */
    if (NULL == zfgets(line, BLK_SIZE, fp) || *line != DELIM_FASTQ) {
	verror(ERR_WARN, "fastaq_next",
	       "Error: Expected '%c' got '%c' reading fastq entry %.1000s",
	       DELIM_FASTQ, *line, e->name);
	return -1;
    }
    e->line++;
    while (line[BLK_SIZE - 1] != '*' && line[BLK_SIZE - 2] != '\n') {
	line[BLK_SIZE - 1] = '*';
	if (NULL == zfgets(line, BLK_SIZE, fp))
	    return -1; /* eof */
    }

    /* Read quality, no more than e->seq_len chars */
    pos = 0;
    more = 0;
    while (pos < e->seq_len && !zfeof(fp)) {
	char *src = line, *dest;
	if (NULL == zfgets(line, BLK_SIZE, fp))
	    break;
	if (!more) e->line++;
	more = (line[BLK_SIZE - 1] != '*' && line[BLK_SIZE - 2] != '\n');
	l = strlen(line);
	if (e->max_qual_len < pos + l + 1) {
	    if (0 != grow_char_string(&e->qual, &e->max_qual_len, pos + l + 1))
		return -1;
	}
	
	for (dest = e->qual + pos; *src; src++) {
	    if (!isspace(*src)) *dest++ = *src;
	}
	pos = dest - e->qual;
    }
    if (e->qual)
	e->qual[pos] = '\0';
    if (pos != e->seq_len) {
	verror(ERR_WARN, "fastaq_next",
	       "Error: differing number of sequence and quality "
	       "characters for entry '%.1000s' at %s line %lu",
	       e->name, e->fn, e->line);
	return -1;
    }

    return 0;
}

/* ----------------------------------------------------------------------
 * tg_index interface below.
 */
int parse_fasta_or_fastq(GapIO *io, char *fn, tg_args *a, int format) {
    int nseqs = 0, ret = 0;
    zfp *fp;
    struct stat sb;
    fastq_entry_t ent = { fn, 0, NULL, NULL, NULL, 0, 0, 0, 0, '\0', '\0' };
    //fastq_entry_t *(*next_seq)(zfp *fp);
    contig_t *c = NULL;
    int last_perd = 1;
    int res;

    vmessage("Loading %s...\n", fn);
    if (-1 == stat(fn, &sb) ||
	NULL == (fp = zfopen(fn, "r"))) {
	perror(fn);
	return -1;
    }

    /* Fetch sequences */
    while ((res = fastaq_next(fp, &ent)) == 0) {
	seq_t seq;
	static int dummy_qual_len = 0;
	static int8_t *dummy_qual = NULL;

	// printf("@%s\n%s\n+\n%s\n", ent->name, ent->seq, ent->qual);
	// printf("%d\tSeq %s len %d / %d\n",
	// nseqs, ent->name, (int)strlen(ent->seq), ent->seq_len);

	if (ent.seq_len <= 0) {
	    verror(ERR_WARN, "parse_fasta_or_fastq",
		   "Sequence named '%.1000s' appears to be blank", ent.name);
	    continue;
	}

	/* Create 1 read contig */
	create_new_contig(io, &c, ent.name, 0);
	
	seq.rec         = 0;
	seq.parent_rec  = 0;
	seq.parent_type = 0;

	seq.pos      = 1;
	seq.flags    = 0;
	seq.seq_tech = STECH_UNKNOWN;
	seq.format   = SEQ_FORMAT_CNF1;
	seq.left     = 1;
	seq.right    = ent.seq_len;

	seq.name_len = strlen(ent.name);
	seq.name     = strdup(ent.name);
	seq.template_name_len = seq.name_len;

	seq.seq      = ent.seq;
	seq.len      = ent.seq_len;

	seq.mapping_qual = 0;

	if (dummy_qual_len < ent.seq_len) {
	    dummy_qual_len = ent.seq_len;
	    dummy_qual = realloc(dummy_qual, dummy_qual_len);
	    if (!dummy_qual) {
		res = -1;
		break;
	    }
	}
	
	seq.conf = dummy_qual;
	assert(seq.conf);

	if (ent.qual) {
	    int i;
	    for (i = 0; i < ent.seq_len; i++) {
		int q = ent.qual[i] - '!';
		if (q < 0)
		    q = 0;
		if (q > 100)
		    q = 100;
		seq.conf[i] = q;
	    }
	} else {
	    assert(dummy_qual);
	    memset(dummy_qual, 0, dummy_qual_len);
	}

	seq.trace_name     = NULL;
	seq.trace_name_len = 0;
	seq.alignment      = NULL;
	seq.alignment_len  = 0;
	seq.sam_aux        = NULL;
	seq.aux_len        = 0;
	seq.anno           = NULL;

	save_range_sequence(io,
			    &seq,
			    0,    /* mapping qual */
			    NULL, /* pair array */
			    0,    /* is_pair */
			    NULL, /* template name */
			    c,    /* contig */
			    a,    /* args */
			    GRANGE_FLAG_TYPE_SINGLE,
			    NULL, /* library */
			    NULL  /* bin return rec */
			    );
			 
	if ((++nseqs & 0xff) == 0) {
	    int perc = 0;
	    off_t pos = zftello(fp);

	    perc = 100.0 * pos / sb.st_size;
	    if (perc > last_perd * 10) {
		vmessage("%c%d%%\n", (nseqs & 0xfff) ? '.' : '*', perc);
		last_perd = perc / 10 + 1;
	    } else {
		vmessage("%c", (nseqs & 0xfff) ? '.' : '*');
	    }
	    UpdateTextOutput();

	    if ((nseqs & 0xfff) == 0)
		cache_flush(io);
	}
    }
    vmessage("100%%\n");

    if (NULL != ent.name) free(ent.name);
    if (NULL != ent.seq)  free(ent.seq);
    if (NULL != ent.qual) free(ent.qual);

    if (res != 1) {
	ret = -1;
    }

    vmessage("Loaded %d sequences\n", nseqs);

    zfclose(fp);
    cache_flush(io);

    return ret;
}
