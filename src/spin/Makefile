#
# Makefile for spin routines
#

LIBS = spin
PROGS = lib$(LIBS) 
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(TK_INC)  $(TKUTILS_INC) $(SEQUTILS_INC) $(IOLIB_INC)
CFLAGS += $(SHLIB_CFLAGS)

OBJS = \
	seq_reg.o \
	seq_reg_cmds.o \
	raster_cmds.o \
	seq_raster.o \
	sequtils.o \
	seq_results.o \
	seq_sendto.o \
	spin_globals.o \
	seq_plot_funcs.o\
	seqed.o\
	tkSeqed.o \
	tkSeqedUtils.o \
	tkSeqedNames.o \
	seqed_write.o \
	seqedInterface.o\
	seqed_translate.o\
	seqed_restriction_enzymes.o\
	seqed_search.o\
	sequence_pair_display.o\
	emboss_input_funcs.o \
	readpam.o\
	compare_spans.o \
	sip_hash.o \
	sip_similar_spans.o\
	sip_find_identity.o\
	sip_quick_scan.o\
	sip_align.o \
	sip_globals.o \
	sip_cmds.o \
	sip_results.o \
	rescan_matches.o \
	probs.o \
	sim.o \
	sip_sim.o\
	init.o\
	nip_globals.o \
	nip_cmds.o \
	nip_base_comp.o \
	nip_gene_search.o \
	nip_trna_search.o\
	codon_content.o \
	nip_stop_codon.o \
	trna_search.o \
	splice_search.o \
	nip_splice_search.o \
	nip_string_search.o \
	nip_wtmatrix_search.o \
	nip_restriction_enzymes.o\
	nip_canvas_box.o\
	dinuc_freqs.o

#
# Main dependency
#
lib$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

DEPS = $(TK_LIB) \
       $(IOUTILS_LIB) \
       $(TKUTILS_LIB) \
       $(SEQUTILS_LIB) \
       $(MISC_LIB)


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/tclIndex $(DIRNAME)
	cp $(S)/Makefile $(S)/spin.in $(S)/spin.bat $(DIRNAME)
	cp $(S)/*.wts $(S)/niprc $(S)/siprc $(S)/spinrc $(DIRNAME)

install:
	$(INSTALL) spin $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/spin
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/spin
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/*.wts $(S)/niprc $(S)/siprc $(S)/spinrc $(INSTALLETC)

# DO NOT DELETE THIS LINE -- make depend depends on it.

codon_content.o: $(PWD)/staden_config.h
compare_spans.o: $(PWD)/staden_config.h
emboss_input_funcs.o: $(PWD)/staden_config.h
nip_base_comp.o: $(PWD)/staden_config.h
nip_canvas_box.o: $(PWD)/staden_config.h
nip_cmds.o: $(PWD)/staden_config.h
nip_gene_search.o: $(PWD)/staden_config.h
nip_restriction_enzymes.o: $(PWD)/staden_config.h
nip_splice_search.o: $(PWD)/staden_config.h
nip_stop_codon.o: $(PWD)/staden_config.h
nip_string_search.o: $(PWD)/staden_config.h
nip_trna_search.o: $(PWD)/staden_config.h
nip_wtmatrix_search.o: $(PWD)/staden_config.h
probs.o: $(PWD)/staden_config.h
raster_cmds.o: $(PWD)/staden_config.h
readpam.o: $(PWD)/staden_config.h
rescan_matches.o: $(PWD)/staden_config.h
seq_plot_funcs.o: $(PWD)/staden_config.h
seq_raster.o: $(PWD)/staden_config.h
seq_reg.o: $(PWD)/staden_config.h
seq_reg_cmds.o: $(PWD)/staden_config.h
seq_results.o: $(PWD)/staden_config.h
seq_sendto.o: $(PWD)/staden_config.h
seqed.o: $(PWD)/staden_config.h
seqed_restriction_enzymes.o: $(PWD)/staden_config.h
seqed_search.o: $(PWD)/staden_config.h
sequence_pair_display.o: $(PWD)/staden_config.h
sequtils.o: $(PWD)/staden_config.h
sim.o: $(PWD)/staden_config.h
sip_align.o: $(PWD)/staden_config.h
sip_cmds.o: $(PWD)/staden_config.h
sip_find_identity.o: $(PWD)/staden_config.h
sip_globals.o: $(PWD)/staden_config.h
sip_hash.o: $(PWD)/staden_config.h
sip_quick_scan.o: $(PWD)/staden_config.h
sip_results.o: $(PWD)/staden_config.h
sip_sim.o: $(PWD)/staden_config.h
sip_similar_spans.o: $(PWD)/staden_config.h
splice_search.o: $(PWD)/staden_config.h
tkSeqed.o: $(PWD)/staden_config.h
tkSeqedNames.o: $(PWD)/staden_config.h
tkSeqedUtils.o: $(PWD)/staden_config.h
trna_search.o: $(PWD)/staden_config.h
