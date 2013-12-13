LIBS = gap5
# tg_view.bin deprecated
PROGS = tg_index.bin
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(TK_INC) \
	      $(IOLIB_INC) $(TKUTILS_INC) $(SEQUTILS_INC) \
	      $(P3_INC) $(CURSES_INC) $(ZLIB_INC) $(LZMA_INC)
	      
CFLAGS += $(SHLIB_CFLAGS)
FFLAGS += $(SHLIB_FFLAGS)

# For gprof profiling
#CFLAGS   += -pg
#CLDFLAGS += -pg

# For debugging g-alloc
#CFLAGS += -DVALGRIND
#CFLAGS += -D_FORTIFY_SOURCE=2 -fstack-protector

GAP5_LIB=-lgap5

#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.

# "g" library components, originally from Gap4's src/g dir
GOBJS = \
	g-files.o \
	g-db.o \
	g-struct.o \
	g-error.o \
	g-alloc.o \
	g-io.o \
	g-connect.o \
	g-request.o

# "tg" layer. This is low level I/O - caching, formatting, etc.
TG_IO = \
	tg_iface_g.o \
	b+tree2.o \
	hache_table.o \
	tg_gio.o \
	tg_cache.o \
	tg_contig.o \
	tg_sequence.o \
	tg_track.o \
	tg_bin.o \
	tg_utils.o \
	tg_register.o \
	tg_tracks.o \
	tg_anno.o \
	tg_library.o \
	tg_scaffold.o \
	tg_check.o \
        actf.o \
	zfio.o \
	gap_cli_arg.o \
	tg_tcl.o

CEDITOR = \
	editor_view.o \
	editor_join.o \
	editor_search.o \
	tkEdNames.o\
	tkEditor.o\
	tman_interface.o\
	tman_display.o\
	primlib.o\
	editor_oligo.o

CEDITOR_SCRAP=\
	contigEditor.o\
	tkEditor.o\
	tkEdUtils.o\
	edInterface.o\
	edUtils2.o \
	tagU1.o\
	undo.o\
	edExtend.o\
	edCommands.o\
	edMutations.o\
	tagEditor.o\
	searchUtils.o\
	tman_interface.o\
	tman_display.o\
	tman_cons.o\
	tman_diff.o\
	join.o\
	oligo.o

TG_INDEX=\
	tg_index_common.o\
	baf.o \
	maq.o \
	maqmap.o \
	sam_index.o \
	sam_pileup.o \
	ace.o \
	caf.o \
	fasta.o \
	afg.o \
	depad_seq_tree.o

GAP5=\
	newgap5_cmds.o \
	gap4_compat.o \
	gap_hash.o\
	hash_lib.o \
	qual.o \
	qualIO.o \
	gap_globals.o \
	tagdb.o\
	init.o \
	notedb.o\
	active_tags.o\
	list_proc.o\
	gap-error.o\
	stack_dump.o\
	tk-io-reg.o\
	tkAppInit.o \
	consen.o\
	contig_selector.o\
	cs-object.o\
	find_repeats.o\
	readpair.o\
	check_assembly.o\
	consensus.o\
	read_depth.o\
	gap_range.o \
	template_display.o \
	depth_track.o \
	template_draw.o \
	quality_plot.o\
	tag_plot.o \
	fij.o\
	do_fij.o\
	break_contig.o \
	dis_readings.o \
	export_contigs.o \
	find_oligo.o \
	import_gff.o \
	shuffle_pads.o \
	restriction_enzymes.o \
	gap_canvas_box.o \
	contig_extend.o \
	auto_break.o \
	$(GOBJS) \
	$(TG_IO) \
	$(CEDITOR) \
	$(TG_INDEX)

GAP5_LIBS=\
	$(IOLIB_LIB) \
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(ZLIB_LIB) \
	$(LZMA_LIB) \
	$(TK_LIB) \
	$(SEQUTILS_LIB) \
	$(P3_LIB) \
	$(MISC_LIB)

$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	-@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(GAP5) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(LINK_PATHFLAG)$(L) $(GAP5) $(GAP5_LIBS) $(SHLIB_DEP)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(GAP5)
	$(MKDEFL) $@ $(GAP5)


TG_IND_OBJ = \
	$(TG_IO) \
	tg_index.o \
#	mapass2.o

TGILIBS = \
	$(GAP5_LIB) \
	$(P3_LIB) \
	$(ZLIB_LIB) \
	$(LZMA_LIB) \
	$(TKUTILS_LIB) \
	$(SEQUTILS_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB)

tg_index.bin: $(TG_IND_OBJ) $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_IND_OBJ) $(TGILIBS) $(LIBSC)

TG_VIEW_OBJ = \
	$(TG_IO) \
	tg_view.o

TGVLIBS = \
	$(TGILIBS) \
	$(CURSES_LIB)

tg_view.bin: $(TG_VIEW_OBJ) $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TG_VIEW_OBJ) $(TGVLIBS) $(LIBSC)

DEPEND_OBJ = $(GAP5) $(TG_IND_OBJ) $(TG_VIEW_OBJ)

install:
	$(INSTALL) gap5 $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/gap5
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/gap5
	cp $(S)/../gap4/acd2tag.tcl $(INSTALLTCL)/gap5
	cp $(S)/gap5rc $(S)/gap5rc_menu_full $(INSTALLETC)
	cp $(PROGLIBS) $(INSTALLLIB)
	@#cp tg_view.bin $(INSTALLBIN)
	cp tg_index.bin $(INSTALLBIN)
	$(INSTALL) tg_index $(INSTALLBIN)
	@#$(INSTALL) tg_view $(INSTALLBIN)
	$(INSTALL) gap5_cmd $(INSTALLBIN)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/*.in $(S)/tclIndex $(S)/Makefile \
	   $(S)/gap5rc $(S)/gap5rc_menu_full $(DIRNAME)

# DO NOT DELETE THIS LINE -- make depend depends on it.

ace.o: $(PWD)/staden_config.h
actf.o: $(PWD)/staden_config.h
active_tags.o: $(PWD)/staden_config.h
afg.o: $(PWD)/staden_config.h
auto_break.o: $(PWD)/staden_config.h
b+tree2.o: $(PWD)/staden_config.h
baf.o: $(PWD)/staden_config.h
break_contig.o: $(PWD)/staden_config.h
caf.o: $(PWD)/staden_config.h
check_assembly.o: $(PWD)/staden_config.h
consen.o: $(PWD)/staden_config.h
consensus.o: $(PWD)/staden_config.h
contig_extend.o: $(PWD)/staden_config.h
contig_selector.o: $(PWD)/staden_config.h
cs-object.o: $(PWD)/staden_config.h
depad_seq_tree.o: $(PWD)/staden_config.h
depth_track.o: $(PWD)/staden_config.h
dis_readings.o: $(PWD)/staden_config.h
do_fij.o: $(PWD)/staden_config.h
editor_join.o: $(PWD)/staden_config.h
editor_oligo.o: $(PWD)/staden_config.h
editor_search.o: $(PWD)/staden_config.h
editor_view.o: $(PWD)/staden_config.h
export_contigs.o: $(PWD)/staden_config.h
fasta.o: $(PWD)/staden_config.h
fij.o: $(PWD)/staden_config.h
find_oligo.o: $(PWD)/staden_config.h
find_repeats.o: $(PWD)/staden_config.h
g-alloc.o: $(PWD)/staden_config.h
g-connect.o: $(PWD)/staden_config.h
g-db.o: $(PWD)/staden_config.h
g-files.o: $(PWD)/staden_config.h
g-io.o: $(PWD)/staden_config.h
g-request.o: $(PWD)/staden_config.h
g-struct.o: $(PWD)/staden_config.h
gap-error.o: $(PWD)/staden_config.h
gap4_compat.o: $(PWD)/staden_config.h
gap_canvas_box.o: $(PWD)/staden_config.h
gap_cli_arg.o: $(PWD)/staden_config.h
gap_globals.o: $(PWD)/staden_config.h
gap_hash.o: $(PWD)/staden_config.h
gap_range.o: $(PWD)/staden_config.h
hache_table.o: $(PWD)/staden_config.h
hash_lib.o: $(PWD)/staden_config.h
import_gff.o: $(PWD)/staden_config.h
init.o: $(PWD)/staden_config.h
list_proc.o: $(PWD)/staden_config.h
maq.o: $(PWD)/staden_config.h
newgap5_cmds.o: $(PWD)/staden_config.h
notedb.o: $(PWD)/staden_config.h
primlib.o: $(PWD)/staden_config.h
qual.o: $(PWD)/staden_config.h
qualIO.o: $(PWD)/staden_config.h
quality_plot.o: $(PWD)/staden_config.h
read_depth.o: $(PWD)/staden_config.h
readpair.o: $(PWD)/staden_config.h
restriction_enzymes.o: $(PWD)/staden_config.h
sam_index.o: $(PWD)/staden_config.h
sam_pileup.o: $(PWD)/staden_config.h
shuffle_pads.o: $(PWD)/staden_config.h
stack_dump.o: $(PWD)/staden_config.h
tag_plot.o: $(PWD)/staden_config.h
tagdb.o: $(PWD)/staden_config.h
template_display.o: $(PWD)/staden_config.h
tg_anno.o: $(PWD)/staden_config.h
tg_bin.o: $(PWD)/staden_config.h
tg_cache.o: $(PWD)/staden_config.h
tg_check.o: $(PWD)/staden_config.h
tg_contig.o: $(PWD)/staden_config.h
tg_gio.o: $(PWD)/staden_config.h
tg_iface_g.o: $(PWD)/staden_config.h
tg_index.o: $(PWD)/staden_config.h
tg_index_common.o: $(PWD)/staden_config.h
tg_library.o: $(PWD)/staden_config.h
tg_register.o: $(PWD)/staden_config.h
tg_scaffold.o: $(PWD)/staden_config.h
tg_sequence.o: $(PWD)/staden_config.h
tg_tcl.o: $(PWD)/staden_config.h
tg_track.o: $(PWD)/staden_config.h
tg_tracks.o: $(PWD)/staden_config.h
tg_utils.o: $(PWD)/staden_config.h
tg_view.o: $(PWD)/staden_config.h
tk-io-reg.o: $(PWD)/staden_config.h
tkAppInit.o: $(PWD)/staden_config.h
tkEdNames.o: $(PWD)/staden_config.h
tkEditor.o: $(PWD)/staden_config.h
tman_display.o: $(PWD)/staden_config.h
tman_interface.o: $(PWD)/staden_config.h
zfio.o: $(PWD)/staden_config.h
