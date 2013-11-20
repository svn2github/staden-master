#
# Makefile for seq_utils routines
#

LIBS 	= seq_utils
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CFLAGS += $(SHLIB_CFLAGS) $(SEQ_UTILS_DLL)

INCLUDES_E += $(MISC_INC) $(TEXTUTILS_INC)

#
# Objects
#
OBJS = \
	align.o \
	align_lib_old.o \
	align_ss.o \
	align_ss2.o \
	align_sv.o \
	dna_utils.o \
	genetic_code.o \
	renz_utils.o \
	sequence_formats.o \
	scramble.o \
	base_comp.o \
	open_reading_frames.o\
	edge.o\
	search_utils.o\
	align_lib.o\
	read_matrix.o\
	filter_words.o


SU_LIBS = \
	$(TKUTILS_LIB) \
	$(MISC_LIB) \
	$(TK_LIB)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@


$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(SU_LIBS)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.gbl $(S)/Makefile $(DIRNAME)
	cp $(S)/align_lib_nuc_matrix $(S)/nuc_matrix $(DIRNAME)

S=$(SRCROOT)/seq_utils
install:
	cp $(PROGLIBS) $(INSTALLLIB)
	cp $(S)/align_lib_nuc_matrix $(S)/nuc_matrix $(INSTALLETC)

# DO NOT DELETE THIS LINE -- make depend depends on it.

align.o: $(PWD)/staden_config.h
align_lib.o: $(PWD)/staden_config.h
align_lib_old.o: $(PWD)/staden_config.h
align_ss.o: $(PWD)/staden_config.h
align_ss2.o: $(PWD)/staden_config.h
align_sv.o: $(PWD)/staden_config.h
base_comp.o: $(PWD)/staden_config.h
dna_utils.o: $(PWD)/staden_config.h
edge.o: $(PWD)/staden_config.h
genetic_code.o: $(PWD)/staden_config.h
open_reading_frames.o: $(PWD)/staden_config.h
read_matrix.o: $(PWD)/staden_config.h
renz_utils.o: $(PWD)/staden_config.h
sequence_formats.o: $(PWD)/staden_config.h
