#
# Makefile for Mutation Library
#
LIBS = mutlib
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CXXFLAGS += $(SHLIB_CXXFLAGS) $(MISC_DLL)

INCLUDES_E += $(MISC_INC) $(MUT_INC) $(SEQUTILS_INC) $(IOLIB_INC)

DEFINES += -DMUTLIB_EXPORTS -DDYNMAT

#
# Objects
#
OBJS = \
	align.o \
	basecall.o \
	caller.o \
	caller_snr.o \
	caller_base.o \
	caller_level.o \
	dnatable.o \
	muttag.o \
	mutationtag.o \
	mutationtag_utils.o \
	pathutil.o \
	peakcall.o \
	stringlist.o \
	tagarray.o \
	trace.o \
	validate.o \
	read_matrix.o \
	sp_align_lib.o \
	sp_alignment.o \
	sp_hash_lib.o \
	sp_matrix.o \
	tracediff.o \
	tracediff_helper.o \
	tracediff_mutations.o \
	tracediff_parameters.o \
	tracediff_validate.o \
	tracealign.o \
	tracealign_cache.o \
	tracealign_helper.o \
	tracealign_interpolate.o \
	tracealign_preprocess.o \
	tracealign_validate.o \
	mutscan.o \
	mutscan_helper.o \
	mutscan_validate.o \
	mutscan_analysis.o \
	mutscan_preprocess.o \
	mutscan_parameters.o


MUT_DEP=\
	$(SEQUTILS_LIB) \
	$(TKUTILS_LIB) \
	$(IOLIB_LIB) \
	$(MISC_LIB) \
	$(MATH_LIB) \
	$(TK_LIB)

#ifeq ($(MACHINE),alpha)
#MUT_DEP += -Bstatic -lgcc
#endif

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(L)/.dir
	-cp $(MUTSRC)/$(LIBS).def $(DEF_FILE) 2>&-
	$(SHLIB_LDXX) $(SHLIB_LDXXFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(MUT_DEP)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.[ch]pp $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)

# DO NOT DELETE THIS LINE -- make depend depends on it

align.o: $(PWD)/staden_config.h
caller.o: $(PWD)/staden_config.h
caller_base.o: $(PWD)/staden_config.h
caller_level.o: $(PWD)/staden_config.h
caller_snr.o: $(PWD)/staden_config.h
mutationtag.o: $(PWD)/staden_config.h
mutationtag_utils.o: $(PWD)/staden_config.h
mutscan.o: $(PWD)/staden_config.h
mutscan_analysis.o: $(PWD)/staden_config.h
mutscan_helper.o: $(PWD)/staden_config.h
mutscan_parameters.o: $(PWD)/staden_config.h
mutscan_preprocess.o: $(PWD)/staden_config.h
mutscan_validate.o: $(PWD)/staden_config.h
muttag.o: $(PWD)/staden_config.h
read_matrix.o: $(PWD)/staden_config.h
sp_align_lib.o: $(PWD)/staden_config.h
sp_alignment.o: $(PWD)/staden_config.h
sp_hash_lib.o: $(PWD)/staden_config.h
sp_matrix.o: $(PWD)/staden_config.h
tagarray.o: $(PWD)/staden_config.h
trace.o: $(PWD)/staden_config.h
tracealign.o: $(PWD)/staden_config.h
tracealign_cache.o: $(PWD)/staden_config.h
tracealign_helper.o: $(PWD)/staden_config.h
tracealign_interpolate.o: $(PWD)/staden_config.h
tracealign_preprocess.o: $(PWD)/staden_config.h
tracealign_validate.o: $(PWD)/staden_config.h
tracediff.o: $(PWD)/staden_config.h
tracediff_helper.o: $(PWD)/staden_config.h
tracediff_mutations.o: $(PWD)/staden_config.h
tracediff_parameters.o: $(PWD)/staden_config.h
tracediff_validate.o: $(PWD)/staden_config.h
validate.o: $(PWD)/staden_config.h
