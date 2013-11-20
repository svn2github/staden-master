#
# Makefile for convert (Alpha version)
#

PROGS = staden_convert

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(GAPDB_EXT_INC) $(TK_INC) $(IOLIB_INC) $(TKUTILS_INC) $(P3_INC) $(TCL_INC) $(GAP4_INC) $(G_INC)

OBJS = \
	list.o \
	main.o \
	dapIO.o \
	dapDB.o \
	bapIO.o \
	bapDB.o \
	gapDB.o \
	process.o \
	flat_sd.o

DEPS =  $(GAP_LIB) \
	$(G_LIB) \
	$(P3_LIB) \
	$(MUT_LIB) \
	$(TKUTILS_LIB) \
	$(SEQUTILS_LIB) \
	$(MISC_LIB) \
	$(IOLIB_LIB)

staden_convert: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(LIBS) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(S)/convert.doc $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

# DO NOT DELETE THIS LINE -- make depend depends on it.

bapDB.o: $(PWD)/staden_config.h
bapIO.o: $(PWD)/staden_config.h
dapDB.o: $(PWD)/staden_config.h
dapIO.o: $(PWD)/staden_config.h
flat_sd.o: $(PWD)/staden_config.h
gapDB.o: $(PWD)/staden_config.h
main.o: $(PWD)/staden_config.h
process.o: $(PWD)/staden_config.h
