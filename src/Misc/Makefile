#
# Makefile for Miscellaneous routines
#
LIBS = misc
PROGS= $(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

CFLAGS += $(SHLIB_CFLAGS) $(MISC_DLL)

INCLUDES_S += $(CORBA_INC)

# pread/pwrite implementations for MinGW/Windows
ifeq ($(findstring MINGW,$(OS)),MINGW)
WIN_FUNCS=win_funcs.o
endif

#
# Objects
#
OBJS = \
	crash.o \
	date.o \
	filenames.o \
	files.o \
	find.o \
	shell.o \
	strings.o \
	usleep.o \
	error.o \
	xerror.o \
	xalloc.o \
	vlen.o \
	FtoC.o \
	getfile.o \
	bitmap.o \
	array.o \
	locks.o \
	array_arith.o \
	parse_db.o\
	dstring.o \
	string_alloc.o\
	strtol64.o\
	$(WIN_FUNCS)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(CORBA_OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(CORBA_OBJS) $(CORBA_DEP) $(BIOLIMS_LIB) $(MATH_LIB)

$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

#find.c:
#	ln -s ../io_lib/utils/find.c

DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)


# DO NOT DELETE THIS LINE -- make depend depends on it.

FtoC.o: $(PWD)/staden_config.h
array_arith.o: $(PWD)/staden_config.h
bitmap.o: $(PWD)/staden_config.h
crash.o: $(PWD)/staden_config.h
date.o: $(PWD)/staden_config.h
dstring.o: $(PWD)/staden_config.h
error.o: $(PWD)/staden_config.h
filenames.o: $(PWD)/staden_config.h
files.o: $(PWD)/staden_config.h
find.o: $(PWD)/staden_config.h
getfile.o: $(PWD)/staden_config.h
parse_db.o: $(PWD)/staden_config.h
shell.o: $(PWD)/staden_config.h
strings.o: $(PWD)/staden_config.h
strtol64.o: $(PWD)/staden_config.h
vlen.o: $(PWD)/staden_config.h
xalloc.o: $(PWD)/staden_config.h
xerror.o: $(PWD)/staden_config.h
