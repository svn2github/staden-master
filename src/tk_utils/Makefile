#
# Makefile for tk_utils routines
#

LIBS = tk_utils
PROGS=$(LIBS)
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_E += $(MISC_INC) $(TKUTILS_INC) $(IOLIB_INC) $(TK_INC) \
	      $(SEQUTILS_INC) $(BIOLIMS_INC) $(PNG_INC)

CFLAGS += $(SHLIB_CFLAGS) $(TK_UTILS_DLL)

#DEFINES += -DCHECK_LICENCE -DTCL_LOCAL_APPINIT=Stash_AppInit
#DEFINES += -DTCL_LOCAL_APPINIT=Stash_AppInit

# Windows requires tkWinX.o too to replace some of the missing X11 functions
ifeq ($(findstring MINGW,$(OS)),MINGW)
TKUTILS_EXTRAS=tkWinX.o
endif

#
# Objects
#
OBJS = \
	cli_arg.o \
	tclXkeylist.o \
	tclXutil.o \
	tcl_utils.o \
	tcl_debug.o \
	misc.o \
	init.o \
	text_output.o \
	tkRaster.o \
	tkRasterBuiltIn.o \
	sheet.o \
	tkSheet.o \
	tkSheet_common.o \
	trace_print.o \
	postscript.o \
	split.o \
	tkTrace.o \
	tkTraceComp.o \
	tkTraceIO.o \
	tkTraceDisp.o \
	capture.o \
	canvas_box.o \
	ruler_tick.o \
	restriction_enzyme_map.o \
	element_canvas.o\
	container.o\
	container_ruler.o\
	tcl_io_lib.o\
	$(TKUTILS_EXTRAS)

# Spin2 only objects, which currently isn't compiled and has compilation
# issues on some systems.
#
#	tkCanvGraph.o\
#	container_cmds.o\
#	tclCanvGraph.o\
#	matrix.o

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

TKU_LIBS = \
	$(IOLIB_LIB) \
	$(TK_LIB) \
	$(MISC_LIB) \
	$(PNG_LIB) \
	$(ZLIB_LIB) \
	$(SOCKET)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE) $(L)/.dir
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(TKU_LIBS)


$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)

# Stash interpreter executable
# Uses the main tcl interpreter, but overrides the appinit to use a local
# function.
ifneq ($(MACHINE),windows)
STASH_OBJS = \
	stash.o \
	tclAppInit.o

STASH_LIBS = \
	$(TKUTILS_LIB) \
	$(TK_LIB)

stash: $(STASH_OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMWIN) $(STASH_OBJS) $(STASH_LIBS) $(LIBSC)
endif

# Copy the tk_utils source
distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/*.tcl $(S)/tclIndex $(S)/Makefile $(DIRNAME)
	cp $(S)/*.bmp $(S)/help_config $(S)/tk_utilsrc $(DIRNAME)
	cp $(S)/stash $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)
	-mkdir $(INSTALLTCL)/tk_utils
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/tk_utils
	cp $(S)/stash $(INSTALLBIN)
ifneq ($(MACHINE),windows)
endif
	cp $(S)/2bars_*.bmp $(S)/help_config $(S)/tk_utilsrc $(INSTALLETC)

DEPEND_OBJ = $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

canvas_box.o: $(PWD)/staden_config.h
capture.o: $(PWD)/staden_config.h
cli_arg.o: $(PWD)/staden_config.h
container.o: $(PWD)/staden_config.h
container_ruler.o: $(PWD)/staden_config.h
element_canvas.o: $(PWD)/staden_config.h
init.o: $(PWD)/staden_config.h
misc.o: $(PWD)/staden_config.h
postscript.o: $(PWD)/staden_config.h
restriction_enzyme_map.o: $(PWD)/staden_config.h
ruler_tick.o: $(PWD)/staden_config.h
sheet.o: $(PWD)/staden_config.h
split.o: $(PWD)/staden_config.h
tcl_debug.o: $(PWD)/staden_config.h
tcl_utils.o: $(PWD)/staden_config.h
text_output.o: $(PWD)/staden_config.h
tkRaster.o: $(PWD)/staden_config.h
tkSheet.o: $(PWD)/staden_config.h
tkTrace.o: $(PWD)/staden_config.h
tkTraceComp.o: $(PWD)/staden_config.h
tkTraceDisp.o: $(PWD)/staden_config.h
tkTraceIO.o: $(PWD)/staden_config.h
trace_print.o: $(PWD)/staden_config.h
