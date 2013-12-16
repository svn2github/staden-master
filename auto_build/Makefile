#
# Targets:
#
# dist	Make both staden and io_lib tarballs.
#       Manually remove .dist_io_lib and .dist_staden to force this to occur.
#
# dist_staden / dist_io_lib
#	Make just one of the above.
#
# all (default)
#	A full rebuild from scratch on all platforms.
#	This also forces a rebuild of the tarballs.
#
# macosx / linux-x86_64 / ...
#	Specific platforms only.
#	This only rebuilds from the existing tarball. If you've modified the
#	source in your local SVN checkout then manually remove
#	.dist_io_lib and/or .dist_staden first or use the all: target.
#

all: .all-recurse
.all-recurse: dist_clean
	$(MAKE) .all-hosts
.all-hosts: macosx linux-x86_64 linux-i386 windows-x86_64 windows-i386

dist_clean:
#	-rm .dist_staden
#	-rm .dist_io_lib

.dist_staden:
	cd $(HOME)/staden/trunk/build.seq3; $(MAKE) distsrc
	touch .dist_staden
dist_staden: .dist_staden

.dist_io_lib:
	cd $(HOME)/io_lib/build.seq3; $(MAKE) dist
	touch .dist_io_lib
dist_io_lib: .dist_io_lib

dist: dist_staden dist_io_lib

macosx: dist
	@echo Starting $@ X build. See make.$@ for output
	bash -x ./local.macosx > make.$@ 2>&1
	@echo $@ build complete


linux: linux-i386 linux-x86_64

linux-i386: dist
	@echo Starting $@ X build. See make.$@ for output
	sh -x ./local.lucid-i386 > make.$@ 2>&1
	@echo $@ build complete

linux-x86_64: dist
	@echo Starting $@ X build. See make.$@ for output
	sh -x ./local.precise-x86_64 > make.$@ 2>&1
	@echo $@ build complete

windows: windows-i386 windows-x86_64

windows-i386: dist
	@echo Starting $@ X build. See make.$@ for output
	bash -x ./local.$@ > make.$@ 2>&1
	@echo $@ build complete

windows-x86_64: dist
	@echo Starting $@ X build. See make.$@ for output
	bash -x ./local.$@ > make.$@ 2>&1
	@echo $@ build complete

clean:
	-rm staden-2.*-MacOSX.dmg
	-rm staden-2.*-linux-*.tar.gz
	-rm staden-2.*-windows.msi
	-rm staden-2.*-windows.exe
