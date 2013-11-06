YEAR=$(shell date +%Y)

DIST=/tmp/staden-2.0.0b10-MacOSX

all: gap4 pregap4 gap5 trev spin 

wrapper:
	echo year=$(YEAR)
	$(CC) -arch i386 -arch x86_64 -O2 -o wrapper wrapper.c

gap4: VERS=$(shell sed -n 's/set GAP_VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/gap/gap.tcl)
gap4: wrapper
	@echo
	@echo === Gap4 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Gap4 $(VERS) Copyright 2011-$(YEAR)</string/%" Gap4.app/Contents/Info.plist > _tmp && mv _tmp Gap4.app/Contents/Info.plist
	cp wrapper Gap4.app/Contents/MacOS/gap4

pregap4: VERS=$(shell sed -n 's/.*Pregap4 version \([^$$]*\).*/\1/p' $(DIST)/share/staden/tcl/pregap4/gui.tcl)
pregap4: wrapper
	@echo
	@echo === Pregap4 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Pregap4 $(VERS) Copyright 2011-$(YEAR)</string/%" Pregap4.app/Contents/Info.plist > _tmp && mv _tmp Pregap4.app/Contents/Info.plist
	cp wrapper Pregap4.app/Contents/MacOS/pregap4

gap5: VERS=$(shell sed -n 's/set GAP_VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/gap5/gap5.tcl)
gap5: wrapper
	@echo
	@echo === Gap5 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Gap5 $(VERS) Copyright 2011-$(YEAR)</string/%" Gap5.app/Contents/Info.plist > _tmp && mv _tmp Gap5.app/Contents/Info.plist
	cp wrapper Gap5.app/Contents/MacOS/gap5

trev: VERS=$(shell sed -n 's/set VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/trev/trev.tcl)
trev: wrapper
	@echo
	@echo === Trev version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Trev $(VERS) Copyright 2011-$(YEAR)</string/%" Trev.app/Contents/Info.plist > _tmp && mv _tmp Trev.app/Contents/Info.plist
	cp wrapper Trev.app/Contents/MacOS/trev

spin: VERS=$(shell grep SPIN_VERSION $(DIST)/share/staden/tcl/spin/spin.tcl |head -1|awk '{print $$NF}')
spin: wrapper 
	@echo
	@echo === Spin version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Spin $(VERS) Copyright 2011-$(YEAR)</string/%" Spin.app/Contents/Info.plist > _tmp && mv _tmp Spin.app/Contents/Info.plist
	cp wrapper Spin.app/Contents/MacOS/spin
