YEAR=$(shell date +%Y)

DIST=/tmp/staden-2.0.0b10-MacOSX

all: gap4 pregap4 gap5 trev spin 

# # -arch x86_64 (or both) seems to crash wrapper.c when run on 10.5.8. I've
# # no idea why. It also outputs warnings about implicit declaration of
# # RunApplicationEventLoop yet this is defined in a 32-bit environment.
# wrapper: wrapper.c
# 	echo year=$(YEAR)
# 	$(CC) -arch i386 -framework Carbon -O2 -o wrapper wrapper.c

gap4: VERS=$(shell sed -n 's/set GAP_VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/gap/gap.tcl)
gap4:
	@echo
	@echo === Gap4 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Gap4 $(VERS) Copyright 2011-$(YEAR)</string/%" Gap4.app/Contents/Info.plist > _tmp.g4 && mv _tmp.g4 Gap4.app/Contents/Info.plist

pregap4: VERS=$(shell sed -n 's/.*Pregap4 version \([^$$]*\).*/\1/p' $(DIST)/share/staden/tcl/pregap4/gui.tcl)
pregap4:
	@echo
	@echo === Pregap4 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Pregap4 $(VERS) Copyright 2011-$(YEAR)</string/%" Pregap4.app/Contents/Info.plist > _tmp.p4 && mv _tmp.p4 Pregap4.app/Contents/Info.plist

gap5: VERS=$(shell sed -n 's/set GAP_VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/gap5/gap5.tcl)
gap5:
	@echo
	@echo === Gap5 version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Gap5 $(VERS) Copyright 2011-$(YEAR)</string/%" Gap5.app/Contents/Info.plist > _tmp.g5 && mv _tmp.g5 Gap5.app/Contents/Info.plist

trev: VERS=$(shell sed -n 's/set VERSION "\([^$$]*\).*"/\1/p' $(DIST)/share/staden/tcl/trev/trev.tcl)
trev:
	@echo
	@echo === Trev version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Trev $(VERS) Copyright 2011-$(YEAR)</string/%" Trev.app/Contents/Info.plist > _tmp.t && mv _tmp.t Trev.app/Contents/Info.plist

spin: VERS=$(shell grep SPIN_VERSION $(DIST)/share/staden/tcl/spin/spin.tcl |head -1|awk '{print $$NF}')
spin: 
	@echo
	@echo === Spin version $(VERS) ===
	sed "s%.* Copyright 2011.*%  <string>Spin $(VERS) Copyright 2011-$(YEAR)</string/%" Spin.app/Contents/Info.plist > _tmp.s && mv _tmp.s Spin.app/Contents/Info.plist
