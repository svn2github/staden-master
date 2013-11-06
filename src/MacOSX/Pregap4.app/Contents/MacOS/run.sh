#!/bin/sh

STADENROOT=`echo $0 | sed 's:[^/]*/Contents/MacOS/run.sh::'`
export STADENROOT

STADEN_PREPEND=1; export STADEN_PREPEND
. "$STADENROOT"/share/staden/staden.profile

export TCL_LIBRARY=$STADENROOT/share/tcltk/tcl8.5
export TK_LIBRARY=$STADENROOT/share/tcltk/tk8.5

# Find tclsh vs tclsh8.5
exec tclsh8.5 "$STADTCL/pregap4/pregap4.tcl"
