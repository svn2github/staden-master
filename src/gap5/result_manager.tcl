#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Creates a results list window
# 
proc result_list_create {io} {
    global .results.Reg
    if {[xtoplevel [set t .results]] == ""} return
    wm title .results "Results manager"

    listbox $t.l -yscrollcommand "$t.s set" -width 37
    scrollbar $t.s -orient vertical -command "$t.l yview"
    frame $t.bb
    button $t.bb.q -text "Quit" -command "result_list_destroy $io"
    button $t.bb.h -text "Help" -command "show_help gap5 {Results}"

    bind $t.l <Motion> {
	%W selection clear 0 end
	%W selection set [%W index @%x,%y]
	%W selection anchor [%W index @%x,%y]
    }

    pack $t.bb -side bottom -fill x
    pack $t.bb.q -side left -expand 1
    pack $t.bb.h -side right -expand 1
    pack $t.s -side right -fill y
    pack $t.l -side left -expand 1 -fill both

    # Register ourselves for updates
    set .results.Reg [contig_register -io $io -contig 0 \
			  -command "result_list_callback $io" \
			  -flags [list REGISTER DEREGISTER QUERY_NAME]]

    wm protocol .results WM_DELETE_WINDOW "result_list_destroy $io"
}

#
# Callback to update the list on a new result (de)register event.
#
proc result_list_callback {io type id cdata args} {
    if {$type == "QUERY_NAME"} {
	return "Result Manager"
    } else {
	result_list_update $io
    }
}

#
# Tidy up on window destroy
#
proc result_list_destroy {io} {
    global .results.Reg
    contig_deregister -io $io -id [set .results.Reg]
    destroy .results
}

#
# Updates the results list window (always .results)
# 
proc result_list_update {io} {
    global gap5_defs
    
    set f [keylget gap5_defs CONTIG_SEL.WIN]

    # Contig selector Results and File->Save menu
    set csmenu  $f.menubar.results
    set csmenu2 $f.menubar.file.save
    if {[winfo exists $csmenu]} {
	set cse 1; # Contig Selector Exists
    } else {
	set cse 0
    }

    set t .results; # toplevel window name

    # Grab the registration list. Of format:
    # {contig id string} ?{contig id string}? ...
    set results [result_names -io $io]

    # Sort list by id to a format of:
    # {id name contig contig ...}
    set list ""
    set done "0"; #Ignore number 1
    foreach i $results {
	set id [lindex $i 1]
	if {[lsearch -exact $done $id] == -1} {
	    lappend done $id
	    set clist ""
       	    foreach j $results {
	        if {[lindex $j 1] == $id} {
		    lappend clist [lindex $j 0]
	        }
	    }

	    lappend list [list $id [lindex $i 2] $clist]
	}
    }

    # Clear and remove old binding
    set x [winfo exists $t]
    if $x {
        $t.l delete 0 end
        bind $t.l <<menu>> {}
    }

    if {$cse} {
        destroy $csmenu
        menu $csmenu
        destroy $csmenu2
        menu $csmenu2
    }

    set count 0
    # Add items to the list
    foreach i $list {
	set cl ""
	foreach cn [lindex $i 0] {
	    lappend cl "\#[lindex $cn 0]"
	}
        set ti [result_time -io $io -id [lindex $i 0]]
        set n "$ti : [lindex $i 1] ($cl)"
        if $x {$t.l insert end $n}

         if {$cse} {
             if {[result_is_2d -io $io -id [lindex $i 0]]} {
                 $csmenu add cascade -label $n -menu $csmenu.m$count
                 menu $csmenu.m$count -tearoff 0 
                 result_list_popup_single $io [lindex $i 0] \
     	            [reg_get_ops -io $io -id [lindex $i 0]] \
     	            $csmenu.m$count

		 $csmenu2 add command -label $n \
		     -command [list CSSavePlot $io $f [lindex $i 0]]

                 incr count
             }
 	}
    }

    # Reenable our binding with our new list.
    if $x {
        bind $t.l <<menu>> [concat \
    	    "result_list_popup $io {$list}" {%W [%W index @%x,%y] \
    	    [expr [winfo rootx %W]+%x-20] \
    	    [expr [winfo rooty %W]+%y-10]}]
    }
}

#
# Displays a menu from a given result list item
#
proc result_list_popup_single {io id ops m} {
    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $m add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
	        $m add command \
	            -label $i \
	            -command [list reg_invoke_op -io $io -id $id -option $count]
	    }
	    incr count
	}
    }
}

proc result_list_popup {io list lbox index x y} {
    if {[winfo exists $lbox.m]} {destroy $lbox.m}

    if {$list == ""} {return}

    set id [lindex [lindex $list $index] 0]
    set m [create_popup $lbox.m [$lbox get $index]]

    set ops [reg_get_ops -io $io -id $id]

    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $m add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
	        $m add command \
	            -label $i \
	            -command "destroy $m; \
			      [list reg_invoke_op -io $io -id $id -option $count]"
	    }
	    incr count
	}
    }

    $lbox selection clear 0 end
    $lbox selection set $index
    tk_popup $m $x $y
}


#
# The "list results" gap command.
#
# Collates information from the registration scheme and then pops up
# the window.
#
proc ListResults {io} {
    # Bring up the window
    result_list_create $io

    # And update it
    result_list_update $io
}

#consistency results menu
proc consistency_result_list_update {io c_win c_id} {
    global gap5_defs

    set cons_menu $c_win.menubar.results

    if {![winfo exists $cons_menu]} {
	return
    }

    # Grab the registration list. Of format:
    # {contig regnum id string} ?{contig regnum id string}? ...
    set results [result_names -io $io]

    # Sort list by id to a format of:
    # {id name {{contig reg} ?{contig reg}? ...}}
    set list ""
    set done "0"; #Ignore number 1
    foreach i $results {
	set id [lindex $i 2]
	if {[lsearch -exact $done $id] == -1} {
	    lappend done $id
	    set clist ""
       	    foreach j $results {
	        if {[lindex $j 2] == $id} {
		    lappend clist "[lindex $j 0] [lindex $j 1]"
	        }
	    }

	    lappend list "$id {[lindex $i 3]} {$clist}"
	}
    }

    destroy $cons_menu
    menu $cons_menu

    set count 0
    # Add items to the list
    foreach i $list {
        set ti [lindex [lindex [lindex $i 2] 0] 0]
        #set ti [result_time -io $io -contig $ti -id [lindex $i 0]]
        set ti [result_time -io $io -id [lindex $i 0]]
        set n "$ti : [lindex $i 1] (#[lindex $i 0])"

	if {[result_is_consistency -io $io -id [lindex $i 0] -cons_id $c_id]} {
	    $cons_menu add cascade -label $n -menu $cons_menu.m$count
	    menu $cons_menu.m$count -tearoff 0 
	    result_list_popup_single $io [lindex $i 0] \
    	            [reg_get_ops -io $io -id [lindex $i 0]] \
    	            $cons_menu.m$count
	    
	    incr count
        }
    }
}
