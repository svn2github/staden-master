#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc Cons_SelMask { parent f } {
    global gap5_defs

    keylset sm SELMASK [keylget gap5_defs CONSENSUS.SELMASK]
    set b1 [keylget sm SELMASK.BUTTON.1]
    set b2 [keylget sm SELMASK.BUTTON.2]
    set b3 [keylget sm SELMASK.BUTTON.3]

    frame $f -bd 2 -relief groove
    button $f.but \
	    -text "Select tags" \
	    -command "TagDialog CONSENSUS.TAGS \
			$parent[keylget gap5_defs SELECT_TAGS.WIN] {}"

    radiolist $f.rl \
	    -title [keylget sm SELMASK.NAME] \
	    -default [keylget sm SELMASK.VALUE]\
	    -buttons [format { \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    { %s -command { %s configure -state disabled}} } \
	    [list $b1] [list $f.but] CONSENSUS.TAGS \
	    [list $b2] [list $f.but] CONSENSUS.TAGS \
	    [list $b3] [list $f.but] ]
    pack $f.rl -side left
    pack $f.but -side right
}

proc NormalDialog { io } {
    global gap5_defs

    set f [keylget gap5_defs CONSENSUS.1.WIN]
    global $f.format.expt.Radio
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: normal"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap5_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

 
    ###########################################################################
    #select masking
    Cons_SelMask $f $f.sel_mask

    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap5_defs CONSENSUS.NORMAL.STRIP_PADS]

    yes_no $f.ambig \
	    -title "Output ambiguity codes" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap5_defs CONSENSUS.NORMAL.AMBIGUITY_CODES]

    ###########################################################################
    #Reading annotations
    radiolist $f.annos \
	-title "Output reading annotations  " \
	-default [keylget gap5_defs CONSENSUS.READ_ANNOTATIONS] \
	-orient horizontal \
	-buttons {{all} {{non-cutoff}} {none}}

    ###########################################################################
    #Reading notes
    yes_no $f.notes \
	    -title "Output reading notes" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap5_defs CONSENSUS.READ_NOTES]

    radiolist $f.template \
	-title "Name consensus by" \
	-default [keylget gap5_defs CONSENSUS.NAME_BY] \
	-orient horizontal \
	-buttons {{{left-most reading}} {{left-most template}}}

    ###########################################################################
    #format

    frame $f.format

    keylset ex EXPT [keylget gap5_defs CONSENSUS.EXPT]
    set $f.format.expt.Radio [keylget ex EXPT.VALUE]
    keylset fo FORMAT [keylget gap5_defs CONSENSUS.NORMAL.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
    set b3 [keylget fo FORMAT.BUTTON.3]
 
    frame $f.format.dummy
    radiolist $f.format.main \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s -command {radiolist_configure %s -state disabled;\
			   yes_no_configure %s -state disabled} } \
	    { %s -command {radiolist_configure %s -state disabled;\
			   yes_no_configure %s -state disabled} } \
	    { %s -command {radiolist_configure %s -state normal;\
			   yes_no_configure %s -state normal} } }\
	    [list $b1] [list $f.annos] [list $f.notes] \
	    [list $b2] [list $f.annos] [list $f.notes] \
	    [list $b3] [list $f.annos] [list $f.notes] ] 

    pack $f.format.main -fill x
    pack $f.format.dummy -fill x

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap5_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Normal_OK_Pressed $f $io $f.infile $f.id \
	    $f.sel_mask.rl $f.pads $f.ambig $f.notes $f.template $f.annos \
	    $f.format.main $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {Con-Normal}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.sel_mask -fill x
    pack $f.pads -fill x
    pack $f.ambig -fill x
    pack $f.template -fill x
    pack $f.format -fill x
    pack $f.annos -fill x
    pack $f.notes -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x

}

proc Normal_OK_Pressed {f io infile id sel_mask strippads ambig notes template annos format output} {
    global gap5_defs

    set gel_anno 0; #no gel annotations with expt file
    set truncate 1; #no gel annotations in hidden data with expt file
    set note_val 0
    set active_tags {}

    #single or list or file
    if {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set masking [radiolist_get $sel_mask]

    #if masking mode is 1 or 2 (mark or mask active tags)
    if {($masking == 1) || ($masking == 2)} {
	set active_tags [GetDefaultTags CONSENSUS.TAGS]
    }
    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]
    set ambig [yes_no_get $ambig]

    #expt format chosen
    if { $out_format == 3 } {
	set expt [radiolist_get $annos]
	if {$expt == 1} {
	    set gel_anno 1
	    set truncate 0
	} elseif {$expt == 2 } {
	    set gel_anno 1
	    set truncate 1
	}
	set note_val [yes_no_get $notes]
    }
    #set out_file [entrybox_get $output.entry]
    set out_file [getFname_in_name $output]

    if {$out_file == ""} return

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type normal\
	    -mask [lindex {mask mark none} [expr $masking-1]] \
	    -format $out_format\
	    -annotations $gel_anno \
            -truncate $truncate \
	    -notes $note_val \
	    -outfile $out_file \
	    -tag_types $active_tags \
	    -strip_pads $strip \
	    -hets $ambig \
	    -name_format [radiolist_get $template]
    ClearBusy
    destroy $f
}

proc ExtendedDialog { io } {
    global gap5_defs

    set f [keylget gap5_defs CONSENSUS.2.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: extended"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap5_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

    
    ###########################################################################
    #hidden data
    frame $f.hidden -bd 2 -relief groove

    HiddenParametersDialogInit $f.ops
    button $f.hidden.options \
	-text "Edit hidden data paramaters" \
	-command "HiddenParametersDialog $f $f.ops"

    pack $f.hidden.options -side top -anchor w


    ###########################################################################
    #select masking
    Cons_SelMask $f $f.sel_mask

    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap5_defs CONSENSUS.EXTENDED.STRIP_PADS]

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap5_defs CONSENSUS.EXTENDED.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s } {%s } }\
	    [list $b1] [list $b2] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap5_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Extended_OK_Pressed $f $io $f.infile $f.id \
	    $f.sel_mask.rl $f.pads $f.format $f.output $f.ops" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {Con-Extended}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.hidden -fill x
    pack $f.sel_mask -fill x
    pack $f.pads -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x
}

proc Extended_OK_Pressed {f io infile id sel_mask strippads format output \
			  hidden_ops} {
    global gap5_defs

    set gel_anno 0; #no gel annotations with expt file
    set truncate 0; #no gel annotations in hidden data with expt file
    set active_tags {}

    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    #Hidden data parameters
    upvar #0 $hidden_ops data
    set win_size $data($data(which_mode)_win_size)
    set max_dash $data(base_max_dash)
    set min_conf $data(conf_min_conf)
    set use_conf [lsearch {base conf} $data(which_mode)]

    set masking [radiolist_get $sel_mask]
    #if masking mode is 1 or 2 (mark or mask active tags)
    if {($masking == 1) || ($masking == 2)} {
	set active_tags [GetDefaultTags CONSENSUS.TAGS]
    }
    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]

    #set out_file [entrybox_get $output]
    set out_file [getFname_in_name $output]

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type extended\
	    -mask [lindex {mask mark none} [expr $masking-1]] \
	    -format $out_format\
	    -win_size $win_size \
	    -max_dashes $max_dash \
	    -min_conf $min_conf \
	    -use_conf $use_conf \
	    -outfile $out_file \
	    -tag_types $active_tags \
	    -strip_pads $strip
    ClearBusy
    destroy $f
}

proc UnfinishedDialog { io } {
    global gap5_defs

    set f [keylget gap5_defs CONSENSUS.3.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: unfinished"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap5_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

  
    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap5_defs CONSENSUS.UNFINISHED.STRIP_PADS]

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap5_defs CONSENSUS.UNFINISHED.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s } {%s } }\
	    [list $b1] [list $b2] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap5_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Unfinished_OK_Pressed $f $io $f.infile $f.id \
	    $f.pads $f.format $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {Con-Unfinished}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.pads -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x

}
proc Unfinished_OK_Pressed {f io infile id strippads format output} {
    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]

    #set out_file [entrybox_get $output]
    set out_file [getFname_in_name $output]

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type unfinished\
	    -format $out_format\
	    -outfile $out_file \
	    -strip_pads $strip
    ClearBusy
    destroy $f
}

proc QualityDialog { io } {
    global gap5_defs

    set f [keylget gap5_defs CONSENSUS.4.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: quality"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap5_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap5_defs CONSENSUS.QUALITY.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { { %s } } [list $b1] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap5_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Quality_OK_Pressed $f $io $f.id $f.infile $f.format\
	    $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {Con-Quality}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x


}

proc Quality_OK_Pressed {f io id infile format output} {
    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set out_format [radiolist_get $format]
    set out_file [getFname_in_name $output]

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type quality\
	    -format $out_format\
	    -outfile $out_file
    ClearBusy
    destroy $f
}

proc Strip_Pads {cons qual new_cons_var new_qual_var} {
    upvar $new_cons_var new_cons
    upvar $new_qual_var new_qual

    set new_cons ""
    set new_qual ""
    set pos 0
    foreach c [split $cons *] {
        set length [string length $c]
        append new_cons $c
        append new_qual [string range $qual $pos [expr {$pos+$length-1}]]
        incr pos [expr {$length+1}]
    }
}

# Applies map to $str over the (start,end) tuples in $pos.
# Assumes the positions are sorted by start coordinate.
#
# This has linear complexity, unlike simply modifying and returning the
# full str every time.
proc map_str {str map pos} {
    set last 0
    set new_str ""
    foreach {start end} $pos {
	if {$end < $last} continue
	if {$start > $last} {
	    append new_str [string range $str $last [expr {$start-1}]]
	} else {
	    set start $last
	}
	set sub [string range $str $start $end]
	append new_str [string map $map $sub]
	set last [expr {$end+1}]
    }
    append new_str [string range $str $last [string length $str]]
    return $new_str
}

# Masks (defi) or marks (acgt) consensus using the tags array
proc mask_consensus {io crec start end cons tags mask} {
    if {$tags != "*"} {
	foreach t $tags {
	    if {[scan $t "%c%c%c%c" a b c d] != 4} {
		puts stderr "Badly formatted tag type '$t'"
		continue
	    }
	    set id [expr {($a<<24)+($b<<16)+($c<<8)+$d}]
	    set filter($id) 1
	}
    }

    set c [$io get_contig $crec]

    if {$mask == "mask"} {
	set map {A d C e G f T i N n}
    } elseif {$mask == "mark"} {
	set map {A a C c G g T t N n}
    } elseif {$mask == "" || $mask == "none"} {
	return $cons
    } else {
	if {[string length $mask] == 1} {
	    set m $mask$mask$mask$mask$mask
	} else {
	    set m ${mask}nnnnn
	}
	set map [list \
		     A [string index $m 0] \
		     C [string index $m 1] \
		     G [string index $m 2] \
		     T [string index $m 3] \
		     N [string index $m 4]]
    }

    set pos_list {}
    foreach anno [$c anno_in_range $start $end] {	
	set t [lindex $anno 3]
	set t [format "%c%c%c%c" \
		   [expr {($t>>24)&0xff}] \
		   [expr {($t>>16)&0xff}] \
		   [expr {($t>> 8)&0xff}] \
		   [expr {($t>> 0)&0xff}]]
	if {$tags != "*" && ![info exists filter([lindex $anno 3])]} continue

	foreach {tag_st tag_en} $anno break
	if {$tag_st < $start} {set tag_st $start}
	if {$tag_en > $end}   {set tag_en $end}
	incr tag_st [expr {-($start)}]
	incr tag_en [expr {-($start)}]
	lappend pos_list $tag_st $tag_en
    }

    $c delete

    return [map_str $cons $map $pos_list]
}

proc get_consensus {args} {
    foreach {key value} $args {
	set opt($key) $value
    }
    set io $opt(-io)

    if {[catch {set fd [open $opt(-outfile) w]} err]} {
	tk_messageBox \
	    -message "Failed to write to $opt(-outfile)" \
	    -icon error \
	    -type ok
	return
    }
    fconfigure $fd -translation binary

    array set contig_done ""
    foreach contig $opt(-contigs) {
	foreach {id start end} $contig {
	    set crec [cname2crec $io $id]
	    if {[info exists contig_done($crec)]} break
	    set contig_done($crec) 1

	    set c [$io get_contig $crec]
	    set id [$c get_name]
	    if {$start == ""} {set start [$c get_visible_start]}
	    if {$end   == ""} {set end   [$c get_visible_end]}
	    $c delete

	    set cons [calc_consensus -io $io \
			  -contigs "{=$crec $start $end}" \
			  -hets $opt(-hets)]
	    if {[info exists opt(-mask)] && $opt(-mask) != "none" && \
		    $opt(-mask) != "" && [info exists opt(-tag_types)] && \
		    $opt(-tag_types) != ""} {
		set cons [mask_consensus $io $crec $start $end \
			      $cons $opt(-tag_types) $opt(-mask)]
	    }
	    switch $opt(-format) {
		1 {
		    # Fastq
		    set qual [calc_quality -io $io \
				  -contigs "{=$crec $start $end}" \
				  -hets $opt(-hets)]
		    if {$opt(-strip_pads)} {
			Strip_Pads $cons $qual new_cons new_qual
			set cons $new_cons; unset new_cons
			set qual $new_qual; unset new_qual
		    }
		    
		    puts $fd @$id
		    set c60 [reformat_sequence -fold 60 -str $cons]
		    puts $fd $c60
		    puts $fd +
		    set q60 [reformat_sequence \
				 -fold 60 \
				 -shift 33 \
				 -min 33 \
				 -max 126 \
				 -str $qual]
		    puts $fd $q60
		}

		2 {
		    # Fasta
		    puts $fd ">$id"
		    if {$opt(-strip_pads)} {
			set cons [join [split $cons *] {}]
		    }
		    set c60 [reformat_sequence -fold 60 -str $cons]
		    puts $fd $c60
		}

		3 {
		    # Experiment File
		    puts "Experiment file format not supported yet"
		}
	    }
	}
    }
    close $fd;
}
