proc exprand { median } {
    # Returns exponentially-distributed random variables
    if { $median <= 0 } { return 0 }
    return [expr { int(-log(rand()) / (log(2) / $median)) }]
}

proc make_tmp { prefix name_out } {
    upvar $name_out name
    set p [pid]
    for { set c 0 } { $c < 1000000 } { incr c } {
	set name "/tmp/${prefix}_$p.$c"
	if {![catch {open $name {CREAT EXCL RDWR} 0600} fd]} {
	    return $fd
	}
	if {[lindex $::errorCode 1] ne "EEXIST" \
		|| [lindex $::errorCode 0] ne "POSIX"} {
	    break
	}
    }
    return -code error -errorcode $::errorCode
}

proc gen_dna { len } {
    set bases { A C G T }
    set dna ""
    # puts stderr "Generating random sequence of $len bp"
    for { set pos 0 } { $pos < $len } { incr pos } {
	append dna [lindex $bases [expr { int(rand() * 4) } ]]
    }
    return $dna
}

proc write_sam_reads { samfd dna reads } {
    foreach rd [lsort -integer -index 1 $reads] {
	foreach { rnum pos l dirn paired clip_left clip_right tag } $rd break
	# puts stderr "$rnum $pos $l $dirn $paired $clip_left $clip_right"
	set lc [gen_dna $clip_left]
	set rc [gen_dna $clip_right]
	set seq [string range $dna $pos [expr { $pos + $l - 1 }]]
	set rdna "$lc$seq$rc"
	set lcigar [expr { $clip_left > 0 ? "${clip_left}S" : "" }]
	set rcigar [expr { $clip_right > 0 ? "${clip_right}S" : "" }]
	set cigar [format "%s%dM%s" $lcigar $l $rcigar]
	if {$tag} {
	    set el [expr {$l + $clip_left + $clip_right}]
	    set tleft [expr {int(rand() * ($el - 1))}]
	    set tright [expr { $tleft + int(rand() * ($el - $tleft))}]
	    set annostr "\tPT:Z:[expr {$tleft + 1}];[expr {$tright + 1}];+;COMM;Note=[expr {$tleft + $pos}]..[expr {$tright + $pos}]"
	} else {
	    set annostr ""
	}
	puts $samfd [format \
			 "r%06d\t%d\tc%04d\t%d\t255\t%s\t*\t0\t0\t%s\t*%s" \
			 $rnum \
			 [ expr { ($paired \
				       ? ($dirn ? 0x63 : 0x93) \
				       : ($dirn ? 0x10 : 0x00))} ] 0 \
			 [ expr { $pos + 1 } ] $cigar $rdna $annostr]
    }
}

proc gen_reads_pattern { dna samfd len spc dirn lclip rclip {tag 0} } {
    set rnum 1
    set clen [string length $dna]
    set reads {}
    for {set pos 0} {$pos < $clen} {incr pos $spc} {
	set l [expr {$pos + $len < $clen ? $len : $clen - $pos }]
	lappend reads [list $rnum $pos $l $dirn 0 $lclip $rclip $tag]
	incr rnum
    }
    write_sam_reads $samfd $dna $reads
}

proc gen_reads { dna samfd rnum_out {tag 0} } {
    upvar $rnum_out rnum
    set rnum 1
    set median_tlen [exprand 500]
    set median_len  [exprand 100]
    puts stderr "Median template length = [expr {$median_tlen + 100}]"
    puts stderr "Median read length     = [expr {$median_len + 50}]"
    
    set tpos [expr {-[exprand $median_tlen]}]
    set rnum 0
    set clen [string length $dna]
    set reads {}
    while { $tpos < $clen } {
	set tl [expr { [exprand $median_tlen] + 100 } ]
	if { $tpos + $tl < 0 } continue
	set trev  [expr { rand() > 0.5 ? 1 : 0 }]
	set fl [expr { [exprand $median_tlen] + 50 } ]
	if { $fl > $tl } { set fl $tl }
	if { $tpos + $fl < 0 } { set fl 0 }
	set rl [expr { [exprand $median_tlen] + 50 } ]
	if { $rl > $tl } { set rl $tl }
	if { $tpos + $tl - $rl >= $clen } { set rl 0 }
	set fpos $tpos
	set rpos [expr {$tpos + $tl - $rl}]
	if { $fl && $fpos < 0 } {
	    incr fl $fpos
	    set fpos 0
	}
	if { $rl && $rpos < 0 } {
	    incr rl $rpos
	    set rpos 0
	}
	if { $fl && $fpos + $fl > $clen } { set fl [expr { $clen - $fpos }]}
	if { $rl && $rpos + $rl > $clen } { set rl [expr { $clen - $rpos }]}
	if {$fl > 0} {
	    lappend reads [list $rnum $fpos $fl 0 [expr { $rl > 0 }] \
			       [exprand 0.4] [exprand 0.4] $tag]
	}
	if {$rl > 0} {
	    lappend reads [list $rnum $rpos $rl 1 [expr { $fl > 0 }] \
			       [exprand 0.4] [exprand 0.4] $tag]
	}
	while { 1 } {
	    set dist [ exprand 20 ]
	    if { $dist < $fl + 1 } { break }
	}
	set tpos [expr {$tpos + $dist}]
	incr rnum
    }
    write_sam_reads $samfd $dna $reads
}

proc add_contig_tags { io cnum clen } {
    for { set p 1 } { $p < $clen } { incr p [expr { int(rand() * 30) - 5 }]} {
	if { $p < 1 } { set p 1 }
	set tlen [expr { int(rand() * 10) + 1 }]
	if { $p + $tlen > $clen } {
	    set tlen [expr { $clen - $p }]
	}
	set rec [$io new_anno_ele 17 $cnum $p [expr { $p + $tlen} ]]
	set t [$io get_anno_ele $rec]
	$t set_comment "$p..[expr {$p + $tlen}]"
	$t delete
    }
}

proc basic_add_tags { tags_out io clen crec srec taglist } {
    upvar $tags_out tags
    
    set s [$io get_sequence $srec]
    set st [$s get_position]
    $s delete

    foreach t $taglist {
	foreach { rec start end } $t break
	puts stderr "rec = $rec start = $start end = $end"
	if {$rec == $crec} {
	    set arec [$io new_anno_ele 17 $rec $start $end]
	    puts stderr "Tag $arec (17 $rec $start $end)"
	} else {
	    set arec [$io new_anno_ele 18 $rec \
			  [expr {$start - $st}] [expr {$end - $st}]]
	    puts stderr "Tag $arec (18 $rec $start $end)"
	}
	set t [$io get_anno_ele $arec]
	$t set_comment "$rec $start..$end"
	$t delete
	set tags($arec) [list $rec $start $end]
    }
}

proc check_contig { io crec { expected "" } } {
    set c [$io get_contig $crec]
    set checked [$c check]
    if {[lindex $checked 0] != 0} {
	puts stderr "Errors found ($checked), quitting."
	return 1
    }
    set vs [$c get_visible_start]
    set ve [$c get_visible_end]
    $c delete
    if { $expected ne "" } {
	# puts stderr "calc_consensus -io $io -contigs {=$crec $vs $ve}"
	set cons [calc_consensus -io $io -contigs "{=$crec $vs $ve}"]
	if { $cons ne $expected } {
	    puts stderr "Consensus mismatch, quitting."
	    set el [string length $expected]
	    if {$el > 100} {
		puts stderr "Expected [string range $expected 0 4] .. [expr {$el - 10}] bases .. [string range $expected [expr {$el - 5}] $el]]"
	    } else {
		puts stderr "Expected $expected"
	    }
	    set cl [string length $cons]
	    if {$cl > 100} {
		puts stderr "Got      [string range $cons 0 4] .. [expr {$cl - 10}] bases .. [string range $cons [expr {$cl - 5}] $cl]]"
	    } else {
		puts stderr "Got      $cons"
	    }
	    return 1
	}
    }
    return 0
}

proc check_sequence { io srec expected comp } {
    set s [$io get_sequence $srec]
    set dna [$s get_seq]
    $s delete
    if { $comp } {
	set dna [ string_reverse [ string map { A T C G G C T A } $dna ] ]
    }
    if { $dna ne $expected } {
	puts stderr "Sequence mismatch, quitting."
	set el [string length $expected]
	if {$el > 100} {
	    puts stderr "Expected [string range $expected 0 4] .. [expr {$el - 10}] bases .. [string range $expected [expr {$el - 5}] $el]]"
	} else {
	    puts stderr "Expected $expected"
	}
	set dl [string length $dna]
	if {$dl > 100} {
	    puts stderr "Got      [string range $dna 0 4] .. [expr {$dl - 10}] bases .. [string range $dna [expr {$dl - 5}] $dl]]"
	} else {
	    puts stderr "Got      $dna"
	}
	return 1
    }
    return 0
}

proc check_tags { io crec clen etags } {
    upvar $etags tags
    set contig [$io get_contig $crec]
    set annos [$contig anno_in_range 1 $clen]
    $contig delete

    array set seen {}
    foreach anno $annos {
	foreach {astart aend arec} $anno break
	set rec [lindex $anno 8]
	set seen($arec) 1
	if { $tags($arec) ne {}} {
	    foreach { trec tstart tend } $tags($arec) break
	    if { $trec != $rec || $tstart != $astart || $tend != $aend } {
		puts stderr "Tag mismatch for $arec: expected ($trec, $tstart, $tend) got ($rec $astart $aend)"
		return 1
	    }
	} else {
	    puts stderr "Tag mismatch for $arec: expected (deleted) got ($rec $astart $aend)"
	    return 1
	}
    }
    foreach arec [array names tags] {
	if { $tags($arec) eq {} } continue
	if { ![info exists seen($arec)] } {
	    foreach { trec tstart tend } $tags($arec) break
	    puts stderr "Tag mismatch for $arec: expected ($trec, $tstart, $tend) got (absent)"
	    return 1
	}
    }
    return 0
}

proc sim_ctrl_delete { undo exp pos {tags {}}} {
    upvar $undo un
    upvar $exp expected
    upvar $pos p
    if {$tags ne ""} {
	upvar $tags t
    }

    if {[string length $expected] <= 1} return

    set u [list $p [string range $expected $p $p]]
    set expected [string replace $expected $p $p]
    
    if {$tags ne ""} {
	set cp [expr {$p + 1}]
	foreach { arec loc } [array get t] {
	    if { $loc eq "" } continue
	    foreach { rec start end } $loc break
	    if { $end >= $cp } {
		lappend u $arec $rec $start $end
		if { $start == $cp && $end == $cp } {
		    set t($arec) {}
		} else {
		    set t($arec) \
			[list $rec [expr {$start > $cp ? $start - 1 : $start}] \
			     [expr {$end - 1}]]
		}
	    }
	}
    }
    lappend un $u
}

proc sim_ctrl_backspace { undo exp pos {tags {}}} {
    upvar $undo u
    upvar $exp expected
    upvar $pos p
    if {$tags ne ""} {
	upvar $tags t
    }

    if { $p > 0 } {
	incr p -1
	if {$tags ne ""} {
	    sim_ctrl_delete u expected p t
	} else {
	    sim_ctrl_delete u expected p
	}
    }

}


proc sim_undo { undo uip exp {tags {}} } {
    upvar $exp expected
    upvar $uip ui
    if {$tags ne ""} {
	upvar $tags t
    }

    if { $ui >= 0 } {
	foreach {p s} [lindex $undo $ui] break
	# puts stderr "p = $p s = $s"
	set expected "[string range $expected 0 [expr {$p - 1}]]$s[string range $expected $p end]"
	if { $tags ne "" } {
	    foreach { arec rec start end } [lrange [lindex $undo $ui] 2 end] {
		set t($arec) [list $rec $start $end]
	    }
	}
	incr ui -1
    } 
}

proc basic_tests { base_io ed clen crec srec dna tags_array pos comp } {
    global $ed
    upvar \#0 $ed opt
    upvar $tags_array tags
    set w $opt(curr_editor)
    set $opt(Cutoffs) [lindex [$w configure -display_cutoffs 1] 4]
    update

    if {$pos < $clen} {
	set u {}
	set expected $dna
	array set etags [array get tags]
	set strpos [expr {$pos - 1}]
	for { set i $pos } { $i < $clen } { incr i } {
	    $w set_cursor 17 $crec $pos
	    editor_delete_base $w [list 17 $crec $pos] 1 1 1
	    update idletasks
	    # after 100
	    sim_ctrl_delete u expected strpos etags
	    if {[check_contig [$w io] $crec $expected] != 0
		|| [check_sequence [$w io] $srec $expected $comp] != 0
		|| [check_tags [$w io] $crec $clen etags] != 0} {
		editor_save $ed
		$base_io flush
		return 1
	    }
	}
	editor_save $ed
	$base_io flush
	if {[check_contig $base_io $crec $expected] != 0} {
	    return 1
	}

	set ui [expr {[llength $u] - 1}]
	for {set i $pos} {$i < $clen} {incr i} {
	    editor_undo $ed
	    update idletasks
	    # after 100
	    sim_undo $u ui expected etags
	    if {[check_contig [$w io] $crec $expected] != 0
		|| [check_sequence [$w io] $srec $expected $comp] != 0
		|| [check_tags [$w io] $crec $clen etags] != 0} {
		editor_save $ed
		$base_io flush
		return 1
	    }
	}
	editor_save $ed
	$base_io flush
	if {[check_contig $base_io $crec $expected] != 0} {
	    return 1
	}
    }

    if {$pos > 1} {
	set u {}
	set expected $dna
	set strpos [expr {$pos}]
	array set etags [array get tags]
	for { set i [expr {$pos + 1}] } { $i > 0 } { incr i -1 } {
	    $w set_cursor 17 $crec $i
	    editor_delete_base $w [list 17 $crec $i] 1 0 1
	    update idletasks
	    sim_ctrl_backspace u expected strpos etags
	    if {[check_contig [$w io] $crec $expected] != 0
		|| [check_sequence [$w io] $srec $expected $comp] != 0
		|| [check_tags [$w io] $crec $clen etags] != 0} {
		editor_save $ed
		$base_io flush
		return 1
	    }
	}
	editor_save $ed
	$base_io flush
	if {[check_contig $base_io $crec $expected] != 0} {
	    return 1
	}
	
	set ui [expr {[llength $u] - 1}]
	for { set i $pos } { $i > 0 } { incr i -1 } {
	    editor_undo $ed
	    update idletasks
	    sim_undo $u ui expected etags
	    if {[check_contig [$w io] $crec $expected] != 0
		|| [check_sequence [$w io] $srec $expected $comp] != 0
		|| [check_tags [$w io] $crec $clen etags] != 0} {
		editor_save $ed
		$base_io flush
		return 1
	    }
	}
	editor_save $ed
	$base_io flush
	if {[check_contig $base_io $crec $expected] != 0} {
	    return 1
	}
    }
    return 0
}

proc del_undo_test { base_io ed clen crec dna } {
    global $ed
    upvar \#0 $ed opt
    set w $opt(curr_editor)
    set $opt(Cutoffs) [lindex [$w configure -display_cutoffs 1] 4]
    update

    set u {}
    set expected $dna
    set strpos 0
    for {set i 0} {$i <= $clen} {incr i} {
	$w set_cursor 17 $crec 1
	editor_delete_base $w [list 17 $crec 1] 1 1 1
	sim_ctrl_delete u expected strpos
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
    }
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }

    set ui [expr {[llength $u] - 1}]
    for {set i 0} {$i <= $clen} {incr i} {
	editor_undo $ed
	sim_undo $u ui expected
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
    }
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }

    set u {}
    set expected $dna
    set strpos [string length $expected]
    for {set i [expr { $clen + 1 }]} {$i >= 0} {incr i -1} {
	$w set_cursor 17 $crec $i
	editor_delete_base $w [list 17 $crec $i] 1 0 1
	sim_ctrl_backspace u expected strpos
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
    }
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }

    set ui [expr {[llength $u] - 1}]
    for {set i 0} {$i <= $clen} {incr i} {
	editor_undo $ed
	sim_undo $u ui expected
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
    }
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }

    set mid [expr { int($clen / 2 + 1) }]
    set u {}
    set expected $dna
    set strpos [expr {$mid - 1}]
    for {set i 0} {$i <= $mid} {incr i} {
	set edpos [expr {$mid - $i}]
	$w set_cursor 17 $crec $edpos
	editor_delete_base $w [list 17 $crec $edpos] 1 1 1
	sim_ctrl_delete u expected strpos
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
	editor_delete_base $w [list 17 $crec $edpos] 1 0 1
	sim_ctrl_backspace u expected strpos
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	# puts stderr $expected
	update idletasks
    }
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }

    set ui [expr {[llength $u] - 1}]
    for {set i 0} {$i <= $mid} {incr i} {
	editor_undo $ed
	sim_undo $u ui expected
	# puts stderr $expected
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
	editor_undo $ed
	sim_undo $u ui expected
	if {[check_contig [$w io] $crec $expected] != 0} {
	    editor_save $ed
	    $base_io flush
	    return 1
	}
	update idletasks
    }

    puts stderr "Saving."
    editor_save $ed
    $base_io flush
    if {[check_contig $base_io $crec $expected] != 0} {
	return 1
    }
    return 0
}

proc close_editor { w } {
    global $w

    set detach ""
    foreach ed [set ${w}(all_editors)] {
	global $ed
	set id [set ${ed}(reg)]
	contig_deregister -io [set ${w}(io_base)] -id $id
	set id [set ${ed}(reg_all)]
	contig_deregister -io [set ${w}(io_base)] -id $id
	lappend detach [$ed contig_rec]
    }

    destroy $w

    foreach crec $detach {
	io_detach $crec
    }
}

proc get_curr_editor {} {
    return [lsearch -inline -regexp [winfo children .] {^\.e\d+$}]
}

proc do_basic_test { clen comp } {
    puts stderr "Generating test data..."
    set dna [gen_dna $clen]
    set sam_name ""
    set samfd [ make_tmp "deltest" sam_name ]
    puts $samfd [ format "@SQ\tSN:c%04d\tLN:%d" 0 $clen ]
    set reads [list [list 1 0 $clen 0 0 0 0 0]]
    write_sam_reads $samfd $dna $reads
    close $samfd

    if {[catch { exec tg_index -o $sam_name $sam_name >@stdout 2>@stderr } msg]} {
	puts stderr "Error running tg_index -o $sam_name $sam_name : $msg"
	return 1
    }
    if {[catch \
     	     {set io [g5::open_database -name $sam_name -access rw]} \
     	     err]} {
     	puts stderr "Couldn't open database '$sam_name': $err"
     	return 1
    }
    $io debug_level 1

    set cnum [cname2crec $io "c0000"]

    if { $comp } {
	complement_contig -io $io -contigs =$cnum
	set dna [ string_reverse [ string map { A T C G G C T A } $dna ] ]
    }

    set snum [$io seq_name2rec "r000001"]
    set half [expr { int($clen / 2) }]
    set quart [expr { int($clen / 4) }]
    set tq [expr { int($clen * 3 / 4) }]
    array set tags {}
    basic_add_tags tags $io $clen $cnum $snum \
	[list [list $cnum 1 2] \
	     [list $snum 1 2] \
	     [list $cnum $quart $quart] [list $snum $quart $quart] \
	     [list $cnum [expr { $half - 5 }] [expr { $half + 5 }]] \
	     [list $snum [expr { $half - 5 }] [expr { $half + 5 }]] \
	     [list $cnum $tq $tq] [list $snum $tq $tq] \
	     [list $cnum [expr { $clen - 2 }] $clen] \
	     [list $snum [expr { $clen - 2 }] $clen]]
    edit_contig -io $io -contig $cnum
    set ed [get_curr_editor]
    if {[basic_tests $io $ed $clen $cnum $snum $dna tags 1 $comp] != 0} {
	tkwait window $ed
	exit 1
    }
    if {[basic_tests $io $ed $clen $cnum $snum $dna tags [expr { int($clen / 2) }] $comp] != 0 } {
	tkwait window $ed
	exit 1
    }
    if {[basic_tests $io $ed $clen $cnum $snum $dna tags $clen $comp] != 0} {
	tkwait window $ed
	exit 1
    }
    close_editor $ed
    $io close
    unset io
    return 0
}

proc do_cutoff_test { clean } {
    puts stderr "Generating test data..."
    set clen 32
    set dna [gen_dna $clen]
    set sam_name ""
    set samfd [ make_tmp "deltest" sam_name ]
    puts $samfd [ format "@SQ\tSN:c%04d\tLN:%d" 0 $clen ]
    set reads [list [list 1 0 [expr {$clen / 2 - 2}] 0 0 0 [expr {$clen  / 2 + 2} ] 0]\
		   [list 2 [expr {$clen / 2}] [expr {$clen / 2}] 0 0 [expr {$clen / 2}] 0 0]]
    for { set p 0 } { $p < $clen } { incr p } {
	lappend reads [list [expr { $p + 3 } ] $p 1 0 0 0 0 0]
    }
    write_sam_reads $samfd $dna $reads
    close $samfd

    if {[catch { exec tg_index -z 1 -o $sam_name $sam_name >@stdout 2>@stderr } msg]} {
	puts stderr "Error running tg_index -o $sam_name $sam_name : $msg"
	return 1
    }
    if {[catch \
     	     {set io [g5::open_database -name $sam_name -access rw]} \
     	     err]} {
     	puts stderr "Couldn't open database '$sam_name': $err"
     	return 1
    }
    $io debug_level 1

    set cname "c0000"
    set cnum [cname2crec $io $cname]

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.before"
    edit_contig -io $io -contig $cnum
    set ed [get_curr_editor]

    if {[del_undo_test $io $ed $clen $cnum $dna] != 0} {
	tkwait window $ed
	exit 1
    }

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.after"

    if [catch {exec cmp "$sam_name.before" "$sam_name.after" } res] {
	if {[lindex $::errorCode 0] eq "CHILDSTATUS"} {
	    puts stderr "cmp failed."
	    tkwait window $ed
	    return 1
	} else {
	    puts stderr "Error running cmp : $res"
	}
    }

    close_editor $ed

    complement_contig -io $io -contigs =$cnum
    set dna [ string_reverse [ string map { A T C G G C T A } $dna ] ]

    edit_contig -io $io -contig $cnum
    set ed [get_curr_editor]

    if {[del_undo_test $io $ed $clen $cnum $dna] != 0} {
	tkwait window $ed
	exit 1
    }

    close_editor $ed

    $io close
    unset io
    if { $clean } {
	puts stderr "Tidying up..."
	file delete "$sam_name"
	file delete "${sam_name}.g5d"
	file delete "${sam_name}.g5x"
	file delete "${sam_name}.log"
    }
    return 0
}

proc do_test { clean mode clen args } {
    puts stderr "Generating test data..."
    # set clen [ expr { 100 + [exprand 1000] } ]
    # set clen 10
    set dna [gen_dna $clen]
    set sam_name ""
    set samfd [ make_tmp "deltest" sam_name ]

    puts $samfd [ format "@SQ\tSN:c%04d\tLN:%d" 0 $clen ]
    # gen_reads $dna $samfd rnum 1
    # gen_reads_pattern $dna $samfd 100 99 0 1 1 1
    if { $mode eq "pattern" } {
	eval [list gen_reads_pattern $dna $samfd] $args
    } else {
	gen_reads $dna $samfd rnum 1
    }
    close $samfd

    if {[catch { exec tg_index -o $sam_name $sam_name >@stdout 2>@stderr } msg]} {
	puts stderr "Error running tg_index -o $sam_name $sam_name : $msg"
	return 1
    }
    #global io
    if {[catch \
     	     {set io [g5::open_database -name $sam_name -access rw]} \
     	     err]} {
     	puts stderr "Couldn't open database '$sam_name': $err"
     	return 1
    }
    #if {[catch {DB_Load $sam_name} err]} {
#	puts stderr "Couldn't open database '$sam_name': $err"
#	return 1
#    }
    $io debug_level 1

    set cname "c0000"
    set cnum [cname2crec $io $cname]

    add_contig_tags $io $cnum $clen

    if { $mode eq "pattern" } {
	return 0
    }

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.before"

    edit_contig -io $io -contig $cnum
    set ed [get_curr_editor]

    if {[del_undo_test $io $ed $clen $cnum $dna] != 0} {
    	tkwait window $ed
    	exit 1
    }

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.after"

    if [catch {exec cmp "$sam_name.before" "$sam_name.after" } res] {
	if {[lindex $::errorCode 0] eq "CHILDSTATUS"} {
	    puts stderr "cmp failed."
	    tkwait window $ed
	    return 1
	} else {
	    puts stderr "Error running cmp : $res"
	}
    }

    close_editor $ed

    complement_contig -io $io -contigs =$cnum
    set dna [ string_reverse [ string map { A T C G G C T A } $dna ] ]

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.cbefore"

    edit_contig -io $io -contig $cnum
    set ed [get_curr_editor]

    if {[del_undo_test $io $ed $clen $cnum $dna] != 0} {
	tkwait window $ed
	exit 1
    }

    export_contigs -io $io -contigs "{=$cnum 1 $clen}" -format sam \
	-outfile "$sam_name.cafter"

    if [catch {exec cmp "$sam_name.cbefore" "$sam_name.cafter" } res] {
	if {[lindex $::errorCode 0] eq "CHILDSTATUS"} {
	    puts stderr "cmp failed."
	    tkwait window $ed
	    return 1
	} else {
	    puts stderr "Error running cmp : $res"
	}
    }

    close_editor $ed

    $io close
    unset io
    if { $clean } {
	puts stderr "Tidying up..."
	file delete "$sam_name"
	file delete "$sam_name.g5d"
	file delete "$sam_name.g5x"
	file delete "$sam_name.log"
	file delete "$sam_name.before"
	file delete "$sam_name.after"
	file delete "$sam_name.cbefore"
	file delete "$sam_name.cafter"
    }

    return 0
}

proc run_child { seed clean } {
    load_package tk_utils
    tk_utils_init
    load_package gap5
    
    InitLists
    wm withdraw .

    puts stderr "Using seed $seed"
    expr srand($seed)

    if {[do_basic_test 80 0]} {
	puts stderr "Basic test failed."
	exit 1
    }
    if {[do_basic_test 80 1]} {
	puts stderr "Complement basic test failed."
	exit 1
    }

    if {[do_cutoff_test $clean]} {
    	puts stderr "Cutoff test failed."
    	exit 1
    }

    set patterns [list [list 100 1 1 0 1 1 1] \
		      [ list 100 1 1 1 1 1 1] \
		      [ list 8300 200 190 0 1 1 1] \
		     ]

    foreach pat $patterns {
	if { [eval [list do_test $clean pattern] $pat] != 0 } {
	    puts stderr "Pattern test $pat failed"
	    exit 1
	}
    }

    if { [do_test $clean sim [ expr { 100 + [exprand 1000] } ] ] } {
	puts stderr "Simulated data test failed"
	exit 1
    }
}

#-----------------------------------------------------------------------------
# MAIN
# Startup code
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}

# source $env(STADTCL)/gap5/gap5.tcl
package require Tk

set seed [expr {int(rand() * 2147483647)}]
if {[lindex $argv 0] != ""} {
    set seed [lindex $argv 0]
}
set clean 0

run_child $seed $clean

exit 0

set case1 {@RG	ID:XX	SM:unknown	LB:XX
@SQ	SN:c0000	LN:24
r000001	0	c0000	1	255	1S21M1S	*	0	0	CAGTAGGAAACAAGTCAGAGGGG	!!!!!!!!!!!!!!!!!!!!!!!	RG:Z:XX	PT:Z:19;22;+;COMM;Note=a
*	768	c0000	1	255	11M	*	0	0	*	*CT:Z:.;COMM;Note=b
*	768	c0000	21	255	1M	*	0	0	*	*CT:Z:.;COMM;Note=c
r000003	0	c0000	22	255	1S1M1S	*	0	0	CGC	!!!	RG:Z:XX
*	768	c0000	22	255	1M	*	0	0	*	*CT:Z:.;COMM;Note=d
}
