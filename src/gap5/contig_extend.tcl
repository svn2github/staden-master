proc ContigExtend {io} {
    global gap5_defs

    set f .contig_extend
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Contig Trim / Extend"

    #--- contig identifier widget
    contig_id $f.id -io $io

    lorf_in $f.infile [keylget gap5_defs EXTEND_CONTIGS.INFILE] \
        "{contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state normal}
        " -bd 2 -relief groove

    #--- Scoring extension parameters
    labelframe $f.extend -text "Contig Extend"

    xyn $f.extend.do_extend \
	-label "Extend contigs" \
	-orient horizontal \
	-default 1

    xentry $f.extend.min_depth \
	-label "Minimum extension depth" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MIN_DEPTH] \
	-type int

    xentry $f.extend.match_score \
	-label "Score for a match" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MATCH_SCORE] \
	-type int

    xentry $f.extend.mismatch_score \
	-label "Score for a mismatch" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MISMATCH_SCORE] \
	-type int

    pack $f.extend.do_extend $f.extend.min_depth $f.extend.match_score \
	$f.extend.mismatch_score -side top -fill both

    #--- Contig trimming parameters
    labelframe $f.trim -text "Contig Trim"
    
    xyn $f.trim.do_trim \
	-label "Trim contigs" \
	-orient horizontal \
	-default 1

    xentry $f.trim.trim_depth \
	-label "Minimum trim depth" \
	-default [keylget gap5_defs EXTEND_CONTIGS.TRIM_DEPTH] \
	-type int

    pack $f.trim.do_trim $f.trim.trim_depth -side top -fill both

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ContigExtend2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ContigExtend" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.id $f.trim $f.extend $f.ok \
	-side top -fill both
}

proc ContigExtend2 {io f} {
    if {[lorf_in_get $f.infile] == 4} {
	set gel_name [contig_id_gel $f.id]
	set lreg [contig_id_lreg $f.id]
	set rreg [contig_id_rreg $f.id]
	set list [list [list $gel_name $lreg $rreg]]
    } elseif {[lorf_in_get $f.infile] == 3} {
	set list [CreateAllContigList=Numbers $io]
    } else {
	set list [lorf_get_list $f.infile]
    }

    set do_extend      [$f.extend.do_extend get]
    set min_depth      [$f.extend.min_depth get]
    set match_score    [$f.extend.match_score get]
    set mismatch_score [$f.extend.mismatch_score get]

    set do_trim        [$f.trim.do_trim get]
    set trim_depth     [$f.trim.trim_depth get]

    if {!$do_trim && !$do_extend} {
	bell
	return
    }

    SetBusy
    log_call contig_extend \
	-io $io \
	-contigs $list \
	-extend $do_extend \
	-min_depth $min_depth \
	-match_score $match_score \
	-mismatch_score $mismatch_score \
	-trim $do_trim \
	-trim_depth $trim_depth
    ClearBusy

    destroy $f
}
