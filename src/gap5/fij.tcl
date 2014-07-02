#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FIJDialog { w io } {
    global gap5_defs
    global LREG
    global $w.ops

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Find internal joins"

    if {[info exists $w.ops]} {
	unset $w.ops
    }

    # 3 Panes
    # 
    # Contig inputs + display data
    #    Contig set 1
    #    Contig set 2
    #    Maximum alignment list to display (remove?)
    #
    # Search params
    #    Algorithm selection (sensitive, quick, fastest)
    #    Use hidden data?
    #    Word length (1)
    #
    #    Sensitive
    #        Diagonal Threshold
    #        Alignment band size (percent) (2)
    #
    #     Quick/Fast
    #        Uses banded + band size (width) (2)
    #        Minimum initial match length
    #        Auto-filter repeats
    #        Maximum word over-representation
    #
    # Alignment filtering
    #    Maximum percent mismatch
    #    Min/max overlap length
    #    Min/Max sequence depth
    #    Read-pair filtering
    #        Size of contig ends
    #        Minimum mapping qual
    #        Minimum spanning pair count
    #        Minimum spanning pair percentage
    #        Library selection (separate dialogue?)

    set b [ttk::notebook $w.book]
    pack $b -side top -fill both -expand 1

    bind $w <Alt-Left>  "$b select \[expr {(\[$b select\]-1)%%3}\]"
    bind $w <Alt-Right> "$b select \[expr {(\[$b select\]+1)%%3}\]"
    bind $w <Alt-c> "$b select 0"
    bind $w <Alt-s> "$b select 1"
    bind $w <Alt-f> "$b select 2"

    ###########################################################################
    # Contig identification
    #
    # contig_id boxes for selecting single contigs
    # need both first as FIJ_config_contig_ids gets called by lorf_in

    set f [frame $b.f_i -padx 5 -pady 5]
    set f_i $f

    $b add $f -text "Contigs"

    labelframe $f.l1 -text "List 1"
    contig_id $f.l1.id1 -io $io -range 0 -trace 2 -frame_relief flat

    frame $f.padding -height 10

    labelframe $f.l2 -text "List 2"
    contig_id $f.l2.id2 -io $io -range 0 -trace 0 -frame_relief flat

    #--- input set 1
    lorf_in $f.l1.infile1 [keylget gap5_defs FIJ.INFILE1] \
	"{FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 normal}" \
	-bd 2 -relief flat

    pack $f.l1.infile1 -fill x
    pack $f.l1.id1 -fill x

    #--- input set 2
    lorf_in $f.l2.infile2 [keylget gap5_defs FIJ.INFILE2] \
	"{FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 normal}" \
	-bd 2 -relief flat

    pack $f.l2.infile2 -fill x
    pack $f.l2.id2 -fill x

    #--- selecting a contig_id box makes it the next to be updated
    bind [entrybox_path $f.l1.id1.ent] <<select>> "FIJ_config_contig_ids $f id1"
    bind [entrybox_path $f.l2.id2.ent] <<select>> "FIJ_config_contig_ids $f id2"
    
    FIJ_config_contig_ids $f id1

    pack $f.l1 $f.padding $f.l2 -fill both

    ###########################################################################
    # Contig searching options
    #
    set f [frame $b.f_s  -padx 5 -pady 5]
    set f_s $f
    $b add $f -text "Searching"

    frame $f.padding1 -height 10
    labelframe $f.s -text "Sensitive algorithm"

    frame $f.padding2 -height 10
    labelframe $f.f -text "Quick algorithm"

    #--- hidden data
    labelframe $f.hidden -text "Hidden/cutoff data"

    HiddenParametersDialogInit $w.ops
    button $f.hidden.options \
	-text "Edit hidden data paramaters" \
	-command "HiddenParametersDialog $w $w.ops"

    set hd [keylget gap5_defs FIJ.HIDDEN]
    xyn $f.hidden.yn \
	-label [keylget hd NAME] \
	-orient horizontal \
	-ycommand "$f.hidden.options configure -state normal" \
	-ncommand "$f.hidden.options configure -state disabled" \
	-variable $w.ops(use_hidden)
    set $w.ops(use_hidden) [keylget hd VALUE]

    pack $f.hidden.yn -side top -fill x
    pack $f.hidden.options -side top -anchor w


    #--- select mode
    SetDefaultTags FIJ.TAGS

    labelframe $f.sel_mode -text "Consensus masking"
    button $f.sel_mode.but \
	    -text "Select tags" \
	    -command "TagDialog FIJ.TAGS $f[keylget gap5_defs SELECT_TAGS.WIN] {}"

    radiolist $f.sel_mode.l \
	    -bd 0 \
            -title "Remove tagged sequence" \
	    -orient horizontal \
	    -default [keylget gap5_defs FIJ.SELMODE.VALUE]\
	    -buttons [format {
		{None -command { %s configure -state disabled}} \
		{Mark -command { %s configure -state normal; \
		                SetDefaultTags %s } } \
		{Mask -command { %s configure -state normal; \
				SetDefaultTags %s } } } \
	    [list $f.sel_mode.but] \
	    [list $f.sel_mode.but] FIJ.TAGS \
	    [list $f.sel_mode.but] FIJ.TAGS ]
    pack $f.sel_mode.l -side top -fill x
    pack $f.sel_mode.but -side top -anchor w

    #---- select word length
    set st [keylget gap5_defs FIJ.WORDLENGTH]
    set b1 [keylget st BUTTON.1]
    set b2 [keylget st BUTTON.2]
    set b3 [keylget st BUTTON.3]

    radiolist $f.word_length \
	    -title [keylget st NAME]\
	    -bd 0 \
	    -relief groove \
	    -orient horizontal \
	    -default [keylget st VALUE] \
	    -buttons [format { {%s} {%s} {%s} } \
			  [list $b1] [list $b2] [list $b3]]

    pack $f.word_length -fill x
    #frame $f.padding -relief groove -bd 2 -height 2
    #pack $f.padding -fill x -pady 5

    #--- Sensitive
    set mm [keylget gap5_defs FIJ.MAXDIAG]
    xentry $f.s.max_prob \
	-label "[keylget mm NAME] ([keylget mm MIN] to [keylget mm MAX])"\
	-default [keylget mm VALUE]\
	-type "check_float [keylget mm MIN] [keylget mm MAX]" \
	-textvariable $w.ops(max_prob)

    set mm [keylget gap5_defs FIJ.BANDSIZE]
    xentry $f.s.band_size \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-type "check_int [keylget mm MIN] [keylget mm MAX]" \
	-textvariable $w.ops(band_size)

    #pack $f.s.label -anchor w -padx 50
    pack $f.s.max_prob -fill x
    pack $f.s.band_size -fill x


    #--- Quick
    set mm [keylget gap5_defs FIJ.USEBAND] 
    xyn $f.f.use_band \
        -label [keylget mm NAME] \
        -orient horiz \
	-variable $w.ops(use_band)
    set $w.ops(use_band) [keylget mm VALUE]

    set mm [keylget gap5_defs FIJ.MINMATCH]
    xentry $f.f.min_match \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(min_match)

    set mm [keylget gap5_defs FIJ.USEFILTERWORDS]
    xyn $f.f.use_filter \
	-label [keylget mm NAME] \
	-orient horizontal \
	-variable $w.ops(use_filter)
    set $w.ops(use_filter) [keylget mm VALUE]

    set mm [keylget gap5_defs FIJ.FILTERWORDS]
    xentry $f.f.filter_cutoff \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(filter_words)

    #--- Algorithm choice
    radiolist $f.blocks \
	-title "Alignment algorithm" \
	-bd 0 \
	-orient horizontal \
	-default [keylget gap5_defs FIJ.ALIGN_MODE]\
	-buttons [format { \
            {fastest   -command {%s configure -state disabled; \
				 %s configure -state disabled;
				 %s configure -state normal;
				 %s configure -state normal;
				 %s configure -state normal;
		                 %s configure -state normal;
		                 set %s(min_match) 25}} \
            {quick     -command {%s configure -state disabled; \
				 %s configure -state disabled;
				 %s configure -state normal;
				 %s configure -state normal;
				 %s configure -state normal;
		                 %s configure -state normal;
	                         set %s(min_match) 20}} \
            {sensitive -command {%s configure -state normal; \
				 %s configure -state normal;
				 %s configure -state disabled;
				 %s configure -state disabled;
				 %s configure -state disabled;
		                 %s configure -state disabled}}} \
            $f.s.band_size $f.s.max_prob $f.f.min_match $f.f.use_band \
		      $f.f.use_filter $f.f.filter_cutoff $w.ops \
            $f.s.band_size $f.s.max_prob $f.f.min_match $f.f.use_band \
		      $f.f.use_filter $f.f.filter_cutoff $w.ops \
            $f.s.band_size $f.s.max_prob $f.f.min_match $f.f.use_band \
		      $f.f.use_filter $f.f.filter_cutoff]


    #--- Match location
    radiolist $f.location \
	-title "Alignment location" \
	-bd 0 \
	-orient horizontal \
	-default [keylget gap5_defs FIJ.ALIGN_LOCATION] \
	-buttons {end/end containment both}

    pack $f.blocks -fill x
    pack $f.location -fill x
    pack $f.hidden -fill x
    pack $f.sel_mode -fill x
    pack $f.padding1
    pack $f.s -fill x
    pack $f.padding2
    pack $f.f -fill x
    pack $f.f.use_band -fill x
    pack $f.f.min_match -fill x
    pack $f.f.use_filter -fill x
    pack $f.f.filter_cutoff -fill x

    ###########################################################################
    # Result filtering

    set f [frame $b.f_f  -padx 5 -pady 5]
    set f_f $f
    $b add $f -text "Filtering"

    #--- Percentage identity
    set mm [keylget gap5_defs FIJ.MAXMIS]
    xentry $f.max_mis \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(max_mismatch)

    #--- Min/max overlap size
    set mm [keylget gap5_defs FIJ.MINOVERLAP]
    xentry $f.min_overlap \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(min_overlap)

    set mm [keylget gap5_defs FIJ.MAXOVERLAP]
    xentry $f.max_overlap \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(max_overlap)


    #--- Min/max depth
    set mm [keylget gap5_defs FIJ.MINDEPTH]
    xentry $f.min_depth \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(min_depth)

    set mm [keylget gap5_defs FIJ.MAXDEPTH]
    xentry $f.max_depth \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE] \
	-textvariable $w.ops(max_depth)


    #--- Read pair filtering
    labelframe $f.rp -text "Read pairs"

    ReadPairParametersDialogInit $io $w.ops

    set rp $f.rp
    set sub ""

    set end_size [keylget gap5_defs FIJ.READPAIR.END_SIZE.VALUE]
#    if { $data(max_end_size) > $end_size } {
#	set end_size $data(max_end_size)
#    }
    xentry $f.rp.end_size \
	-label [keylget gap5_defs FIJ.READPAIR.END_SIZE.NAME] \
	-default $end_size \
	-textvariable $w.ops(rp_end_size)
    lappend sub $f.rp.end_size

    set mm [keylget gap5_defs FIJ.READPAIR.MIN_MAP_QUAL]
    xentry $f.rp.min_mq \
	-label [keylget mm NAME] \
        -default [keylget mm VALUE] \
        -type "check_int [keylget mm MIN] [keylget mm MAX]" \
	-textvariable $w.ops(rp_min_mq)
    lappend sub $f.rp.min_mq

    set mm [keylget gap5_defs FIJ.READPAIR.MIN_FREQ]
    xentry $f.rp.min_freq \
	-label [keylget mm NAME] \
        -default [keylget mm VALUE] \
        -type "check_int [keylget mm MIN] [keylget mm MAX]" \
	-textvariable $w.ops(rp_min_freq)
    lappend sub $f.rp.min_freq

    set mm [keylget gap5_defs FIJ.READPAIR.MIN_PERC]
    xentry $f.rp.min_perc \
	-label [keylget mm NAME] \
        -default [keylget mm VALUE] \
        -type "check_int [keylget mm MIN] [keylget mm MAX]" \
	-textvariable $w.ops(rp_min_perc)
    lappend sub $f.rp.min_perc

    button $f.rp.libs \
	-text "Select libraries" \
	-command [list ReadPairParametersDialog $f $io $w.ops]
    lappend sub $f.rp.libs

    set fij_rp_sc_d "foreach i {$sub} {\$i configure -state disabled}"
    set fij_rp_sc_e "foreach i {$sub} {\$i configure -state normal}"

    radiolist $f.rp.mode \
	-title [keylget gap5_defs FIJ.READPAIR.MODE.NAME] \
	-default 4 \
	-orient horizontal \
	-buttons [list [list [keylget gap5_defs FIJ.READPAIR.MODE.BUTTON1] \
			    -command $fij_rp_sc_e ] \
		      [list [keylget gap5_defs FIJ.READPAIR.MODE.BUTTON2] \
			   -command $fij_rp_sc_e ] \
		      [list [keylget gap5_defs FIJ.READPAIR.MODE.BUTTON3] \
			   -command $fij_rp_sc_e ] \
		      [list [keylget gap5_defs FIJ.READPAIR.MODE.BUTTON4] \
			   -command $fij_rp_sc_d ] \
		 ]
    pack $f.rp.mode $f.rp.end_size $f.rp.min_mq $f.rp.min_freq \
	$f.rp.min_perc -side top -fill x
    pack $f.rp.libs -side top -anchor w

 
    #--- Filter to unique ends only
    xyn $f.unique_ends \
	-orient horizontal \
	-label "Filter non-unique contig matches" \
	-variable $w.ops(unique_ends);
    set $w.ops(unique_ends) [keylget gap5_defs FIJ.UNIQUE_ENDS]

    #--- Alignment display size
    xentry $f.align_length \
	-label "Maximum alignment length to list (bp)" \
	-default [keylget gap5_defs FIJ.MAX_ALIGNMENT] \
	-type "check_int 0 10000000" \
	-textvariable $w.ops(align_length)


    pack $f.min_overlap -fill x
    pack $f.max_overlap -fill x
    pack $f.min_depth -fill x
    pack $f.max_depth -fill x
    pack $f.max_mis -fill x
    pack $f.align_length -fill x
    pack $f.unique_ends -fill x
    pack $f.rp -fill x


    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $w.ok_cancel \
	-ok_command "FIJ_OK_Pressed $w $io $w.ops \
                    $f_i.l1.infile1 $f_i.l1.id1 $f_i.l2.infile2 $f_i.l2.id2 \
                    $f_s.word_length $f_s.blocks $f_s.location $f_s.sel_mode.l\
                    $f_f.rp.mode" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {FIJ-Dialogue}" \
	-bd 2 \
	-relief groove
    ###########################################################################

#    pack $f.sel_task -fill x
#    pack $f.sc -fill x

    $w.book select 0
    pack $w.ok_cancel -fill x

}

proc FIJ_config_contig_ids { f id { state unchanged } } {
    # Set which of the contig_id boxes get updated when clicking in the
    # contig list or contig selector, and in which order.
    if { $state != "unchanged" } {
	if {$id == "id1"} {
	    contig_id_configure "$f.l1.$id" -state $state
	} else {
	    contig_id_configure "$f.l2.$id" -state $state
	}
    }
    set id1_state [eval [entrybox_path $f.l1.id1.ent] cget -state]
    set id2_state [eval [entrybox_path $f.l2.id2.ent] cget -state]

    set boxen {}
    if { "$id" == "id1" } {
	if { $id1_state == "normal" } { lappend boxen "$f.l1.id1" }
	if { $id2_state == "normal" } { lappend boxen "$f.l2.id2" }
    } else {
	if { $id2_state == "normal" } { lappend boxen "$f.l2.id2" }
	if { $id1_state == "normal" } { lappend boxen "$f.l1.id1" }	
    }
    # Need at least one item or we get errors...
    if { [llength $boxen] == 0 } { lappend boxen "$f.l1.id1" }
    SetCurFrame $f.l1 $boxen
}

###########################################################################
proc FIJ_OK_Pressed { f io aname infile1 id1 infile2 id2 
		      word_length blocks location sel_mode rp_mode } {
    
    global CurContig
    global LREG
    global RREG
    global NGRec
    global gap5_defs

    upvar #0 $aname data

    if {[lorf_in_get $infile1] == 3} {
	set list1 [CreateAllContigList $io]
    } elseif {[lorf_in_get $infile1] == 4} {
	set gel_name [contig_id_gel $id1]
	set list1 "{$gel_name}"
    } else {
	set list1 [lorf_get_list $infile1]
    }

    if {[lorf_in_get $infile2] == 3} {
	set list2 [CreateAllContigList $io]
    } elseif {[lorf_in_get $infile2] == 4} {
	set gel_name [contig_id_gel $id2]
	set list2 "{$gel_name}"
    } else {
	set list2 [lorf_get_list $infile2]
    }

    if {$data(use_hidden)} {
	set win_size $data($data(which_mode)_win_size)
	set max_dash $data(base_max_dash)
	set min_conf $data(conf_min_conf)
	set use_conf [lsearch {base conf} $data(which_mode)]
	set use_hidden 1
    } else {
	set win_size 0
	set max_dash 0
	set min_conf 0
	set use_conf 0
	set use_hidden 0
    }

    set word_length [lindex {? 12 8 4} [radiolist_get $word_length]]
    
    set rp_mode_str [lindex {all_all end_all end_end off } \
		     [expr {[radiolist_get $rp_mode] - 1} ]]

    set fast_mode 0
    if {[radiolist_get $blocks] <= 2} {
        # Quick method
	set min_match $data(min_match)
	set band_size $data(use_band)

	if {[radiolist_get $blocks] == 1} {
	    set fast_mode 1
	}
    } else {
	# Sensitive method
	set min_match 0
	set band_size $data(band_size)
    }

    set ends         [expr {[radiolist_get $location]&1}]
    set containments [expr {([radiolist_get $location]&2)==2}]
    set masking      [radiolist_get $sel_mode]
    if {($masking == 2) || ($masking == 3)} {
        set active_tags [GetDefaultTags FIJ.TAGS]
    } else {
	set active_tags {}
    }

    if {$data(use_filter) == 0} {
	set filter_words 0
    } else {
	set filter_words $data(filter_words)
    }


    # Destroy dialog before showing plot to avoid window activation problems
    destroy $f


    #draw the contig comparator and dot plot
    ContigComparator $io 


    SetBusy
    set id [log_call find_internal_joins \
		-mask [lindex {"" none mark mask} $masking] \
		-tag_types $active_tags \
		-io            $io \
		-min_overlap   $data(min_overlap) \
		-max_overlap   $data(max_overlap) \
		-min_depth     $data(min_depth) \
		-max_depth     $data(max_depth) \
		-max_pmismatch $data(max_mismatch) \
		-word_length   $word_length \
		-max_prob      $data(max_prob)\
		-min_match     $min_match \
		-band          $band_size \
		-win_size      $win_size \
		-max_dashes    $max_dash \
		-min_conf      $min_conf \
		-use_conf      $use_conf \
		-use_hidden    $data(use_hidden) \
		-max_display   $data(align_length) \
		-fast_mode     $fast_mode \
		-filter_words  $filter_words \
		-rp_mode       $rp_mode_str \
		-rp_end_size   $data(rp_end_size) \
		-rp_min_mq     $data(rp_min_mq) \
		-rp_min_freq   $data(rp_min_freq) \
		-rp_min_perc   $data(rp_min_perc) \
		-rp_libraries  $data(rp_libs) \
		-unique_ends   $data(unique_ends) \
		-containments  $containments \
		-ends          $ends \
		-contigs1      $list1 \
		-contigs2      $list2]

    if {$id > 0} {
	# Draw it too
	result_notify \
	    -io $io \
	    -id $id \
	    -type GENERIC \
	    -args "{task TASK_CS_PLOT}"
    }

    ClearBusy
}

# Hidde data parameters box, used by Find Internal Joins and Calculate
# Consensus
# $w is the parent window
# $aname is a global array name in which to fill out the results
proc HiddenParametersDialog {w aname} {
    global gap5_defs
    upvar #0 $aname data

    set f $w.hiddenparam
    if {![winfo exists $f]} {
	xtoplevel $f -resizable 0
	wm title $f "Hidden data parameters"
    } else {
	raise $f
	wm deiconify $f
    }

    ###################################################
    # extend by confidence values

    set w [frame $f.conf -bd 2 -relief groove]

    label $w.label -text [keylget gap5_defs HIDDEN.CONF.LABEL]

    set ws [keylget gap5_defs HIDDEN.CONF.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(conf_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap5_defs HIDDEN.CONF.MINCONF]
    scalebox $w.min_conf \
	    -title [keylget md NAME] \
	    -orient horizontal \
	    -to [keylget md MAX] \
	    -from [keylget md MIN]\
	    -default $data(conf_min_conf) \
	    -width 5 \
	    -type CheckInt

    pack $w.label -side top -anchor c
    pack $w.win_size -side top -fill x
    pack $w.min_conf -side top -fill x

    ###################################################
    # extend by base calls

    set w [frame $f.unc -bd 2 -relief groove]

    label $w.label -text [keylget gap5_defs HIDDEN.UNC.LABEL]

    set ws [keylget gap5_defs HIDDEN.UNC.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(base_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap5_defs HIDDEN.UNC.MAXDASH]
    scalebox $w.max_dash \
	    -title [keylget md NAME] \
	    -orient horizontal \
	    -to [keylget md MAX] \
	    -from [keylget md MIN]\
	    -default $data(base_max_dash) \
	    -width 5 \
	    -type CheckInt

    pack $w.label -side top -anchor c
    pack $w.win_size -side top -fill x
    pack $w.max_dash -side top -fill x

    ###################################################
    # Main frame and question
    set yn [keylget gap5_defs HIDDEN.MODE]
    yes_no $f.which \
	-title [keylget yn NAME] \
	-bd 2 \
	-relief groove \
	-orient horizontal \
	-ycommand "scalebox_configure $f.conf.win_size -state normal;\
                   scalebox_configure $f.conf.min_conf -state normal;\
                   scalebox_configure $f.unc.win_size  -state disabled;\
                   scalebox_configure $f.unc.max_dash  -state disabled" \
	-ncommand "scalebox_configure $f.conf.win_size -state disabled;\
                   scalebox_configure $f.conf.min_conf -state disabled;\
                   scalebox_configure $f.unc.win_size  -state normal;\
                   scalebox_configure $f.unc.max_dash  -state normal" \
	-default [lsearch {base conf} $data(which_mode)]

    ###################################################
    # Final bits

    okcancelhelp $f.ok_cancel \
	-ok_command "HiddenParametersDialog_OK $f $aname" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {FIJ-Dialogue}" \

    pack $f.which $f.conf $f.unc $f.ok_cancel \
	-side top -fill both -expand 1
}

proc HiddenParametersDialogInit {aname} {
    global gap5_defs
    upvar #0 $aname data

    set data(which_mode)    [lindex {base conf} [keylget gap5_defs HIDDEN.MODE.VALUE]]
    set data(base_win_size) [keylget gap5_defs HIDDEN.UNC.WINSIZE.VALUE]
    set data(base_max_dash) [keylget gap5_defs HIDDEN.UNC.MAXDASH.VALUE]
    set data(conf_win_size) [keylget gap5_defs HIDDEN.CONF.WINSIZE.VALUE]
    set data(conf_min_conf) [keylget gap5_defs HIDDEN.CONF.MINCONF.VALUE]
}

proc HiddenParametersDialog_OK {w aname} {
    upvar #0 $aname data

    set data(which_mode)    [lindex {base conf} [yes_no_get   $w.which]]
    set data(base_win_size) [scalebox_get $w.unc.win_size]
    set data(base_max_dash) [scalebox_get $w.unc.max_dash]
    set data(conf_win_size) [scalebox_get $w.conf.win_size]
    set data(conf_min_conf) [scalebox_get $w.conf.min_conf]

    destroy $w
}

proc FIJRPUpdateLibList { io libinfo } {
    upvar $libinfo li
    set db [$io get_database]
    set nl [$db get_num_libraries]
    set li(num_libs) $nl
    set max_end_size 0
    for {set i 0} {$i < $nl} {incr i} {
	set rec [$db get_library_rec $i]
        set lib [$io get_library $rec]
        $lib update_stats
	
	set count [$lib get_count]
        set size  [$lib get_insert_size]
	set sd    [$lib get_insert_sd]
	for {set max 0; set k 0; set j 0} {$j < 3} {incr j} {
	    if {$max < [lindex $count $j]} {
                set max [lindex $count $j]
                set k $j
            }
        }
	set count [lindex $count $k]
        set size  [lindex $size  $k]
	set sd    [lindex $sd    $k]
	set end_size [expr { int(($size + 3 * $sd + 9) / 10) * 10 }]
	if { $max_end_size < $end_size } { set max_end_size $end_size }
	set li($i) [ list $rec [$lib get_name] $count $size ]
	$lib delete
    }
    set li(max_end_size) $max_end_size
}

proc ReadPairParametersDialog { w io aname } {
    global gap5_defs
    upvar #0 $aname data

    set f $w.readpair

    if {![winfo exists $f]} {
	xtoplevel $f -resizable 0
	wm title $f [keylget gap5_defs FIJ.READPAIR.DIALOG.NAME]
    } else {
	wm deiconify $f
	raise $f
	return
    }
    
    frame $f.libs
    label $f.libs.label -text "By default templates in all libraries are used to spot read pairs. To restrict to specific libraries, select them from the list below." -wrap 400 -justify left

    tablelist $f.libs.tl \
        -height 5 \
        -columns {15 Name 10 "Pair count" 10 "Insert size"} \
        -selectmode extended \
        -exportselection 0 \
        -stretch 0 \
        -yscrollcommand [list $f.libs.ys set]
    scrollbar $f.libs.ys -command "$f.libs.tl yview" -orient vertical
    pack $f.libs.label -side top -fill both
    pack $f.libs.tl -side left -expand 1 -fill both
    pack $f.libs.ys -side right -fill both

    upvar #0 $f.libs.tl l_rec

    FIJRPUpdateLibList $io data
    set nl $data(num_libs)
    for {set i 0} {$i < $nl} {incr i} {
	set li $data($i)
	foreach { rec name count size } $li {
	    $f.libs.tl insert end [list $name $count $size]
	    set l_rec($i) $rec
	}
    }

    okcancelhelp $f.ok_cancel \
	-ok_command "ReadPairParametersDialogOK $f $aname" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {FIJ-Dialogue}" \
	-bd 2 \
	-relief groove
    
    pack $f.libs $f.ok_cancel \
	-side top -fill both -expand 1
}

proc ReadPairParametersDialogInit {io aname} {
    global gap5_defs
    upvar #0 $aname data

    FIJRPUpdateLibList $io data
    set data(rp_libs)     {}
}

proc ReadPairParametersDialogOK {w aname} {
    upvar #0 $aname data
    upvar #0 $w.libs.tl l_rec

    set libs {}
    foreach idx [$w.libs.tl curselection] {
	lappend libs $l_rec($idx)
    }
    set data(rp_libs) $libs

    destroy $w
}
