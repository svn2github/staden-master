#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

set seq_techs [list Sanger Illumina SOLiD 454 Helico IonTorrent PacBio ONT]

#-----------------------------------------------------------------------------
# Sets the consensus and quality cutoff values
#
proc ConfigureCutoffs {io} {
    global gap5_defs consensus_iub

    set l [keylget gap5_defs CONFIGURE]
    set t [keylget l WIN]

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Consensus algorithm"

    scalebox $t.consensus \
	-title [keylget l CONSENSUS_NAME] \
	-orient horizontal \
	-from 1 \
	-to 100 \
	-width 5 \
	-default [expr [keylget gap5_defs CONSENSUS_CUTOFF]*100]\
	-type CheckInt

    scalebox $t.quality \
	-title [keylget l QUALITY_NAME] \
	-orient horizontal \
	-from -1 \
	-to 100 \
	-width 5 \
	-default [keylget gap5_defs QUALITY_CUTOFF] \
	-type CheckInt

    checkbutton $t.chem \
	-text [keylget l CHEMISTRY_NAME] \
	-variable $t.chem.v
    global $t.chem.v chem_as_double
    set $t.chem.v $chem_as_double

    checkbutton $t.iub \
	-text "Produce IUB codes in consensus" \
	-variable $t.iub.v
    global $t.iub.v
    set $t.iub.v $consensus_iub

    global seq_techs default_seq_tech
    xcombobox $t.seq_tech \
	-label "Seq. tech. for unknown datasets" \
	-fixed_list 1 \
	-values $seq_techs
    $t.seq_tech set [lindex $seq_techs [expr {$default_seq_tech-1}]]

    radiolist $t.mode \
	-title [keylget l CMODE_NAME] \
	-bd 2 -relief groove \
	-default [expr [keylget gap5_defs CONSENSUS_MODE]+1] \
	-buttons [format {{%s -command %s} {%s -command %s} {%s -command %s}} \
		[list [keylget l CMODE_BUTTON1]] \
		[list "scalebox_configure $t.quality -state disabled"] \
		[list [keylget l CMODE_BUTTON2]] \
		[list "scalebox_configure $t.quality -state normal"] \
		[list [keylget l CMODE_BUTTON3]] \
		[list "scalebox_configure $t.quality -state disabled"]]

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "ConfigureCutoffs2 $io 0 $t $t.consensus $t.quality \
		$t.chem.v $t.mode $t.iub.v $t.seq_tech" \
	-perm_command "ConfigureCutoffs2 $io 1 $t $t.consensus $t.quality \
		$t.chem.v $t.mode $t.iub.v" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Con-Calculation}"

    pack $t.mode $t.consensus $t.quality -side top -fill both
    pack $t.chem -pady 7 -side top -anchor w
    pack $t.iub -pady 7 -side top -anchor w
    pack $t.seq_tech -side top -anchor w
    pack $t.ok -side top -fill both
}

proc ConfigureCutoffs2 {io perm t consensus quality chem mode iub seq_tech} {
    global gap5_defs consensus_cutoff quality_cutoff chem_as_double $chem env
    global consensus_mode consensus_iub $iub default_seq_tech seq_techs

    set consensus_mode [expr [radiolist_get $mode]-1]
    if {$consensus_mode == 0} {
	set quality_cutoff -1
    } else {
        set quality_cutoff [scalebox_get $quality]
	if {$quality_cutoff == -1} {
	    set quality_cutoff 0
	}
    }
    set consensus_cutoff [expr [scalebox_get $consensus]/100.0]
    set chem_as_double [set $chem]
    set consensus_iub [set $iub]
    set seq_tech [lsearch $seq_techs [$seq_tech get]]
    if {$seq_tech >= 0} {
	incr seq_tech
	if {$seq_tech != $default_seq_tech} {
	    # Warning: changing this invalidates cached consensus. We should
	    # warn the user and clear it.
	    set yn [tk_messageBox -type yesno -title "Clear cached consensus" -message "Changing the default sequencing technology may invalidate any precomputed consensus.  Hence the consensus cache will need to be purged, which may take a few minutes.  Do you wish to make this configuration change?"]

	    if {$yn == "yes"} {
		set default_seq_tech $seq_tech
		invalidate_consensus_cache $io
		set_database_param $io default_seq_tech $default_seq_tech
	    }
	}
    }

    keylset gap5_defs CONSENSUS_MODE   $consensus_mode
    keylset gap5_defs CONSENSUS_CUTOFF $consensus_cutoff
    keylset gap5_defs QUALITY_CUTOFF   $quality_cutoff
    keylset gap5_defs CHEM_AS_DOUBLE   $chem_as_double
    keylset gap5_defs CONSENSUS_IUB    $consensus_iub

    if {$perm} {
        update_defs gap5_defs $env(HOME)/.gap5rc \
	    CONSENSUS_MODE CONSENSUS_CUTOFF QUALITY_CUTOFF CHEM_AS_DOUBLE \
	    CONSENSUS_IUB
    }

    destroy $t
}

# Stores variable 'var'/'val' in the global database CNFG annotation.
proc set_database_param {io var val} {
    set db [$io get_database]
    if {[set arec [$db get_config_anno]] == 0} {
	set arec [$io new_anno_ele 16 0 0 0]
	if {$arec == 0} {
	    return -1
	}

	set ae [$io get_anno_ele $arec]
	$ae set_type CNFG

	$db set_config_anno $arec
    } else {
	set ae [$io get_anno_ele $arec]
    }

    set comment [$ae get_comment]
    array set vars ""

    foreach line [split $comment "\n"] {
	if {$line == "\n"} continue
	if {[regexp {set ([^ ]*) (.*)} $line _ a b] != 1} continue

	set vars($a) $b
    }
    set vars($var) $val

    set comment ""
    foreach key [array names vars] {
	append comment "set $key $vars($key)\n"
    }
    
    $ae set_comment $comment

    $io flush

    return 0
}

# Invalidates the consensus cache for all contigs
proc invalidate_consensus_cache {io} {
    set db [$io get_database]
    set nc [$db get_num_contigs]

    vfuncheader "Purging consensus cache"
    for {set i 0} {$i < $nc} {incr i} {
	set crec [$io contig_order $i]
	set c [$io get_contig $crec]
	vmessage "Purging consensus for contig [expr {$i+1}] of $nc (=$crec)"
	$c invalidate_consensus
	$c delete
	$io flush
    }
}

#-----------------------------------------------------------------------------
# Sets the 'maxseq' variable
proc SetMaxseq {io} {
    global gap5_defs maxseq maxdb

    set l [keylget gap5_defs CONFMAXSEQ]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set maxseq"

    if {[catch {set db [io_read_database $io]}]} {
	set min1 10000
	set min2 800
    } else {
	set nc [keylget db num_contigs]
	set nr [keylget db num_readings]
	set tl [db_info t_contig_length $io]
	set min1 [expr round($tl+2+21*$nc)]
	set min2 [expr {$nc+$nr+2}]
    }
	
    entrybox $t.maxseq_val \
	-title "[keylget l MAXSEQ.NAME] (min $min1)"\
	-default $maxseq \
	-width 10 \
	-type "CheckIntMin $min1"

    entrybox $t.maxdb_val \
	-title "[keylget l MAXDB.NAME] (min $min2)"\
	-default $maxdb \
	-width 10 \
	-type "CheckIntMin $min2"

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "SetMaxseq2 $t $t.maxseq_val $t.maxdb_val" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Conf-Set Maxseq}"

    pack $t.maxseq_val $t.maxdb_val $t.ok -side top -fill both
}

proc SetMaxseq2 {t maxseq_w maxdb_w} {
    global maxseq maxdb
    if {[set new [entrybox_get $maxseq_w]] == ""} {bell; return}
    set maxseq $new
    if {[set new [entrybox_get $maxdb_w]] == ""} {bell; return}
    set maxdb $new
    verror ERR_WARN gap "setting maxseq parameter to $maxseq"
    verror ERR_WARN gap "setting maxdb parameter to $maxdb"
    destroy $t
}



# -----------------------------------------------------------------------------
# Menu specifications (eg 'user levels')

proc ConfigureMenus {} {
    global env gap5_defs

    set w .configure_menus

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Configure menus"

    label $w.label \
	    -text "Please select user-level for menus"
    pack $w.label -side top
    global $w.Menu
    foreach m [keylget gap5_defs MENU_LEVELS] {
	foreach {l r} $m {}
	radiobutton $w.menu_$r \
		-text $l \
		-variable $w.Menu \
		-value $r
#		-command "ConfigureMenus2 $w $r"
	pack $w.menu_$r -side top -anchor w
    }
    set $w.Menu [keylget gap5_defs MENU_LEVEL]

    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "ConfigureMenus2 0 $w \[set $w.Menu\]" \
	-perm_command "ConfigureMenus2 1 $w \[set $w.Menu\]" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap4 {Conf-Configure Menus}"
    pack $w.ok -side top -fill both
}

proc ConfigureMenus2 {perm w menu} {
    global gap5_defs env gap_menu

    destroy $w

    keylset gap5_defs MENU_LEVEL $menu

    set mpath .mainwin.menus
    $mpath delete 0 end
    foreach child [winfo children $mpath] {
	destroy $child
    }

    create_menus $gap_menu $mpath $menu
    reset_menu_state gap_menu .mainwin.menus

    if {$perm} {
        update_defs gap5_defs $env(HOME)/.gap5rc MENU_LEVEL
    }
}

# -----------------------------------------------------------------------------
# Selection of genetic code

proc SetGeneticCode {} {
    global gap5_defs
    set w .genetic_code

    if {[xtoplevel $w] == ""} return
    wm title $w "Set genetic code"

    frame $w.codes
    pack $w.codes -side top -fill both -expand 1

    listbox $w.codes.l \
	    -xscrollcommand "$w.codes.xs set" \
	    -yscrollcommand "$w.codes.ys set" \
	    -width 25 \
	    -height 15
    scrollbar $w.codes.xs \
	    -orient horizontal \
	    -command "$w.codes.l xview"
    scrollbar $w.codes.ys \
	    -orient vertical \
	    -command "$w.codes.l yview"
    grid rowconfigure $w.codes 0 -weight 1
    grid columnconfigure $w.codes 0 -weight 1
    grid $w.codes.l $w.codes.ys -sticky nsew
    grid $w.codes.xs -sticky nsew

    set fd [open [keylget gap5_defs GENETIC_CODE_DIR]/code_index r]
    set i 0
    set index ""
    while {[gets $fd line] != -1} {
	regexp {([^ ]*) *(.*)} $line tmp file name
	$w.codes.l insert end "$name"
	lappend index $file
	incr i
    }
    close $fd
    bind $w.codes.l <<use>> "SetGeneticCode2 0 $w $w.codes.l [list $index]"

    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "SetGeneticCode2 0 $w $w.codes.l [list $index]" \
	-perm_command "SetGeneticCode2 1 $w $w.codes.l [list $index]" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Conf-Set Genetic Code}"
    pack $w.ok -side top -fill both
}

proc SetGeneticCode2 {perm w listbox index} {
    global gap5_defs env

    if {[set cur [$listbox curselection]] == ""} {
	bell
	return
    }
    set file [lindex $index [lindex $cur 0]]

    if {[load_genetic_code \
	    -filename [keylget gap5_defs GENETIC_CODE_DIR]/$file] != -1} {
	vfuncheader "Set genetic code"
	vmessage "[$listbox get [lindex $cur 0]] genetic code:\n"
	set fd [open [keylget gap5_defs GENETIC_CODE_DIR]/$file r]
	vmessage [read $fd]
	close $fd
    }

    keylset gap5_defs GENETIC_CODE [keylget gap5_defs GENETIC_CODE_DIR]/$file

    if {$perm} {
        update_defs gap5_defs $env(HOME)/.gap5rc GENETIC_CODE
    }

    destroy $w
}

# -----------------------------------------------------------------------------
# Adjustment of alignment scores (although not matrix).

proc SetAlignmentScores {} {
    global gap5_defs align_open_cost align_extend_cost
    set w .alignment_scores

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Alignment scores"

    scalebox $w.open \
	    -orient horizontal \
	    -width 3 \
	    -title "Open penalty" \
	    -from [keylget gap5_defs ALIGNMENT.OPEN.MIN] \
	    -to [keylget gap5_defs ALIGNMENT.OPEN.MAX] \
	    -default $align_open_cost
    scalebox $w.extend \
	    -orient horizontal \
	    -width 3 \
	    -title "Extend penalty" \
	    -from [keylget gap5_defs ALIGNMENT.EXTEND.MIN] \
	    -to [keylget gap5_defs ALIGNMENT.EXTEND.MAX] \
	    -default $align_extend_cost
    pack $w.open $w.extend -side top -fill both

    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "SetAlignmentScores2 0 $w $w.open $w.extend" \
	-perm_command "SetAlignmentScores2 1 $w $w.open $w.extend" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Conf-Alignment Scores}"
    pack $w.ok -side top -fill both
}

proc SetAlignmentScores2 {perm w open extend} {
    global env gap5_defs align_open_cost align_extend_cost

    set align_open_cost [scalebox_get $open]
    set align_extend_cost [scalebox_get $extend]

    keylset gap5_defs ALIGNMENT.OPEN.COST $align_open_cost
    keylset gap5_defs ALIGNMENT.EXTEND.COST $align_extend_cost

    if {$perm} {
        update_defs gap5_defs $env(HOME)/.gap5rc \
		ALIGNMENT.OPEN.COST ALIGNMENT.EXTEND.COST
    }

    destroy $w
}

# -----------------------------------------------------------------------------
# Setting of RAWDATA environment variable

proc SetRawData {io} {
    global gap5_defs env
    set w .raw_data_path

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Trace file location"

    if {[set nn [FindRawDataNote $io]] != 0} {
	set n [io_read_note $io $nn]
	if {[keylget n annotation] != 0} {
	    set rawdata [io_read_text $io [keylget n annotation]]
	} else {
	    set rawdata ""
	}
    } elseif {[info exists env(RAWDATA)]} {
	set rawdata $env(RAWDATA)
    } else {
	set rawdata ""
    }

    entrybox $w.val \
	-title "Trace file directories"\
	-default "$rawdata" \
	-width 40 \
	-command "SetRawData2 $io $w"

    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "SetRawData2 $io $w \[entrybox_get $w.val\]" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Conf-Trace File Location}"

    pack $w.val $w.ok -side top -fill both
}

# Looks for the RAWD note.
# Returns note number if found, or zero if not.
proc FindRawDataNote {io} {
    set db [io_read_database $io]
    if {[keylget db notes] == 0} {
	return 0
    }
    for {set i [keylget db notes]} {$i != 0} {set i [keylget n next]} {
	set n [io_read_note $io $i]
	if {[string compare [keylget n type] "RAWD"] == 0} {
	    return $i
	}
    }
    return 0
}

# Sets the RAWD note.
# If there is not an existing RAWD note, set note_num to zero, otherwise
# it should be set to the note number.
proc WriteRawDataNote {io note_num path} {
    global env

    # Set the environment variable
    set env(RAWDATA) $path

    if {$note_num == 0} {
	set note_num [new_note -io $io -type RAWD -to "database" -number 0]
    }
    edit_note -io $io -note $note_num -type RAWD -comment $path
}

proc SetRawData2 {io w rawdata} {
    global read_only

    destroy $w

    # If database is writable, also set the RAWD database note
    if {$read_only == 0} {
	set nn [FindRawDataNote $io]
	WriteRawDataNote $io $nn $rawdata
    }
}

#appends to an existing RAWD note
proc AppendRawDataNote {io path} {
    global env
    
    set nn 0
    if {[set nn [FindRawDataNote $io]] != 0} {
	set n [io_read_note $io $nn]
	if {[keylget n annotation] != 0} {
	    set rawdata [io_read_text $io [keylget n annotation]]
	} else {
	    set rawdata ""
	}
    } elseif {[info exists env(RAWDATA)]} {
	#NB this will add the environment variable into the notebook as
	#otherwise this variable will be overwritten by the appended note 
	#which is added into the notebook.
	set rawdata $env(RAWDATA)
    } else {
	set rawdata ""
    }

    #FIXME: need to deal with : in windows

    #remove \n from rawdata
    regsub "\n$" $rawdata {} rawdata 
    set rawdata $rawdata:$path

    WriteRawDataNote $io $nn $rawdata
}

proc SetTemplateStatusConfig {} {
    global template_size_tolerance min_vector_len
    global ignore_all_ptype ignore_custom_ptype

    set w .template_size_tolerance

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Template Status Configurations"

    xentry $w.vector \
	-label "Minimum valid vector tag length" \
	-default $min_vector_len

    xentry $w.size \
	-label "Size limit scale factor" \
	-default $template_size_tolerance

    xyn $w.ignore_all_ptype \
	-label "Ignore all primer-type values" \
	-default $ignore_all_ptype \

    xyn $w.ignore_custom_ptype \
	-label "Ignore custom primer-type values" \
	-default $ignore_custom_ptype \

    okcancelhelp $w.ok \
	-ok_command "SetTemplateStatusConfig2 $w 0" \
	-perm_command "SetTemplateStatusConfig2 $w 1" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Conf-Template Status}"

    pack $w.size $w.vector -side top -fill both
    pack $w.ignore_all_ptype $w.ignore_custom_ptype -side top -fill both -expand 1
    pack $w.ok -side top -fill both
}

proc SetTemplateStatusConfig2 {w perm} {
    global template_size_tolerance min_vector_len gap5_defs env
    global ignore_all_ptype ignore_custom_ptype
   
    set new [$w.size get]
    if {$new <= 0} {
	bell
	return
    }

    set veclen [$w.vector get]
    if {$veclen < 0} {
	bell
	return
    }

    set ignore_all_ptype    [$w.ignore_all_ptype get]
    set ignore_custom_ptype [$w.ignore_custom_ptype get]

    destroy $w

    set template_size_tolerance $new
    set min_vector_len $veclen

    keylset gap5_defs TEMPLATE_TOLERANCE  $template_size_tolerance
    keylset gap5_defs MIN_VECTOR_LENGTH   $min_vector_len
    keylget gap5_defs IGNORE_ALL_PTYPE    $ignore_all_ptype
    keylget gap5_defs IGNORE_CUSTOM_PTYPE $ignore_custom_ptype

    global template_check_flags
    set template_check_flags \
	[expr {$ignore_all_ptype*4 + $ignore_custom_ptype*8}]

    if {$perm} {
	update_defs gap5_defs $env(HOME)/.gap5rc TEMPLATE_TOLERANCE
	update_defs gap5_defs $env(HOME)/.gap5rc MIN_VECTOR_LENGTH
	update_defs gap5_defs $env(HOME)/.gap5rc IGNORE_ALL_PTYPE
	update_defs gap5_defs $env(HOME)/.gap5rc IGNORE_CUSTOM_PTYPE
    }
}
