proc set_global_defaults {} {
    uplevel #0 {
	global gap5_defs
	set consensus_mode          [keylget gap5_defs CONSENSUS_MODE]
	set consensus_cutoff        [keylget gap5_defs CONSENSUS_CUTOFF]
	set quality_cutoff          [keylget gap5_defs QUALITY_CUTOFF]
	set chem_as_double          [keylget gap5_defs CHEM_AS_DOUBLE]
	set consensus_iub           [keylget gap5_defs CONSENSUS_IUB]
	set template_size_tolerance [keylget gap5_defs TEMPLATE_TOLERANCE]
	set min_vector_len          [keylget gap5_defs MIN_VECTOR_LENGTH]
	set align_open_cost         [keylget gap5_defs ALIGNMENT.OPEN.COST]
	set align_extend_cost       [keylget gap5_defs ALIGNMENT.EXTEND.COST]
	load_alignment_matrix       [keylget gap5_defs ALIGNMENT.MATRIX_FILE]
	set ignore_all_ptype        [keylget gap5_defs IGNORE_ALL_PTYPE]
	set ignore_custom_ptype     [keylget gap5_defs IGNORE_CUSTOM_PTYPE]
    }
}

proc set_database_defaults {io} {
    global default_seq_tech
    set default_seq_tech 2

    set db [$io get_database]
    set arec [$db get_config_anno]
    if {$arec <= 0} return

    set ae [$io get_anno_ele $arec]
    foreach line [split [$ae get_comment] "\n"] {
	if {$line == "\n"} continue
	if {[regexp {set ([^ ]*) (.*)} $line _ a b] != 1} continue

	switch $a {
	    default_seq_tech {
		global $a
		set $a $b
	    }

	    default {
		puts stderr "Unknown configuration variable \"$a\" = \"$b\""
	    }
	}
    }
}