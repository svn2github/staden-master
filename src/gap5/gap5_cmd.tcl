#!/bin/sh

# A generic interface to running a selection of gap5 fucntions as command
# line scripts.
#
# The general syntax is gap5_cmd sub-command [command-options] database.vers
#
#\
exec tclsh $0 ${@+"$@"}

#-----------------------------------------------------------------------------
# Startup code
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package gap5

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
	
# Error reporting, override the tk route
catch {rename tk_messageBox {}}
proc tk_messageBox {args} {
    foreach {a b} $args {
	set opt($a) $b
    }

    if {[info exists opt(-icon)] && $opt(-icon) == "error"} {
	global errorCode errorInfo
	puts stderr "ERROR: $opt(-message)"

	puts "\nError code: $errorCode"
	puts "\nError message:\n$errorInfo"
    } else {
	puts $opt(-message)
    }
}

# Generates usage text
proc usage {cmd e opts} {
    puts "Usage: gap5_cmd $cmd \[options\] DBNAME.VERS"
    puts ""
    puts {Where [options] are any of the following:}
    
    foreach {opt arg def param help} $opts {
	if {$opt == ""} {
	    puts ""
	    continue
	}

	if {$arg} {
	    append help " \[$def\]"
	}

	set h ""
	set l 0
	foreach word [split $help " \n\t"] {
	    set sl [string length $word]
	    if {$l + $sl > 52} {
		append h "\n                          "
		set l 0
	    }
	    incr l $sl
	    incr l
	    append h " $word"
	}

	regsub {\|} $opt {|-} opt

	if {$arg} {
	    puts [format "  -%-20s %s" "$opt '$param'" $h]
	} else {
	    puts [format "  -%-20s %s" $opt $h]
	}
    }
    puts ""
    exit $e
}

# Fills out the global opt() array
proc parse_opts {argv cmd opts} {
    global opt

    foreach {tag val def param help} $opts {
	regsub {.*\|} $tag {} tag
	set opt($tag) $def
    }

    while {[string match "-*" [lindex $argv]]} {
	set opcode [string range [lindex $argv 0] 1 end]

	set skip 0
	foreach {_tag val def param help} $opts {
	    set tlist [split $_tag "|"]
	    foreach tag $tlist {
		if {$opcode != $tag} continue

		if {$val == 0} {
		    set opt([lindex $tlist end]) 1
		    set skip 1
		    catch {::cmd::${cmd}::_[lindex $tlist 0] opt}
		    break
		} else {
		    if {[llength $argv] < 2} {
			puts stderr "Option '$opcode' requires an argument"
			exit 1
		    }
		    set opt([lindex $tlist end]) [lindex $argv 1]
		    set skip 2
		    catch {::cmd::${cmd}::_[lindex $tlist 0] opt}
		    break
		}
	    }
	}

	if {$skip == 0} {
	    puts stderr "Unknown option '$opcode'"
	    usage $cmd 1 $opts
	}

	set argv [lrange $argv $skip end]
    }

    return $argv
}

# General DB opening
proc db_open {db access} {
    if {[catch {set io [g5::open_database -name $db -access $access]} err]} {
	puts stderr "Couldn't open database '$db': $err"
	exit 1
    }

    if {$access == "rw" && [$io read_only]} {
	puts stderr "Exiting as failed to open with read/write access."
	exit 1
    }

    return $io
}

#-----------------------------------------------------------------------------
# COMMAND: export

namespace eval cmd::export {
    set name "Export sequences"
}

set ::cmd::export::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}

    contigs 1 {*}  {list} {Output only specific contigs. 'list' is a space separated list of contig names.}
    f|format  1 sam  {fmt}  {Controls the output format. 'fmt' should be one of sam, ace, baf, fasta or fastq.}
    o|out     1 "out.*" {filename} {Where to write output. Suffix defaults to 'fmt'. Output "-" indicates stdout.}
    fixmates 0 0   {}     {Attempts to fix up SAM mate pairing flags prior to output.}
}

proc ::cmd::export::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname ro]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set opt(out) [regsub {\*} $opt(out) $opt(format)]
    
    if {[catch {export_contigs \
		    -io $io \
		    -contigs $opt(contigs) \
		    -format $opt(format) \
		    -fixmates $opt(fixmates) \
		    -outfile $opt(out)} err]} {
	puts stderr "Failed in export_contigs call: $err"
    }

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: consensus
namespace eval cmd::consensus {
    set name "Save consensus"
}

set ::cmd::consensus::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}

    contigs    1 {*}      list     {Output only specific contigs. 'list' is a space separated list of contig names}
    f|format     1 fastq    {fmt}    {Controls the output format. 'fmt' should fasta or fastq.}
    strip_pads 0 0        {}       {Removes padding characters.}
    o|out        1 "cons.*" filename {Where to write output. Suffix defaults to 'fmt'.}
}

proc ::cmd::consensus::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname ro]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList $io]
    }

    set opt(out) [regsub {\*} $opt(out) $opt(format)]
    set format [lsearch {- fastq fasta exp} $opt(format)]

    if {[catch {get_consensus \
		    -io $io \
		    -contigs $opt(contigs) \
		    -format $format \
		    -strip_pads $opt(strip_pads) \
		    -outfile $opt(out)} err]} {
	puts stderr "Failed in get_consensus call: $err"
	$io close
	exit 1
    }

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: check
namespace eval cmd::check {
    set name "Run Check Database function"
}

set ::cmd::check::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}

    contigs    1 {*}      list     {Output only specific contigs. 'list' is a space separated list of contig names}
    l|level    1 2        num      {Set check level to '1' or '2'.}
    f|fix      0 0        {}       {Attempts to fix the database.                  *PLEASE BACK UP THE DB FIRST*}
    }

proc ::cmd::check::run {dbname _options} {
    upvar $_options opt

    if {$opt(fix)} {
	set acc rw
    } else {
	set acc ro
    }

    set io [db_open $dbname $acc]

    if {$opt(contigs) == "*"} {
	puts "=== checking entire DB ==="
	set err [$io check $opt(fix) $opt(level)]
    } else {
	foreach crec $opt(contigs) {
	    set c [$io get_contig $crec]
	    puts "=== checking contig #$crec ==="
	    incr err [$c check $opt(fix) $opt(level)]
	    $c delete
	}
    }

    if {$acc == "rw"} {
	puts "Saving changes"
	$io flush
	$io close
    }

    exit [expr {$err == 0 ? 0 : 1}]
}


#-----------------------------------------------------------------------------
# COMMAND: shuffle_pads
namespace eval cmd::shuffle_pads {
    set name "Locally realigns data"
}

set ::cmd::shuffle_pads::opts "
    h       0 0    {}     {Shows this help.}
    help    0 0    {}     {Shows this help.}
    {} {} {} {} {}

    contigs    1 {*}      list     {Output only specific contigs. 'list' is a space separated list of contig names}
    band_size  1 [keylget gap5_defs SHUFFLE_PADS.BAND_SIZE] {num} {Set the alignment band size.}
    soft_clips 1 [keylget gap5_defs SHUFFLE_PADS.SOFT_CLIPS] {num} {Adjust and align soft-clips/cutoff positions.}
"

proc ::cmd::shuffle_pads::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname rw]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    if {[catch {shuffle_pads \
		    -io $io \
		    -contigs $opt(contigs) \
		    -band $opt(band_size) \
	            -soft_clips $opt(soft_clips)} err]} {
	puts stderr "Failed in shuffle_pads call: $err"
	#$io close
	exit 1
    }

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: auto_break
namespace eval cmd::auto_break {
    set name "Find and automatically break apart suspect joins"
}

set ::cmd::auto_break::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Compare only specific contigs against each other. 'list' is a space separated list of contig names}
    n|dry_run     0 0   {}   {Find breaks only, but do not make or tag them}
    b|no_break    0 0   {}   {Find and tag breaks, but do not make them}
    {} {} {} {} {}    	  
    min_mqual     1 0   int  {Minimum mapping quality to consider in assembly. Used during problem finding. 0 => use all readings.}
    filter_consensus 1 1 int {Analyses the consensus sequence to find repeated sequence words. When disabled only single and dinucleotide repeats are filtered.}
    repeat_score  1 10  int  {Only used if filter_consensus is enabled. Words more than N times more frequently than expected as filtered out.}
    unique_mqual  1 20  int  {Minimum mapping quality used to define a pair as being "uniquely mapped"}
    good_unique   1 30  int  {Weight for all uniquely mapping good read pairs; >0}
    good          1 10  int  {Weight for all other good read pairs; >0}
    bad_unique    1 -25 int  {Weight for all uniquely mapping bad read pairs; <0}
    bad           1 -5  int  {Weight for all other bad read pairs; <0}
    large_unique  1 -20 int  {Weight for slightly poor uniquely mapping sized read pairs}
    large         1 -3  int   {Weight for slightly poor other sized read pairs}
    spanning_unique 1 -5 int  {Weight for contig-spanning unique read pairs; <= 0}
    spanning      1 0   int   {Weight for contig-spanning non-unique read pairs; <= 0}
    singleton_unique 1 -1 int {Weight for unique singletons that should be pairs; <=0}
    singleton     1 0    int  {Weight for non-unique singletons that should be pairs; <=0}
    min_score     1 -100 int  {Minimum combined score after applying the weights above; may be negative.}
    end_skip      1 1000 int {Skip problems within X bases of the contig ends}
}

proc lreverse l {set r "";foreach i $l {set r [linsert $r 0 $i]}; return $r}

proc ::cmd::auto_break::run {dbname _options} {
    upvar $_options opt

    if {$opt(dry_run)} {
	set io [db_open $dbname ro]
    } else {
	set io [db_open $dbname rw]
    }

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    # Auto_break does the analysis and returns a list of things to do,
    # rather than making the breaks itself.
    puts "\n>>>\n>>> Stage 1: finding regions to break apart\n>>>"

    set r [auto_break \
	       -io                      $io \
	       -contigs                 $opt(contigs) \
	       -filter_consensus        $opt(filter_consensus) \
	       -repeat_score	        $opt(repeat_score) \
	       -unique_mqual	        $opt(unique_mqual) \
	       -good_weight             $opt(good) \
	       -good_unique_weight      $opt(good_unique) \
	       -bad_unique_weight       $opt(bad_unique) \
	       -bad_weight              $opt(bad) \
	       -large_unique_weight     $opt(large_unique) \
	       -large_weight            $opt(large) \
	       -spanning_unique_weight  $opt(spanning_unique) \
	       -spanning_weight         $opt(spanning) \
	       -singleton_unique_weight $opt(singleton_unique) \
	       -singleton_weight        $opt(singleton) \
	       -min_score               $opt(min_score) \
	       -min_mqual               $opt(min_mqual) \
	       -end_skip                $opt(end_skip)]

    if {$opt(dry_run) || $opt(no_break)} {
	$io close
	return
    }

    puts "\n===\n=== Stage 2: Breaking\n==="

    # Process them in reverse order, as multiple holes in the same
    # contig need to be broken right to left.
    set r [lreverse $r]

    foreach hole $r {
	set crec [lindex $hole 0]
	set pos  [lindex $hole 1]
	set seqs [lrange $hole 2 end]

	puts "> Processing contig $crec at position $pos"

	if {$seqs != ""} {
	    # Disassemble first. Move 2 => move group to new contig
	    puts "Disassembling [llength $seqs] sequence(s) to new contig"
	    if {[catch {log_call disassemble_readings \
			    -io             $io \
			    -readings       $seqs \
			    -move           2 \
			    -remove_holes   0 \
			    -duplicate_tags 1} err]} {
		puts "Disassemble_readings produced error message: $err"
	    }
	}

	puts "Breaking contig #$crec at position $pos"
	if {[catch {log_call break_contig \
			-io          $io \
			-contig      $crec \
			-pos         $pos \
			-break_holes 0} err]} {
	    puts "Break_contig produced error message: $err"
	}
    }


    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: auto_join
namespace eval cmd::auto_join {
    set name "Find and automatically make joins"
}

set ::cmd::auto_join::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Compare only specific contigs against each other. 'list' is a space separated list of contig names}
    contigs1      1 {*} list {Compare set 1 against set 2. 'list' is a space separated list of contig names}
    contigs2      1 {*} list {Compare set 1 against set 2. 'list' is a space separated list of contig names}
    strict        0 0   {}   {Strict parameter set; filtered >=50bp @ 0% mismatch}
    lenient       0 0   {}   {Lenient parameter set; >=30bp @ 5% mismatch}
    n|dry_run     0 0   {}   {Find joins only, but do not make them}
    {} {} {} {} {}    	  
    min_overlap   1 50  val  {Minimum length of alignment}
    max_overlap   1 0   val  {Maximum length of alignment (0 => no maximum)}
    word_length   1 12  val  {Word length for hashing}
    min_match     1 50  val  {Minimum length of initial exact match}
    max_pmismatch 1 0.0 val  {Maximum percentage mismatch}
    rp_mode       1 off mode {Read-pair checking mode, one of off, end_end, end_all, all_all}
    rp_end_size   1 5000 val {Size of portion to treat as contig end}
    rp_min_mq     1 10  val  {Minimum mapping quality for read-pair}
    rp_min_freq   1 2   val  {Minimum number of spanning read-pairs}
    rp_min_perc   1 70  val  {Minimum percentage of spanning read-pairs}
    rp_libs       1 {*} list {Libraries to use for read-pair scanning, specify a blank list to see the choices available.}
    filter_words  1 5   N    {Filter repeat words occuring more than N times}
    fastest       0 0   {}   {Enable "fastest" mode}
    min_depth     1 0   val  {Minimum depth (0 => no minimum)}
    max_depth     1 0   val  {Maximum depth (0 => no maximum)}
    no_containments 0 0 {}   {Disallow containment joins}
    no_ends         0 0 {}   {Disallow overlap joins between contig ends}
    unique_ends   0 0   {}   {Only make joins between contig ends that match uniquely to one place.}
}

proc ::cmd::auto_join::_contigs {_options} {
    upvar $_options opt
    set opt(contigs1) $opt(contigs)
    set opt(contigs2) $opt(contigs)
}

proc ::cmd::auto_join::_lenient {_options} {
    upvar $_options opt
    set opt(min_overlap)   30
    set opt(min_match)     20
    set opt(max_pmismatch) 5.0
    set opt(rp_mode)       off
    set opt(filter_words)  10
}

proc ::cmd::auto_join::_strict {_options} {
    upvar $_options opt
    set opt(min_overlap)   50
    set opt(min_match)     50
    set opt(max_pmismatch) 0.0
    set opt(rp_mode)       end_end
    set opt(rp_min_freq)   2
    set opt(filter_words)  5
}

proc ::cmd::auto_join::run {dbname _options} {
    upvar $_options opt

    if {$opt(dry_run)} {
	set io [db_open $dbname ro]
    } else {
	set io [db_open $dbname rw]
    }

    if {$opt(rp_libs) == "*"} {
	set opt(rp_libs) ""; #Gap5 interprets that as all

    } elseif {$opt(rp_libs) == ""} {
	puts "List of available libraries, as #num or name:"
	set db [$io get_database]
	set nl [$db get_num_libraries]
	for {set i 0} {$i < $nl} {incr i} {
	    set rec [$db get_library_rec $i]
	    set lib [$io get_library $rec]
	    puts "\t#$rec\t[$lib get_name]"
	    $lib delete
	}
	exit 0

    } else {
	set libs ""

	set db [$io get_database]
	set nl [$db get_num_libraries]
	foreach lib_id $opt(rp_libs) {
	    if {[regexp {^#([0-9]+)} $lib_id _ rec]} {
		lappend libs $rec
	    } else {
		for {set i 0} {$i < $nl} {incr i} {
		    set rec [$db get_library_rec $i]
		    set lib [$io get_library $rec]
		    if {[string match [$lib get_name] $lib_id]} {
			lappend libs $rec
			$lib delete
			break;
		    }
		    $lib delete
		}
		if {$i == $nl} {
		    puts stderr "Library name \"$lib_id\" not found."
		    $io close
		    exit 1
		}
	    }
	}
	set opt(rp_libs) $libs
    }

    if {$opt(contigs1) == "*"} {
	set opt(contigs1) [CreateAllContigList=Numbers $io]
    }
    if {$opt(contigs2) == "*"} {
	set opt(contigs2) [CreateAllContigList=Numbers $io]
    }

    set fast_mode 0
    if {$opt(fastest)} {incr fast_mode}

    set id [find_internal_joins \
		-io            $io \
		-contigs1      $opt(contigs1) \
		-contigs2      $opt(contigs2) \
		-min_overlap   $opt(min_overlap) \
		-max_overlap   $opt(max_overlap) \
		-max_pmismatch $opt(max_pmismatch) \
		-word_length   $opt(word_length) \
		-min_match     $opt(min_match) \
		-filter_words  $opt(filter_words) \
		-fast_mode     $fast_mode \
		-rp_mode       $opt(rp_mode) \
		-rp_end_size   $opt(rp_end_size) \
		-rp_min_mq     $opt(rp_min_mq) \
		-rp_min_freq   $opt(rp_min_freq) \
		-rp_min_perc   $opt(rp_min_perc) \
		-rp_libraries  $opt(rp_libs) \
	        -min_depth     $opt(min_depth) \
	        -max_depth     $opt(max_depth) \
	        -containments  [expr {!$opt(no_containments)}] \
	        -ends          [expr {!$opt(no_ends)}] \
	        -unique_ends   $opt(unique_ends) \
		-use_conf      0 \
		-min_conf      0 \
		-use_hidden    0]

    if {!$opt(dry_run)} {
	result_notify \
	    -io $io \
	    -id $id \
	    -type GENERIC \
	    -args "{task TASK_AUTO_JOIN}"
    }

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: fij
namespace eval cmd::fij {
    set name "Find Internal Joins"
}

set ::cmd::fij::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Compare only specific contigs against each other. 'list' is a space separated list of contig names}
    contigs1      1 {*} list {Compare set 1 against set 2. 'list' is a space separated list of contig names}
    contigs2      1 {*} list {Compare set 1 against set 2. 'list' is a space separated list of contig names}
    strict        0 0   {}   {Strict parameter set; filtered >=50bp @ 0% mismatch}
    lenient       0 0   {}   {Lenient parameter set; >=30bp @ 5% mismatch}
    o|out         1 fij.out {fn} {Output filename for Contig Comparator}
    {} {} {} {} {}    	  
    min_overlap   1 50  val  {Minimum length of alignment}
    max_overlap   1 0   val  {Maximum length of alignment (0 => no maximum)}
    word_length   1 12  val  {Word length for hashing}
    min_match     1 50  val  {Minimum length of initial exact match}
    max_pmismatch 1 0.0 val  {Maximum percentage mismatch}
    rp_mode       1 off mode {Read-pair checking mode, one of off, end_end, end_all, all_all}
    rp_end_size   1 5000 val {Size of portion to treat as contig end}
    rp_min_mq     1 10  val  {Minimum mapping quality for read-pair}
    rp_min_freq   1 2   val  {Minimum number of spanning read-pairs}
    rp_min_perc   1 70  val  {Minimum percentage of spanning read-pairs}
    rp_libs       1 {*} list {Libraries to use for read-pair scanning, specify a blank list to see the choices available.}
    filter_words  1 5   N    {Filter repeat words occuring more than N times}
    fastest       0 0   {}   {Enable "fastest" mode}
    min_depth     1 -1  val  {Minimum median depth of combined contig (-1 => no minimum)}
    max_depth     1 -1  val  {Maximum median depth of combined contig (-1 => no maximum)}
    containments  1 1   bool {Allow containment joins}
    ends          1 1   bool {Allow overlap joins between contig ends}
}

proc ::cmd::fij::_contigs {_options} {
    upvar $_options opt
    set opt(contigs1) $opt(contigs)
    set opt(contigs2) $opt(contigs)
}

proc ::cmd::fij::_lenient {_options} {
    upvar $_options opt
    set opt(min_overlap)   30
    set opt(min_match)     20
    set opt(max_pmismatch) 5.0
    set opt(rp_mode)       off
    set opt(filter_words)  10
}

proc ::cmd::fij::_strict {_options} {
    upvar $_options opt
    set opt(min_overlap)   50
    set opt(min_match)     50
    set opt(max_pmismatch) 0.0
    set opt(rp_mode)       end_end
    set opt(rp_min_freq)   2
    set opt(filter_words)  5
}

proc ::cmd::fij::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname ro]

    if {$opt(rp_libs) == "*"} {
	set opt(rp_libs) ""; #Gap5 interprets that as all

    } elseif {$opt(rp_libs) == ""} {
	puts "List of available libraries, as #num or name:"
	set db [$io get_database]
	set nl [$db get_num_libraries]
	for {set i 0} {$i < $nl} {incr i} {
	    set rec [$db get_library_rec $i]
	    set lib [$io get_library $rec]
	    puts "\t#$rec\t[$lib get_name]"
	    $lib delete
	}
	exit 0

    } else {
	set libs ""

	set db [$io get_database]
	set nl [$db get_num_libraries]
	foreach lib_id $opt(rp_libs) {
	    if {[regexp {^#([0-9]+)} $lib_id _ rec]} {
		lappend libs $rec
	    } else {
		for {set i 0} {$i < $nl} {incr i} {
		    set rec [$db get_library_rec $i]
		    set lib [$io get_library $rec]
		    if {[string match [$lib get_name] $lib_id]} {
			lappend libs $rec
			$lib delete
			break;
		    }
		    $lib delete
		}
		if {$i == $nl} {
		    puts stderr "Library name \"$lib_id\" not found."
		    $io close
		    exit 1
		}
	    }
	}
	set opt(rp_libs) $libs
    }

    if {$opt(contigs1) == "*"} {
	set opt(contigs1) [CreateAllContigList=Numbers $io]
    }

    if {$opt(contigs2) == "*"} {
	set opt(contigs2) [CreateAllContigList=Numbers $io]
    }

    set fast_mode 0
    if {$opt(fastest)} {incr fast_mode}

    set id [find_internal_joins \
		-io            $io \
		-contigs1      $opt(contigs1) \
		-contigs2      $opt(contigs2) \
		-min_overlap   $opt(min_overlap) \
		-max_overlap   $opt(max_overlap) \
		-max_pmismatch $opt(max_pmismatch) \
		-word_length   $opt(word_length) \
		-min_match     $opt(min_match) \
		-filter_words  $opt(filter_words) \
		-fast_mode     $fast_mode \
		-rp_mode       $opt(rp_mode) \
		-rp_end_size   $opt(rp_end_size) \
		-rp_min_mq     $opt(rp_min_mq) \
		-rp_min_freq   $opt(rp_min_freq) \
		-rp_min_perc   $opt(rp_min_perc) \
		-rp_libraries  $opt(rp_libs) \
	        -min_depth     $opt(min_depth) \
	        -max_depth     $opt(max_depth) \
	        -containments  $opt(containments) \
	        -ends          $opt(ends) \
		-use_conf      0 \
		-min_conf      0 \
		-use_hidden    0]

    result_notify \
	-io $io \
	-id $id \
	-type GENERIC \
	-args "{task TASK_CS_SAVE data [list $opt(out)]}"

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: find_repeats
namespace eval cmd::find_repeats {
    set name "Find Repeats"
}

set ::cmd::find_repeats::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Output only specific contigs. 'list' is a space separated list of contig names}
    o|out         1 repeats.out {fn} {Output filename for Contig Comparator}
    min_match     1 50  val  {Minimum length repeat}
    direction     1 1   dir  {Direction; 1=fwd, 2=rev, 3=both}
    tag_file      1 {}  fn   {Create file of repeat tags}
}

proc ::cmd::find_repeats::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname rw]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set id [find_repeats \
		-io            $io \
		-contigs       $opt(contigs) \
		-min_match     $opt(min_match) \
		-outfile       $opt(tag_file) \
		-direction     $opt(direction)]

    result_notify \
	-io $io \
	-id $id \
	-type GENERIC \
	-args "{task TASK_CS_SAVE data [list $opt(out)]}"

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: find_read_pairs
namespace eval cmd::find_read_pairs {
    set name "Find Read Pairs"
}

set ::cmd::find_read_pairs::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Output only specific contigs. 'list' is a space separated list of contig names}
    o|out         1 read_pairs.out {fn} {Output filename for Contig Comparator}
    {} {} {} {} {}    	  
    mode         1 end_end mode {Read-pair checking mode, one of end_end, end_all, all_all}
    end_size     1 2000 val {Size of portion to treat as contig end}
    min_mq       1 10  val  {Minimum mapping quality for read-pair}
    min_freq     1 2   val  {Minimum number of spanning read-pairs}
    libraries    1 {}  list {Space separated list of library record numbers}
}

proc ::cmd::find_read_pairs::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname ro]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set id [find_read_pairs \
		-io            $io \
		-contigs       $opt(contigs) \
		-mode          $opt(mode) \
		-end_size      $opt(end_size) \
		-min_map_qual  $opt(min_mq) \
		-min_freq      $opt(min_freq) \
		-libraries     $opt(libraries)]

    result_notify \
	-io $io \
	-id $id \
	-type GENERIC \
	-args "{task TASK_CS_SAVE data [list $opt(out)]}"

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: find_oligos
namespace eval cmd::find_oligos {
    set name "Find Oligos (Sequence Search)"
}

set ::cmd::find_oligos::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Output only specific contigs. 'list' is a space separated list of contig names}
    o|out         1 oligos.out  {fn} {Output filename for Contig Comparator}
    min_pmatch 1  100.0 val  {Minimum percentage match}
    seq           1 {}  str  {DNA sequence to search for}
    c|consensus_only 0 0 {} {Only search for matches in the consensus}
    cutoffs       1 0   val  {Also look for matches in clipped/cutoff sequence}
    rp_min_freq   1 2   val  {Minimum number of spanning read-pairs}
    seq_file      1 {}  fn   {Use 'fn' as a file of sequences to search for}
}

proc ::cmd::find_oligos::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname rw]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set id [find_oligo \
		-io             $io \
		-contigs        $opt(contigs) \
		-seq            $opt(seq) \
		-consensus_only $opt(consensus_only) \
		-cutoffs        $opt(cutoffs) \
		-file           $opt(seq_file) \
		-min_pmatch     $opt(min_pmatch)]

    result_notify \
	-io $io \
	-id $id \
	-type GENERIC \
	-args "{task TASK_CS_SAVE data [list $opt(out)]}"

    $io close
}


#-----------------------------------------------------------------------------
# COMMAND: check_assembly
namespace eval cmd::check_assembly {
    set name "Check Assembly"
}

set ::cmd::check_assembly::opts {
    h|help  0 0    {}     {Shows this help.}
    {} {} {} {} {}
    
    contigs       1 {*} list {Output only specific contigs. 'list' is a space separated list of contig names}
    o|out         1 errs.out {fn} {Output filename for Contig Comparator}
    {} {} {} {} {}    	  
    max_pmismatch 1 15.0 val  {Maximum percentage mismatch}
    win_size      1 30   val  {Window size over which to look for mismatches}
    N|ignore_N    0 0    {}   {Treat N as a matching base}
}

proc ::cmd::check_assembly::run {dbname _options} {
    upvar $_options opt
    set io [db_open $dbname rw]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set id [check_assembly \
		-io            $io \
		-contigs       $opt(contigs) \
		-max_pmismatch $opt(max_pmismatch) \
		-win_size      $opt(win_size) \
		-ignore_N      $opt(ignore_N)]

    result_notify \
	-io $io \
	-id $id \
	-type GENERIC \
	-args "{task TASK_CS_SAVE data [list $opt(out)]}"

    $io close
}


#-----------------------------------------------------------------------------
# Main entry

set cmd [lindex $argv 0]
set argv [lrange $argv 1 end]

if {$cmd == "" || ![namespace exists ::cmd::$cmd]} {
    if {$cmd != ""} {
	puts "Unknown sub-command '$cmd'\n"
    }
    puts "Usage: gap5_cmd <command> \[options\]\n"

    set n 0
    foreach ns [lsort [namespace children ::cmd]] {
	set h [set ${ns}::name]
	regsub {.*::} $ns {} ns
	puts [format "%8s %-15s %s" [expr {$n?"":"Command:"}] $ns $h]
	incr n
    }

    exit 1
}

set argv [parse_opts $argv $cmd [set ::cmd::${cmd}::opts]]

if {$opt(help)} {
    usage $cmd 0 [set ::cmd::${cmd}::opts]
}

if {[llength $argv] != 1} {
    usage $cmd 1 [set ::cmd::${cmd}::opts]
}

set dbname [lindex $argv 0]

::cmd::${cmd}::run [lindex $argv 0] opt
