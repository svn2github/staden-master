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
		    -band $opt(band_size)} err]} {
	puts stderr "Failed in shuffle_pads call: $err"
	#$io close
	exit 1
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
    
    contigs       1 {*} list {Output only specific contigs. 'list' is a space separated list of contig names}
    strict        0 0   {}   {Strict parameter set; filtered >=50bp @ 0% mismatch}
    lenient       0 0   {}   {Lenient parameter set; >=30bp @ 5% mismatch}
    n|dry_run     0 0   {}   {Find joins only, but do not make them}
    {} {} {} {} {}    	  
    min_overlap   1 50  val  {Minimum length of alignment}
    word_length   1 12  val  {Word length for hashing}
    min_match     1 50  val  {Minimum length of initial exact match}
    max_pmismatch 1 0.0 val  {Maximum percentage mismatch}
    rp_mode       1 off mode {Read-pair checking mode, one of off, end_end, end_all, all_all}
    rp_end_size   1 5000 val {Size of portion to treat as contig end}
    rp_min_mq     1 10  val  {Minimum mapping quality for read-pair}
    rp_min_freq   1 2   val  {Minimum number of spanning read-pairs}
    filter_words  1 5   N    {Filter repeat words occuring more than N times}
    fastest       0 0   {}   {Enable "fastest" mode}
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
    set io [db_open $dbname rw]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList=Numbers $io]
    }

    set fast_mode 1
    if {$opt(fastest)} {incr fast_mode}

    set id [find_internal_joins \
		-io            $io \
		-contigs1      $opt(contigs) \
		-contigs2      $opt(contigs) \
		-min_overlap   $opt(min_overlap) \
		-max_pmismatch $opt(max_pmismatch) \
		-word_length   $opt(word_length) \
		-min_match     $opt(min_match) \
		-filter_words  $opt(filter_words) \
		-fast_mode     $fast_mode \
		-rp_mode       $opt(rp_mode) \
		-rp_end_size   $opt(rp_end_size) \
		-rp_min_mq     $opt(rp_min_mq) \
		-rp_min_freq   $opt(rp_min_freq) \
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
