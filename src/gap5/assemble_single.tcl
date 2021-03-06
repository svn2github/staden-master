# GUI to import fasta/fastq reads as single-read contigs.
# This is only desireable for *small* numbers of reads, but in that case
# it offers an easy mechanism for getting consensus or finishing reads
# assembled so we can then run find internal joins.
#
# More generally users should be using the mapped reads assembly method
# instead.
proc AssemblySingle {io} {
    global gap5_defs

    set w .assembly_single
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Assemble Fasta/Fastq"

    # Some basic help
    label $w.label -wraplength 300 -justify left \
	-text "This function should only be used sparingly, for small numbers of reads.  The reads will not be aligned against existing contigs, instead creating one contig per read in the fasta/fastq file.  Consider using Map Reads instead if you desire a more usual assembly process."

    frame $w.sep -height 2 -bd 1 -relief groove

    # Input format
    radiolist $w.format \
	-title "Input readings from" \
	-orient horizontal \
	-buttons {fofn fasta fastq} \
	-default [keylget gap5_defs ASSEMBLE_SINGLE.FORMAT]

    # fasta/q input file name
    getFname $w.file "File name" load
    
    okcancelhelp $w.ok \
	-ok_command "AssemblySingle2 $io $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Assembly-Fasta}"

    pack $w.label -side top -fill both -expand 1
    pack $w.sep -side top -fill both -expand 1 -padx 10 -pady 10
    pack $w.format $w.file $w.ok -side top -fill both -expand 1
    
}

proc AssemblySingle2 {io w} {
    set fn [getFname_in_name $w.file]
    if {$fn == ""} {
	bell
	return
    }

    set prefix [tmpnam]

    if {[radiolist_get $w.format] == 1} {
	generate_fastq $fn $prefix.fastq
	set fn $prefix.fastq
	set format Q
    } else {
	set format [lindex {- - F Q} [radiolist_get $w.format]]
    }

    # Import
    if {![quit_displays -io $io -msg "assemble_fasta"]} {
	# Someone's too busy to shutdown?
	bell
	return
    }

    destroy $w
    SetBusy

    if { [ catch { log_call import_reads \
		       -io $io \
		       -append 1 \
		       -file $fn \
		       -format $format \
		       -index_names 1 } ]
     } {
	verror ERR_WARN "AssemblySingle" "Import failed"
	tk_messageBox -icon error -type ok -title "Import failed" \
		-message "Error detected while importing reads."
	$io flush
	ClearBusy
	AssemblySinglePostLoad $io
	return
    }

    catch {glob file delete $prefix.fastq}

    $io flush
    ClearBusy
    AssemblySinglePostLoad $io
}

proc AssemblySinglePostLoad {io} {
    global gap5_defs

    # Display contig selector if appropriate.
    set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
    global do_csel
    if {![winfo exists $cs_win] && $do_csel} {
	if {$do_csel == 2 || [db_info num_contigs $io] <= 1000} {
	    ContigSelector $io
	} else {
	    vmessage "\nSkipping contig selector due to large number of contigs."
	}
    } else {
	ContigInitReg $io
	catch {raise $cs_win}
    }

    # We may be going from 0 contigs to N contigs, so "un-grey" menus if so.
    ActivateMenu_Open 
    Menu_Check_RO $io

    # Check CurContig/LREG/RREG if not already done
    global CurContig
    if {![info exists CurContig] || $CurContig == ""} {
	InitContigGlobals $io
    }
}