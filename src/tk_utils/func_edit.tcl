# Allows editing of Tcl functions from within the console
proc func_edit {func} {
    if {[catch {info arg $func} err]} {
	tk_messageBox -message $err
	return
    }

    set t "[tmpnam].tcl"
    set fd [open $t w]
    puts $fd "proc $func [list [info arg $func]] {[info body $func]}"
    close $fd

    #exec xterm -e emacs -nw $t
    exec xemacs $t

    set fd [open $t r]
    set code [read $fd]
    close $fd
    file delete $t

    if {[catch {uplevel #0 $code} err]} {
	tk_messageBox -message $err
	return
    }
}