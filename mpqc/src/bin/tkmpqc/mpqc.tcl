#!/home/midway/SUN4/bin/wish -f
# Program: mpqc
# Tcl version: 7.0 (Tcl/Tk/XF)
# Tk version: 3.3
# XF version: 2.2
#

# module inclusion
global env
global xfLoadPath
global xfLoadInfo
set xfLoadInfo 0
if {[info exists env(XF_LOAD_PATH)]} {
  if {[string first $env(XF_LOAD_PATH) /usr/local/lib/] == -1} {
    set xfLoadPath $env(XF_LOAD_PATH):/usr/local/lib/
  } {
    set xfLoadPath /usr/local/lib/
  }
} {
  set xfLoadPath /usr/local/lib/
}

global argc
global argv
global tkVersion
set tmpArgv ""
for {set counter 0} {$counter < $argc} {incr counter 1} {
  case [string tolower [lindex $argv $counter]] in {
    {-xfloadpath} {
      incr counter 1
      set xfLoadPath "[lindex $argv $counter]:$xfLoadPath"
    }
    {-xfstartup} {
      incr counter 1
      source [lindex $argv $counter]
    }
    {-xfbindfile} {
      incr counter 1
      set env(XF_BIND_FILE) "[lindex $argv $counter]"
    }
    {-xfcolorfile} {
      incr counter 1
      set env(XF_COLOR_FILE) "[lindex $argv $counter]"
    }
    {-xfcursorfile} {
      incr counter 1
      set env(XF_CURSOR_FILE) "[lindex $argv $counter]"
    }
    {-xffontfile} {
      incr counter 1
      set env(XF_FONT_FILE) "[lindex $argv $counter]"
    }
    {-xfmodelmono} {
      if {$tkVersion >= 3.0} {
        tk colormodel . monochrome
      }
    }
    {-xfmodelcolor} {
      if {$tkVersion >= 3.0} {
        tk colormodel . color
      }
    }
    {-xfloading} {
      set xfLoadInfo 1
    }
    {-xfnoloading} {
      set xfLoadInfo 0
    }
    {default} {
      lappend tmpArgv [lindex $argv $counter]
    }
  }
}
set argv $tmpArgv
set argc [llength $tmpArgv]
unset counter
unset tmpArgv


# procedure to show window ShowWindow.input
proc ShowWindow.input { args} {
# xf ignore me 7

  # build widget .input
  if {"[info procs XFEdit]" != ""} {
    catch "XFDestroy .input"
  } {
    catch "destroy .input"
  }
  toplevel .input 

  # Window manager configurations
  global tkVersion
  wm positionfrom .input program
  wm sizefrom .input program
  wm maxsize .input 1000 900
  wm minsize .input 10 10
  wm title .input {mpqc input}


  # build widget .input.basis
  button .input.basis  -text {basis}

  # build widget .input.done
  button .input.done  -command {DestroyWindow.input}  -state {active}  -text {done}

  # build widget .input.force
  button .input.force  -text {force}

  # build widget .input.geometry
  button .input.geometry  -text {geometry}

  # build widget .input.intco
  button .input.intco  -text {intco}

  # build widget .input.misc
  button .input.misc  -text {misc}

  # build widget .input.scf
  button .input.scf  -text {scf}

  # build widget .input.write
  button .input.write  -text {write}

  # pack widget .input
  pack append .input  .input.scf {left frame center}  .input.force {left frame center}  .input.intco {left frame center}  .input.basis {left frame center}  .input.geometry {left frame center}  .input.misc {left frame center}  .input.write {left frame center}  .input.done {left frame center}

  if {"[info procs XFEdit]" != ""} {
    catch "XFMiscBindWidgetTree .input"
    after 2 "catch {XFEditSetShowWindows}"
  }
}

proc DestroyWindow.input {} {# xf ignore me 7
  if {"[info procs XFEdit]" != ""} {
    if {"[info commands .input]" != ""} {
      global xfShowWindow.input
      set xfShowWindow.input 0
      XFEditSetPath .
      after 2 "XFSaveAsProc .input; XFEditSetShowWindows"
    }
  } {
    catch "destroy .input"
    update
  }
}


# procedure to show window .
proc ShowWindow. {args} {# xf ignore me 7

  # Window manager configurations
  global tkVersion
  wm positionfrom . program
  wm sizefrom . program
  wm maxsize . 1000 900
  wm minsize . 683 408
  wm title . {tkmpqc 0.1}


  # build widget .inframe
  frame .inframe \
    -background {#b0d0ff} \
    -borderwidth {2}

  # build widget .inframe.lscroll
  scrollbar .inframe.lscroll \
    -background {#b0d0ff} \
    -command {.inframe.tframe.text yview} \
    -foreground {#7f7f7f} \
    -relief {sunken}

  # build widget .inframe.label
  label .inframe.label \
    -background {#b0d0ff} \
    -text {mpqc.in}

  # build widget .inframe.tframe
  frame .inframe.tframe \
    -borderwidth {2} \
    -relief {sunken}

  # build widget .inframe.tframe.text
  text .inframe.tframe.text \
    -relief {sunken} \
    -width {37} \
    -yscrollcommand {.inframe.lscroll set}

  # pack widget .inframe.tframe
  pack append .inframe.tframe \
    .inframe.tframe.text {top frame center expand filly}

  # pack widget .inframe
  pack append .inframe \
    .inframe.label {top frame s} \
    .inframe.lscroll {left frame center filly} \
    .inframe.tframe {top frame center expand filly}

  # build widget .mainops
  frame .mainops \
    -borderwidth {2}

  # build widget .mainops.input
  button .mainops.input \
    -background {#7f7f7f} \
    -command {ShowWindow.input} \
    -text {Input} \
    -width {14}

  # build widget .mainops.quit
  button .mainops.quit \
    -background {#7f7f7f} \
    -command {destroy .} \
    -text {Quit} \
    -width {14}

  # build widget .mainops.run
  button .mainops.run \
    -background {#7f7f7f} \
    -text {Run} \
    -width {14}

  # build widget .mainops.spreadsheet
  button .mainops.spreadsheet \
    -background {#7f7f7f} \
    -text {Spreadsheet} \
    -width {14}

  # build widget .mainops.submit
  button .mainops.submit \
    -background {#7f7f7f} \
    -text {Submit} \
    -width {14}

  # pack widget .mainops
  pack append .mainops \
    .mainops.input {top frame s pady 22} \
    .mainops.run {top frame center} \
    .mainops.submit {top frame center} \
    .mainops.spreadsheet {top frame center} \
    .mainops.quit {top frame center}

  # build widget .outframe
  frame .outframe \
    -borderwidth {2}

  # build widget .outframe.lscroll
  scrollbar .outframe.lscroll \
    -command {.outframe.tframe.box yview} \
    -relief {sunken}

  # build widget .outframe.label
  label .outframe.label \
    -text {mpqc output}

  # build widget .outframe.tframe
  frame .outframe.tframe \
    -borderwidth {2} \
    -relief {sunken}

  # build widget .outframe.tframe.bscroll
  scrollbar .outframe.tframe.bscroll \
    -command {.outframe.tframe.box xview} \
    -orient {horizontal} \
    -relief {sunken}

  # build widget .outframe.tframe.box
  listbox .outframe.tframe.box \
    -font {-adobe-times-medium-r-normal--14-140-75-75-p-74-iso8859-1} \
    -geometry {10x2} \
    -xscrollcommand {.outframe.tframe.bscroll set} \
    -yscrollcommand {.outframe.lscroll set}

  # pack widget .outframe.tframe
  pack append .outframe.tframe \
    .outframe.tframe.box {top frame center expand fill} \
    .outframe.tframe.bscroll {top frame center fillx}

  # pack widget .outframe
  pack append .outframe \
    .outframe.label {top frame center} \
    .outframe.lscroll {left frame center filly} \
    .outframe.tframe {top frame center expand fill}

  # build widget .paths
  frame .paths \
    -borderwidth {2}

  # build widget .paths.expath
  frame .paths.expath \
    -borderwidth {2}

  # build widget .paths.expath.entry
  entry .paths.expath.entry \
    -relief {sunken} \
    -textvariable {mpqc_path} \
    -width {60}

  # build widget .paths.expath.label
  label .paths.expath.label \
    -text {Executable:}

  # pack widget .paths.expath
  pack append .paths.expath \
    .paths.expath.label {left frame e padx 44} \
    .paths.expath.entry {left frame center expand fillx}

  # build widget .paths.mpqcin
  frame .paths.mpqcin \
    -borderwidth {2}

  # build widget .paths.mpqcin.entry
  entry .paths.mpqcin.entry \
    -relief {sunken} \
    -textvariable {input_dir} \
    -width {60}

  # build widget .paths.mpqcin.label
  label .paths.mpqcin.label \
    -text {Input Directory:}

  # pack widget .paths.mpqcin
  pack append .paths.mpqcin \
    .paths.mpqcin.label {left frame e padx 19} \
    .paths.mpqcin.entry {right frame nw expand fillx}

  # build widget .paths.output
  frame .paths.output \
    -borderwidth {2}

  # build widget .paths.output.entry
  entry .paths.output.entry \
    -relief {sunken} \
    -textvariable {outfile}

  # build widget .paths.output.label
  label .paths.output.label \
    -text {Output File:}

  # pack widget .paths.output
  pack append .paths.output \
    .paths.output.label {left frame e padx 44} \
    .paths.output.entry {top frame center fillx}

  # pack widget .paths
  pack append .paths \
    .paths.mpqcin {top frame e expand fillx} \
    .paths.output {top frame e expand fillx} \
    .paths.expath {top frame e expand fillx}

  # build widget .title
  frame .title \
    -borderwidth {2}

  # build widget .title.label
  label .title.label \
    -background {#b0d0ff} \
    -relief {sunken} \
    -text {TkMPQC 0.1}

  # build widget .title.status
  message .title.status \
    -aspect {1500} \
    -padx {5} \
    -pady {2} \
    -text {Welcome to MPQC (or is it SC?)} \
    -textvariable {statusmsg}

  # pack widget .title
  pack append .title \
    .title.label {left frame e padx 22 fillx} \
    .title.status {right frame w expand}

  # pack widget .
  pack append . \
    .title {top frame w pady 11 fillx} \
    .paths {top frame w pady 20 fillx} \
    .mainops {left frame n padx 13} \
    .inframe {left frame center filly} \
    .outframe {right frame center expand fill}

  .inframe.tframe.text insert end {}



  if {"[info procs XFEdit]" != ""} {
    catch "XFMiscBindWidgetTree .xfEdit"
    after 2 "catch {XFEditSetShowWindows}"
  }
}


# User defined procedures


# Internal procedures



# application parsing procedure
proc XFLocalParseAppDefs {xfAppDefFile} {
  global xfAppDefaults

  # basically from: Michael Moore
  if {[file exists $xfAppDefFile] &&
      [file readable $xfAppDefFile] &&
      "[file type $xfAppDefFile]" == "link"} {
    catch "file type $xfAppDefFile" xfType
    while {"$xfType" == "link"} {
      if {[catch "file readlink $xfAppDefFile" xfAppDefFile]} {
        return
      }
      catch "file type $xfAppDefFile" xfType
    }
  }
  if {!("$xfAppDefFile" != "" &&
        [file exists $xfAppDefFile] &&
        [file readable $xfAppDefFile] &&
        "[file type $xfAppDefFile]" == "file")} {
    return
  }
  if {![catch "open $xfAppDefFile r" xfResult]} {
    set xfAppFileContents [read $xfResult]
    close $xfResult
    foreach line [split $xfAppFileContents "\n"] {
      # backup indicates how far to backup.  It applies to the
      # situation where a resource name ends in . and when it
      # ends in *.  In the second case you want to keep the *
      # in the widget name for pattern matching, but you want
      # to get rid of the . if it is the end of the name. 
      set backup -2  
      set line [string trim $line]
      if {[string index $line 0] == "#" || "$line" == ""} {
        # skip comments and empty lines
        continue
      }
      set list [split $line ":"]
      set resource [string trim [lindex $list 0]]
      set i [string last "." $resource]
      set j [string last "*" $resource]
      if {$j > $i} { 
        set i $j
        set backup -1
      }
      incr i
      set name [string range $resource $i end]
      incr i $backup
      set widname [string range $resource 0 $i]
      set value [string trim [lindex $list 1]]
      if {"$widname" != "" && "$widname" != "*"} {
        # insert the widget and resourcename to the application
        # defaults list.
        set xfAppDefaults($widname:[string tolower $name]) $value
      }
    }
  }
}

# application loading procedure
proc XFLocalLoadAppDefs {xfClasses {xfPriority "startupFile"} {xfAppDefFile ""}} {
  global env

  if {"$xfAppDefFile" == ""} {
    set xfFileList ""
    if {[info exists env(XUSERFILESEARCHPATH)]} {
      append xfFileList [split $env(XUSERFILESEARCHPATH) :]
    }
    if {[info exists env(XAPPLRESDIR)]} {
      append xfFileList [split $env(XAPPLRESDIR) :]
    }
    if {[info exists env(XFILESEARCHPATH)]} {
      append xfFileList [split $env(XFILESEARCHPATH) :]
    }
    append xfFileList " /usr/lib/X11/app-defaults"
    append xfFileList " /usr/X11/lib/X11/app-defaults"

    foreach xfCounter1 $xfClasses {
      foreach xfCounter2 $xfFileList {
        set xfPathName $xfCounter2
        if {[regsub -all "%N" "$xfPathName" "$xfCounter1" xfResult]} {
          set xfPathName $xfResult
        }
        if {[regsub -all "%T" "$xfPathName" "app-defaults" xfResult]} {
          set xfPathName $xfResult
        }
        if {[regsub -all "%S" "$xfPathName" "" xfResult]} {
          set xfPathName $xfResult
        }
        if {[regsub -all "%C" "$xfPathName" "" xfResult]} {
          set xfPathName $xfResult
        }
        if {[file exists $xfPathName] &&
            [file readable $xfPathName] &&
            ("[file type $xfPathName]" == "file" ||
             "[file type $xfPathName]" == "link")} {
          catch "option readfile $xfPathName $xfPriority"
          if {"[info commands XFParseAppDefs]" != ""} {
            XFParseAppDefs $xfPathName
          } {
            if {"[info commands XFLocalParseAppDefs]" != ""} {
              XFLocalParseAppDefs $xfPathName
            }
          }
        } {
          if {[file exists $xfCounter2/$xfCounter1] &&
              [file readable $xfCounter2/$xfCounter1] &&
              ("[file type $xfCounter2/$xfCounter1]" == "file" ||
               "[file type $xfCounter2/$xfCounter1]" == "link")} {
            catch "option readfile $xfCounter2/$xfCounter1 $xfPriority"
            if {"[info commands XFParseAppDefs]" != ""} {
              XFParseAppDefs $xfCounter2/$xfCounter1
            } {
              if {"[info commands XFLocalParseAppDefs]" != ""} {
                XFLocalParseAppDefs $xfCounter2/$xfCounter1
              }
            }
          }
        }
      }
    }
  } {
    # load a specific application defaults file
    if {[file exists $xfAppDefFile] &&
        [file readable $xfAppDefFile] &&
        ("[file type $xfAppDefFile]" == "file" ||
         "[file type $xfAppDefFile]" == "link")} {
      catch "option readfile $xfAppDefFile $xfPriority"
      if {"[info commands XFParseAppDefs]" != ""} {
        XFParseAppDefs $xfAppDefFile
      } {
        if {"[info commands XFLocalParseAppDefs]" != ""} {
          XFLocalParseAppDefs $xfAppDefFile
        }
      }
    }
  }
}

# application setting procedure
proc XFLocalSetAppDefs {{xfWidgetPath "."}} {
  global xfAppDefaults

  if {![info exists xfAppDefaults]} {
    return
  }
  foreach xfCounter [array names xfAppDefaults] {
    if {[string match "${xfWidgetPath}*" $xfCounter]} {
      set widname [string range $xfCounter 0 [expr [string first : $xfCounter]-1]]
      set name [string range $xfCounter [expr [string first : $xfCounter]+1] end]
      # Now lets see how many tcl commands match the name
      # pattern specified.
      set widlist [info command $widname]
      if {"$widlist" != ""} {
        foreach widget $widlist {
          # make sure this command is a widget.
          if {![catch "winfo id $widget"]} {
            catch "$widget configure -[string tolower $name] $xfAppDefaults($xfCounter)" 
          }
        }
      }
    }
  }
}



# end source
proc EndSrc {} {
  global input_dir mpqc_path statusmsg

  set input_dir [exec pwd]
  set mpqc_path "mpqcic"

  set statusmsg "Welcome to MPQC (or is it SC?)"
}

# prepare auto loading
global auto_path
global tk_library
global xfLoadPath
foreach xfElement [eval list [split $xfLoadPath :] $auto_path] {
  if {[file exists $xfElement/tclIndex]} {
    lappend auto_path $xfElement
  }
}
catch "unset auto_index"

catch "unset auto_oldpath"

catch "unset auto_execs"


# initialize global variables
proc InitGlobals {} {
  global {fontInFile}
  set {fontInFile} {file4}
  global {fontLine}
  set {fontLine} {}
  global {fontReadList}
  set {fontReadList} {--symbol-medium-r-normal--0-0-0-0-p-0--symbol
-adobe-courier-bold-o-normal--10-100-75-75-m-60-iso8859-1
-adobe-courier-bold-o-normal--12-120-75-75-m-70-iso8859-1
-adobe-courier-bold-o-normal--14-140-75-75-m-90-iso8859-1
-adobe-courier-bold-o-normal--18-180-75-75-m-110-iso8859-1
-adobe-courier-bold-o-normal--24-240-75-75-m-150-iso8859-1
-adobe-courier-bold-o-normal--8-80-75-75-m-50-iso8859-1
-adobe-courier-bold-r-normal--10-100-75-75-m-60-iso8859-1
-adobe-courier-bold-r-normal--12-120-75-75-m-70-iso8859-1
-adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1
-adobe-courier-bold-r-normal--18-180-75-75-m-110-iso8859-1
-adobe-courier-bold-r-normal--24-240-75-75-m-150-iso8859-1
-adobe-courier-bold-r-normal--8-80-75-75-m-50-iso8859-1
-adobe-courier-medium-o-normal--10-100-75-75-m-60-iso8859-1
-adobe-courier-medium-o-normal--12-120-75-75-m-70-iso8859-1
-adobe-courier-medium-o-normal--14-140-75-75-m-90-iso8859-1
-adobe-courier-medium-o-normal--18-180-75-75-m-110-iso8859-1
-adobe-courier-medium-o-normal--24-240-75-75-m-150-iso8859-1
-adobe-courier-medium-o-normal--8-80-75-75-m-50-iso8859-1
-adobe-courier-medium-r-normal--10-100-75-75-m-60-iso8859-1
-adobe-courier-medium-r-normal--12-120-75-75-m-70-iso8859-1
-adobe-courier-medium-r-normal--14-140-75-75-m-90-iso8859-1
-adobe-courier-medium-r-normal--18-180-75-75-m-110-iso8859-1
-adobe-courier-medium-r-normal--24-240-75-75-m-150-iso8859-1
-adobe-courier-medium-r-normal--8-80-75-75-m-50-iso8859-1
-adobe-helvetica-bold-o-normal--10-100-75-75-p-60-iso8859-1
-adobe-helvetica-bold-o-normal--12-120-75-75-p-69-iso8859-1
-adobe-helvetica-bold-o-normal--14-140-75-75-p-82-iso8859-1
-adobe-helvetica-bold-o-normal--18-180-75-75-p-104-iso8859-1
-adobe-helvetica-bold-o-normal--24-240-75-75-p-138-iso8859-1
-adobe-helvetica-bold-o-normal--8-80-75-75-p-50-iso8859-1
-adobe-helvetica-bold-r-normal--10-100-75-75-p-60-iso8859-1
-adobe-helvetica-bold-r-normal--12-120-75-75-p-70-iso8859-1
-adobe-helvetica-bold-r-normal--14-140-75-75-p-82-iso8859-1
-adobe-helvetica-bold-r-normal--18-180-75-75-p-103-iso8859-1
-adobe-helvetica-bold-r-normal--24-240-75-75-p-138-iso8859-1
-adobe-helvetica-bold-r-normal--8-80-75-75-p-50-iso8859-1
-adobe-helvetica-medium-o-normal--10-100-75-75-p-57-iso8859-1
-adobe-helvetica-medium-o-normal--12-120-75-75-p-67-iso8859-1
-adobe-helvetica-medium-o-normal--14-140-75-75-p-78-iso8859-1
-adobe-helvetica-medium-o-normal--18-180-75-75-p-98-iso8859-1
-adobe-helvetica-medium-o-normal--24-240-75-75-p-130-iso8859-1
-adobe-helvetica-medium-o-normal--8-80-75-75-p-47-iso8859-1
-adobe-helvetica-medium-r-normal--10-100-75-75-p-56-iso8859-1
-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1
-adobe-helvetica-medium-r-normal--14-140-75-75-p-77-iso8859-1
-adobe-helvetica-medium-r-normal--18-180-75-75-p-98-iso8859-1
-adobe-helvetica-medium-r-normal--24-240-75-75-p-130-iso8859-1
-adobe-helvetica-medium-r-normal--8-80-75-75-p-46-iso8859-1
-adobe-new century schoolbook-bold-i-normal--10-100-75-75-p-66-iso8859-1
-adobe-new century schoolbook-bold-i-normal--12-120-75-75-p-76-iso8859-1
-adobe-new century schoolbook-bold-i-normal--14-140-75-75-p-88-iso8859-1
-adobe-new century schoolbook-bold-i-normal--18-180-75-75-p-111-iso8859-1
-adobe-new century schoolbook-bold-i-normal--24-240-75-75-p-148-iso8859-1
-adobe-new century schoolbook-bold-i-normal--8-80-75-75-p-56-iso8859-1
-adobe-new century schoolbook-bold-r-normal--10-100-75-75-p-66-iso8859-1
-adobe-new century schoolbook-bold-r-normal--12-120-75-75-p-77-iso8859-1
-adobe-new century schoolbook-bold-r-normal--14-140-75-75-p-87-iso8859-1
-adobe-new century schoolbook-bold-r-normal--18-180-75-75-p-113-iso8859-1
-adobe-new century schoolbook-bold-r-normal--24-240-75-75-p-149-iso8859-1
-adobe-new century schoolbook-bold-r-normal--8-80-75-75-p-56-iso8859-1
-adobe-new century schoolbook-medium-i-normal--10-100-75-75-p-60-iso8859-1
-adobe-new century schoolbook-medium-i-normal--12-120-75-75-p-70-iso8859-1
-adobe-new century schoolbook-medium-i-normal--14-140-75-75-p-81-iso8859-1
-adobe-new century schoolbook-medium-i-normal--18-180-75-75-p-104-iso8859-1
-adobe-new century schoolbook-medium-i-normal--24-240-75-75-p-136-iso8859-1
-adobe-new century schoolbook-medium-i-normal--8-80-75-75-p-50-iso8859-1
-adobe-new century schoolbook-medium-r-normal--10-100-75-75-p-60-iso8859-1
-adobe-new century schoolbook-medium-r-normal--12-120-75-75-p-70-iso8859-1
-adobe-new century schoolbook-medium-r-normal--14-140-75-75-p-82-iso8859-1
-adobe-new century schoolbook-medium-r-normal--18-180-75-75-p-103-iso8859-1
-adobe-new century schoolbook-medium-r-normal--24-240-75-75-p-137-iso8859-1
-adobe-new century schoolbook-medium-r-normal--8-80-75-75-p-50-iso8859-1
-adobe-symbol-medium-r-normal--10-100-75-75-p-61-adobe-fontspecific
-adobe-symbol-medium-r-normal--12-120-75-75-p-74-adobe-fontspecific
-adobe-symbol-medium-r-normal--14-140-75-75-p-85-adobe-fontspecific
-adobe-symbol-medium-r-normal--18-180-75-75-p-107-adobe-fontspecific
-adobe-symbol-medium-r-normal--24-240-75-75-p-142-adobe-fontspecific
-adobe-symbol-medium-r-normal--8-80-75-75-p-51-adobe-fontspecific
-adobe-times-bold-i-normal--10-100-75-75-p-57-iso8859-1
-adobe-times-bold-i-normal--12-120-75-75-p-68-iso8859-1
-adobe-times-bold-i-normal--14-140-75-75-p-77-iso8859-1
-adobe-times-bold-i-normal--18-180-75-75-p-98-iso8859-1
-adobe-times-bold-i-normal--24-240-75-75-p-128-iso8859-1
-adobe-times-bold-i-normal--8-80-75-75-p-47-iso8859-1
-adobe-times-bold-r-normal--10-100-75-75-p-57-iso8859-1
-adobe-times-bold-r-normal--12-120-75-75-p-67-iso8859-1
-adobe-times-bold-r-normal--14-140-75-75-p-77-iso8859-1
-adobe-times-bold-r-normal--18-180-75-75-p-99-iso8859-1
-adobe-times-bold-r-normal--24-240-75-75-p-132-iso8859-1
-adobe-times-bold-r-normal--8-80-75-75-p-47-iso8859-1
-adobe-times-medium-i-normal--10-100-75-75-p-52-iso8859-1
-adobe-times-medium-i-normal--12-120-75-75-p-63-iso8859-1
-adobe-times-medium-i-normal--14-140-75-75-p-73-iso8859-1
-adobe-times-medium-i-normal--18-180-75-75-p-94-iso8859-1
-adobe-times-medium-i-normal--24-240-75-75-p-125-iso8859-1
-adobe-times-medium-i-normal--8-80-75-75-p-42-iso8859-1
-adobe-times-medium-r-normal--10-100-75-75-p-54-iso8859-1
-adobe-times-medium-r-normal--12-120-75-75-p-64-iso8859-1
-adobe-times-medium-r-normal--14-140-75-75-p-74-iso8859-1
-adobe-times-medium-r-normal--18-180-75-75-p-94-iso8859-1
-adobe-times-medium-r-normal--24-240-75-75-p-124-iso8859-1
-adobe-times-medium-r-normal--8-80-75-75-p-44-iso8859-1
-b&h-lucida bright-demibold-i-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucida bright-demibold-i-normal--10-100-72-72-p-60-iso8859-1
-b&h-lucida bright-demibold-i-normal--12-120-72-72-p-73-iso8859-1
-b&h-lucida bright-demibold-i-normal--14-140-72-72-p-84-iso8859-1
-b&h-lucida bright-demibold-i-normal--6-60-72-72-p-37-iso8859-1
-b&h-lucida bright-demibold-i-normal--8-80-72-72-p-48-iso8859-1
-b&h-lucida bright-demibold-r-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucida bright-demibold-r-normal--10-100-72-72-p-60-iso8859-1
-b&h-lucida bright-demibold-r-normal--12-120-72-72-p-71-iso8859-1
-b&h-lucida bright-demibold-r-normal--14-140-72-72-p-84-iso8859-1
-b&h-lucida bright-demibold-r-normal--6-60-72-72-p-37-iso8859-1
-b&h-lucida bright-demibold-r-normal--8-80-72-72-p-48-iso8859-1
-b&h-lucida bright-medium-i-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucida bright-medium-i-normal--10-100-72-72-p-57-iso8859-1
-b&h-lucida bright-medium-i-normal--12-120-72-72-p-68-iso8859-1
-b&h-lucida bright-medium-i-normal--14-140-72-72-p-81-iso8859-1
-b&h-lucida bright-medium-i-normal--6-60-72-72-p-35-iso8859-1
-b&h-lucida bright-medium-i-normal--8-80-72-72-p-45-iso8859-1
-b&h-lucida bright-medium-r-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucida bright-medium-r-normal--10-100-72-72-p-57-iso8859-1
-b&h-lucida bright-medium-r-normal--12-120-72-72-p-68-iso8859-1
-b&h-lucida bright-medium-r-normal--14-140-72-72-p-81-iso8859-1
-b&h-lucida bright-medium-r-normal--6-60-72-72-p-35-iso8859-1
-b&h-lucida bright-medium-r-normal--8-80-72-72-p-46-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-0-0-0-0-m-0-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-10-100-72-72-m-60-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-12-120-72-72-m-70-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-14-140-72-72-m-90-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-18-180-72-72-m-110-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-6-60-72-72-m-39-iso8859-1
-b&h-lucida sans typewriter-bold-r-normal-sans-8-80-72-72-m-50-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-0-0-0-0-m-0-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-10-100-72-72-m-60-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-12-120-72-72-m-70-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-14-140-72-72-m-90-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-18-180-72-72-m-110-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-6-60-72-72-m-40-iso8859-1
-b&h-lucida sans typewriter-medium-r-normal-sans-8-80-72-72-m-50-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-10-100-72-72-p-67-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-12-120-72-72-p-79-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-14-140-72-72-p-93-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-18-180-72-72-p-120-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-6-60-72-72-p-38-iso8859-1
-b&h-lucida sans-bold-i-normal-sans-8-80-72-72-p-49-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-10-100-72-72-p-67-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-12-120-72-72-p-79-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-14-140-72-72-p-92-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-18-180-72-72-p-120-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-6-60-72-72-p-38-iso8859-1
-b&h-lucida sans-bold-r-normal-sans-8-80-72-72-p-50-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-10-100-72-72-p-60-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-12-120-72-72-p-71-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-14-140-72-72-p-82-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-18-180-72-72-p-106-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-6-60-72-72-p-34-iso8859-1
-b&h-lucida sans-medium-i-normal-sans-8-80-72-72-p-46-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-10-100-72-72-p-59-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-12-120-72-72-p-71-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-14-140-72-72-p-82-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-18-180-72-72-p-106-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-6-60-72-72-p-35-iso8859-1
-b&h-lucida sans-medium-r-normal-sans-8-80-72-72-p-45-iso8859-1
-b&h-lucida-bold-i-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida-bold-r-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida-medium-i-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucida-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-b&h-lucidabright-demibold-i-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucidabright-demibold-r-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucidabright-medium-i-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucidabright-medium-r-normal--0-0-0-0-p-0-iso8859-1
-b&h-lucidatypewriter-bold-r-normal-sans-0-0-0-0-m-0-iso8859-1
-b&h-lucidatypewriter-medium-r-normal-sans-0-0-0-0-m-0-iso8859-1
-b&h-menu-medium-r-normal--12-120-75-75-p-75-iso8859-1
-bitstream-charter-bold-i-normal--10-100-75-75-p-62-iso8859-1
-bitstream-charter-bold-i-normal--12-120-75-75-p-74-iso8859-1
-bitstream-charter-bold-i-normal--15-140-75-75-p-93-iso8859-1
-bitstream-charter-bold-i-normal--19-180-75-75-p-117-iso8859-1
-bitstream-charter-bold-i-normal--25-240-75-75-p-154-iso8859-1
-bitstream-charter-bold-i-normal--8-80-75-75-p-50-iso8859-1
-bitstream-charter-bold-r-normal--10-100-75-75-p-63-iso8859-1
-bitstream-charter-bold-r-normal--12-120-75-75-p-75-iso8859-1
-bitstream-charter-bold-r-normal--15-140-75-75-p-94-iso8859-1
-bitstream-charter-bold-r-normal--19-180-75-75-p-119-iso8859-1
-bitstream-charter-bold-r-normal--25-240-75-75-p-157-iso8859-1
-bitstream-charter-bold-r-normal--8-80-75-75-p-50-iso8859-1
-bitstream-charter-medium-i-normal--10-100-75-75-p-55-iso8859-1
-bitstream-charter-medium-i-normal--12-120-75-75-p-65-iso8859-1
-bitstream-charter-medium-i-normal--15-140-75-75-p-82-iso8859-1
-bitstream-charter-medium-i-normal--19-180-75-75-p-103-iso8859-1
-bitstream-charter-medium-i-normal--25-240-75-75-p-136-iso8859-1
-bitstream-charter-medium-i-normal--8-80-75-75-p-44-iso8859-1
-bitstream-charter-medium-r-normal--10-100-75-75-p-56-iso8859-1
-bitstream-charter-medium-r-normal--12-120-75-75-p-67-iso8859-1
-bitstream-charter-medium-r-normal--15-140-75-75-p-84-iso8859-1
-bitstream-charter-medium-r-normal--19-180-75-75-p-106-iso8859-1
-bitstream-charter-medium-r-normal--25-240-75-75-p-139-iso8859-1
-bitstream-charter-medium-r-normal--8-80-75-75-p-45-iso8859-1
-dec-menu-medium-r-normal--12-100-75-75-p-75-iso8859-1
-dec-menu-medium-r-normal--12-120-75-75-p-75-iso8859-1
-dec-terminal-bold-r-normal--14-140-75-75-c-80-dec-dectech
-dec-terminal-bold-r-normal--14-140-75-75-c-80-iso8859-1
-dec-terminal-medium-r-narrow--12-140-75-75-c-80-iso8859-1
-dec-terminal-medium-r-normal--14-140-75-75-c-80-dec-dectech
-dec-terminal-medium-r-normal--14-140-75-75-c-80-iso8859-1
-itc-avantgarde-book-o-normal--0-0-0-0-p-0-iso8859-1
-itc-avantgarde-book-r-normal--0-0-0-0-p-0-iso8859-1
-itc-avantgarde-demi-o-normal--0-0-0-0-p-0-iso8859-1
-itc-avantgarde-demi-r-normal--0-0-0-0-p-0-iso8859-1
-itc-bookman-demi-i-normal--0-0-0-0-p-0-iso8859-1
-itc-bookman-demi-r-normal--0-0-0-0-p-0-iso8859-1
-itc-bookman-light-i-normal--0-0-0-0-p-0-iso8859-1
-itc-bookman-light-r-normal--0-0-0-0-p-0-iso8859-1
-itc-courier-bold-o-normal--0-0-0-0-m-0-iso8859-1
-itc-courier-bold-r-normal--0-0-0-0-m-0-iso8859-1
-itc-courier-medium-o-normal--0-0-0-0-m-0-iso8859-1
-itc-courier-medium-r-normal--0-0-0-0-m-0-iso8859-1
-itc-zapfchancery-medium-i-normal--0-0-0-0-p-0-iso8859-1
-itc-zapfdingbats-medium-r-normal--0-0-0-0-p-0--dingbats
-linotype-helvetica-bold-o-narrow--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-10-100-72-72-p-46-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-12-120-72-72-p-54-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-14-140-72-72-p-63-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-6-60-72-72-p-28-iso8859-1
-linotype-helvetica-bold-o-narrow-sans-8-80-72-72-p-37-iso8859-1
-linotype-helvetica-bold-o-normal--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-o-normal-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-o-normal-sans-10-100-72-72-p-57-iso8859-1
-linotype-helvetica-bold-o-normal-sans-12-120-72-72-p-67-iso8859-1
-linotype-helvetica-bold-o-normal-sans-14-140-72-72-p-79-iso8859-1
-linotype-helvetica-bold-o-normal-sans-6-60-72-72-p-34-iso8859-1
-linotype-helvetica-bold-o-normal-sans-8-80-72-72-p-44-iso8859-1
-linotype-helvetica-bold-r-narrow--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-10-100-72-72-p-46-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-12-120-72-72-p-54-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-14-140-72-72-p-63-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-6-60-72-72-p-28-iso8859-1
-linotype-helvetica-bold-r-narrow-sans-8-80-72-72-p-37-iso8859-1
-linotype-helvetica-bold-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-r-normal-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-bold-r-normal-sans-10-100-72-72-p-57-iso8859-1
-linotype-helvetica-bold-r-normal-sans-12-120-72-72-p-67-iso8859-1
-linotype-helvetica-bold-r-normal-sans-14-140-72-72-p-79-iso8859-1
-linotype-helvetica-bold-r-normal-sans-6-60-72-72-p-34-iso8859-1
-linotype-helvetica-bold-r-normal-sans-8-80-72-72-p-44-iso8859-1
-linotype-helvetica-medium-o-narrow--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-10-100-72-72-p-45-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-12-120-72-72-p-52-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-14-140-72-72-p-61-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-6-60-72-72-p-27-iso8859-1
-linotype-helvetica-medium-o-narrow-sans-8-80-72-72-p-36-iso8859-1
-linotype-helvetica-medium-o-normal--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-o-normal-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-o-normal-sans-10-100-72-72-p-55-iso8859-1
-linotype-helvetica-medium-o-normal-sans-12-120-72-72-p-65-iso8859-1
-linotype-helvetica-medium-o-normal-sans-14-140-72-72-p-76-iso8859-1
-linotype-helvetica-medium-o-normal-sans-6-60-72-72-p-32-iso8859-1
-linotype-helvetica-medium-o-normal-sans-8-80-72-72-p-42-iso8859-1
-linotype-helvetica-medium-r-narrow--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-10-100-72-72-p-45-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-12-120-72-72-p-52-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-14-140-72-72-p-61-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-6-60-72-72-p-27-iso8859-1
-linotype-helvetica-medium-r-narrow-sans-8-80-72-72-p-36-iso8859-1
-linotype-helvetica-medium-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-linotype-helvetica-medium-r-normal-sans-10-100-72-72-p-55-iso8859-1
-linotype-helvetica-medium-r-normal-sans-12-120-72-72-p-65-iso8859-1
-linotype-helvetica-medium-r-normal-sans-14-140-72-72-p-76-iso8859-1
-linotype-helvetica-medium-r-normal-sans-6-60-72-72-p-32-iso8859-1
-linotype-helvetica-medium-r-normal-sans-8-80-72-72-p-42-iso8859-1
-linotype-new century schoolbook-bold-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-new century schoolbook-bold-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-new century schoolbook-bold-r-normal--10-100-72-72-p-60-iso8859-1
-linotype-new century schoolbook-bold-r-normal--12-120-72-72-p-71-iso8859-1
-linotype-new century schoolbook-bold-r-normal--14-140-72-72-p-84-iso8859-1
-linotype-new century schoolbook-bold-r-normal--6-60-72-72-p-36-iso8859-1
-linotype-new century schoolbook-bold-r-normal--8-80-72-72-p-49-iso8859-1
-linotype-new century schoolbook-medium-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-new century schoolbook-medium-i-normal--10-100-72-72-p-56-iso8859-1
-linotype-new century schoolbook-medium-i-normal--12-120-72-72-p-67-iso8859-1
-linotype-new century schoolbook-medium-i-normal--14-140-72-72-p-79-iso8859-1
-linotype-new century schoolbook-medium-i-normal--6-60-72-72-p-34-iso8859-1
-linotype-new century schoolbook-medium-i-normal--8-80-72-72-p-46-iso8859-1
-linotype-new century schoolbook-medium-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-new century schoolbook-medium-r-normal--10-100-72-72-p-56-iso8859-1
-linotype-new century schoolbook-medium-r-normal--12-120-72-72-p-67-iso8859-1
-linotype-new century schoolbook-medium-r-normal--14-140-72-72-p-78-iso8859-1
-linotype-new century schoolbook-medium-r-normal--6-60-72-72-p-34-iso8859-1
-linotype-new century schoolbook-medium-r-normal--8-80-72-72-p-45-iso8859-1
-linotype-palatino-bold-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-palatino-bold-i-normal--10-100-72-72-p-55-iso8859-1
-linotype-palatino-bold-i-normal--12-120-72-72-p-66-iso8859-1
-linotype-palatino-bold-i-normal--14-140-72-72-p-78-iso8859-1
-linotype-palatino-bold-i-normal--6-60-72-72-p-33-iso8859-1
-linotype-palatino-bold-i-normal--8-80-72-72-p-44-iso8859-1
-linotype-palatino-bold-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-palatino-bold-r-normal--10-100-72-72-p-56-iso8859-1
-linotype-palatino-bold-r-normal--12-120-72-72-p-67-iso8859-1
-linotype-palatino-bold-r-normal--14-140-72-72-p-79-iso8859-1
-linotype-palatino-bold-r-normal--6-60-72-72-p-34-iso8859-1
-linotype-palatino-bold-r-normal--8-80-72-72-p-45-iso8859-1
-linotype-palatino-medium-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-palatino-medium-i-normal--10-100-72-72-p-52-iso8859-1
-linotype-palatino-medium-i-normal--12-120-72-72-p-63-iso8859-1
-linotype-palatino-medium-i-normal--14-140-72-72-p-74-iso8859-1
-linotype-palatino-medium-i-normal--6-60-72-72-p-32-iso8859-1
-linotype-palatino-medium-i-normal--8-80-72-72-p-42-iso8859-1
-linotype-palatino-medium-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-palatino-medium-r-normal--10-100-72-72-p-55-iso8859-1
-linotype-palatino-medium-r-normal--12-120-72-72-p-65-iso8859-1
-linotype-palatino-medium-r-normal--14-140-72-72-p-77-iso8859-1
-linotype-palatino-medium-r-normal--6-60-72-72-p-34-iso8859-1
-linotype-palatino-medium-r-normal--8-80-72-72-p-44-iso8859-1
-linotype-times-bold-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-times-bold-i-normal--10-100-72-72-p-52-iso8859-1
-linotype-times-bold-i-normal--12-120-72-72-p-63-iso8859-1
-linotype-times-bold-i-normal--14-140-72-72-p-73-iso8859-1
-linotype-times-bold-i-normal--6-60-72-72-p-31-iso8859-1
-linotype-times-bold-i-normal--8-80-72-72-p-42-iso8859-1
-linotype-times-bold-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-times-bold-r-normal--10-100-72-72-p-54-iso8859-1
-linotype-times-bold-r-normal--12-120-72-72-p-65-iso8859-1
-linotype-times-bold-r-normal--14-140-72-72-p-75-iso8859-1
-linotype-times-bold-r-normal--6-60-72-72-p-32-iso8859-1
-linotype-times-bold-r-normal--8-80-72-72-p-43-iso8859-1
-linotype-times-medium-i-normal--0-0-0-0-p-0-iso8859-1
-linotype-times-medium-i-normal--10-100-72-72-p-51-iso8859-1
-linotype-times-medium-i-normal--12-120-72-72-p-61-iso8859-1
-linotype-times-medium-i-normal--14-140-72-72-p-72-iso8859-1
-linotype-times-medium-i-normal--6-60-72-72-p-31-iso8859-1
-linotype-times-medium-i-normal--8-80-72-72-p-41-iso8859-1
-linotype-times-medium-r-normal--0-0-0-0-p-0-iso8859-1
-linotype-times-medium-r-normal--10-100-72-72-p-51-iso8859-1
-linotype-times-medium-r-normal--12-120-72-72-p-62-iso8859-1
-linotype-times-medium-r-normal--14-140-72-72-p-73-iso8859-1
-linotype-times-medium-r-normal--6-60-72-72-p-31-iso8859-1
-linotype-times-medium-r-normal--8-80-72-72-p-42-iso8859-1
-misc-fixed-bold-r-normal--13-100-100-100-c-70-iso8859-1
-misc-fixed-bold-r-normal--13-100-100-100-c-80-iso8859-1
-misc-fixed-bold-r-normal--13-120-75-75-c-70-iso8859-1
-misc-fixed-bold-r-normal--13-120-75-75-c-80-iso8859-1
-misc-fixed-bold-r-normal--15-120-100-100-c-90-iso8859-1
-misc-fixed-bold-r-normal--15-140-75-75-c-90-iso8859-1
-misc-fixed-bold-r-semicondensed--13-100-100-100-c-60-iso8859-1
-misc-fixed-bold-r-semicondensed--13-120-75-75-c-60-iso8859-1
-misc-fixed-medium-r-normal--10-100-75-75-c-60-iso8859-1
-misc-fixed-medium-r-normal--10-70-100-100-c-60-iso8859-1
-misc-fixed-medium-r-normal--13-100-100-100-c-70-iso8859-1
-misc-fixed-medium-r-normal--13-100-100-100-c-80-iso8859-1
-misc-fixed-medium-r-normal--13-120-75-75-c-70-iso8859-1
-misc-fixed-medium-r-normal--13-120-75-75-c-80-iso8859-1
-misc-fixed-medium-r-normal--14-110-100-100-c-70-iso8859-1
-misc-fixed-medium-r-normal--14-130-75-75-c-140-jisx0208.1983-0
-misc-fixed-medium-r-normal--14-130-75-75-c-70-iso8859-1
-misc-fixed-medium-r-normal--14-130-75-75-c-70-jisx0201.1976-0
-misc-fixed-medium-r-normal--15-120-100-100-c-90-iso8859-1
-misc-fixed-medium-r-normal--15-140-75-75-c-90-iso8859-1
-misc-fixed-medium-r-normal--20-140-100-100-c-100-iso8859-1
-misc-fixed-medium-r-normal--20-200-75-75-c-100-iso8859-1
-misc-fixed-medium-r-normal--8-60-100-100-c-50-iso8859-1
-misc-fixed-medium-r-normal--8-80-75-75-c-50-iso8859-1
-misc-fixed-medium-r-normal--9-80-100-100-c-60-iso8859-1
-misc-fixed-medium-r-normal--9-90-75-75-c-60-iso8859-1
-misc-fixed-medium-r-semicondensed--12-110-75-75-c-60-iso8859-1
-misc-fixed-medium-r-semicondensed--12-90-100-100-c-60-iso8859-1
-misc-fixed-medium-r-semicondensed--13-100-100-100-c-60-iso8859-1
-misc-fixed-medium-r-semicondensed--13-120-75-75-c-60-iso8859-1
-monotype-bembo-bold-i-normal--0-0-0-0-p-0-iso8859-1
-monotype-bembo-bold-i-normal--10-100-75-75-p-48-iso8859-1
-monotype-bembo-bold-i-normal--12-120-75-75-p-58-iso8859-1
-monotype-bembo-bold-i-normal--14-140-75-75-p-67-iso8859-1
-monotype-bembo-bold-i-normal--6-60-72-72-p-32-iso8859-1
-monotype-bembo-bold-i-normal--8-80-72-72-p-43-iso8859-1
-monotype-bembo-bold-r-normal--0-0-0-0-p-0-iso8859-1
-monotype-bembo-bold-r-normal--10-100-75-75-p-51-iso8859-1
-monotype-bembo-bold-r-normal--12-120-75-75-p-59-iso8859-1
-monotype-bembo-bold-r-normal--14-140-75-75-p-69-iso8859-1
-monotype-bembo-bold-r-normal--6-60-72-72-p-33-iso8859-1
-monotype-bembo-bold-r-normal--8-80-72-72-p-45-iso8859-1
-monotype-bembo-medium-i-normal--0-0-0-0-p-0-iso8859-1
-monotype-bembo-medium-i-normal--10-100-75-75-p-45-iso8859-1
-monotype-bembo-medium-i-normal--12-120-75-75-p-55-iso8859-1
-monotype-bembo-medium-i-normal--14-140-75-75-p-64-iso8859-1
-monotype-bembo-medium-i-normal--6-60-72-72-p-30-iso8859-1
-monotype-bembo-medium-i-normal--8-80-72-72-p-41-iso8859-1
-monotype-bembo-medium-r-normal--0-0-0-0-p-0-iso8859-1
-monotype-bembo-medium-r-normal--10-100-75-75-p-47-iso8859-1
-monotype-bembo-medium-r-normal--12-120-75-75-p-55-iso8859-1
-monotype-bembo-medium-r-normal--14-140-75-75-p-63-iso8859-1
-monotype-bembo-medium-r-normal--6-60-72-72-p-32-iso8859-1
-monotype-bembo-medium-r-normal--8-80-72-72-p-41-iso8859-1
-monotype-gill sans-bold-i-normal--10-100-75-75-p-50-iso8859-1
-monotype-gill sans-bold-i-normal--12-120-75-75-p-60-iso8859-1
-monotype-gill sans-bold-i-normal--14-140-75-75-p-69-iso8859-1
-monotype-gill sans-bold-i-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill sans-bold-i-normal-sans-6-60-72-72-p-32-iso8859-1
-monotype-gill sans-bold-i-normal-sans-8-80-72-72-p-43-iso8859-1
-monotype-gill sans-bold-r-normal--10-100-75-75-p-51-iso8859-1
-monotype-gill sans-bold-r-normal--12-120-75-75-p-63-iso8859-1
-monotype-gill sans-bold-r-normal--14-140-75-75-p-73-iso8859-1
-monotype-gill sans-bold-r-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill sans-bold-r-normal-sans-6-60-72-72-p-34-iso8859-1
-monotype-gill sans-bold-r-normal-sans-8-80-72-72-p-45-iso8859-1
-monotype-gill sans-medium-i-normal--10-100-75-75-p-44-iso8859-1
-monotype-gill sans-medium-i-normal--12-120-75-75-p-52-iso8859-1
-monotype-gill sans-medium-i-normal--14-140-75-75-p-61-iso8859-1
-monotype-gill sans-medium-i-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill sans-medium-i-normal-sans-6-60-72-72-p-28-iso8859-1
-monotype-gill sans-medium-i-normal-sans-8-80-72-72-p-39-iso8859-1
-monotype-gill sans-medium-r-normal--10-100-75-75-p-43-iso8859-1
-monotype-gill sans-medium-r-normal--12-120-75-75-p-56-iso8859-1
-monotype-gill sans-medium-r-normal--14-140-75-75-p-64-iso8859-1
-monotype-gill sans-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill sans-medium-r-normal-sans-6-60-72-72-p-31-iso8859-1
-monotype-gill sans-medium-r-normal-sans-8-80-72-72-p-41-iso8859-1
-monotype-gill-bold-i-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill-bold-r-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-gill-normal-i-normal-sans-0-0-0-0-p-0-iso8859-1
-monotype-rockwell-bold-i-normal--0-0-0-0-p-0-iso8859-1
-monotype-rockwell-bold-i-normal--10-100-75-75-p-53-iso8859-1
-monotype-rockwell-bold-i-normal--12-120-75-75-p-63-iso8859-1
-monotype-rockwell-bold-i-normal--14-140-75-75-p-74-iso8859-1
-monotype-rockwell-bold-i-normal--6-60-72-72-p-34-iso8859-1
-monotype-rockwell-bold-i-normal--8-80-72-72-p-46-iso8859-1
-monotype-rockwell-bold-r-normal--0-0-0-0-p-0-iso8859-1
-monotype-rockwell-bold-r-normal--10-100-75-75-p-52-iso8859-1
-monotype-rockwell-bold-r-normal--12-120-75-75-p-61-iso8859-1
-monotype-rockwell-bold-r-normal--14-140-75-75-p-74-iso8859-1
-monotype-rockwell-bold-r-normal--6-60-72-72-p-35-iso8859-1
-monotype-rockwell-bold-r-normal--8-80-72-72-p-47-iso8859-1
-monotype-rockwell-medium-i-normal--0-0-0-0-p-0-iso8859-1
-monotype-rockwell-medium-i-normal--10-100-75-75-p-49-iso8859-1
-monotype-rockwell-medium-i-normal--12-120-75-75-p-59-iso8859-1
-monotype-rockwell-medium-i-normal--14-140-75-75-p-68-iso8859-1
-monotype-rockwell-medium-i-normal--6-60-72-72-p-32-iso8859-1
-monotype-rockwell-medium-i-normal--8-80-72-72-p-43-iso8859-1
-monotype-rockwell-medium-r-normal--0-0-0-0-p-0-iso8859-1
-monotype-rockwell-medium-r-normal--10-100-75-75-p-47-iso8859-1
-monotype-rockwell-medium-r-normal--12-120-75-75-p-56-iso8859-1
-monotype-rockwell-medium-r-normal--14-140-75-75-p-67-iso8859-1
-monotype-rockwell-medium-r-normal--6-60-72-72-p-33-iso8859-1
-monotype-rockwell-medium-r-normal--8-80-72-72-p-44-iso8859-1
-schumacher-clean-bold-r-normal--10-100-75-75-c-60-iso8859-1
-schumacher-clean-bold-r-normal--10-100-75-75-c-80-iso8859-1
-schumacher-clean-bold-r-normal--12-120-75-75-c-60-iso8859-1
-schumacher-clean-bold-r-normal--12-120-75-75-c-80-iso8859-1
-schumacher-clean-bold-r-normal--13-130-75-75-c-80-iso8859-1
-schumacher-clean-bold-r-normal--14-140-75-75-c-80-iso8859-1
-schumacher-clean-bold-r-normal--15-150-75-75-c-90-iso8859-1
-schumacher-clean-bold-r-normal--16-160-75-75-c-80-iso8859-1
-schumacher-clean-bold-r-normal--8-80-75-75-c-80-iso8859-1
-schumacher-clean-medium-i-normal--12-120-75-75-c-60-iso8859-1
-schumacher-clean-medium-i-normal--8-80-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--10-100-75-75-c-50-iso8859-1
-schumacher-clean-medium-r-normal--10-100-75-75-c-60-iso8859-1
-schumacher-clean-medium-r-normal--10-100-75-75-c-70-iso8859-1
-schumacher-clean-medium-r-normal--10-100-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--12-120-75-75-c-60-iso8859-1
-schumacher-clean-medium-r-normal--12-120-75-75-c-70-iso8859-1
-schumacher-clean-medium-r-normal--12-120-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--13-130-75-75-c-60-iso8859-1
-schumacher-clean-medium-r-normal--13-130-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--14-140-75-75-c-70-iso8859-1
-schumacher-clean-medium-r-normal--14-140-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--15-150-75-75-c-90-iso8859-1
-schumacher-clean-medium-r-normal--16-160-75-75-c-80-iso8859-1
-schumacher-clean-medium-r-normal--6-60-75-75-c-40-iso8859-1
-schumacher-clean-medium-r-normal--6-60-75-75-c-50-iso8859-1
-schumacher-clean-medium-r-normal--6-60-75-75-c-60-iso8859-1
-schumacher-clean-medium-r-normal--8-80-75-75-c-50-iso8859-1
-schumacher-clean-medium-r-normal--8-80-75-75-c-60-iso8859-1
-schumacher-clean-medium-r-normal--8-80-75-75-c-70-iso8859-1
-schumacher-clean-medium-r-normal--8-80-75-75-c-80-iso8859-1
-sony-fixed-medium-r-normal--16-120-100-100-c-80-iso8859-1
-sony-fixed-medium-r-normal--16-120-100-100-c-80-jisx0201.1976-0
-sony-fixed-medium-r-normal--16-150-75-75-c-80-iso8859-1
-sony-fixed-medium-r-normal--16-150-75-75-c-80-jisx0201.1976-0
-sony-fixed-medium-r-normal--24-170-100-100-c-120-iso8859-1
-sony-fixed-medium-r-normal--24-170-100-100-c-120-jisx0201.1976-0
-sony-fixed-medium-r-normal--24-230-75-75-c-120-iso8859-1
-sony-fixed-medium-r-normal--24-230-75-75-c-120-jisx0201.1976-0
-sun-open look cursor-----12-120-75-75-p-455-sunolcursor-1
-sun-open look glyph-----10-100-75-75-p-106-sunolglyph-1
-sun-open look glyph-----12-120-75-75-p-116-sunolglyph-1
-sun-open look glyph-----14-140-75-75-p-136-sunolglyph-1
-sun-open look glyph-----19-190-75-75-p-163-sunolglyph-1
-urw-itc avant garde-demi-o-normal-sans-0-0-0-0-p-0-iso8859-1
-urw-itc avant garde-demi-o-normal-sans-10-100-72-72-p-56-iso8859-1
-urw-itc avant garde-demi-o-normal-sans-12-120-72-72-p-68-iso8859-1
-urw-itc avant garde-demi-o-normal-sans-14-140-72-72-p-79-iso8859-1
-urw-itc avant garde-demi-o-normal-sans-6-60-72-72-p-34-iso8859-1
-urw-itc avant garde-demi-o-normal-sans-8-80-72-72-p-45-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-0-0-0-0-p-0-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-10-100-72-72-p-56-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-12-120-72-72-p-68-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-14-140-72-72-p-79-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-6-60-72-72-p-34-iso8859-1
-urw-itc avant garde-demi-r-normal-sans-8-80-72-72-p-45-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-0-0-0-0-p-0-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-10-100-72-72-p-57-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-12-120-72-72-p-68-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-14-140-72-72-p-79-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-6-60-72-72-p-34-iso8859-1
-urw-itc avant garde-medium-o-normal-sans-8-80-72-72-p-45-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-0-0-0-0-p-0-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-10-100-72-72-p-57-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-12-120-72-72-p-68-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-14-140-72-72-p-79-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-6-60-72-72-p-34-iso8859-1
-urw-itc avant garde-medium-r-normal-sans-8-80-72-72-p-45-iso8859-1
-urw-itc bookman-demi-i-normal--10-100-72-72-p-61-iso8859-1
-urw-itc bookman-demi-i-normal--12-120-72-72-p-72-iso8859-1
-urw-itc bookman-demi-i-normal--14-140-72-72-p-85-iso8859-1
-urw-itc bookman-demi-i-normal--6-60-72-72-p-36-iso8859-1
-urw-itc bookman-demi-i-normal--8-80-72-72-p-48-iso8859-1
-urw-itc bookman-demi-r-normal--10-100-72-72-p-60-iso8859-1
-urw-itc bookman-demi-r-normal--12-120-72-72-p-73-iso8859-1
-urw-itc bookman-demi-r-normal--14-140-72-72-p-84-iso8859-1
-urw-itc bookman-demi-r-normal--6-60-72-72-p-36-iso8859-1
-urw-itc bookman-demi-r-normal--8-80-72-72-p-48-iso8859-1
-urw-itc bookman-light-i-normal--10-100-72-72-p-56-iso8859-1
-urw-itc bookman-light-i-normal--12-120-72-72-p-66-iso8859-1
-urw-itc bookman-light-i-normal--14-140-72-72-p-80-iso8859-1
-urw-itc bookman-light-i-normal--6-60-72-72-p-35-iso8859-1
-urw-itc bookman-light-i-normal--8-80-72-72-p-45-iso8859-1
-urw-itc bookman-light-r-normal--10-100-72-72-p-57-iso8859-1
-urw-itc bookman-light-r-normal--12-120-72-72-p-69-iso8859-1
-urw-itc bookman-light-r-normal--14-140-72-72-p-80-iso8859-1
-urw-itc bookman-light-r-normal--6-60-72-72-p-35-iso8859-1
-urw-itc bookman-light-r-normal--8-80-72-72-p-46-iso8859-1
10x20
12x24
12x24kana
12x24romankana
5x8
6x10
6x12
6x13
6x13bold
6x9
7x13
7x13bold
7x14
7x14romankana
8x13
8x13bold
8x16
8x16kana
8x16romankana
9x15
9x15bold
a12biluc
a12bluci
a12butto
a12iluci
a12lucid
a12sbarh
a12sbarv
a12sldrh
a12sldrv
a14
adobe-helvetica
adobe-helvetica-bold
adobe-helvetica-boldoblique
adobe-helvetica-oblique
adobe-newcenturyschlbk-bold
adobe-newcenturyschlbk-italic
adobe-newcenturyschlbk-roman
adobe-times-bold
adobe-times-bolditalic
adobe-times-italic
adobe-times-roman
avantgarde-book
avantgarde-bookoblique
avantgarde-demi
avantgarde-demioblique
b12bluci
b12butto
b12lucid
b12sbarh
b12sbarv
b12sldrh
b12sldrv
bembo
bembo-bold
bembo-bolditalic
bembo-italic
bookman-demi
bookman-demiitalic
bookman-light
bookman-lightitalic
c12bluci
c12butto
c12iluci
c12lucid
c12sbarh
c12sbarv
c12sldrh
c12sldrv
charb
charbi
chari
charr
charter-black
charter-black-italic
charter-italic
charter-roman
circles
cmr.b.
cmr.r.
courb
courbo
courier
courier-bold
courier-boldoblique
courier-oblique
couro
courr
cursor
d12bluci
d12butto
d12iluci
d12lucid
d12sbarh
d12sbarv
d12sldrh
d12sldrv
decw$cursor
decw$session
e12bluci
e12butto
e12lucid
e12sbarh
e12sbarv
e12sldrh
e12sldrv
f12bluci
f12butto
f12lucid
f12sbarh
f12sbarv
f12sldrh
f12sldrv
fixed
g12bluci
g12butto
g12lucid
g12sbarh
g12sbarv
g12sldrh
g12sldrv
gallant.r.
gillsans
gillsans-bold
gillsans-bolditalic
gillsans-italic
h12bluci
h12butto
h12iluci
h12lucid
h12sbarh
h12sbarv
h12sldrh
h12sldrv
helvb
helvbo
helvetica
helvetica-bold
helvetica-boldoblique
helvetica-narrow
helvetica-narrow-bold
helvetica-narrow-boldoblique
helvetica-narrow-oblique
helvetica-oblique
helvo
helvr
icon
k14
kana14
kanji
lucida-bright
lucida-brightdemibold
lucida-brightdemibolditalic
lucida-brightitalic
lucidabright
lucidabright-demi
lucidabright-demiitalic
lucidabright-italic
lucidasans
lucidasans-bold
lucidasans-bolditalic
lucidasans-italic
lucidasans-typewriter
lucidasans-typewriterbold
lucidasanstypewriter
lucidasanstypewriter-bold
ncenb
ncenbi
nceni
ncenr
newcenturyschlbk-bold
newcenturyschlbk-bolditalic
newcenturyschlbk-italic
newcenturyschlbk-roman
newscursor
nil2
olcursor
olglyph
palatino-bold
palatino-bolditalic
palatino-italic
palatino-roman
r14
r16
r24
rk14
rk16
rk24
rockwell
rockwell-bold
rockwell-bolditalic
rockwell-italic
screen
screen-bold
screen.b.
screen.r.
symbb
symbol
terminal
terminal-bold
terminal-bold-normal
terminal-normal
timb
timbi
times-bold
times-bolditalic
times-italic
times-roman
timi
timr
variable
vshd
vtbold
vtsingle
zapfchancery-mediumitalic
zapfdingbats
}
  global {input_dir}
  set {input_dir} {/home/midway/seidl/src/SC/src.bin/tkmpqc}
  global {mpqc_path}
  set {mpqc_path} {mpqcic}
  global {outfile}
  set {outfile} {mpqc.out}
  global {statusmsg}
  set {statusmsg} {Welcome to MPQC (or is it SC?)}

  # please don't modify the following
  # variables. They are needed by xf.
  global {autoLoadList}
  set {autoLoadList(main.tcl)} {0}
  set {autoLoadList(mpqc.tcl)} {0}
  global {internalAliasList}
  set {internalAliasList} {}
  global {moduleList}
  set {moduleList(mpqc.tcl)} {}
  global {preloadList}
  set {preloadList(xfInternal)} {}
  global {symbolicName}
  set {symbolicName(root)} {.}
  global {xfWmSetPosition}
  set {xfWmSetPosition} {}
  global {xfWmSetSize}
  set {xfWmSetSize} {}
  global {xfAppDefToplevels}
  set {xfAppDefToplevels} {}
}

# initialize global variables
InitGlobals

# display/remove toplevel windows.
ShowWindow.

global xfShowWindow.input
set xfShowWindow.input 0

# load default bindings.
if {[info exists env(XF_BIND_FILE)] &&
    "[info procs XFShowHelp]" == ""} {
  source $env(XF_BIND_FILE)
}

# parse and apply application defaults.
XFLocalLoadAppDefs Mpqc
XFLocalSetAppDefs

# end source
EndSrc

# eof
#

