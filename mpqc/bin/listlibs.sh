#!/bin/csh -f

# Given a list of cpp options like -Dxxx and -Iyyy and the name
# of an include file this will run cpp -P on the file with the
# given options.  The resulting file is interpreting as a list
# if libraries that must be linked in.  Redundant names are removed
# from the list.

set path = ($path /usr/lib)

#echo listlibs: start > ll.trace

#echo $* >> ll.trace
#cpp -P $* >> ll.trace

set redundant_libs = `cpp -P $*`
set backwards_red_libs =
foreach i ($redundant_libs)
  set backwards_red_libs = ( $i $backwards_red_libs )
  end

#echo $redundant_libs >> ll.trace

set nonredundant_libs =
foreach i ($backwards_red_libs)
  set found = 0
  foreach j ($nonredundant_libs)
    if (XX$i == XX$j) then
      set found = 1
      continue;
      endif
    end
  if ($found == 0) then
    set nonredundant_libs = (`echo $i` $nonredundant_libs)
    endif
  end


# The tr is needed to get rid of double quotes which are needed
# for certain library names.
echo $nonredundant_libs | tr \" " "

#echo listlibs: end >> ll.trace
