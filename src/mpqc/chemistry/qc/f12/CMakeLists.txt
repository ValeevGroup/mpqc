set(sources
cabs_singles.h
ccsd_f12.cpp
ccsd_f12.h
ccsd_t_f12.cpp
ccsd_t_f12.h
db_f12_intermediates.h
dbccsd_f12.cpp
dbccsd_f12.h
dbmp2f12.cpp
dbmp2f12.h
f12_intermediates.h
gf2f12.cpp
gf2f12.h
linkage.h
mp2f12.cpp
mp2f12.h
)

add_mpqc_library(f12 sources sources "MPQCscf;MPQCmbpt;MPQCcc" "mpqc/chemistry/qc/f12")

