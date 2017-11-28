set(sources
cabs_singles.h
ccsd_f12.cpp
ccsd_f12.h
ccsd_t_f12.cpp
ccsd_t_f12.h
db_f12_intermediates.h
dbgf2f12.h
dbmp2f12.cpp
dbmp2f12.h
dbmp2f12_impl.h
f12_intermediates.h
gf2f12.cpp
gf2f12.h
linkage.h
mp2f12.cpp
mp2f12.h
mp2f12_impl.h
)

add_mpqc_library(lcao_f12 sources sources "MPQClcao_scf;MPQClcao_mbpt;MPQClcao_cc" "mpqc/chemistry/qc/lcao/f12")

