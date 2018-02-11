set(public_headers
    gen_seq.h
    get_type.h
    make_array.h
    predicates.h
    sfinae.h
    store_variadic_pack.h
    type_traits.h
   )
   
add_mpqc_hdr_library(util_meta public_headers "" "mpqc/util/meta")
