
#include <scconfig.h>
#include <util/container/carray.h>

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void delete_c_array2<int>(int**);
template int** new_c_array2<int>(int,int,int);
template int** new_zero_c_array2<int>(int,int,int);
template void delete_c_array3<int>(int***);
template int*** new_c_array3<int>(int,int,int,int);
template int*** new_zero_c_array3<int>(int,int,int,int);

template void delete_c_array2<double>(double**);
template double** new_c_array2<double>(int,int,double);
template double** new_zero_c_array2<double>(int,int,double);
template void delete_c_array3<double>(double***);
template double*** new_c_array3<double>(int,int,int,double);
template double*** new_zero_c_array3<double>(int,int,int,double);
#endif
