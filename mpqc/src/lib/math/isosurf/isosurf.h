
extern int pg_isosurface_ts( int type, float *data, float *valdata, 
                             int nx_in, int ny_in, int nz_in, 
                             double value, P_Point *corner1, P_Point *corner2, 
                             int show_inside, int ftn_order,
                             TriangulatedSurface& surface,
                             UniformLattice& lattice);
