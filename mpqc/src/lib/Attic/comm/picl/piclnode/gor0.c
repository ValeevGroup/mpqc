/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  PICL source code                                               *
 *                                                                 *
 *  We welcome questions, comments, and bug reports, and request   *
 *  that you share any modifications with us.                      *
 *                                                                 *
 *  Patrick Worley                                                 *
 *  Oak Ridge National Laboratory                                  *
 *  worley@msr.epm.ornl.gov                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


gor0 ( buf, items, datatype, msgtype, root )
char *buf ;
int items, datatype, msgtype, root ;
/*
 *  Componentwise logical or of a vector over all processors,
 *  using given topology.
 *
 *  buf        array of data to be combined by logical or.
 *  items      number of items in array buf.
 *  datatype   code number for type of data:
 *  datatype = 0 : char
 *             1 : short   (int*2 in Fortran)
 *             2 : int
 *             3 : long    (int*4 in Fortran)
 *  msgtype    user-defined id to distinguish messages.
 *  root       processor in which final result will reside.
 *
 *  Caution: buf is overwritten.
 */
{
	void or0() ;

	gcomb0(buf, items, datatype, msgtype, root, or0) ;
}

void or0 ( p1, p2, items, datatype)
	char *p1, *p2 ;
	int items, datatype ;
{
	char *c1, *c2 ;
	short *s1, *s2 ;
	int *i1, *i2 ;
	long *l1, *l2 ;
	int i, bytes, size, datasize0() ;

	size = datasize0(datatype) ;
	bytes = items*size ;
	switch (datatype) {
	case 0:
		for (i = 0 ; i < bytes ; i += size) {
			c1 = (char *)(p1 + i) ;
			c2 = (char *)(p2 + i) ;
			*c1 |= *c2 ;
		}
		break ;
	case 1:
		for (i = 0 ; i < bytes ; i += size) {
			s1 = (short *)(p1 + i) ;
			s2 = (short *)(p2 + i) ;
			*s1 |= *s2 ;
		}
		break ;
	case 2:
		for (i = 0 ; i < bytes ; i += size) {
			i1 = (int *)(p1 + i) ;
			i2 = (int *)(p2 + i) ;
			*i1 |= *i2 ;
		}
		break ;
	case 3:
		for (i = 0 ; i < bytes ; i += size) {
			l1 = (long *)(p1 + i) ;
			l2 = (long *)(p2 + i) ;
			*l1 |= *l2 ;
		}
		break ;
	}
}

