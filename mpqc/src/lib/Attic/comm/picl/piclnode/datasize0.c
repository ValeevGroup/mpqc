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


int datasize0 ( datatype )
	int datatype ;
{
	int size ;

	switch (datatype) {
	case 0:
		size = sizeof(char) ;
		break ;
	case 1:
		size = sizeof(short) ;
		break ;
	case 2:
		size = sizeof(int) ;
		break ;
	case 3:
		size = sizeof(long) ;
		break ;
	case 4:
		size = sizeof(float) ;
		break ;
	case 5:
		size = sizeof(double) ;
		break ;
	}
	return (size) ;
}

