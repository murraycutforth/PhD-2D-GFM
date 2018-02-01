#ifndef BBRANGE_H
#define BBRANGE_H

class BBrange {
	
	public:
	
	int imin, imax, jmin, jmax;
	
	BBrange (int imin, int imax, int jmin, int jmax)
	:
		imin	(imin),
		imax	(imax),
		jmin	(jmin),
		jmax	(jmax)
	{}
};


#endif
