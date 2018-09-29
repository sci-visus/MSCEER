#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#define INDEX_TYPE long long //int	// type to index elements of a mesh - need to support big meshes! - keep this signed to make arithmetic consistent
//#define FLOATTYPE float				// type used for computing numerical integration
#define FLOATTYPE float				// type used for computing numerical integration
#define INT_TYPE int					// regular old ints
#define DIM_TYPE unsigned char			// used to query the dimension of a cell - we usually have values between 0-3
#define BYTE_TYPE unsigned char			// just a regular old byte
#define BOUNDARY_TYPE unsigned char		// used to query the boundary classification of a cell - usually a small integer
#define ASSIGNED_TYPE unsigned char		// used to query if a cell has been assigned - usually 0, 1 values

#ifdef WIN32
enum DestType : BYTE_TYPE {
#else
enum DestType {
#endif	
	BACKGROUND,
	UNASSIGNED,
	ASSIGNED,
	CERTAIN_TERMINAL,
	CERTAIN_NONTERMINAL
};

  
#endif
