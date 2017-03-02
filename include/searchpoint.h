#ifndef SGSEARCHPOINT_H
#define SGSEARCHPOINT_H

#include <iostream>

namespace SG{

struct searchpoint
{
	double val;
	unsigned int indexGrid;
	unsigned int indexBin;

	std::ostream& outputToStream(std::ostream& raus) const;

	searchpoint();
	searchpoint(double valIn, unsigned int indexGridIn, unsigned int indexBinIn);
};

bool wayToSortSearchpoints(searchpoint const& i, searchpoint const& j);

std::ostream& operator<<(std::ostream& lhs, searchpoint const& rhs);

bool operator<(searchpoint const& a, searchpoint const& b);

bool operator<=(searchpoint const& a, searchpoint const& b);

bool operator>(searchpoint const& a, searchpoint const& b);

bool operator>=(searchpoint const& a, searchpoint const& b);

} // end of namespace SG

#endif 
