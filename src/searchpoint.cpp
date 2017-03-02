/**
 * sgeisselsoeder
 *
 * June 2016
 *
 */

#include "searchpoint.h"
#include <iostream>

namespace SG{

std::ostream& searchpoint::outputToStream(std::ostream& raus) const
{
	raus << indexGrid << " " << indexBin << " " << val;
	return raus;
}

searchpoint::searchpoint()
{
	val = 0.0;
	indexGrid = 0;
	indexBin = 0;
}

searchpoint::searchpoint(double valIn, unsigned int indexGridIn, unsigned int indexBinIn)
{
	val = valIn;
	indexGrid = indexGridIn;
	indexBin = indexBinIn;
}


bool wayToSortSearchpoints(searchpoint const& i, searchpoint const& j)
{
	return i.val > j.val;
}

std::ostream& operator<<(std::ostream& lhs, searchpoint const& rhs)
{
	return rhs.outputToStream(lhs);
}

bool operator<(searchpoint const& a, searchpoint const& b)
{
	return a.val < b.val;
}

bool operator<=(searchpoint const& a, searchpoint const& b)
{
	return a.val <= b.val;
}

bool operator>(searchpoint const& a, searchpoint const& b)
{
	return a.val > b.val;
}

bool operator>=(searchpoint const& a, searchpoint const& b)
{
	return a.val >= b.val;
}

} // end of namespace SG

