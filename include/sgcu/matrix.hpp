/*
 * sgcu.hpp
 * Contains a collection of sometimes useful common code sniplets
 *
 *  Created on: 2016-05-25
 *  Based on utils.hpp, 27.2.2012
 *  Author: sgeisselsoeder
 */

#ifndef SGCOMMONUTILS_MATRIX_H_
#define SGCOMMONUTILS_MATRIX_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>    // std::sort

namespace sgcu
{

template <typename T>
std::ostream& operator<<(std::ostream& raus, std::vector<std::vector<T> > const & data)
{
	for (unsigned int i = 0; i < data.size(); ++i)
	{
		raus << data[i] << std::endl;
	}
	return raus;
}

template <typename T>
void saveMatrixToFile(std::vector<std::vector<T> > const & data, std::string fileName)
{
	std::ofstream raus;
	raus.open(fileName.c_str());
	if ( (! raus.is_open()) || (! raus.good()) )
	{
		std::cerr << "Failed to open file " << fileName << " for matrix output." << std::endl;
		throw (-1);
	}
	raus << data;
	if (! raus.good() )
	{
		std::cerr << "Failed to output matrix to file " << fileName << std::endl;
		throw(-2);
	}
	raus.close();
}

//TODO
//template <typename T>
//void loadMatrixFromFile(std::string fileName, std::vector<T> & data)
//{
//	std::ifstream rein(fileName.c_str());
//	if (! rein.good())
//	{
//		std::cerr << "Failed to open file " << fileName << std::endl;
//		throw -1;
//	}
//	data.clear();
//	while (rein.good())
//	{
//		T val;
//		rein >> val;
//		if (rein.good())
//		{
//			data.push_back(val);
//		}
//	}
//	std::cout << "Loaded " << data.size() << " entries from file " << fileName << std::endl;
//	rein.close();
//}


template <typename T>
void transpose(std::vector<std::vector<T> > const & data, std::vector<std::vector<T> >& dataTransposed)
{
	unsigned int numRows = data.size();
	unsigned int numColumns = 0;
	if (numRows > 0) numColumns = data[0].size();
	if (numRows*numColumns == 0) return;	// check if data actually contains entries, nothing to do for empty data

	dataTransposed.clear();
	dataTransposed.resize(numColumns);
//	std::cerr << "asd" << std::endl;
	for (unsigned int j = 0; j < numColumns; ++j)
	{
//		std::cerr << "j=" << j << std::endl;
		dataTransposed[j].resize(numRows);
		for (unsigned int i = 0; i < numRows; ++i)
		{
//			std::cerr << "i=" << i << std::endl;
			dataTransposed[j][i] = data[i][j];
		}
	}
}

template <typename T>
void normalizeRowwise(std::vector<std::vector<T> >& data, unsigned int skipFirstNRows = 0)
{
	for (unsigned int i = skipFirstNRows; i < data.size(); ++i)
		sgcu::normalize(data[i]);
}

template <typename T>
void normalizeColumnwise(std::vector<std::vector<T> >& data, unsigned int skipFirstNColumns = 0)
{
	std::vector<std::vector<T> > dataTransposed;
	transpose(data, dataTransposed);
	normalizeRowwise(dataTransposed, skipFirstNColumns);
	transpose(dataTransposed, data);
}

template <typename T>
void removeColumns(std::vector<std::vector<T> > & data, std::vector<unsigned int> const & columnsToRemove)
{
	// sort columnsToRemove descending to start removing the latest index first (otherwise you mess up the data)
	std::vector<unsigned int> columnsToRemoveSorted = columnsToRemove;
	std::sort(columnsToRemoveSorted.begin(), columnsToRemoveSorted.end(), std::greater<unsigned int>());
//	std::cout << columnsToRemoveSorted << std::endl;

	std::vector<std::vector<T> > dataTransposed;
	transpose(data, dataTransposed);
	for (unsigned int i = 0; i < columnsToRemoveSorted.size(); ++i)
		dataTransposed.erase(dataTransposed.begin()+columnsToRemoveSorted[i]);
	transpose(dataTransposed, data);
}

template <typename T>
void getColumns(std::vector<std::vector<T> > const & data, std::vector<unsigned int> const & columnsToGet, std::vector<std::vector<T> > & dataExtracted)
{
	std::vector<std::vector<T> > dataTransposed;
	transpose(data, dataTransposed);
	std::vector<std::vector<T> > dataExtractedTransposed;
	for (unsigned int i = 0; i < columnsToGet.size(); ++i)
		dataExtractedTransposed.push_back(dataTransposed[columnsToGet[i]]);
	transpose(dataExtractedTransposed, dataExtracted);
}

} // end of namespace sgcu

#endif /* SGCOMMONUTILS_VECTORS_H_ */
