/*
 * vtkWriter.h
 *
 *  Created on: 09.11.2012
 *  Author: sgeisselsoeder
 */
 

#ifndef SGVTKWRITER_H_
#define SGVTKWRITER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "candidate.hpp"
#include <boost/lexical_cast.hpp>

namespace SG {

typedef float vtkDataType_t;

class vtkDataSet{
private:
	int check(std::vector<std::vector<vtkDataType_t> > const& data, unsigned int entrySize, unsigned int gridSize);
	
public:
	std::vector<std::vector<vtkDataType_t> > data_;
	std::string dataName_;
	
	// set the data with a check for apropriate dataformat, gridSize is zero for a new grid and the size of the already set grid for additional data
	int set(std::vector<std::vector<vtkDataType_t> > const& data, std::string name, unsigned int entrySize, unsigned int gridSize);
	
	inline unsigned int size(){
		return data_.size();
	}
	
	std::ostream& print(std::ostream& raus) const;
	
	std::ostream& printData(std::ostream& raus) const;
};


inline std::ostream& operator<< (std::ostream& raus, vtkDataSet const& data){
	return data.print(raus);
}


class vtkWriter{
private:
	// allData_ holds all dataSets; the first one is considered to be the geometry defining grid and therefor must be a 3D vector
	std::vector<vtkDataSet> allData_;

public:
	vtkWriter(){
		// nothing to do here
	}

	~vtkWriter(){
		// nothing to do here
	}
	
	// See if data with the name dataSet are already present.
	// Returns an iterator to the dataSet if found or to allData_.end() otherwise
	std::vector<vtkDataSet>::iterator find(std::string dataName);
	
	// sets the vtk writer to the initial state
	inline void clear(){
		allData_.clear();
	}
	
	// add a dataset of type T and name dataName; n defines how many entries per row; overwrites preexisting data of the name dataName
	// returns 0 if successful, negative int otherwise
	int add(std::vector<std::vector<vtkDataType_t> > const& newData, std::string dataName);

	// add a vector of candidates as 2 to 3 dimensional azimuth-zenith-proximity space
	// returns 0 if successful, negative int otherwise
	int add2D(std::vector<candidate_t> const& newData, std::string dataName);
	
	// add a vector of candidates displayed as a 3 dimensional sphere
	// returns 0 if successful, negative int otherwise
	int add3D(std::vector<candidate_t> const& newData, std::string dataName);
	
	// erases the data of name dataName; does nothing if dataName was not present
	void erase(std::string dataName);
	
	// writes all dataSets to the vtk file fileName
	// return 0 if everything works, negative numbers otherwise
	int write(std::string fileName);
	
	inline void write(std::vector<candidate_t> const& newData, std::string dataName){
		//std::cerr << "DEBUG write in .h started" << std::endl;
		clear();
		//std::cerr << "DEBUG write in .h 1" << std::endl;
		std::string dataNameWithVtk = dataName;
		//std::cerr << "DEBUG write in .h 2 started" << std::endl;
		std::string dataNameWithoutVtk = dataName;
		//std::cerr << "DEBUG write in .h 3 started" << std::endl;
		if (dataName.find(".vtk") == std::string::npos){
		//std::cerr << "DEBUG write in .h 4 started" << std::endl;
			dataNameWithVtk += ".vtk";
		//std::cerr << "DEBUG write in .h 5 started" << std::endl;
		} else {
		//std::cerr << "DEBUG write in .h 6started" << std::endl;
			dataNameWithoutVtk = dataName.substr(0, dataName.find(".vtk"));
		//std::cerr << "DEBUG write in .h 7started" << std::endl;
		}
		//std::cerr << "DEBUG write in .h 8started" << std::endl;
		add3D(newData, dataNameWithoutVtk);
		//std::cerr << "DEBUG write in .h 9started" << std::endl;
		write(dataNameWithVtk);
		//std::cerr << "DEBUG write in .h 0started" << std::endl;
	}
	
	inline void write(std::vector<candidate_t> const& newData, std::vector<std::vector<vtkDataType_t> > ratios, std::string dataName){
		clear();
		
		if (ratios.size() == 0){
			std::cerr << "Ratios contains no data! Not writing anything!" << std::endl;
			return;
		}
		unsigned int numBins = ratios[0].size();
		for (unsigned int i = 0; i < numBins; ++i){
			std::string dataNameWithVtk = dataName;
			std::string dataNameWithoutVtk = dataName;
			if (dataName.find(".vtk") == std::string::npos){
				dataNameWithVtk += ".vtk";
			} else {
				dataNameWithoutVtk = dataName.substr(0, dataName.find(".vtk"));
			}
			
			add3D(newData, dataNameWithoutVtk);
			std::string temp = "Ratios";
			std::vector<std::vector<vtkDataType_t> > valuesI;
			for (unsigned int j = 0; j < ratios.size(); ++j){
				std::vector<vtkDataType_t> entry;
				entry.push_back(ratios[j][i]);
				valuesI.push_back(entry);
			}
			add(valuesI, temp);
			dataNameWithVtk = dataNameWithoutVtk + boost::lexical_cast<std::string>(i) + ".vtk";
			write(dataNameWithVtk);
		}
	}
};

/*
// HOWTO USE:
	vtkWriter w;

	std::vector< std::vector<double> > data;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			for (int k = 0; k < 2; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i));
				d.push_back(static_cast<double>(j));
				d.push_back(static_cast<double>(k));
				data.push_back(d);
			}
	w.add(data, "cool_grid");
	
	std::vector< std::vector<double> > data2;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			for (int k = 0; k < 2; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i%2));
				d.push_back(static_cast<double>(j%1));
				d.push_back(static_cast<double>(k%2));
				data2.push_back(d);
			}
	w.add(data2, "testVectors");
			
	std::vector< std::vector<double> > data3;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			for (int k = 0; k < 2; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i%2*j*k));
				data3.push_back(d);
			}
	w.add(data3, "coolScalars");

	
	w.write("test.vtk");
	
	w.clear();
	
	std::vector< std::vector<double> > data222;
	for (int i = 0; i < 10; ++i)
		for (int j = 0; j < 10; ++j)
			for (int k = 0; k < 10; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i));
				d.push_back(static_cast<double>(j));
				d.push_back(static_cast<double>(k));
				data222.push_back(d);
			}
	w.add(data222, "cool_grid2");
	
	std::vector< std::vector<double> > data22;
	for (int i = 0; i < 10; ++i)
		for (int j = 0; j < 10; ++j)
			for (int k = 0; k < 10; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i%2));
				d.push_back(static_cast<double>(j%1));
				d.push_back(static_cast<double>(k%2));
				data22.push_back(d);
			}
	w.add(data22, "testVectors2");
			
	std::vector< std::vector<double> > data32;
	for (int i = 0; i < 10; ++i)
		for (int j = 0; j < 10; ++j)
			for (int k = 0; k < 10; ++k){
				std::vector<double> d;
				d.push_back(static_cast<double>(i%2*j*k));
				data32.push_back(d);
			}
	w.add(data32, "coolScalars2");

	
	w.write("test2.vtk");
//	*/

} // end of namespace SG

#endif /* SGVTKWRITER_H_ */
