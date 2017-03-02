/*
 * vtkWriter.cpp
 *
 *  Created on: 09.11.2012
 *  Author: sgeisselsoeder
 */


#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vtkWriter.h"

namespace SG {

int vtkDataSet::check(std::vector<std::vector<vtkDataType_t> > const& data, unsigned int entrySize, unsigned int gridSize){
	if (gridSize != 0){
		if (data.size() != gridSize){
			std::cerr << "Failed to set the data because the number of vectors (" << data.size() << ") and the grid size (" << gridSize << ") do not agree." << std::endl;
			return -1;
		}
	}
	for (std::vector<std::vector<vtkDataType_t> >::const_iterator it = data.begin(); it != data.end(); ++it){
		if (it->size() != entrySize){
			std::cerr << "Failed to set the data because one vector has " << it->size() << " entries instead of " << entrySize << "." << std::endl;
			return -1;
		}
	}
		
	return 0;
}
	
// set the data with a check for apropriate dataformat, gridSize is zero for a new grid and the size of the already set grid for additional data
int vtkDataSet::set(std::vector<std::vector<vtkDataType_t> > const& data, std::string name, unsigned int entrySize, unsigned int gridSize){
	int ret = 0;
		
	// check if the data matches the format they are supposed to be saved
	ret = check(data, entrySize, gridSize);
	if (ret) {
		std::cerr << "Failed to set the data " << name << "." << std::endl;
		return ret;
	}
		
	// store the data
	for (std::vector<std::vector<vtkDataType_t> >::const_iterator it = data.begin(); it != data.end(); ++it){
		data_.push_back(*it);
	}
	dataName_ = name;
	return 0;
}
	
std::ostream& vtkDataSet::print(std::ostream& raus) const{
	if (data_.begin()->size() == 1){
		// output scalar subheader
		raus << "SCALARS " << dataName_ << " float" << std::endl;
		raus << "LOOKUP_TABLE default" << std::endl;
	} else if (data_.begin()->size() == 3){
		// output vector subheader
		raus << "VECTORS " << dataName_ << " float" << std::endl;
	} else {
		std::cerr << "Failed to print the data //TODO" << std::endl;
	}
		
	// output of the actual data
	return printData(raus);
}
	
std::ostream& vtkDataSet::printData(std::ostream& raus) const {
	raus.precision(10);
	// simply output the vector data in two dimensions
	for (std::vector<std::vector<vtkDataType_t> >::const_iterator itData = data_.begin(); itData != data_.end(); ++itData){
		for (std::vector<vtkDataType_t>::const_iterator thisPointIt = itData->begin(); thisPointIt != itData->end(); ++thisPointIt){
			raus << *thisPointIt << " ";
		}
		raus << std::endl;
	}
	raus << std::endl;
	
	return raus;
}

	
// See if data with the name dataSet are already present.
// Returns an iterator to the dataSet if found or to allData_.end() otherwise
std::vector<vtkDataSet>::iterator vtkWriter::find(std::string dataName){
	for (std::vector<vtkDataSet>::iterator it = allData_.begin(); it != allData_.end(); ++it)
		if (it->dataName_ == dataName)
			return it;
	
	return allData_.end();
}

std::vector<std::vector<vtkDataType_t> > gridAzZen;

// add a vector of candidates as 2 to 3 dimensional azimuth-zenith-proximity space
// returns 0 if successful, negative int otherwise
int vtkWriter::add2D(std::vector<candidate_t> const& newData, std::string dataName){
	if (newData.size() == 0) {
		std::cout << "vtkWriter: Not adding empty 2D DataSet." << std::endl;
		return 0;
	}
	std::vector<std::vector<vtkDataType_t> > grid;
	std::vector<std::vector<vtkDataType_t> > gridP;
	std::vector<std::vector<vtkDataType_t> > time;
	std::vector<std::vector<vtkDataType_t> > energy;
	std::vector<std::vector<vtkDataType_t> > proximity;
// 	std::vector<std::vector<vtkDataType_t> > proximityAlt;
	std::vector<std::vector<vtkDataType_t> > multiPurpose;
	
	for (std::vector<candidate_t>::const_iterator it = newData.begin(); it != newData.end(); ++it){
		std::vector<vtkDataType_t> temp;
		temp.clear(); temp.push_back(it->azimuth_); temp.push_back(it->zenith_); temp.push_back(0.0);
		grid.push_back(temp);
		temp.clear(); temp.push_back(it->azimuth_); temp.push_back(it->zenith_); temp.push_back(it->proximity_);
		gridP.push_back(temp);
		temp.clear(); temp.push_back(it->time_);
		time.push_back(temp);
		temp.clear(); temp.push_back(it->energy_);
		energy.push_back(temp);
		temp.clear(); temp.push_back(it->proximity_);
		proximity.push_back(temp);
// 		temp.clear(); temp.push_back(it->proximityAlt_);
// 		proximityAlt.push_back(temp);
		temp.clear(); temp.push_back(it->multiPurposeValue_);
		multiPurpose.push_back(temp);
		/*
		// I want this syntax to work!
		grid.push_back(std::vector<vtkDataType_t>(it->azimuth_, it->zenith_, 0.0));
		gridP.push_back(std::vector<vtkDataType_t>(it->azimuth_, it->zenith_, it->proximity_));
		time.push_back(std::vector<vtkDataType_t>(it->time_));
		energy.push_back(std::vector<vtkDataType_t>(it->energy_));
		proximity.push_back(std::vector<vtkDataType_t>(it->proximity_));
		*/
	}
	
	int ret = 0;
	ret += add(grid, dataName);
	ret += add(gridP, "proximityScaledGrid");
	ret += add(time, "Time");
	ret += add(energy, "Energy");
	ret += add(proximity, "Proximity");
// 	ret += add(proximity, "ProximityAlternat");
	ret += add(multiPurpose, "MultiPurposeValue");
	return ret;
}


// add a vector of candidates displayed as a 3 dimensional sphere
// returns 0 if successful, negative int otherwise
int vtkWriter::add3D(std::vector<candidate_t> const& newData, std::string dataName){
//std::cerr << "DEBUG add3D start" << std::endl;
	if (newData.size() == 0) {
//std::cerr << "DEBUG add3D 1start" << std::endl;
		std::cout << "vtkWriter: Not adding empty 3D DataSet." << std::endl;
		return 0;
	}
//std::cerr << "DEBUG add3D 2start" << std::endl;
	std::vector<std::vector<vtkDataType_t> > grid;
	std::vector<std::vector<vtkDataType_t> > gridP;
	std::vector<std::vector<vtkDataType_t> > time;
	std::vector<std::vector<vtkDataType_t> > time2;
	std::vector<std::vector<vtkDataType_t> > energy;
	std::vector<std::vector<vtkDataType_t> > proximity;
	std::vector<std::vector<vtkDataType_t> > multiPurpose;
	
//std::cerr << "DEBUG add3D 3start" << std::endl;
	for (std::vector<candidate_t>::const_iterator it = newData.begin(); it != newData.end(); ++it){
		std::vector<vtkDataType_t> tgrid;
		std::vector<vtkDataType_t> tgridP;
		std::vector<vtkDataType_t> ttime;
		std::vector<vtkDataType_t> ttime2;
		std::vector<vtkDataType_t> tenergy;
		std::vector<vtkDataType_t> tproximity;
		std::vector<vtkDataType_t> tmultiPurpose;

//std::cerr << "DEBUG add3D 4start" << std::endl;
		tgrid.push_back(it->x()); tgrid.push_back(it->y()); tgrid.push_back(it->z());
//std::cerr << "DEBUG add3D 5start" << std::endl;
		grid.push_back(tgrid);
		tgridP.push_back(it->xP(1000)); tgridP.push_back(it->yP(1000)); tgridP.push_back(it->zP(1000));
//std::cerr << "DEBUG add3D 6start" << std::endl;
		gridP.push_back(tgridP);
		ttime.push_back(it->time_);
		// unsigned int dist = it - newData.begin();
		// if (dist == 90000 || dist == 0 || dist == 60000)
		// 	std::cout << dist << " DEBUG vtwriter.cpp add3D time " << it->time_ << std::endl;
//std::cerr << "DEBUG add3D 7start" << std::endl;
		time.push_back(ttime);
		tenergy.push_back(it->energy_);
//std::cerr << "DEBUG add3D 8start" << std::endl;
		energy.push_back(tenergy);

		ttime2.push_back(it->time_);
		time2.push_back(ttime2);

		tproximity.push_back(it->proximity_);
//std::cerr << "DEBUG add3D 9start" << std::endl;
		proximity.push_back(tproximity);
//std::cerr << "multipurpose = " << it->multiPurposeValue_ << std::endl;
		tmultiPurpose.push_back(it->multiPurposeValue_);
//std::cerr << "DEBUG add3D 10start" << std::endl;
		multiPurpose.push_back(tmultiPurpose);
//std::cerr << "DEBUG add3D 11start" << std::endl;
	}
	
//std::cerr << "DEBUG add3D 12start" << std::endl;
	int ret = 0;
	ret += add(grid, dataName);
//std::cerr << "DEBUG add3D 13start" << std::endl;
	ret += add(gridP, "proximityScaledGrid");
//std::cerr << "DEBUG add3D 14start" << std::endl;
	ret += add(time, "Time");
//std::cerr << "DEBUG add3D 15start" << std::endl;
	ret += add(energy, "Energy");
	ret += add(time, "Time2");
//std::cerr << "DEBUG add3D 16start" << std::endl;
	ret += add(proximity, "Proximity");
//std::cerr << "DEBUG add3D 17start" << std::endl;
	ret += add(multiPurpose, "MultiPurposeValue");
//std::cerr << "DEBUG add3D 18start" << std::endl;
//std::cerr << "DEBUG add3D finished" << std::endl;
	return ret;
}
	
// add a dataset of type T and name dataName; n defines how many entries per row; overwrites preexisting data of the name dataName
// returns 0 if successful, negative int otherwise
int vtkWriter::add(std::vector<std::vector<vtkDataType_t> > const& newData, std::string dataName){
	int ret = 0;
		
	if ((newData.begin()->size() != 1) && (newData.begin()->size() != 3)){
		std::cerr << "The first entry of " << dataName << " has " << newData.begin()->size() << " entries, but only data with all vectors of sizes 1 and 3 are supported at the moment - sorry." << std::endl;
		return -2;
	}

	vtkDataSet d;
					
	// Do we already have data named like dataName?
	std::vector<vtkDataSet>::iterator it = find(dataName);
	if ( it != allData_.end()){	// yes: overwrite it
	
		// is it the grid?
		if (it == allData_.begin()){	// a grid must be defined!
			ret = d.set(newData, dataName, 3, 0);
			if (ret){
				std::cerr << "A three dimensional grid must be added first. Failed to interpret " << dataName << " as 3D grid." << std::endl; 
				return -1;
			}
				
		} else {	// not the grid: then every data possible is accepted
			ret = d.set(newData, dataName, newData.begin()->size(), allData_.begin()->size());
			if (ret){
				std::cerr << "The interpretation of " << dataName << " with " << newData.begin()->size() << " failed." << std::endl; 
				return -1;
			}
		}
		*it = d;
		
		
	} else { // not yet present: add it
		
		if (allData_.size() == 0){	// a grid must be define first!
			ret = d.set(newData, dataName, 3, 0);
			if (ret){
				std::cerr << "A three dimensional grid must be added first. Failed to interpret " << dataName << " as 3D grid." << std::endl; 
				return -1;
			}
		} else {	// not the grid: then every data possible is accepted
			ret = d.set(newData, dataName, newData.begin()->size(), allData_.begin()->size());
			if (ret){
				std::cerr << "The interpretation of " << dataName << " with " << newData.begin()->size() << " failed." << std::endl; 
				return -1;
			}
		}
		allData_.push_back(d);
	}
	return 0;
}
	
// erases the data of name dataName; does nothing if dataName was not present
void vtkWriter::erase(std::string dataName){
	std::vector<vtkDataSet>::iterator it = find(dataName);
	if (it == allData_.begin()){
		std::cout << "Erasing the grid has also erased all other data!" << std::endl;
		clear();
	}
	else if (it != allData_.end())
		allData_.erase(it);
}
	
// writes all dataSets to the vtk file fileName
// return 0 if everything works, negative numbers otherwise
int vtkWriter::write(std::string fileName){
//std::cerr << "DEBUG write start" << std::endl;
	// open the file
	std::ofstream raus;
	raus.open(fileName.c_str());
	if (! raus.is_open()){
		std::cerr << "Failed to open " << fileName << " for writing data in vtk format." << std::endl;
		return -1;
	}
	
//std::cerr << "DEBUG write is open, alldats size = " << allData_.size() << std::endl;
	raus.precision(10);
	if (allData_.size() >= 1){

//std::cerr << "DEBUG write 1" << std::endl;
		// output header
		raus << "# vtk DataFile Version 2.0" << std::endl;
		raus << allData_.begin()->dataName_ << std::endl;
		raus << "ASCII" << std::endl << std::endl;
	
//std::cerr << "DEBUG write 2" << std::endl;
		// output grid geometry subheader
		raus << "DATASET UNSTRUCTURED_GRID" << std::endl;
		raus << "POINTS " << allData_.begin()->size() << " float" << std::endl;
			
//std::cerr << "DEBUG write 3" << std::endl;
		// output the actual grid data without a vector header
		allData_.begin()->printData(raus);
			
//std::cerr << "DEBUG write 4" << std::endl;
		// output to say: "hey, the upcoming data belong the already set grid!"
		raus << "POINT_DATA " << allData_.begin()->size() << std::endl;
			
//std::cerr << "DEBUG write 5" << std::endl;
		// output the dataSet
		for (std::vector<vtkDataSet>::const_iterator it = allData_.begin()+1; it != allData_.end(); ++it){

//std::cerr << "DEBUG write 6" << std::endl;
			raus << *it << std::endl;
		}
		raus << std::endl;
//std::cerr << "DEBUG write 7" << std::endl;
	}
	
	// close the file(stream) - it has done it's job :)
	raus.close();
//std::cerr << "DEBUG write 8" << std::endl;
//std::cerr << "DEBUG write finished" << std::endl;
	
	return 0;
}
	
} // end of namespace SG
