/*
 * sgcu.hpp
 * Contains a collection of sometimes useful common code sniplets
 *
 *  Created on: 2016-05-25
 *  Based on utils.hpp, 27.2.2012
 *  Author: sgeisselsoeder
 */

#ifndef SGCOMMONUTILS_ARGUMENTS_H_
#define SGCOMMONUTILS_ARGUMENTS_H_

#include "boost/program_options.hpp"
#include <iostream>

namespace sgcu
{

template <typename T>
class argumentParser
{
private:
boost::program_options::variables_map vm_;
T temp_;
bool alreadyParsed_;

public:
boost::program_options::options_description desc_;

void defaultParse(int argc, char** argv)
{
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc_), vm_);
	boost::program_options::notify(vm_);
	alreadyParsed_ = true;
	if (vm_.count("help"))
	{
		std::cout << desc_ << "\n";
		std::cout << "Terminating now ..." << std::endl;
		throw -1;
	}
}

template <typename T2>
T2 parseVariable(std::string const & variableName, T2 defaultValue)
{
	if (! alreadyParsed_)
	{
		std::cerr << "You need to call void defaultParse(int argc, char** argv) before parsing variables!" << std::endl;
		std::cerr << "Example: myArgumentParser.defaultParse(argc, argv);" << std::endl;
		throw -2;
	}
	T2 variable = defaultValue;
	try
	{
		if (vm_.count(variableName))
		{
			variable = vm_[variableName].as<T2>();
			std::cout << variableName + " was set to " << variable << std::endl;
		}
		else
		{
			std::cout << variableName + " was not set. Using the default: " << variable << std::endl;
		}
	}
	catch (...)
	{
		std::cerr << "parseVariable: Failed to use boost::program_options::variables_map to parse variable " << variableName << std::endl;
		std::cerr << vm_[variableName].as<std::string>() << std::endl;
		throw (-1);
	}
	return variable;
}

argumentParser()
{
	alreadyParsed_ = false;
}

};

} // end of namespace sgcu

#endif /* SGCOMMONUTILS_ARGUMENTS_H_ */
