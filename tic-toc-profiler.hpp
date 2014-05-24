/*
The MIT License (MIT)

Copyright (c) 2013 Romain Futrzynski

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef TIC_TOC_PROFILER_HPP
#define TIC_TOC_PROFILER_HPP


//#define USE_PROFILER // Should be at the program level, not library level

#include <map>
#include <vector>
#include <string>
#include <fstream>

#include "yaml-cpp/yaml.h"


#ifdef USE_BOOST_CHRONO
#define BOOST_CHRONO_HEADER_ONLY
#define BOOST_ERROR_CODE_HEADER_ONLY
#include <boost/chrono.hpp>
namespace thischrono = boost::chrono;

#else

#include <chrono>
namespace thischrono = std::chrono;
#endif

using namespace std;

class event;
class eventlist;

class eventlist
{
public:
	map<string, event> eventlist_;
	eventlist *parentlist_;

	eventlist() {};
	/*eventlist(eventlist *creator) : parentlist_(creator)
	{};*/
};

class event
{
public:
	thischrono::high_resolution_clock::time_point tic_time, toc_time;
	vector<double> samples;
	eventlist *subEvents;

	event()
	{
		subEvents = new eventlist;
		samples.reserve(256);
	};
	/*event(eventlist *creator)
	{
	subEvents = new eventlist(creator);
	samples.reserve(256);
	}*/

	~event()
	{
		delete subEvents;
	};
};


/**		   \class DoProfiler
 *		   \brief The actual tic-toc-profiler class.
 *         The class profiler refers to this class if the program is compiled with USE_PROFILER
 *
 *  Detailed description starts here.
 */
class DoProfiler
{
	eventlist *current;
	string out_file;
	eventlist root;

public:
	/**
	 * Contains the resolution in seconds of the clock used.
	 */
	const double resolution;

	/**
	* Default constructor.
	* If the default constructor is used and no filename is given when using dump(), the timings will be printed in the default file "tic-toc-profiler.yml" in the execution directory.
	*/
	DoProfiler() : resolution((double)thischrono::high_resolution_clock::period::num / thischrono::high_resolution_clock::period::den), out_file("tic-toc-profiler.yml"), root(), current(&root)
	{}

	/**
	* Constructor with default filename where to save timings.
	*@param fname default filename where to save timings.
	*/
	DoProfiler(const string& fname) : resolution((double)thischrono::high_resolution_clock::period::num / thischrono::high_resolution_clock::period::den), out_file(fname), root(), current(&root)
	{}

	/**
	* Starts timing a section of code.
	*@param key name of the section being timed.
	*/
	inline void tic(const string& key)
	{
		(*(*current).eventlist_[key].subEvents).parentlist_ = current;
		current = (*current).eventlist_[key].subEvents;
		(*(*current).parentlist_).eventlist_[key].tic_time = thischrono::high_resolution_clock::now();
	}

	/**
	 * Stops timing a section of code.
	 * If a section with the same name on the same level has already been timed, all the samples are kept in a vector.
	 *@param key name of the section being timed.
	 */
	inline double toc(const string& key)
	{
		(*current).parentlist_->eventlist_[key].toc_time = thischrono::high_resolution_clock::now();
		current = (*current).parentlist_;
		const double sample = ((*current).eventlist_[key].toc_time - (*current).eventlist_[key].tic_time).count() * resolution;
		(*current).eventlist_[key].samples.push_back(sample);
		//cout << ((*parent)[key].toc_time - (*parent)[key].tic_time).count() << " " << (double)boost::chrono::steady_clock::period::num/boost::chrono::steady_clock::period::den  << endl;

		return sample;
	}

	/**
	 * Stops timing a section of code and prints the duration to std::cout.
	 * If a section with the same name on the same level has already been timed, all the samples are kept in a vector.
 	 *@param key name of the section being timed.
	 *@param text text to print in std::cout before the sample value.
     */
	inline double toc(const string& key, const string& text)
	{
		(*current).parentlist_->eventlist_[key].toc_time = thischrono::high_resolution_clock::now();
		current = (*current).parentlist_;
		const double sample = ((*current).eventlist_[key].toc_time - (*current).eventlist_[key].tic_time).count() * resolution;
		(*current).eventlist_[key].samples.push_back(sample);

		cout << text << sample << endl;
		return sample;
	}

	/**
	 * Prints samples to the file specified on construction.
	 */
	void dump();

	/**
	* Prints samples to the file specified by fname.
	*@param fname name of the sfile where to print the samples.
	*/
	void dump(const string& fname);

	/**
	* Clears all the samples recorded.
	*/
	inline void clear()
	{
		root.eventlist_.clear();
		current = &root;
	}
};

/**		   \class NoProfiler
 *		   \brief A dummy tic-toc-profiler class.
 *         The class profiler refers to this dummy class if the program is not compiled with USE_PROFILER
 *
 *  This class has the same members as the DoProfiler class, but they don't do anything. This is intended for use as an easier alternative to removing all calls to the profiler when no timing information about the code is required. 
 * @see DoProfiler for a description of the members.
 */
class NoProfiler
{
	eventlist *current = nullptr;
	string out_file = "";

public:
	const double factor = 0;
	eventlist root;
	NoProfiler() {};
	NoProfiler(const string&) {};

	inline void tic(const string&) {};
	inline double toc(const string&) { return -1; };
	inline double toc(const string&, const string&) { return -1; };
	void dump() {};
	void dump(string) {};
	inline void clear(){};
};

/**		\def USE_PROFILER
 *	Indicates if the tic-toc-profiler should be used. If USE_PROFILER is not defined, the type profiler refers to the dummy, empty class NoProfiler. If it is defined, profiler refers to the class DoProfiler insstead.
*/
#ifdef USE_PROFILER
typedef DoProfiler profiler;
#else
typedef NoProfiler profiler;
#endif


YAML::Emitter& operator << (YAML::Emitter& out, const eventlist& events);
YAML::Emitter& operator << (YAML::Emitter& out, const DoProfiler& prof);

inline YAML::Emitter& operator << (YAML::Emitter& out, const DoProfiler& prof)
{
	/*
	out << YAML::BeginMap;
	out << YAML::Key << "factor";
	out << YAML::Value << prof.factor;
	out << YAML::Comment("for time in seconds");
	out << YAML::EndMap;
	*/
	//std::cout << "unit: " << prof.factor << endl;
	out << prof.root;

	return out;
}

inline YAML::Emitter& operator << (YAML::Emitter& out, const eventlist& events)
{
	out << YAML::BeginMap;

	for (map<string, event>::const_iterator it = events.eventlist_.begin(); it != events.eventlist_.end(); it++)
	{
		out << YAML::Key << (*it).first;
		out << YAML::Value << YAML::BeginMap;
		out << YAML::Key << "samples";
		out << YAML::Value << YAML::Flow << (*it).second.samples;
		if ((*(*it).second.subEvents).eventlist_.size() > 0)
		{
			out << YAML::Key << "subFunctions";
			out << YAML::Value << (*(*it).second.subEvents);
		}
		out << YAML::EndMap;
	}
	out << YAML::EndMap;

	return out;
}

inline void DoProfiler::dump()
{

	YAML::Emitter yout;
	yout.SetDoublePrecision(16);
	yout << *this;
	ofstream fout(out_file.c_str(), ios_base::trunc);
	fout << yout.c_str();
	fout.close();

	//cout << yout.c_str();	

}
inline void DoProfiler::dump(const string& fname)
{

	YAML::Emitter yout;
	yout.SetDoublePrecision(16);
	yout << *this;
	ofstream fout(fname.c_str(), ios_base::trunc);
	fout << yout.c_str();
	fout.close();

}

#endif // %ifdef TIC_TOC_PROFILER_HPP
