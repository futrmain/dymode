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

class event
{
	public:
	//string name;
	thischrono::high_resolution_clock::time_point tic_time, toc_time;
	//vector<unsigned long long int> samples;
	vector<double> samples;
	map<string, event> *subEvents;
	map<string, event> *up;

	event()
	{
		subEvents = new map<string, event>;
		samples.reserve(256);
	}
};
 


class DoProfiler
{
	map<string, event> *child;
	map<string, event> *parent;
	string out_file;

	public:
		const double factor;
		map<string, event> root;
		DoProfiler() : factor((double)thischrono::high_resolution_clock::period::num / thischrono::high_resolution_clock::period::den)
		{
			out_file = "tic-toc-profiler.yml";
			parent = 0;
			child = &root;
		}

		DoProfiler(string fname) : factor((double)thischrono::high_resolution_clock::period::num / thischrono::high_resolution_clock::period::den)
		{
			out_file = fname;
			parent = 0;
			child = &root;
		}

		void tic(string key)
		{

			//(*child)[key].name = key;
			(*child)[key].up = parent;
			//(*child)[key].subEvents = new map<string, event>;
			parent = child;
			child = (*child)[key].subEvents;
			(*parent)[key].tic_time = thischrono::high_resolution_clock::now();
		}
		void toc(string key)
		{

			(*parent)[key].toc_time = thischrono::high_resolution_clock::now();
			(*parent)[key].samples.push_back(((*parent)[key].toc_time - (*parent)[key].tic_time).count() * factor);
			//cout << ((*parent)[key].toc_time - (*parent)[key].tic_time).count() << " " << (double)boost::chrono::steady_clock::period::num/boost::chrono::steady_clock::period::den  << endl;
			child = parent;
			parent = (*parent)[key].up;

		}
		void dump();
		
		void dump(string fname);
			

};

class NoProfiler
{
	map<string, event> *child;
	map<string, event> *parent;

	public:
		double factor;
		map<string, event> root;
		NoProfiler(string) {};

		void tic(string) {};
		void toc(string) {};
		void dump() {};
		void dump(string) {};
};

#ifdef USE_PROFILER
	typedef DoProfiler profiler;
#else
	typedef NoProfiler profiler;
#endif


YAML::Emitter& operator << (YAML::Emitter& out, const map<string, event>& events);
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

inline YAML::Emitter& operator << (YAML::Emitter& out, const map<string, event>& events)
{
	out << YAML::BeginMap;

	for (map<string, event>::const_iterator it = events.begin(); it != events.end(); it++)
	{
		out << YAML::Key << (*it).first;
		out << YAML::Value << YAML::BeginMap;
		out << YAML::Key << "samples";
		out << YAML::Value << YAML::Flow << (*it).second.samples;
		if ((*(*it).second.subEvents).size() > 0)
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

	yout << *this;
	ofstream fout(out_file.c_str());
	fout << yout.c_str();
	fout.close();

	//cout << yout.c_str();	

}
inline void DoProfiler::dump(string fname)
{

	YAML::Emitter yout;

	yout << *this;
	ofstream fout(fname.c_str());
	fout << yout.c_str();
	fout.close();

}

#endif // %ifdef TIC_TOC_PROFILER_HPP
