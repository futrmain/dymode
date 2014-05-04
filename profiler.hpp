#ifndef PROFILERHPP
#define PROFILERHPP

#define BOOST_CHRONO_HEADER_ONLY
#define BOOST_ERROR_CODE_HEADER_ONLY

#define USE_PROFILER

#include <map>
#include <vector>
#include <boost/chrono.hpp>
#include <string>
#include "yaml-cpp/yaml.h"
#include <fstream>

using namespace std;

class event
{
	public:
	//string name;
	boost::chrono::steady_clock::time_point tic_time, toc_time;
	//vector<unsigned long long int> samples;
	vector<double> samples;
	map<string, event> *subEvents;
	map<string, event> *up;

	event();
};
 
event::event()
{
	subEvents = new map<string, event>;
	samples.reserve(256);
}

class DoProfiler
{
	map<string, event> *child;
	map<string, event> *parent;
	string out_file;

	public:
		double factor;
		map<string, event> root;
		DoProfiler();
		DoProfiler(string);

		void tic(string);
		void toc(string);
		void dump();
		void dump(string);



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

DoProfiler::DoProfiler()
{
	out_file = "timing.yml";
	parent = 0;
	child = &root;	
	factor = (double)boost::chrono::steady_clock::period::num/boost::chrono::steady_clock::period::den;
}

DoProfiler::DoProfiler(string fname)
{
	out_file = fname;
	parent = 0;
	child = &root;	
	factor = (double)boost::chrono::steady_clock::period::num/boost::chrono::steady_clock::period::den;
}

void DoProfiler::tic(string key)
{
	
	//(*child)[key].name = key;
	(*child)[key].up = parent;
	//(*child)[key].subEvents = new map<string, event>;
	parent = child;
	child = (*child)[key].subEvents;
	(*parent)[key].tic_time = boost::chrono::steady_clock::now();

}

void DoProfiler::toc(string key)
{

	(*parent)[key].toc_time = boost::chrono::steady_clock::now();
	(*parent)[key].samples.push_back((double)((*parent)[key].toc_time - (*parent)[key].tic_time).count() / 1e+9);
	//cout << ((*parent)[key].toc_time - (*parent)[key].tic_time).count() << " " << (double)boost::chrono::steady_clock::period::num/boost::chrono::steady_clock::period::den  << endl;
	child = parent;
	parent = (*parent)[key].up;	

}

YAML::Emitter& operator << (YAML::Emitter& out, const map<string, event>& events); 
YAML::Emitter& operator << (YAML::Emitter& out, const DoProfiler& prof); 

void DoProfiler::dump()
{

	YAML::Emitter yout;

	yout << *this;
	ofstream fout(out_file.c_str());
	fout << yout.c_str();
	fout.close();

//cout << yout.c_str();	

}

void DoProfiler::dump(string fname)
{

	YAML::Emitter yout;

	yout << *this;
	ofstream fout(fname.c_str());
	fout << yout.c_str();
	fout.close();	

}

YAML::Emitter& operator << (YAML::Emitter& out, const map<string, event>& events) 
{
	out << YAML::BeginMap;
	
	for (map<string, event>::const_iterator it = events.begin(); it !=  events.end(); it++)
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

YAML::Emitter& operator << (YAML::Emitter& out, const DoProfiler& prof) 
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

#endif
