#pragma once
#ifndef MODE_SORT_H
#define MODE_SORT_H

#include <boost/algorithm/string.hpp>

enum sort_type { energy = 0 };
enum sort_order { ascend = 0, descend };

class sort_method
{
private:
	int foo;

public:
	sort_type stype;
	sort_order order;
	bool conjugates;
	int energy_ref;

	// Default constructor
	sort_method() : stype(energy), energy_ref(-11), order(descend), conjugates(false)
	{}

	sort_method(string method)
	{
		set(method);
	}

	void set(string method)
	{
		vector<string> input;

		boost::split(input, method, boost::is_any_of(","));

		// Get sort_type
		if (input[0] == "energy")
		{
			stype = energy;
			energy_ref = -11;
			order = descend;
			conjugates = false;

			// override defaults
			for (int k = 1; k < input.size(); ++k)
			{
				if (!input[k].compare("c"))
					conjugates = true;
				else if (!input[k].compare("i"))
					order = ascend;
				else if (!input[k].compare("mean"))
					energy_ref = -10;
				else if (!input[k].compare("median"))
					energy_ref = -11;
				else if (sscanf(input[k].c_str(), "%i", &foo) > 0)
					energy_ref = foo;
				else
				{
					if (BLACS::myrank == 0)
						cout << "Unknown option " << input[k] << " for the sorting method. Option ignored." << endl;
				}
			}
		}
		else
		{
			if (BLACS::myrank == 0)
				cout << "Sort type " << input[0] << " is not a valid type. Using defaults." << endl;
			stype = energy;
			energy_ref = -11;
			order = descend;
			conjugates = false;
		}
	}
};

class ModeSort
{
public:

};

#endif // MODE_SORT_H