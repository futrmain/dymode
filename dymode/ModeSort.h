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

// General template, implements ascending order by default
template <typename scalar, sort_order order>
struct s
{
	scalar d;
	int index;
	bool operator < (const struct s &other) const 
	{
		return d < other.d;
	}
};

// Specialization for descending order
template <typename scalar>
struct s<scalar, descend>
{
	scalar d;
	int index;
	bool operator < (const struct s &other) const 
	{
		return d > other.d;
	}
};

template<typename MatrixType>
class ModeSort
{
public:
	MatrixXi orderedIdx;

	ModeSort(SharedMatrix<MatrixType>& Modes, MatrixType& eigenvalues, MatrixType& norm, MatrixType& singulars, sort_method& method, int NMAX)
	{

	}

	MatrixXi energy_sort(MatrixType& eigenvalues, MatrixType& norm, sort_method& method, int NMAX)
	{
		vector<s<MatrixType::Scalar>, method.order> v;

		for (int i = 0; i < norm.cols(); ++i) 
		{
			s s_temp;
			if (eigenvalues(0, i).imag() < 0)
			{
				// Should we discard the mode ?
				if (method.conjugates == false)
					s_temp.d = -1;
				else
					s_temp.d = norm(0, i);
			}
			else
			{
				s_temp.d = norm(0, i);
			}

			MatrixType::Scalar correction;
			switch (method.energy_ref)
			{
			case -11:	// median
				correction = std::pow(eigenvalue(0, i).abs(), eigenvalues.cols() / 2);
				break;
			case -10:	// mean
				if (eigenvalue(0, i).abs() == 1)
					correction = 1;
				else
					correction = (1 - std::pow(eigenvalue(0, i).abs(), eigenvalues.cols())) / (eigenvalues.cols() * (1 - eigenvalue(0, i).abs()));
				break;
			default:	// snapshot number
				correction = std::pow(eigenvalue(0, i).abs(), method.energy_ref);
			}
			s_temp.d = s_temp.d * correction;
			s_temp.index = i;
			v.push_back(s);
		}


		std::partial_sort(v.begin(), v.begin() + NMAX, v.end());
		
		MatrixXi indices(1, NMAX);
		for(int i = 0; i < NMAX; ++i)
		{
			indices(0, i) = v[i].index;
		}

		return indices;
	}
};

#endif // MODE_SORT_H