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

template<typename MatrixType>
class ModeSort
{
public:
	MatrixXi orderedIdx;

	ModeSort(){};
	ModeSort(SharedMatrix<MatrixType>& Modes, MatrixType& eigenvalues, Matrix<typename MatrixType::RealScalar, Dynamic, Dynamic>& norm, Matrix<typename MatrixType::RealScalar, Dynamic, Dynamic>& singulars, sort_method& method, int NMAX)
	{
		const int nmax_ = min(NMAX, (int)norm.cols());
		orderedIdx = getOrderedIdx(Modes, eigenvalues, norm, singulars, method, nmax_);
	}

	MatrixXi getOrderedIdx(SharedMatrix<MatrixType>& Modes, MatrixType& eigenvalues, Matrix<typename MatrixType::RealScalar, Dynamic, Dynamic>& norm, Matrix<typename MatrixType::RealScalar, Dynamic, Dynamic>& singulars, sort_method& method, int NMAX)
	{
		if (method.stype == energy)
		{
			return energy_sort(eigenvalues, norm, method, NMAX);
		}
		else //... implement other types of sorting
			return energy_sort(eigenvalues, norm, method, NMAX);
	}

	MatrixXi energy_sort(const MatrixType& eigenvalues, const Matrix<typename MatrixType::RealScalar, Dynamic, Dynamic>& norm, const sort_method& method, const int& NMAX)
	{
		vector<pair<double, int>> v(norm.cols());

		for (int i = 0; i < norm.cols(); ++i)
		{
			v[i].second = i;

			if (eigenvalues(i, 0).imag() < 0)
			{
				// Should we discard the mode ?
				if (method.conjugates == false)
				{
					if (method.order == descend)
						v[i].first = -1;
					else
						v[i].first = std::numeric_limits<double>::infinity();
				}
				else
					v[i].first = norm(0, i);
			}
			else
			{
				v[i].first = norm(0, i);
			}

			typename MatrixType::RealScalar correction;
			switch (method.energy_ref)
			{
			case -11:	// median
				correction = std::pow(abs(eigenvalues(i, 0)), eigenvalues.rows() / 2);
				break;
			case -10:	// mean
				if (abs(eigenvalues(i, 0)) == 1)
					correction = 1;
				else
					correction = (1 - std::pow(abs(eigenvalues(i, 0)), eigenvalues.rows())) / (eigenvalues.rows() * (1 - abs(eigenvalues(i, 0))));
				break;
			default:	// snapshot number
				correction = std::pow(abs(eigenvalues(i, 0)), method.energy_ref);
			}
			v[i].first = v[i].first * correction;
		}

		if (method.order == descend)
			std::partial_sort(v.begin(), v.begin() + NMAX, v.end(), greater<pair<double, int>>());
		else
			std::partial_sort(v.begin(), v.begin() + NMAX, v.end(), less<pair<double, int>>());

		MatrixXi indices(1, NMAX);
		for (int i = 0; i < NMAX; ++i)
		{
			indices(0, i) = v[i].second;
		}

		return indices;
	}
};

#endif // MODE_SORT_H