#include <Eigen/Dense>
#include "BLACS.h"
#include "PBLAS.h"
#include "SharedMatrix.h"

#ifndef PEIGEN_VANDERMONDE_H
#define PEIGEN_VANDERMONDE_H

using namespace peigen;

template <typename MatrixType>
SharedMatrix<MatrixType> vander(MatrixType v, int order, int rblock, int cblock)
{
	SharedMatrix<MatrixType> V(v.rows(), order, rblock, cblock);

	MatrixType v_loc(V.local_matrix.rows(), 1);
	for (int i = 0; i < V.local_matrix.rows(); ++i)
	{
		const int g = BLACS::indxl2g(i, V.rblock(), BLACS::grid_rows, BLACS::myrow);
		v_loc(i,0) = v(g,0);
	}

	for (int j = 0; j < V.local_matrix.cols(); ++j)
	{
		const int g = BLACS::indxl2g(j, V.cblock(), BLACS::grid_cols, BLACS::mycol);
		V.local_matrix.col(j) = v_loc.array().pow(g);
	}

	return V;
}




#endif // PEIGEN_VANDERMONDE_H