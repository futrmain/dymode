#ifndef DYMODE_COLUMNSQUAREDNORM
#define DYMODE_COLUMNSQUAREDNORM


template <typename MatrixType>
Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic> ColumnSquaredNorm(peigen::SharedMatrix<MatrixType> A)
{
	peigen::SharedMatrix<Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic>> SNorms(1, A.cols(), 1, A.cblock());
	Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic> local_snorm;
	if (A.local_matrix.cols() > 0)
	{
		local_snorm = A.local_matrix.colwise().squaredNorm();
	}

	char scope[7] = { 'C', 'O', 'L', 'U', 'M', 'N', '\0' };
	char top[2] = { ' ', '\0' };
	// Sum the norm of columns by process column
	peigen::BLACS::Cdgsum2d(peigen::BLACS::ctxt, scope, top, 1, local_snorm.cols(), local_snorm.data(), local_snorm.rows(), 1, -1);

	if (BLACS::myrow == 0)
	{
		SNorms.local_matrix = local_snorm;
	}
	SNorms.gather(0);

	local_snorm.resize(1, SNorms.cols());
	if (peigen::BLACS::myrank == 0)
	{
		local_snorm = SNorms.global_matrix;
	}
	peigen::BLACS::COMM_ACTIVE.Bcast(local_snorm.data(), local_snorm.cols(), MPI::DOUBLE, 0);

	return local_snorm;
}

#endif // DYMODE_COLUMNSQUAREDNORM
