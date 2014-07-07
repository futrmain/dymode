#ifndef DYMODE_COLUMNSQUAREDNORM
#define DYMODE_COLUMNSQUAREDNORM


template <typename MatrixType>
Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic> ColumnNorm(peigen::SharedMatrix<MatrixType> A)
{
	peigen::SharedMatrix<Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic>> SNorms(1, A.cols(), 1, A.cblock());
	Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic> local_norm(0, 0);
	
	if (A.local_matrix.cols() > 0)
	{
		if (A.local_matrix.rows() > 0)
		{
			local_norm = A.local_matrix.colwise().stableNorm();

			int nrows = 0;
			if (BLACS::myrow == 0)
			{
				Eigen::Matrix<typename MatrixType::RealScalar, Eigen::Dynamic, Eigen::Dynamic, RowMajor> Buffer(BLACS::grid_rows, A.local_matrix.cols());
				Buffer.row(0) = local_norm;
				++nrows;
				for (int r = 1; r < BLACS::grid_rows; ++r)
				{
					if (BLACS::peigen_numroc(A.rows(), A.rblock(), r, 0 /*origin*/, BLACS::grid_rows) > 0)
					{
						BLACS::COMM_ACTIVE.Irecv(Buffer.row(nrows).data(), A.local_matrix.cols(), MPI::DOUBLE, BLACS::Cblacs_pnum(BLACS::ctxt, r, BLACS::mycol), 0/*tag*/);
						++nrows;
					}
				}
				BLACS::COMM_ACTIVE.Barrier();
				local_norm = Buffer.block(0, 0, nrows, A.local_matrix.cols()).colwise().stableNorm();
				SNorms.local_matrix = local_norm;
			}
			else
			{
				BLACS::COMM_ACTIVE.Send(local_norm.data(), local_norm.cols(), MPI::DOUBLE, BLACS::Cblacs_pnum(BLACS::ctxt, 0, BLACS::mycol), 0);
				BLACS::COMM_ACTIVE.Barrier();
			}
		}
		else 
		{
			BLACS::COMM_ACTIVE.Barrier();
		}
	}
	else
	{
		BLACS::COMM_ACTIVE.Barrier();
	}
	
	SNorms.gather(0);

	local_norm.resize(1, SNorms.cols());
	if (peigen::BLACS::myrank == 0)
	{
		local_norm = SNorms.global_matrix;
	}
	peigen::BLACS::COMM_ACTIVE.Bcast(local_norm.data(), local_norm.cols(), MPI::DOUBLE, 0);
	
	return local_norm;
}

#endif // DYMODE_COLUMNSQUAREDNORM
