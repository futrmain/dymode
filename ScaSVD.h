#ifndef PEIGEN_SCASVD_H
#define PEIGEN_SCASVD_H

namespace peigen
{
	template <typename MatrixType>
	class ScaSVD
	{
	private:
	  typedef typename MatrixType::RealScalar RealScalar;
	public:
		MatrixType singularValues;
		SharedMatrix<MatrixType> matrixU, matrixVt;

		ScaSVD(SharedMatrix<MatrixType> M, bool computeU, bool computeVt);

		RealScalar residual(SharedMatrix<MatrixType> original, bool normalized = true)
		{
			SharedMatrix<MatrixType> A(original);
			SharedMatrix<MatrixType> R(matrixU);

			// Multiplication by diagonal elements
			for (int j = 0; j < R.local_matrix.cols(); ++j)
			{
				const int g = BLACS::indxl2g(j, R.cblock(), BLACS::grid_cols, BLACS::mycol);
				R.local_matrix.col(j).noalias() = R.local_matrix.col(j) * singularValues(g, 0);
			}

			A.pgemm(1., R, matrixVt, -1.);

			if (normalized == true)
			{
				// Seems to produce +Inf results
				//A.local_matrix = A.local_matrix.array() / original.local_matrix.array();
			}

			//cout << "Residual MATRIX from SVD RESIDUAL(): " << endl << original << endl << flush;


			// Compute local highest residual
			return (A.localBlock().rows() * A.localBlock().cols()) > 0 ? A.localBlock().cwiseAbs().maxCoeff() : -1;
		}

		RealScalar global_residual(SharedMatrix<MatrixType> original)
		{
			double r_loc = residual(original);
			double r;

			BLACS::COMM_ACTIVE.Allreduce(&r_loc, &r, 1, MPI::DOUBLE, MPI::MAX);

			return r;
		}

	};

	template <typename MatrixType>
	ScaSVD<MatrixType>::ScaSVD(SharedMatrix<MatrixType> M, bool computeU, bool computeVt)
	{
		assert((M.rblock() == M.cblock()) && "SVD: 'M' must be distributed using a square block-cyclic distribution");
		//M.printDetails();

		// Note: MKL only computes the slim matrices or nothing, no full matrix
		char getU, getVt;
		getU = (computeU ? 'V' : 'N');
		getVt = (computeVt ? 'V' : 'N');

		int nvals = min(M.x, M.y);

		//singularValues.resize(nvals, 1);
		singularValues = MatrixType::Zero(nvals, 1);

		if (computeU)
			//matrixU.resize(M.x, nvals, M.rblock(), M.cblock());
			matrixU = SharedMatrix<MatrixType>(M.x, nvals, 'z', M.rblock(), M.cblock());
		else
			matrixU.resize(0, 0);

		if (computeVt)
			//matrixVt.resize(nvals, M.y, M.rblock(), M.cblock());
			matrixVt = SharedMatrix<MatrixType>(nvals, M.y, 'z', M.rblock(), M.cblock());
		else
			matrixVt.resize(0, 0);

		double min_workspace_length = 0;

		// calling pdgesvd with lwork:=-1 only computes the minimum size necessary for the array "work"
		// the computed value is returned as the first element in "dlwork"
		int info = 0;
		PBLAS::pxgesvd(getU, getVt, M.x, M.y, M.localData(), M.i, M.j, M.descriptor(), singularValues.data(), matrixU.local_matrix.data(), matrixU.i, matrixU.j, matrixU.descriptor(), matrixVt.localData(), matrixVt.i, matrixVt.j, matrixVt.descriptor(), &min_workspace_length, /*do a space query*/-1, &info);

		
		int lwork = (int) min_workspace_length;
		//MatrixType workspace(lwork, 1);
		MatrixType workspace = MatrixType::Zero(lwork, 1);

		// calling pdgesvd the allocated workspace "work" and its length in lwork
		// the data in M will be destroyed after this call
		PBLAS::pxgesvd(getU, getVt, M.x, M.y, M.localData(), M.i, M.j, M.descriptor(), singularValues.data(), matrixU.localData(), matrixU.i, matrixU.j, matrixU.descriptor(), matrixVt.localData(), matrixVt.i, matrixVt.j, matrixVt.descriptor(), workspace.data(), lwork, &info);

		
		/*if(BLACS::myrank == 24)
			std::cout << "On (" << BLACS::myrank << "), M equals " << M.local_matrix <<endl <<flush;
		BLACS::Cblacs_barrier(BLACS::ctxt, "All");*/
	}


}	// end namespace peigen
#endif // PEIGEN_SCASVD_H
