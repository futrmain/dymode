#ifndef PEIGEN_SCAEIGENSOLVER_H
#define PEIGEN_SCAEIGENSOLVER_H

#include "ScaHessenberg.h"
#include "ScaSchur.h"

namespace peigen
{
	enum EigenMethod { EigSerial = 0, EigHess, EigSchur };

	template <typename MatrixType>
	class ScaEigenSolver
	{
	private:
		SharedMatrix<MatrixType> S;
		SharedMatrix<Matrix<complex<MatrixType::RealScalar>, Dynamic, Dynamic>> evectors;		
		Matrix<complex<MatrixType::RealScalar>, Dynamic, Dynamic> evalues;
		EigenMethod method_;

	public:
		inline MPI::Datatype ScaEigenSolver<MatrixType>::MPIType()
		{
			return MPI::DOUBLE_COMPLEX;
		}

		ScaEigenSolver<MatrixType>(const SharedMatrix<MatrixType>& A, const bool computeEigenVectors = true, const EigenMethod method = EigSchur);

		inline ScaEigenSolver<MatrixType> eigenvectors() { return evectors; }
		inline MatrixType eigenvalues() { return evalues; }
		inline setMethod(const EigenMethod& method) { method_ = method };

	};

	template <typename MatrixType>
	ScaEigenSolver<MatrixType>::ScaEigenSolver(const SharedMatrix<MatrixType>& A, const bool computeEigenVectors = true, const EigenMethod method = EigSchur) : S(A), method_(method)
	{
		assert(A.rows() == A.cols() && "CALLING EIGEN SOLVER ON NON SQUARE MATRIX");

		//A.gather(0);
		//std::cout << "B as input: " << endl << A.global_matrix << endl;

		//std::cout << "0" << endl;

		EigenSolver<MatrixType> eig;
		if (method_ == Eigen)
		{
			S.gather(0);

			if (BLACS::myrank == 0 )
				eig.compute(S.global_matrix, computeEigenVectors);
		}
		else
		{
			ScaHessenberg<MatrixType> hess(S, computeEigenVectors);

			if (method_ == EigHess)
			{
				SharedMatrix<MatrixType> H = hess.matrixH();
				H.gather(0);
				SharedMatrix<MatrixType> Q = hess.matrixQ();
				Q.gather(0);

				if (BLACS::myrank == 0)
					eig.computeFromHessenberg(H.global_matrix, Q.global_matrix, computeEigenVectors);
			}
			else if (method_ == EigSchur)
			{
				ScaSchur<MatrixType> schur(hess.matrixH(), hess.matrixQ(), computeEigenVectors, computeEigenVectors);

				SharedMatrix<MatrixType> Z = schur.matrixZ();
				Z.gather(0);
				SharedMatrix<MatrixType> T = schur.matrixT();
				T.gather(0);

				if (BLACS::myrank == 0)
					eig.computeFromSchur(T.global_matrix, Z.global_matrix, computeEigenVectors);
			}
			else
			{
				assert(0 && "The method to solve the Eigen problem is incorrect. Possible values are Eigen, EigHess and EigSchur.");
			}
		}

		if (BLACS::myrank == 0)
		{
			evalues = eig.eigenvalues();
			evectors.global_matrix = eig.eigenvectors();
		}

		evectors.dispatch(0, evectors.rblock(), evectors.cblock());
		BLACS::COMM_ACTIVE.Bcast(evalues.data(), evalues.rows(), MPIType(), 0);


		std::cout << "Eigen values from Eigen" << endl << evalues << endl << endl;
		std::cout << "Eigen vectors from Eigen" << endl << evectors << endl << endl;

		//ScaHessenberg<MatrixXd> hess(A);

		//std::cout << "Matrix H " << endl << hess.matrixH() << endl << endl;
		//std::cout << "Matrix Q " << endl << hess.matrixQ() << endl;

		//

		////MPI::COMM_WORLD.Barrier();

		//SharedMatrix<MatrixType> H = hess.matrixH();
		//H.gather(0);
		//SharedMatrix<MatrixType> Q = hess.matrixQ();
		//Q.gather(0);

		//std::cout << "Residuals after Hess: " << endl << Q.global_matrix * H.global_matrix * Q.global_matrix.transpose() - A.global_matrix<< endl;

		////EigenSolver<MatrixXd> eig(H.rows());
		////eig.computeFromHessenberg(H.global_matrix, Q.global_matrix, true);

		///*A.gather(0);

		//EigenSolver<MatrixXd> eig;
		//eig.setMaxIterations(1000000);
		//eig.compute(A.global_matrix);
		//std::cout << "Max number of iterations" << endl << eig.getMaxIterations()  << endl << endl;

		//std::cout << "Eigen values from Eigen" << endl << eig.eigenvalues() << endl << endl;*/


		//ScaSchur<MatrixXd> schur(hess.matrixH(), hess.matrixQ(), true);

		//std::cout << "Matrix T " << endl << schur.matrixT() << endl << endl;
		//std::cout << "Matrix Z " << endl << schur.matrixZ() << endl;

		//std::cout << "Eigen values " << endl << schur.eigenvals << endl << endl;


		//SharedMatrix<MatrixType> reconstruct(A);
		//reconstruct = schur.matrixZ() * schur.matrixT();
		//reconstruct = reconstruct * schur.matrixZ().transpose();
		//reconstruct.gather(0);
		//
		//std::cout << "Residuals after Schur: " << endl << reconstruct.global_matrix - A.global_matrix  << endl;

		//
		//SharedMatrix<MatrixType> Z = schur.matrixZ();
		//Z.gather(0);
		//SharedMatrix<MatrixType> T = schur.matrixT();
		//T.gather(0);
		//eig.computeFromSchur(T.global_matrix, Z.global_matrix, true);
		////eig.compute(A.global_matrix);

		//std::cout << "Eigen values from Eigen: " << endl << eig.eigenvalues() << endl << endl;
		//std::cout << "Eigen vectors from Eigen: " << endl << eig.eigenvectors() << endl << endl;

		//std::cout << "Residuals from Eigen problem: " << endl <<  endl << endl;
		//std::cout << A.global_matrix * eig.eigenvectors() - eig.eigenvectors() * eig.eigenvalues().asDiagonal() << endl;
	}

}	// end namespace peigen
#endif // PEIGEN_SCAEIGENSOLVER_H