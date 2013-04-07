#ifndef PEIGEN_SCAEIGENSOLVE_H
#define PEIGEN_SCAEIGENSOLVE_H

namespace peigen
{

	template <typename MatrixType>
	class ScaEigenSolve
	{
	protected:
		typedef typename MatrixType::Scalar Scalar;

	public:
		SharedMatrix<MatrixType> Hessenberg, eigenvectors, eigenvalues;

		ScaEigenSolve<MatrixType>(SharedMatrix<MatrixType> A);
	};

	template <typename MatrixType>
	ScaEigenSolve<MatrixType>::ScaEigenSolve(SharedMatrix<MatrixType> A) : Hessenberg(A)
	{
		assert(A.rows() == A.cols() && "CALLING EIGEN SOLVER ON NON SQUARE MATRIX");

		int info;
		MatrixType work(1,1);
		MatrixType tau(Hessenberg.rows(),1);

		PBLAS::pxgehrd(Hessenberg.rows(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(), tau.data(), work.data(), /*lwork*/ -1, &info);

		work.resize(work(0,0), 1);

		PBLAS::pxgehrd(Hessenberg.rows(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(), tau.data(), work.data(), /*lwork*/ work.size(), &info);

		
		if (BLACS::myrank==0)
			std::cout << "SIZE FOR EIGEN SOLVER, lwork, ltau:" << work(0,0) << "\t"<<tau(0,0)<< std::endl;

		std::cout << Hessenberg;
	}

}	// end namespace peigen
#endif // PEIGEN_SCAEIGENSOLVE_H