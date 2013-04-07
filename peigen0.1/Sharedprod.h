#ifndef PEIGEN_SHARED_PROD_H
#define PEIGEN_SHARED_PROD_H

namespace peigen
{
	using namespace Eigen;

	template <typename MatrixType>
	class SharedMatrix;

	template <typename MatrixType>
	class Sharedprod
	{
	public:
		typedef typename MatrixType::Scalar Scalar;
		Scalar alpha;
		SharedMatrix<MatrixType> &A;
		SharedMatrix<MatrixType> &B;
		
		Sharedprod<MatrixType>(SharedMatrix<MatrixType> &a, SharedMatrix<MatrixType> &b) : A(a), B(b), alpha(1) {}
		Sharedprod<MatrixType>(Scalar x, Sharedprod<MatrixType> &P) : A(P.A), B(P.B), alpha(x*P.alpha) {}
		Sharedprod<MatrixType>(Sharedprod<MatrixType> &P, SharedMatrix<MatrixType> &b) : A(P.eval()), B(b), alpha(1) {std::cout << A;}

		SharedMatrix<MatrixType> eval();
	};

	template <typename MatrixType>
	SharedMatrix<MatrixType> Sharedprod<MatrixType>::eval()
	{
		// Allocate result matrix
		int m, n, k, rb, cb;
		char opA = (A.use_transpose? 'T' : 'N');
		char opB = (B.use_transpose? 'T' : 'N');

		rb = A.rblock();
		cb = B.cblock();

		assert((A.y==B.x) && "Multiplying matrices of different size");
		SharedMatrix<MatrixType> P(A.x, B.y, rb, cb);

		// Define submatrices
		int iA = A.i;
		int jA = A.j;
		int iB = B.i;
		int jB = B.j;
		int iP = P.i;
		int jP = P.j;

		

		/*PBLAS::pxgemm(opA, opB, m, n, k, (Scalar)1, 
			A.localData(), iA, jA, A.descriptor(),		// matrix A
			B.localData(), iB, jB, B.descriptor(),	// B
			(Scalar)0, P.localData(), iP, jP, P.descriptor());	// resulting temporary product	*/
		PBLAS::pxgemm(opA, opB, A.x, B.y, A.y, alpha, 
			A.localData(), iA, jA, A.descriptor(),		// matrix A
			B.localData(), iB, jB, B.descriptor(),	// B
			(Scalar)0, P.localData(), iP, jP, P.descriptor());	// resulting temporary product

		A.clear();
		B.clear();
		return P;
	}


}	// end namespace peigen
#endif // PEIGEN_SHARED_PROD_H