#include <Eigen/Dense>
#include "Traits.h"
#include "Sharedprod.h"


#ifndef PEIGEN_SHARED_MATRIX_H
#define PEIGEN_SHARED_MATRIX_H

namespace peigen
{
	using namespace Eigen;	

	template <typename MatrixType>
	class SharedMatrix
	{
	protected:
		typedef typename MatrixType::Scalar Scalar;
		typedef MatrixType EigenMType;

		int nrows, ncols;
		
		int desc[9];

	public:
		MatrixType local_matrix;
		MatrixType global_matrix;

		bool use_transpose;
		int i,j,x,y;
		int nrblock, ncblock;

		Scalar *localData() {return local_matrix.data();}
		Scalar *globalData() {return global_matrix.data();}
		inline int rows() const {return nrows;}
		inline int cols() const {return ncols;}
		inline int rblock() const {return nrblock;}
		inline int cblock() const {return ncblock;}
		inline int *descriptor() {return desc;}

		// Constructors
		SharedMatrix<MatrixType>();
		SharedMatrix<MatrixType>(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock);
		SharedMatrix<MatrixType>(int global_M, int global_N, char init, int rBlockSize, int cBlockSize); 
		
		// Copy constructor
		SharedMatrix<MatrixType>(const SharedMatrix<MatrixType> &other);
		SharedMatrix<MatrixType>(Sharedprod<MatrixType> prod);

		// Destructor
		~SharedMatrix<MatrixType>() {}

		// Initializers
		static SharedMatrix<MatrixType> Ones(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'o', rBlockSize, cBlockSize); } 
		static SharedMatrix<MatrixType> Zeros(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'z', rBlockSize, cBlockSize); } 
		static SharedMatrix<MatrixType> Proc(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'p', rBlockSize, cBlockSize); } 
	
		// Operators
		SharedMatrix<MatrixType> &operator=(SharedMatrix<MatrixType> &other)
		{
			nrows = other.nrows;
			ncols = other.ncols; 
			nrblock = other.nrblock; 
			ncblock = other.ncblock;
			local_matrix = other.local_matrix;
			std::copy(other.desc, other.desc + 9, desc);
			other.clear();
			return *this;
		}

		SharedMatrix<MatrixType> &operator=(Sharedprod<MatrixType> prod)
		{
			SharedMatrix<MatrixType>  P = prod.eval();
			nrows = P.nrows;
			ncols = P.ncols; 
			nrblock = P.nrblock; 
			ncblock = P.ncblock;
			local_matrix = P.local_matrix;
			std::copy(P.desc, P.desc + 9, desc);
			prod.A.clear();
			prod.B.clear();
			return *this;
		}

		SharedMatrix<MatrixType> & transpose() {
			use_transpose = !use_transpose; 
			std::swap(x,y); 
			std::swap(nrblock, ncblock);
			return *this;
		}

		SharedMatrix<MatrixType> & block(int bi, int bj, int bx, int by) 
		{
			assert((((bi+bx) <= x)&&((by+bj) <= y)) && "Taking a block that is too large");
			i = bi+i;
			j = bj+j;
			x = bx;
			y = by;
			return *this;
		}

		SharedMatrix<MatrixType> & clear() 
		{
			if (use_transpose)
				transpose();
			i = 1;
			j = 1;
			x = rows();
			y = cols();
			return *this;
		}
	};







	// Constructors
	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix<MatrixType>()
	{
		use_transpose = false;
		i = j = 1;
		x = y = 0;


		// If none is specified, use the default block size of the context
		nrblock = BLACS::rblock;
		ncblock = BLACS::cblock;

		// Global size of the shared matrix
		nrows = 0;
		ncols = 0;

		// Neither .local_matrix nor .desc are set with the default constructor!
	}

	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix<MatrixType>(int global_M, int global_N, int rBlockSize, int cBlockSize)
		: nrows(global_M), ncols(global_N), nrblock(rBlockSize), ncblock(cBlockSize), use_transpose(false), i(1), j(1), x(global_M), y(global_N)
	{
		// Size of the data on the local process
		int m = BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows));
		int n = BLACS::numroc_(&ncols, &ncblock, &(BLACS::mycol), BLACS::iZERO, &(BLACS::grid_cols));

		local_matrix.resize(m, n);

		// Finally fill in the BLACS array descriptor
		int info;
		m = max(1, m);	// necessary because descinit() will throw if llda is 0

		BLACS::descinit_(desc, &nrows, &ncols, &nrblock, &ncblock, BLACS::iZERO, BLACS::iZERO, &(BLACS::ctxt), &m, &info);

		if (info!=0)
		{
			std::cout << "descinit return error code: " << info  << ", on process " << BLACS::myrank << std::endl;
			std::cout << "mn, llda on " << BLACS::myrank << " is " <<m <<", size is " << local_matrix.rows() <<"x"<<local_matrix.cols()<< ", " <<BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows)) <<endl <<flush;
		}
	}





	// pre-filled constructor
	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix<MatrixType>(int global_M, int global_N, char init, int rBlockSize, int cBlockSize)
		: nrows(global_M), ncols(global_N), nrblock(rBlockSize), ncblock(cBlockSize), use_transpose(false), i(1), j(1), x(global_M), y(global_N)
	{
		// Size of the data on the local process
		int m = BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows));
		int n = BLACS::numroc_(&ncols, &ncblock, &(BLACS::mycol), BLACS::iZERO, &(BLACS::grid_cols));

		switch (init)
		{
		case 'o':
			local_matrix = MatrixType::Ones(m, n);
			break;
		case 'z':
			local_matrix = MatrixType::Zero(m, n);
			break;
		case 'p':
			local_matrix = MatrixType::Ones(m, n) * BLACS::myrank;
			break;
		}

		// Finally fill in the BLACS array descriptor
		int info;
		m = max(1, m);	// necessary because descinit() will throw if llda is 0

		BLACS::descinit_(desc, &nrows, &ncols, &nrblock, &ncblock, BLACS::iZERO, BLACS::iZERO, &(BLACS::ctxt), &m, &info);

		if (info!=0)
		{
			std::cout << "descinit return error code: " << info  << ", on process " << BLACS::myrank << std::endl;
			std::cout << "mn, llda on " << BLACS::myrank << " is " <<m <<", size is " << local_matrix.rows() <<"x"<<local_matrix.cols()<< ", " <<BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows)) <<endl <<flush;
		}
	}
	






	// Copy contructor
	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix<MatrixType>(const SharedMatrix<MatrixType> &other)
		: nrows(other.x), ncols(other.y), nrblock(other.nrblock), ncblock(other.ncblock), use_transpose(other.use_transpose), i(1), j(1), x(other.x), y(other.y)
	{
		local_matrix = other.local_matrix;
		std::copy(other.desc, other.desc + 9, desc);
	}

	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix<MatrixType>(Sharedprod<MatrixType> prod) : i(1), j(1)
	{
		use_transpose = false;
		SharedMatrix<MatrixType>  P = prod.eval();
		nrows = x = P.nrows;
		ncols = y = P.ncols; 
		nrblock = P.nrblock; 
		ncblock = P.ncblock;
		local_matrix = P.local_matrix;
		std::copy(P.desc, P.desc + 9, desc);
		prod.A.clear();
		prod.B.clear();
	}







	// cout overload
	template <typename MatrixType>
	std::ostream& operator<<( std::ostream& os, const SharedMatrix<MatrixType>& M )
	{
		/* Printing the global matrix */
		for (int r = 0; r < M.rows(); r++)
		{
			for(int cb = 0; cb < ceil((double)M.cols()/M.cblock()); cb++)
			{
				int rb = floor((double)r/M.rblock()) ;
				if((BLACS::myrow==rb%BLACS::grid_rows) && (BLACS::mycol==cb%BLACS::grid_cols))
				{
					int i,j,l;
					i = M.rblock()*floor(r/(BLACS::grid_rows*M.rblock())) + r%M.rblock();
					j = M.cblock() * floor(cb/BLACS::grid_cols);
					l = min(M.cblock(), M.cols()-cb*M.cblock());
					os << M.local_matrix.block(i, j, 1, l) << "\t" << flush;
				}
				BLACS::Cblacs_barrier(BLACS::ctxt, "All");
			}
			if (BLACS::ROOT)
				os << endl << flush;
			BLACS::Cblacs_barrier(BLACS::ctxt, "All");
		}
		if (BLACS::ROOT)
			os << endl << flush;
		BLACS::Cblacs_barrier(BLACS::ctxt, "All");
		/* global matrix is printed */

		return os;
	}







	template <typename MatrixType>
	Sharedprod<MatrixType> operator*(SharedMatrix<MatrixType> &a, SharedMatrix<MatrixType> &b)
	{return Sharedprod<MatrixType>(a, b);}

	template <typename MatrixType>
	Sharedprod<MatrixType> operator*(Sharedprod<MatrixType> p, SharedMatrix<MatrixType> a)
	{SharedMatrix<MatrixType> partial_result = p.eval();
	return Sharedprod<MatrixType>(partial_result, a);}


	

	


	template<typename MatrixType>
	struct traits<SharedMatrix<MatrixType> >
	{
		enum 
		{
			transpose = false,
		};
	};

}	// end namespace peigen
#endif // PEIGEN_SHARED_MATRIX_H