#include <Eigen/Dense>
#include "BLACS.h"
#include "PBLAS.h"
#include "Sharedprod.h"
#include <typeinfo>


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
		SharedMatrix();
		SharedMatrix(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock);
		SharedMatrix(int global_M, int global_N, char init, int rBlockSize, int cBlockSize); 
		
		// Copy constructor
		SharedMatrix(const SharedMatrix<MatrixType> &other);
		SharedMatrix(Sharedprod<MatrixType> prod);

		// Destructor
		~SharedMatrix<MatrixType>() {}

		// Initializers
		static SharedMatrix<MatrixType> Ones(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'o', rBlockSize, cBlockSize); } 
		static SharedMatrix<MatrixType> Zeros(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'z', rBlockSize, cBlockSize); } 
		static SharedMatrix<MatrixType> Proc(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock) { return SharedMatrix<MatrixType>(global_M, global_N, 'p', rBlockSize, cBlockSize); } 
	
		// Members
		void resize(int global_M, int global_N, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock);
		void printDetails();
		void dispatch(int source, int rBlockSize = BLACS::rblock, int cBlockSize = BLACS::cblock);
		void gather(int sink);
		inline MPI::Datatype MPIType();

		// Classic call to ScaLapack pxgemm
		void pgemm(double alpha, SharedMatrix<MatrixType> A, SharedMatrix<MatrixType> B, double beta)
		{
			// Allocate result matrix
			int m, n, k, rb, cb;

			char opA = (A.use_transpose ? 'T' : 'N');	// FIXME need to distinguish transpose and adjoint
			char opB = (B.use_transpose ? 'T' : 'N');

			rb = A.rblock();
			cb = B.cblock();

			SharedMatrix<MatrixType>& P = *this;

			assert(((A.y == B.x) && (B.y == P.y ) && (A.x == P.x)) && "Multiplying matrices of different size");

			// Define submatrices
			int iA = A.i;
			int jA = A.j;
			int iB = B.i;
			int jB = B.j;
			int iP = P.i;
			int jP = P.j;

			
			PBLAS::pxgemm(opA, opB, A.x, B.y, A.y, alpha,
				A.localData(), iA, jA, A.descriptor(),		// matrix A
				B.localData(), iB, jB, B.descriptor(),	// B
				beta, P.localData(), iP, jP, P.descriptor());	// resulting temporary product

		}


		// FIXME this should account for transpose flags somehow
		MatrixType localBlock()
		{
			int i_loc = 0;
			while (BLACS::indxl2g(i_loc, this->rblock(), BLACS::grid_rows, BLACS::myrow) < this->i - 1)
				++i_loc;

			int j_loc = 0;
			while (BLACS::indxl2g(j_loc, this->cblock(), BLACS::grid_cols, BLACS::mycol) < this->j - 1)
				++j_loc;

			int li_loc = this->local_matrix.rows() - i_loc;
			while (BLACS::indxl2g(i_loc + li_loc, this->rblock(), BLACS::grid_rows, BLACS::myrow) > this->i + this->x - 1)
				--li_loc;

			int lj_loc = this->local_matrix.cols() - j_loc;
			while (BLACS::indxl2g(j_loc + lj_loc, this->cblock(), BLACS::grid_cols, BLACS::mycol) > this->j + this->y - 1)
				--lj_loc;

			return this->local_matrix.block(i_loc, j_loc, li_loc, lj_loc);
		}



		// Operators
		SharedMatrix<MatrixType> &operator=(SharedMatrix<MatrixType> &other)
		{
			nrows = other.nrows;
			ncols = other.ncols; 
			i = other.i;
			j = other.j;
			x = other.x;
			y = other.y;
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
			i = P.i;
			j = P.j;
			x = P.x;
			y = P.y;
			nrblock = P.nrblock; 
			ncblock = P.ncblock;
			local_matrix = P.local_matrix;
			std::copy(P.desc, P.desc + 9, desc);
			//prod.A.clear();
			//prod.B.clear();
			return *this;
		}





		// Transformations
		SharedMatrix<MatrixType> & transpose() {
			use_transpose = !use_transpose; 
			std::swap(x,y); 
			std::swap(nrblock, ncblock);
			return *this;
		}

		SharedMatrix<MatrixType> transposeInPlace() {
			std::swap(x,y); 
			std::swap(nrblock, ncblock);
			std::swap(nrows, ncols);
			local_matrix.transposeInPlace();
			return *this;
		}

		SharedMatrix<MatrixType> & block(int bi, int bj, int bx, int by) 
		{
			//assert((((bi+bx) <= x)&&((by+bj) <= y)) && "Taking a block that is too large");
			i = bi+i;
			//std::cout << " j = " << j << ", bj = " << bj << std::endl;
			j = bj+j;
			//std::cout << " j = " << j << std::endl;
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

		// Cast
		template<typename newScalar>
		SharedMatrix<Matrix<newScalar, Dynamic, Dynamic> > cast()
		{
			SharedMatrix<Matrix<newScalar, Dynamic, Dynamic> > M(nrows, ncols, nrblock, ncblock);
			M.local_matrix = local_matrix.cast<newScalar>();
			M.global_matrix = global_matrix.cast<newScalar>();
			return M;
		}
	};







	// Constructors
	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix()
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
	SharedMatrix<MatrixType>::SharedMatrix(int global_M, int global_N, int rBlockSize, int cBlockSize)
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
	SharedMatrix<MatrixType>::SharedMatrix(int global_M, int global_N, char init, int rBlockSize, int cBlockSize)
		: nrows(global_M), ncols(global_N), nrblock(rBlockSize), ncblock(cBlockSize), use_transpose(false), i(1), j(1), x(global_M), y(global_N)
	{
		// Size of the data on the local process
		int m = BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows));
		int n = BLACS::numroc_(&ncols, &ncblock, &(BLACS::mycol), BLACS::iZERO, &(BLACS::grid_cols));

		switch (init)
		{
		case 'o':	//ones
			local_matrix = MatrixType::Ones(m, n);
			break;
		case 'z':	//zeros
			local_matrix = MatrixType::Zero(m, n);
			break;
		case 'p':	//process id
			local_matrix = MatrixType::Ones(m, n) * BLACS::myrank;
			break;
		case 'r':	//random
			local_matrix = MatrixType::Random(m, n);
			break;
		case 'i':	//identity
			local_matrix = MatrixType::Zero(m, n);
			for (int k=0; k < min(global_M, global_N); ++k)
			{
				if ((BLACS::myrow == BLACS::indxg2p(k, rBlockSize, BLACS::grid_rows)) && (BLACS::mycol == BLACS::indxg2p(k, cBlockSize, BLACS::grid_cols)))
					{
					local_matrix(BLACS::indxg2l(k, rBlockSize, BLACS::grid_rows), BLACS::indxg2l(k, cBlockSize, BLACS::grid_cols)) = 1;
					}				
			}
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
	SharedMatrix<MatrixType>::SharedMatrix(const SharedMatrix<MatrixType> &other)
		: nrows(other.nrows), ncols(other.ncols), nrblock(other.nrblock), ncblock(other.ncblock), use_transpose(other.use_transpose), i(other.i), j(other.j), x(other.x), y(other.y)
	{
		local_matrix = other.local_matrix;
		std::copy(other.desc, other.desc + 9, desc);
	}

	template <typename MatrixType>
	SharedMatrix<MatrixType>::SharedMatrix(Sharedprod<MatrixType> prod) : i(1), j(1)
	{
		use_transpose = false;
		SharedMatrix<MatrixType>  P = prod.eval();
		nrows = x = P.nrows;
		ncols = y = P.ncols; 
		nrblock = P.nrblock; 
		ncblock = P.ncblock;
		local_matrix = P.local_matrix;
		std::copy(P.desc, P.desc + 9, desc);
		//prod.A.clear();
		//prod.B.clear();
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



	// Resize
	template <typename MatrixType>
	void SharedMatrix<MatrixType>::resize(int global_M, int global_N, int rBlockSize, int cBlockSize)
	{
		// Global size of the shared matrix
		nrows = global_M;
		x = nrows;
		ncols = global_N;
		y = ncols;
		nrblock = rBlockSize;
		ncblock = cBlockSize;

		// Size of the data on the local process
		int m = BLACS::numroc_(&nrows, &nrblock, &(BLACS::myrow), BLACS::iZERO, &(BLACS::grid_rows));
		int n = BLACS::numroc_(&ncols, &ncblock, &(BLACS::mycol), BLACS::iZERO, &(BLACS::grid_cols));
		
		local_matrix.resize(m, n);

		// Finally fill in the BLACS array descriptor
		int info;
		m = max(1, m);	// necessary because descinit will throw if llda is 0

		BLACS::descinit_(desc, &nrows, &ncols, &nrblock, &ncblock, BLACS::iZERO, BLACS::iZERO, &(BLACS::ctxt), &m, &info);

		if (info!=0)
		{
			std::cout << "descinit return error code: " << info  << ", on process " << BLACS::myrank << std::endl;
		}
	}


	// Details
	template <typename MatrixType>
	void SharedMatrix<MatrixType>::printDetails()
	{
		if (BLACS::ROOT)
			std::cout << "My global size is " << nrows << " x " << ncols << endl
			<< "My block size is " << i << ","<< j << "\t"<< x <<","<< y << endl
			<< "My blocking factors are " << nrblock << " and " << ncblock << endl
			<< "I am " << use_transpose << " transposed" << endl
			<< "My descriptor values are: " << endl 
			<< "dense matrix " << desc[0] << endl
			<< "BLACS context " << desc[1] << endl
			<< "m " << desc[2] << endl
			<< "n " << desc[3] << endl
			<< "mb " << desc[4] << endl
			<< "nb " << desc[5] << endl
			<< "rsrc " << desc[6] << endl
			<< "csrc " << desc[7] << endl
			<< "lld " << desc[8] << endl
			<< endl << flush;
	}


	template <typename MatrixType>
	Sharedprod<MatrixType> operator*(SharedMatrix<MatrixType> a, SharedMatrix<MatrixType> b)
	{return Sharedprod<MatrixType>(a, b);}

	template <typename MatrixType>
	Sharedprod<MatrixType> operator*(Sharedprod<MatrixType> p, SharedMatrix<MatrixType> a)
	{SharedMatrix<MatrixType> partial_result = p.eval();
	return Sharedprod<MatrixType>(partial_result, a);}

	// MPI gather / dispatch
	template <typename Derived>
	inline MPI::Datatype SharedMatrix<Derived>::MPIType()
	{
		if (typeid(local_matrix.data()) == typeid(double *))
		{
			return MPI::DOUBLE;
		}
		else
		{
			return MPI::DOUBLE_COMPLEX;
		}
	}

	template <typename MatrixType>
	void SharedMatrix<MatrixType>::gather(int sink)
	{
		if (BLACS::myrank==sink)
		{
			//cout << "the sink has started gathering" << endl;
			// Create a matrix of buffer-matrices to receive into
			Matrix<MatrixType, Dynamic, Dynamic> buffer(BLACS::grid_rows, BLACS::grid_cols);
			for (int proc_row = 0; proc_row < BLACS::grid_rows; proc_row++)
			{
				for (int proc_col = 0; proc_col < BLACS::grid_cols; proc_col++)
				{
					int pid = BLACS::Cblacs_pnum((BLACS::ctxt), proc_row, proc_col);
					if (pid != sink)
					{
						int m = BLACS::numroc_(&nrows, &nrblock, &proc_row, BLACS::iZERO, &(BLACS::grid_rows));
						int n = BLACS::numroc_(&ncols, &ncblock, &proc_col, BLACS::iZERO, &(BLACS::grid_cols));

						buffer(proc_row, proc_col).resize(m, n);
						//cout << BLACS::myrank << " is receiving "<< m*n << " values." << endl;
						if (m*n > 0)
							BLACS::COMM_ACTIVE.Recv(buffer(proc_row, proc_col).data(), m*n, MPIType(), pid, 1);
					}
					else
					{
						buffer(proc_row, proc_col) = local_matrix;
					}
				}
			}
			//cout << "is it OK here?"<<endl;
			/* Reorder the block-cycliced buffer into the global matrix */
			global_matrix.resize(nrows, ncols);

			for (int rb = 0; rb < ceil((double)nrows/nrblock); rb++)
			{
				for(int cb = 0; cb < ceil((double)ncols/ncblock); cb++)
				{
					int roffset,coffset,_nrows,_ncols;
					roffset = nrblock * floor(rb / BLACS::grid_rows);
					coffset = ncblock * floor(cb / BLACS::grid_cols);
					_nrows = min(nrblock, nrows-rb*nrblock);
					_ncols = min(ncblock, ncols-cb*ncblock);

					//cout << "rb-cb: " << rb<<"-"<<cb << ", nrblock-ncblock: " << nrblock<<"-"<<ncblock << ", _nrows-_ncols: " << _nrows<<"-"<<_ncols << ", roffset-coffset: " << roffset<<"-"<<coffset << endl;

					global_matrix.block(rb*nrblock, cb*ncblock,  _nrows, _ncols) = 
						buffer(rb%BLACS::grid_rows, cb%BLACS::grid_cols).block(roffset, coffset, _nrows, _ncols);
				}
			}
		}
		else
		{
			//cout << BLACS::myrank << " is sending "<< local_matrix.size() << "values." << endl;
			if (local_matrix.size() > 0)
				BLACS::COMM_ACTIVE.Send(localData(), local_matrix.size(), MPIType(), sink, 1);
			//cout << BLACS::myrank << " has sent"<<endl;
		}
	}

	template <typename MatrixType>
	void SharedMatrix<MatrixType>::dispatch(int source, int rBlockSize, int cBlockSize)
	{
		int size[2];
		if (BLACS::myrank==source)
		{
			size[0] = global_matrix.rows();
			size[1] = global_matrix.cols();
		}
		BLACS::COMM_ACTIVE.Bcast(&size, 2, MPI::INT, source); // FIXME check that MPI and BLACS ranks are identical, it may cause trouble...

		//cout << "size " << size[0] << " " << size[1] << endl;
		resize(size[0], size[1], rBlockSize, cBlockSize);

		if (BLACS::myrank==source)
		{
			// Prepare a send buffer
			Matrix<MatrixType, Dynamic, Dynamic> buffer(BLACS::grid_rows, BLACS::grid_cols);
			for (int proc_row = 0; proc_row < BLACS::grid_rows; proc_row++)
			{
				for (int proc_col = 0; proc_col < BLACS::grid_cols; proc_col++)
				{
					int m = BLACS::numroc_(&nrows, &nrblock, &proc_row, BLACS::iZERO, &(BLACS::grid_rows));
					int n = BLACS::numroc_(&ncols, &ncblock, &proc_col, BLACS::iZERO, &(BLACS::grid_cols));
					buffer(proc_row, proc_col).resize(m, n);
				}
			}

			for (int rb = 0; rb < ceil((double)nrows/nrblock); rb++)
			{
				for(int cb = 0; cb < ceil((double)ncols/ncblock); cb++)
				{
					//cout << "rb cb  " << rb << " " << cb << endl;
					int roffset,coffset,_nrows,_ncols;
					roffset = nrblock * floor(rb / BLACS::grid_rows);
					coffset = ncblock * floor(cb / BLACS::grid_cols);
					_nrows = min(nrblock, nrows-rb*nrblock);
					_ncols = min(ncblock, ncols-cb*ncblock);
					//cout << "roffset, coffset, nrows, ncols " << roffset<< coffset<< nrows<< ncols << endl;
					//cout << "rb*rblock, cb*cblock,  nrows, ncols " << rb*rblock<< " " <<cb*cblock << endl;
					buffer(rb%BLACS::grid_rows, cb%BLACS::grid_cols).block(roffset, coffset, _nrows, _ncols) =
						global_matrix.block(rb*nrblock, cb*ncblock,  _nrows, _ncols);			
				}
			}
			for (int proc_row = 0; proc_row < BLACS::grid_rows; proc_row++)
			{
				for (int proc_col = 0; proc_col < BLACS::grid_cols; proc_col++)
				{
					int pid = BLACS::Cblacs_pnum(BLACS::ctxt, proc_row, proc_col);
					if (pid != source)				
					{
						BLACS::COMM_ACTIVE.Send(buffer(proc_row, proc_col).data(), buffer(proc_row, proc_col).size(), MPIType(), pid, 1);
					}
					else
					{
						local_matrix = buffer(proc_row, proc_col);
					}
				}
			}
		}
		else
		{
			BLACS::COMM_ACTIVE.Recv(localData(), local_matrix.size(), MPIType(), source, 1);
		}

		// free some memory
		global_matrix.resize(0,0);	//technically only source should do that, but whatever
	}



}	// end namespace peigen
#endif // PEIGEN_SHARED_MATRIX_H