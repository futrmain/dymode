#ifndef PEIGEN_SHARED_TRANSP_H
#define PEIGEN_SHARED_TRANSP_H

namespace peigen
{
	template <typename SharedType>
	class SharedTranspose : SharedMatrix<typename SharedType::EigenMType>
	{
	public:
		typedef typename SharedType::Scalar Scalar;
		typedef typename SharedType::EigenMType EigenMType;

		// Constructor
		inline SharedTranspose(SharedType& a_matrix) : matrix_t(a_matrix) {}

		inline int rows() const { return matrix_t.cols(); }
		inline int cols() const { return matrix_t.rows(); }
		inline int rblock() const { return matrix_t.cblock(); }
		inline int cblock() const { return matrix_t.rblock(); }

		SharedType matrix_t;
	};

}	// end namespace peigen
#endif // PEIGEN_SHARED_TRANSP_H