#ifndef PEIGEN_PBLAS_H
#define PEIGEN_PBLAS_H

namespace peigen
{
	namespace PBLAS
	{
		extern "C" {
			void pdgemv_(char*, int*, int*, double*, double*, int*, int*, int*, double* x, int* ix, int* jx, int* descx, int* incx, double* beta, double* y, int* iy, int* jy, int* descy, int* incy);

			void psgemv_(char*, int*, int*, float*, float*, int*, int*, int*, float* x, int* ix, int* jx, int* descx, int* incx, float* beta, float* y, int* iy, int* jy, int* descy, int* incy);

			void pcgemv_(char*, int*, int*, std::complex<float>*, std::complex<float>*, int*, int*, int*, std::complex<float>* x, int* ix, int* jx, int* descx, int* incx, std::complex<float>* beta, std::complex<float>* y, int* iy, int* jy, int* descy, int* incy);

			void pzgemv_(char*, int*, int*, std::complex<double>*, std::complex<double>*, int*, int*, int*, std::complex<double>* x, int* ix, int* jx, int* descx, int* incx, std::complex<double>* beta, std::complex<double>* y, int* iy, int* jy, int* descy, int* incy);




			void pdgemm_(char*transa, char*transb, int*m, int*n, int*k, double*alpha, double*a, int*ia, int*ja, int*desca, double*b, int*ib, int*jb, int*descb, double*beta, double*c, int*ic, int*jc, int*descc);

			void psgemm_(char*transa, char*transb, int*m, int*n, int*k, float*alpha, float*a, int*ia, int*ja, int*desca, float*b, int*ib, int*jb, int*descb, float*beta, float*c, int*ic, int*jc, int*descc);

			void pcgemm_(char*transa, char*transb, int*m, int*n, int*k, std::complex<float>*alpha, std::complex<float>*a, int*ia, int*ja, int*desca, std::complex<float>*b, int*ib, int*jb, int*descb, std::complex<float>*beta, std::complex<float>*c, int*ic, int*jc, int*descc);

			void pzgemm_(char*transa, char*transb, int*m, int*n, int*k, std::complex<double>*alpha, std::complex<double>*a, int*ia, int*ja, int*desca, std::complex<double>*b, int*ib, int*jb, int*descb, std::complex<double>*beta, std::complex<double>*c, int*ic, int*jc, int*descc);



			// FIXME it's not PBLAS, it's from ScaLAPACK
			void pdgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* ia, int* ja, int* desca, double* s, double* u, int* iu, int* ju, int* descu, double* vt, int* ivt, int* jvt, int* descvt, double* work, int* lwork, int* info);

			void pzgesv_(int *, int *, std::complex<double> *, int *, int *, int *, int *, std::complex<double> *, int *, int *, int *, int *);
			void pcgesv_(int *, int *, std::complex<float> *, int *, int *, int *, int *, std::complex<float> *, int *, int *, int *, int *);

			void pdgesv_(int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);
			void psgesv_(int *, int *, float *, int *, int *, int *, int *, float *, int *, int *, int *, int *);
		}

		using namespace std;
		/**************************************************************************************************/
		/*------------------------    OVERLOADING PBLAS/ScaLAPACK FUNCTIONs    ---------------------------*/
		/**************************************************************************************************/
		/////////////////////////				   P*GEMM   				//////////////////////////////
		inline void pxgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, int ia, int ja, int *desca, double *b, int ib, int jb, int *descb, double beta, double *c, int ic, int jc, int *descc)
		{
			pdgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
			return;
		}
		inline void pxgemm(char transa, char transb, int m, int n, int k, float alpha, float *a, int ia, int ja, int *desca, float *b, int ib, int jb, int *descb, float beta, float *c, int ic, int jc, int *descc)
		{
			psgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
			return;
		}
		inline void pxgemm(char transa, char transb, int m, int n, int k, complex<float> alpha, complex<float> *a, int ia, int ja, int *desca, complex<float> *b, int ib, int jb, int *descb, complex<float> beta, complex<float> *c, int ic, int jc, int *descc)
		{
			pcgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
			return;
		}
		inline void pxgemm(char transa, char transb, int m, int n, int k, complex<double> alpha, complex<double> *a, int ia, int ja, int *desca, complex<double> *b, int ib, int jb, int *descb, complex<double> beta, complex<double> *c, int ic, int jc, int *descc)
		{
			pzgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
			return;
		}


	}	// end namespace PBLAS
}	// end namespace peigen

#endif // PEIGEN_PBLAS_H