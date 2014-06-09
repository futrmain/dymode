#ifndef PEIGEN_PBLAS_H
#define PEIGEN_PBLAS_H

namespace peigen
{
	namespace PBLAS
	{
		extern "C" {

			using std::complex;

			void pdgemv_(char*, int*, int*, double*, double*, int*, int*, int*, double* x, int* ix, int* jx, int* descx, int* incx, double* beta, double* y, int* iy, int* jy, int* descy, int* incy);

			void psgemv_(char*, int*, int*, float*, float*, int*, int*, int*, float* x, int* ix, int* jx, int* descx, int* incx, float* beta, float* y, int* iy, int* jy, int* descy, int* incy);

			void pcgemv_(char*, int*, int*, std::complex<float>*, std::complex<float>*, int*, int*, int*, std::complex<float>* x, int* ix, int* jx, int* descx, int* incx, std::complex<float>* beta, std::complex<float>* y, int* iy, int* jy, int* descy, int* incy);

			void pzgemv_(char*, int*, int*, std::complex<double>*, std::complex<double>*, int*, int*, int*, std::complex<double>* x, int* ix, int* jx, int* descx, int* incx, std::complex<double>* beta, std::complex<double>* y, int* iy, int* jy, int* descy, int* incy);




			void pdgemm_(char*transa, char*transb, int*m, int*n, int*k, double*alpha, double*a, int*ia, int*ja, int*desca, double*b, int*ib, int*jb, int*descb, double*beta, double*c, int*ic, int*jc, int*descc);

			void psgemm_(char*transa, char*transb, int*m, int*n, int*k, float*alpha, float*a, int*ia, int*ja, int*desca, float*b, int*ib, int*jb, int*descb, float*beta, float*c, int*ic, int*jc, int*descc);

			void pcgemm_(char*transa, char*transb, int*m, int*n, int*k, std::complex<float>*alpha, std::complex<float>*a, int*ia, int*ja, int*desca, std::complex<float>*b, int*ib, int*jb, int*descb, std::complex<float>*beta, std::complex<float>*c, int*ic, int*jc, int*descc);

			void pzgemm_(char*transa, char*transb, int*m, int*n, int*k, std::complex<double>*alpha, std::complex<double>*a, int*ia, int*ja, int*desca, std::complex<double>*b, int*ib, int*jb, int*descb, std::complex<double>*beta, std::complex<double>*c, int*ic, int*jc, int*descc);



			// FIXME it's not PBLAS, it's from ScaLAPACK
			void psgesvd_(char* jobu, char* jobvt, int* m, int* n, float* a, int* ia, int* ja, int* desca, float* s, float* u, int* iu, int* ju, int* descu, float* vt, int* ivt, int* jvt, int* descvt, float* work, int* lwork, int* info);

			void pdgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* ia, int* ja, int* desca, double* s, double* u, int* iu, int* ju, int* descu, double* vt, int* ivt, int* jvt, int* descvt, double* work, int* lwork, int* info);

			void pcgesvd_(char* jobu, char* jobvt, int* m, int* n, std::complex<float>* a, int* ia, int* ja, int* desca, std::complex<float>* s, std::complex<float>* u, int* iu, int* ju, int* descu, std::complex<float>* vt, int* ivt, int* jvt, int* descvt, std::complex<float>* work, int* lwork, int* info);

			void pzgesvd_(char* jobu, char* jobvt, int* m, int* n, std::complex<double>* a, int* ia, int* ja, int* desca, std::complex<double>* s, std::complex<double>* u, int* iu, int* ju, int* descu, std::complex<double>* vt, int* ivt, int* jvt, int* descvt, std::complex<double>* work, int* lwork, int* info);



			void pzgesv_(int *, int *, std::complex<double> *, int *, int *, int *, int *, std::complex<double> *, int *, int *, int *, int *);

			void pcgesv_(int *, int *, std::complex<float> *, int *, int *, int *, int *, std::complex<float> *, int *, int *, int *, int *);

			void pdgesv_(int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);

			void psgesv_(int *, int *, float *, int *, int *, int *, int *, float *, int *, int *, int *, int *);

			void pzgesvx_(char* fact, char* trans, int* n, int* nrhs, 
				std::complex<double> *a, int* ia, int* ja, int* desca, 
				std::complex<double> *af, int* iaf, int* jaf, int* descaf, 
				int *ipiv, char* equed, double *r, double *c, 
				std::complex<double> *b, int* ib, int* jb, int* descb, 
				std::complex<double> *x, int* ix, int* jx, int* descx, 
				double *rcond, double *ferr, double *berr, 
				std::complex<double> *work, int* lwork, double* rwork, int* lrwork, int *info);


			//////////////////// HESSENBERG REDUCTION
			void psgehrd_(int *n, int *ilo, int *ihi, float *a, int *ia, int *ja, int *desca, float *tau, float *work, int *lwork, int *info);

			void pdgehrd_(int *n, int *ilo, int *ihi, double *a, int *ia, int *ja, int *desca, double *tau, double *work, int *lwork, int *info);

			void pcgehrd_(int *n, int *ilo, int *ihi, complex<float> *a, int *ia, int *ja, int *desca, complex<float> *tau, complex<float> *work, int *lwork, int *info);

			void pzgehrd_(int *n, int *ilo, int *ihi, complex<double> *a, int *ia, int *ja, int *desca, complex<double> *tau, complex<double> *work, int *lwork, int *info);

			////////////////////// Orthogonal Multiplication
			void psormhr_(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, float *a, int *ia, int *ja, int *desca, float *tau, float *c, int *ic, int *jc, int *descc, float *work, int *lwork, int *info);

			void pdormhr_(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, double *a, int *ia, int *ja, int *desca, double *tau, double *c, int *ic, int *jc, int *descc, double *work, int *lwork, int *info);


			///////////////////////
			//////////////////////// Schur
			void pshseqr_(char *job, char* compz, int* n, int* ilo, int* ihi, float* h, int* desch, float* wr, float* wi, float* z, int* descz, float* work, int* lwork, int* iwork, int* liwork, int* info);

			void pdhseqr_(char *job, char* compz, int* n, int* ilo, int* ihi, double* h, int* desch, double* wr, double* wi, double* z, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* info);

			//////////////////////// peigen Schur
			void peigen_pdhseqr_(char *job, char* compz, int* n, int* ilo, int* ihi, double* h, int* desch, double* wr, double* wi, double* z, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* info);

			//auto PEIGEN_PDHSEQR_ = peigen_pdhseqr_;

			//////////////////////// peigen Schur other one
			void pdlahqr_(bool *wantt, bool *wantz, int *n, int *ilo, int *ihi, double *a, int *desca, double *wr, double *wi, int *iloz, int *ihiz, double *z, int *descz, double *work, int *lwork, int *iwork, int *ilwork, int *info);
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


		/////////////////////////				   P*GESVD  				//////////////////////////////
		inline void pxgesvd(char jobu, char jobvt, int m, int n, float* a, int ia, int ja, int* desca, float* s, float* u, int iu, int ju, int* descu, float* vt, int ivt, int jvt, int* descvt, float* work, int lwork, int* info)
		{
			psgesvd_(&jobu, &jobvt, &m, &n, a, &ia, &ja, desca, s, u, &iu, &ju, descu, vt, &ivt, &jvt, descvt, work, &lwork, info);
			return;
		}

		inline void pxgesvd(char jobu, char jobvt, int m, int n, double* a, int ia, int ja, int* desca, double* s, double* u, int iu, int ju, int* descu, double* vt, int ivt, int jvt, int* descvt, double* work, int lwork, int* info)
		{
			pdgesvd_(&jobu, &jobvt, &m, &n, a, &ia, &ja, desca, s, u, &iu, &ju, descu, vt, &ivt, &jvt, descvt, work, &lwork, info);
			return;
		}

		inline void pxgesvd(char jobu, char jobvt, int m, int n, std::complex<float>* a, int ia, int ja, int* desca, std::complex<float>* s, std::complex<float>* u, int iu, int ju, int* descu, std::complex<float>* vt, int ivt, int jvt, int* descvt, std::complex<float>* work, int lwork, int* info)
		{
			pcgesvd_(&jobu, &jobvt, &m, &n, a, &ia, &ja, desca, s, u, &iu, &ju, descu, vt, &ivt, &jvt, descvt, work, &lwork, info);
			return;
		}

		inline void pxgesvd(char jobu, char jobvt, int m, int n, std::complex<double>* a, int ia, int ja, int* desca, std::complex<double>* s, std::complex<double>* u, int iu, int ju, int* descu, std::complex<double>* vt, int ivt, int jvt, int* descvt, std::complex<double>* work, int lwork, int* info)
		{
			pzgesvd_(&jobu, &jobvt, &m, &n, a, &ia, &ja, desca, s, u, &iu, &ju, descu, vt, &ivt, &jvt, descvt, work, &lwork, info);
			return;
		}


		/////////////////////////				   P*GESV  				//////////////////////////////
		inline void pxgesv(int n, int nrhs, float *A, int ia, int ja, int *descA, int *ipiv, float *B, int ib, int jb, int *descB, int *info)
		{
			psgesv_(&n, &nrhs, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, info);
		}

		inline void pxgesv(int n, int nrhs, double *A, int ia, int ja, int *descA, int *ipiv, double *B, int ib, int jb, int *descB, int *info)
		{
			pdgesv_(&n, &nrhs, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, info);
		}

		inline void pxgesv(int n, int nrhs, std::complex<float> *A, int ia, int ja, int *descA, int *ipiv, std::complex<float> *B, int ib, int jb, int *descB, int *info)
		{
			pcgesv_(&n, &nrhs, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, info);
		}

		inline void pxgesv(int n, int nrhs, std::complex<double> *A, int ia, int ja, int *descA, int *ipiv, std::complex<double> *B, int ib, int jb, int *descB, int *info)
		{
			pzgesv_(&n, &nrhs, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, info);
		}

		/////////////////////////				   P*GESVX  				//////////////////////////////
		inline void pxgesvx(char fact, char trans, int n, int nrhs, std::complex<double> *a, int ia, int ja, int *desca, std::complex<double> *af, int iaf, int jaf, int *descaf, int *ipiv, char *equed, double *r, double *c, std::complex<double> *b, int ib, int jb, int *descb, std::complex<double> *x, int ix, int jx, int *descx, double *rcond, double *ferr, double *berr, std::complex<double> *work, int lwork, double *rwork, int lrwork, int *info)
		{
			pzgesvx_(&fact, &trans, &n, &nrhs, a, &ia, &ja, desca, af, &iaf, &jaf, descaf, ipiv, equed, r, c, b, &ib, &jb, descb, x, &ix, &jx, descx, rcond, ferr, berr, work, &lwork, rwork, &lrwork, info);
		}

		




		/////////////////////////				   P*GEMV  				//////////////////////////////
		inline void pxgemv(char opA, int m, int n, float alpha, float *A, int ia, int ja, int *descA, float *x, int ix, int jx, int *descx, int incx, float beta, float *y, int iy, int jy, int *descy, int incy)
		{
			psgemv_(&opA, &m, &n, &alpha, A, &ia, &ja, descA, x, &ix, &jx, descx, &incx, &beta, y, &iy, &jy, descy, &incy);
		}

		inline void pxgemv(char opA, int m, int n, double alpha, double *A, int ia, int ja, int *descA, double *x, int ix, int jx, int *descx, int incx, double beta, double *y, int iy, int jy, int *descy, int incy)
		{
			pdgemv_(&opA, &m, &n, &alpha, A, &ia, &ja, descA, x, &ix, &jx, descx, &incx, &beta, y, &iy, &jy, descy, &incy);
		}

		inline void pxgemv(char opA, int m, int n, std::complex<float> alpha, std::complex<float> *A, int ia, int ja, int *descA, std::complex<float> *x, int ix, int jx, int *descx, int incx, std::complex<float> beta, std::complex<float> *y, int iy, int jy, int *descy, int incy)
		{
			pcgemv_(&opA, &m, &n, &alpha, A, &ia, &ja, descA, x, &ix, &jx, descx, &incx, &beta, y, &iy, &jy, descy, &incy);
		}

		inline void pxgemv(char opA, int m, int n, std::complex<double> alpha, std::complex<double> *A, int ia, int ja, int *descA, std::complex<double> *x, int ix, int jx, int *descx, int incx, std::complex<double> beta, std::complex<double> *y, int iy, int jy, int *descy, int incy)
		{
			pzgemv_(&opA, &m, &n, &alpha, A, &ia, &ja, descA, x, &ix, &jx, descx, &incx, &beta, y, &iy, &jy, descy, &incy);
		}






		//////////////////// HESSENBERG REDUCTION
		inline void pxgehrd(int n, int ilo, int ihi, float *a, int ia, int ja, int *desca, float *tau, float *work, int lwork, int *info)
		{
			psgehrd_(&n, &ilo, &ihi, a, &ia, &ja, desca, tau, work, &lwork, info);
		}

		inline void pxgehrd(int n, int ilo, int ihi, double *a, int ia, int ja, int *desca, double *tau, double *work, int lwork, int *info)
		{
			pdgehrd_(&n, &ilo, &ihi, a, &ia, &ja, desca, tau, work, &lwork, info);
		}

		inline void pxgehrd(int n, int ilo, int ihi, complex<float> *a, int ia, int ja, int *desca, complex<float> *tau, complex<float> *work, int lwork, int *info)
		{
			pcgehrd_(&n, &ilo, &ihi, a, &ia, &ja, desca, tau, work, &lwork, info);
		}

		inline void pxgehrd(int n, int ilo, int ihi, complex<double> *a, int ia, int ja, int *desca, complex<double> *tau, complex<double> *work, int lwork, int *info)
		{
			pzgehrd_(&n, &ilo, &ihi, a, &ia, &ja, desca, tau, work, &lwork, info);
		}

		////////////////////////// Normal multiplication
		inline void pxormhr(char side, char trans, int m, int n, int ilo, int ihi, float *a, int ia, int ja, int *desca, float *tau, float *c, int ic, int jc, int *descc, float *work, int lwork, int *info)
		{
			psormhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &ia, &ja, desca, tau, c, &ic, &jc, descc, work, &lwork, info);
		}

		inline void pxormhr(char side, char trans, int m, int n, int ilo, int ihi, double *a, int ia, int ja, int *desca, double *tau, double *c, int ic, int jc, int *descc, double *work, int lwork, int *info)
		{
			pdormhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &ia, &ja, desca, tau, c, &ic, &jc, descc, work, &lwork, info);
		}

		/////////////////// Schur
		inline void pxhseqr(char job, char compz, int n, int ilo, int ihi, float* h, int* desch, float* wr, float* wi, float* z, int* descz, float* work, int lwork, int* iwork, int liwork, int* info)
		{
			pshseqr_(&job, &compz, &n, &ilo, &ihi, h, desch, wr, wi, z, descz, work, &lwork, iwork, &liwork, info);
		}

		inline void pxhseqr(char job, char compz, int n, int ilo, int ihi, double* h, int* desch, double* wr, double* wi, double* z, int* descz, double* work, int lwork, int* iwork, int liwork, int* info)
		{
			pdhseqr_(&job, &compz, &n, &ilo, &ihi, h, desch, wr, wi, z, descz, work, &lwork, iwork, &liwork, info);
		}


		inline void peigen_pxhseqr(char job, char compz, int n, int ilo, int ihi, double* h, int* desch, double* wr, double* wi, double* z, int* descz, double* work, int lwork, int* iwork, int liwork, int* info)
		{
			peigen_pdhseqr_(&job, &compz, &n, &ilo, &ihi, h, desch, wr, wi, z, descz, work, &lwork, iwork, &liwork, info);
		}

		inline void pxlahqr(bool wantt, bool wantz, int n, int ilo, int ihi, double *a, int *desca, double *wr, double *wi, int iloz, int ihiz, double *z, int *descz, double *work, int lwork, int *iwork, int ilwork, int *info)
		{
			pdlahqr_(&wantt, &wantz, &n, &ilo, &ihi, a, desca, wr, wi, &iloz, &ihiz, z, descz, work, &lwork, iwork, &ilwork, info);
		}
		

	}	// end namespace PBLAS

	namespace ScaLAPACK
	{

	} // end namespace ScaLAPACK
}	// end namespace peigen

#endif // PEIGEN_PBLAS_H