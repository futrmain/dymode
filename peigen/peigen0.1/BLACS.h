#ifndef PEIGEN_BLACS_H
#define PEIGEN_BLACS_H

namespace peigen
{
	namespace BLACS
	{
		extern "C" {
			/* Cblacs declarations */
			void Cblacs_pinfo(int*, int*);
			void Cblacs_get(int, int, int*);
			void Cblacs_gridinit(int*, const char*, int, int);
			void Cblacs_gridinfo(int, int*, int*, int*,int*);
			void Cblacs_pcoord(int, int, int*, int*);
			void Cblacs_gridexit(int);
			void Cblacs_exit(int);
			void Cblacs_barrier(int, const char*);
			void Cdgerv2d(int, int, int, double*, int, int, int);
			void Cdgesd2d(int, int, int, double*, int, int, int);
			int Cblacs_pnum(int , int , int);

			int numroc_(int*, int*, int*, int*, int*);

			void pdpotrf_(char*, int*, double*,
				int*, int*, int*, int*);

			void descinit_( int *, int *, int *, int *, int *, int *, int *,
				int *, int *, int *);

		}

		int ctxt, myrank, myrow, mycol, grid_rows, grid_cols, numproc; 
		int rblock, cblock;
		bool ROOT;
		std::string major;

		// Constants that can be passed to stupid FORTRAN by reference
		int izero = 0;
		int ione = 1;
		double dzero = 0;
		double done = 1;
		int *iZERO = &izero; 
		int *iONE = &ione; 
		double *dONE = &done;
		double *dZERO = &dzero;

		void init(int numtasks)
		{
			// default blocking factors, so that all matrices within this context can fall back to a common default size
			rblock = 3;
			cblock = 2;

			major = "R";	// the process grid will be row-major

			/* Begin Cblas context */
			grid_cols = floor(sqrt(numtasks));
			grid_rows = floor(numtasks / grid_cols);


			Cblacs_pinfo(&myrank, &numproc);
			ROOT = (myrank == 0);
			Cblacs_get(0, 0, &ctxt);	// get the system context
			Cblacs_gridinit(&ctxt, major.c_str(), grid_rows, grid_cols);
			Cblacs_pcoord(ctxt, myrank, &myrow, &mycol);

			if (ROOT)
				std::cout << "We are using " << grid_rows*grid_cols << " processes out of " << numtasks << 
				", in a " << grid_rows << "x" << grid_cols << " fashion. ctxt.id is " << ctxt << std::endl;	
		}
	}	// end namespace BLACS
}	// end namespace peigen

#endif // PEIGEN_BLACS_H