#pragma once
#ifndef PEIGEN_BLACS_H
#define PEIGEN_BLACS_H

#ifdef _WIN32 /* Win32 or Win64 environment */
#define numroc_ NUMROC
#define descinit_ DESCINIT
#endif 

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
			void chk1mat_(const int*, const int*, const int*, const int*, const int*, const int*, int*, const int*, int*);
			void pchk2mat_(const int*, const int*, const int*, const int*, const int*, const int*, int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*, const int*, const int*, int*, const int*,  int*);


			inline int indxg2l(int gidx, int sblock, int nprocs)
			{
				return sblock*floor(gidx / (sblock*nprocs)) + gidx % sblock;
			}
			
			
			inline int indxl2g(int lidx, int sblock, int nprocs, int iproc)
			{
				return nprocs*sblock*(lidx/sblock) + lidx%sblock +  ((nprocs+iproc-0) % nprocs)*sblock ;
			}

			inline int indxg2p(int gidx, int sblock, int nprocs)
			{
				// FIXME: should be "isrc + (gidx/sblock) % nprocs"
				// with isrc The coordinate of the process that possesses the first row / column of the distributed matrix.
				return (gidx/sblock) % nprocs;
			}

			void pdpotrf_(char*, int*, double*,
				int*, int*, int*, int*);

			void descinit_( int *, int *, int *, int *, int *, int *, int *,
				int *, int *, int *);

			void Cdgsum2d( int icontxt, char *scope, char *top, int m, int n, double *A, int lda, int rdest, int cdest );
		}

		int ctxt, myrank, myrow, mycol, grid_rows, grid_cols, numproc; 
		int myrowOrig, mycolOrig, grid_rowsOrig, grid_colsOrig; 
		int rblock, cblock;
		bool ROOT;
		bool active;
		std::string major;

		MPI::Intracomm COMM_ACTIVE;

		// Constants that can be passed to stupid FORTRAN by reference
		extern int izero = 0;
		int ione = 1;
		double dzero = 0.;
		double done = 1;

		int *iZERO ;
		int *iONE ; 
		double *dONE ;
		double *dZERO ;


		inline void init(int numtasks)
		{
			// default blocking factors, so that all matrices within this context can fall back to a common default size
			rblock = 3;
			cblock = 2;

			major = "R";	// the process grid will be row-major

			/* Begin Cblas context for the process grid */
			/* Square pattern */
			
			grid_cols = floor(sqrt(numtasks));
			grid_rows = floor(numtasks / grid_cols);
			
			
			/* Column pattern for the process grid */
			
			//grid_cols = 4;
			//grid_rows = floor(numtasks/4);

			iZERO = new int;
			*iZERO = 0;

			iONE = new int;
			*iONE = 1;

			dZERO = new double;
			*dZERO = 0;

			dONE = new double;
			*dONE = 1;
			
			
			/* Working for statistics */
			//grid_cols = 1;
			//grid_rows = numtasks;
			
			
			
			///cout << "we're here 1" << endl;

			Cblacs_pinfo(&myrank, &numproc);

			//cout << "we're there" << endl;


			ROOT = (myrank == 0);
			Cblacs_get(0, 0, &ctxt);	// get the system context

			

			Cblacs_gridinit(&ctxt, major.c_str(), grid_rows, grid_cols);
				
			

			if (ROOT)
				std::cout << "We are using " << grid_rows*grid_cols << " processes out of " << numtasks << 
				", in a " << grid_rows << "x" << grid_cols << " fashion. Context ID is " << ctxt << std::endl;	

			active = (myrank < grid_rows*grid_cols);
			COMM_ACTIVE = MPI::COMM_WORLD.Split((int)active, BLACS::myrank);
			if (active)
			{
				Cblacs_pcoord(ctxt, myrank, &myrow, &mycol);
			}
			else
			{
				myrow = -1;
				mycol = -1;
			}
			//std::cout << "(" << myrank << ") " << myrow << " x " << mycol << std::endl << std::flush;	
		}
		
		inline void transpGrid()
		{
			int tmp = grid_cols;
			grid_cols = grid_rows;
			grid_rows = tmp;
			
			tmp = myrow;
			myrow = mycol;
			mycol = tmp;
		}
		
		inline void rowify()
		{
			grid_colsOrig = grid_cols;
			grid_rowsOrig = grid_rows;
			
			myrowOrig = myrow;
			mycolOrig = mycol;
			
			grid_rows = 1;
			grid_cols = grid_colsOrig * grid_rowsOrig;
			
			myrow = myrank;
			mycol = 0;
		}
		
		inline void unrowify()
		{
			grid_rows = grid_rowsOrig;
			grid_cols = grid_colsOrig;
			
			myrow = myrowOrig;
			mycol = mycolOrig;
		}

		inline void printGrid()
		{
			using namespace std;

			if (BLACS::ROOT)
					cout << "Here is the process grid: " << endl << flush;
			MPI::COMM_WORLD.Barrier();

			for (int r = 0; r < grid_rows; r++)
			{
				for (int c = 0; c < grid_cols; c++)
				{
					if ( (myrow==r) && (mycol==c) )
					{
						cout << "("<< myrank <<", "<<myrow<<"|"<<mycol<<")\t" << flush;
					}
					MPI::COMM_WORLD.Barrier();
				}
				if (BLACS::ROOT)
					cout << endl << flush;
				MPI::COMM_WORLD.Barrier();
			}

			cout<<flush;
			MPI::COMM_WORLD.Barrier();

			if (BLACS::ROOT)
				cout << endl << "Processes not on the grid: " << endl << flush;

			MPI::COMM_WORLD.Barrier();

			if (mycol == -1)
				cout << "("<< myrank <<", "<<myrow<<"|"<<mycol <<")\t" << flush;

			MPI::COMM_WORLD.Barrier();

			if (BLACS::ROOT)
					cout << endl << flush;

			MPI::COMM_WORLD.Barrier();
		}
	}	// end namespace BLACS
}	// end namespace peigen

#endif // PEIGEN_BLACS_H