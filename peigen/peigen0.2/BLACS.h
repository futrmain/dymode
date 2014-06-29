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

			void Cdgsum2d( int icontxt, char *scope, char *top, int m, int n, double *A, int lda, int rdest, int cdest );
		}

		int ctxt, myrank, myrow, mycol, grid_rows, grid_cols, numproc; 
		int rblock, cblock;
		bool ROOT;
		bool active;
		std::string major;

		MPI::Intracomm COMM_ACTIVE;

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
			/*grid_cols = floor(sqrt(numtasks));
			grid_rows = floor(numtasks / grid_cols);*/
			grid_cols = 1;
			grid_rows = numtasks;

			Cblacs_pinfo(&myrank, &numproc);
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

		void printGrid()
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