#ifndef TEST_H5INDMD_H
#define TEST_H5INDMD_H

#include "H5inDMD.h"

namespace test_H5inDMD
{
	class test_datasetreader
	{
	public:
		phdfp::datasetreader dreader;

		test_datasetreader(int nfiles, string filenameroot) : dreader(phdfp::datasetreader(nfiles, filenameroot))
		{
		}

		void read(string datasetname)
		{
			dreader.read(datasetname);

			Eigen::MatrixXi snapsperproc(1, peigen::BLACS::numproc);
			for (int k = 0; k < snapsperproc.cols(); ++k)
			{
				snapsperproc(0, k) = dreader.snapshots_per_process(k, 0).size();
			}
			//cout << "Process " << peigen::BLACS::myrank << " thinks each process has this many snapshots " << snapsperproc << endl;

			cout << "Process " << peigen::BLACS::myrank << " owns this part of the matrix: " << endl << dreader.input_data << endl;
		}

		void createShared(int rblock, int cblock, int nskip = 1)
		{
			peigen::SharedMatrix<Eigen::MatrixXd> S(dreader.createShared(rblock, cblock, nskip));
			cout << "The shared Matrix on " << peigen::BLACS::myrank << " is " << endl << S.local_matrix << endl;

			//if (peigen::BLACS::myrank==0)
				//cout << "And finally " << endl << S << endl;

		}
	};


	/*
	usage


	test_H5inDMD::test_datasetreader testreader(nfiles, "../../../testmat");
	testreader.read("testgroup_T");
	
	testreader.createShared(5, 5, nskip_step);
	*/
}

#endif