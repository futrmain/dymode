#pragma once

#ifndef ROW2TEXT_H
#define ROW2TEXT_H


template<typename MatrixType>
MatrixXd row2single(SharedMatrix<MatrixType>& TimeCoeff, int row)
{
	MatrixXd exportdata(0, TimeCoeff.cols());

	if (BLACS::myrow == BLACS::indxg2p(row, TimeCoeff.rblock(), BLACS::grid_rows))
	{
		int i_loc = BLACS::indxg2l(row, TimeCoeff.rblock(), BLACS::grid_rows);
		
		// Gather the global row on process column 0
		Matrix<MPI::Request, Dynamic, Dynamic> Irecv_requests;
		Matrix<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> RecvBuffer;
		if (BLACS::mycol == 0)
		{
			exportdata.resize(TimeCoeff.cols(), 1);
			Irecv_requests.resize(BLACS::grid_cols, 1);
			RecvBuffer.resize(BLACS::grid_cols, 1);

			for (int c = 0; c < BLACS::grid_cols; ++c)
			{
				int size = BLACS::peigen_numroc(TimeCoeff.cols(), TimeCoeff.cblock(), c, 0, BLACS::grid_cols);
				RecvBuffer(c, 0).resize(size, 1);
				int powner = BLACS::Cblacs_pnum(BLACS::ctxt, BLACS::myrow, c);
				Irecv_requests(c, 0) = BLACS::COMM_ACTIVE.Irecv(RecvBuffer(c, 0).data(), size, MPI::DOUBLE, powner, 1/*tag*/);
			}
		}

		Matrix<double, Dynamic, Dynamic> SendBuffer = TimeCoeff.local_matrix.row(i_loc);
		
		int row_root = BLACS::Cblacs_pnum(BLACS::ctxt, BLACS::myrow, 0);
		BLACS::COMM_ACTIVE.Send(SendBuffer.data(), SendBuffer.cols(), MPI::DOUBLE, row_root, 1 /*tag*/);

		if (BLACS::mycol == 0)
		{
			MPI::Request::Waitall(BLACS::grid_cols, Irecv_requests.data());

			// Combine the buffers
			for (int cb = 0; cb < ceil((double)exportdata.rows() / TimeCoeff.cblock()); cb++)
			{
				int coffset = TimeCoeff.cblock() * floor(cb / BLACS::grid_cols);
				int _ncols = min(TimeCoeff.cblock(), (int)(exportdata.rows() - cb*TimeCoeff.cblock()));
				int pr_owner = cb % BLACS::grid_cols;
				
				exportdata.block(cb*TimeCoeff.cblock(), 0, _ncols, 1) = RecvBuffer(pr_owner, 0).block(coffset, 0, _ncols, 1);
			}
		}
		
	}

	// Send it to rank 0
	if (BLACS::myrank == 0 && BLACS::mycol == 0 && BLACS::myrow == BLACS::indxg2p(row, TimeCoeff.rblock(), BLACS::grid_rows))
	{
		// noting to do
	}
	else
	{
		if (BLACS::myrank == 0)
		{
			exportdata.resize(TimeCoeff.cols(), 1);

			int row_root = BLACS::Cblacs_pnum(BLACS::ctxt, BLACS::indxg2p(row, TimeCoeff.rblock(), BLACS::grid_rows), 0);
			BLACS::COMM_ACTIVE.Recv(exportdata.data(), exportdata.rows(), MPI::DOUBLE, row_root, 2/*tag*/);

		}
		else if (BLACS::mycol == 0 && BLACS::myrow == BLACS::indxg2p(row, TimeCoeff.rblock(), BLACS::grid_rows))
		{
			BLACS::COMM_ACTIVE.Send(exportdata.data(), exportdata.rows(), MPI::DOUBLE, 0, 2 /*tag*/);
			exportdata.resize(0, TimeCoeff.cols());
		}
	}

	return exportdata;
}

#endif // ROW2TEXT_H