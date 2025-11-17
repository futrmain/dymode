#pragma once

#ifndef COL2ENSIGHT_H
#define COL2ENSIGHT_H

void set80line(string& s)
{
	s.resize(80, ' ');
	s.back() = '\n';
}
void gold_print_header(int mode, string title, string var, FILE *pFile)
{
	stringstream stext;
	string text;

	//stext << complex_part << " of Mode "
	//	<< setfill('0') << setw(6) << mode
	//	<< " for " << var;
	stext << title << " for " << var;
	text = stext.str();
	set80line(text);

	fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
	stext.clear();//clear any bits set
	stext.str(std::string());
}

void gold_print_values(MatrixXf values, geofilereader geo, FILE *pFile)
{
	stringstream stext;
	string text;

	int offset = 0;

	for (auto it_part = geo.parts.begin(); it_part != geo.parts.end(); ++it_part)
	{
		stext << "part";
		text = stext.str();
		set80line(text);

		fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
		stext.clear();//clear any bits set
		stext.str(std::string());

		int part_number = (*it_part).number;
		fwrite(&part_number, 1, 1 * sizeof(int), pFile);

		for (unsigned int k = 0; k < (*it_part).telements.size(); ++k)
		{
			stext << (*it_part).telements[k];
			text = stext.str();
			set80line(text);

			fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
			stext.clear();//clear any bits set
			stext.str(std::string());

			int nelems = (*it_part).nelements[k];
			if (nelems == -1)
				nelems = values.rows();
			fwrite(values.data() + offset, sizeof(float), nelems, pFile);
			offset += (*it_part).nelements[k];
		}
	}
} 

inline Matrix<float, Dynamic, Dynamic> fillBuffer(SharedMatrix<MatrixXd> M, int j)
{
	Matrix<float, Dynamic, Dynamic> buff = M.local_matrix.col(j).cast<float>();

	return buff;
}

inline Matrix<float, Dynamic, Dynamic> fillBuffer(SharedMatrix<MatrixXcd> M, int j)
{
	Matrix<float, Dynamic, Dynamic> buff(M.local_matrix.rows(), 2);

	buff.col(0) = M.local_matrix.col(j).cwiseAbs().cast<float>();
	buff.col(1) = M.local_matrix.col(j).imag().binaryExpr(M.local_matrix.col(j).real(), std::ptr_fun(atan2<double, double>)).cast<float>();

	return buff;
}


template<typename MatrixType>
int col2ensight(SharedMatrix<MatrixType>& Modes, int col, string name, bool complex_values, geofilereader georead, phdfp::datasetreader dreader, options opt)
{
	if (BLACS::mycol == BLACS::indxg2p(col, Modes.cblock(), BLACS::grid_cols))
	{
		int j_loc = BLACS::indxg2l(col, Modes.cblock(), BLACS::grid_cols);
		MatrixXf exportdata;

		const int nc = 1 + complex_values;

		// Gather the global column on row 0
		Matrix<MPI::Request, Dynamic, Dynamic> Irecv_requests;
		Matrix<Matrix<float, Dynamic, Dynamic>, Dynamic, 1> RecvBuffer;
		if (BLACS::myrow == 0)
		{
			exportdata.resize(Modes.rows(), nc);
			Irecv_requests.resize(BLACS::grid_rows, 1);
			RecvBuffer.resize(BLACS::grid_rows, 1);

			for (int r = 0; r < BLACS::grid_rows; ++r)
			{
				int size = BLACS::peigen_numroc(Modes.rows(), Modes.rblock(), r, 0, BLACS::grid_rows);
				RecvBuffer(r, 0).resize(size, nc);
				int powner = BLACS::Cblacs_pnum(BLACS::ctxt, r, BLACS::mycol);
				Irecv_requests(r, 0) = BLACS::COMM_ACTIVE.Irecv(RecvBuffer(r, 0).data(), size * nc, MPI::FLOAT, powner, 1/*tag*/);
			}
		}

		Matrix<float, Dynamic, Dynamic> SendBuffer = fillBuffer(Modes, j_loc);
		
		int col_root = BLACS::Cblacs_pnum(BLACS::ctxt, 0, BLACS::mycol);
		BLACS::COMM_ACTIVE.Send(SendBuffer.data(), SendBuffer.rows() * SendBuffer.cols(), MPI::FLOAT, col_root, 1 /*tag*/);

		if (BLACS::myrow == 0)
		{
			MPI::Request::Waitall(BLACS::grid_rows, Irecv_requests.data());

			// Combine the buffers
			for (int rb = 0; rb < ceil((double)exportdata.rows() / Modes.rblock()); rb++)
			{
				int roffset = Modes.rblock() * floor(rb / BLACS::grid_rows);
				int _nrows = min(Modes.rblock(), (int)(exportdata.rows() - rb*Modes.rblock()));
				int pr_owner = rb % BLACS::grid_rows;
				
				exportdata.block(rb*Modes.rblock(), 0, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 0, _nrows, 1);
				if (complex_values)
				{
					exportdata.block(rb*Modes.rblock(), 1, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 1, _nrows, 1);
				}
			}
		}

		// Print to disk from row 0
		if (BLACS::myrow == 0)
		{
			cout << "(" << BLACS::myrank << ") " << " writing mode " << name << "...";
			FILE *pFile;

			int offset_gold = 0;
			MatrixXf values;
			for (string var : opt.variables)
			{
				if (!(var == "null"))
				{
					stringstream filenameRE;
					//filenameRE << "mode" << setfill('0') << setw(6) << m << "." << var << ".abs";
					filenameRE << name << "." << var;
					if (complex_values)
					{
						filenameRE << ".abs";
					}

					pFile = fopen((opt.outdir + filenameRE.str()).c_str(), "wb");

					gold_print_header(0, name, var, pFile);

					values = exportdata.block(offset_gold, 0, dreader.Np, 1);
					gold_print_values(values, georead, pFile);

					fclose(pFile);

					if (complex_values)
					{
						stringstream filenameIM;
						filenameIM << name << "." << var << ".ang";

						pFile = fopen((opt.outdir + filenameIM.str()).c_str(), "wb");

						gold_print_header(0, name, var, pFile);

						values = exportdata.block(offset_gold, 1, dreader.Np, 1);
						gold_print_values(values, georead, pFile);

						fclose(pFile);
					}

					offset_gold += dreader.Np;
				}
			}
			cout << "\tDONE" << endl;
		}
	}

	return 0;
}

#endif // COL2ENSIGHT_H