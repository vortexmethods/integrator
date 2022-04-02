/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.10   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2021/05/17     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Parallel.cpp                                                     |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Parallel
\author Марчевский Илья Константинович
\version 1.10
\date 17 мая 2021 г.
*/

#include "Parallel.h"

//using namespace VMlib;

MPI_Datatype Parallel::MPI_V3D, Parallel::MPI_I3D, Parallel::MPI_V2D, Parallel::MPI_I2D, Parallel::MPI_PAIRII;

// Распределение задач по процессорам
parProp Parallel::SplitMPIone(size_t n, bool bcastAll) const
{
	parProp par;

	par.comm = comm;
	par.np = np;
	par.rank = rank;

	par.totalLen = static_cast<int>(n);
	MPI_Bcast(&par.totalLen, 1, MPI_INT, 0, comm);

	if (rank == 0)
	{
		par.len.clear();
		par.disp.clear();

		par.len.push_back(static_cast<int>(n));
		par.disp.push_back(0);
		
		for (int s = 1; s < np; ++s)
		{
			par.len.push_back(0);
			par.disp.push_back(static_cast<int>(n-1));
		}
	}

	MPI_Scatter(par.len.data(), 1, MPI_INT, &par.myLen, 1, MPI_INT, 0, comm);
	MPI_Scatter(par.disp.data(), 1, MPI_INT, &par.myDisp, 1, MPI_INT, 0, comm);

	if (bcastAll)
	{
		if (rank != 0)
		{
			par.len.resize(np);
			par.disp.resize(np);
		}
		MPI_Bcast(par.len.data(), np, MPI_INT, 0, comm);
		MPI_Bcast(par.disp.data(), np, MPI_INT, 0, comm);
	}

	return par;

}//SplitMPIone(...)


// Распределение задач по процессорам
parProp Parallel::SplitMPI(size_t n, bool bcastAll) const
{	
	parProp par;
	
	par.comm = comm;
	par.np = np;
	par.rank = rank;

	par.totalLen = static_cast<int>(n);
	MPI_Bcast(&par.totalLen, 1, MPI_INT, 0, comm);
	
	if (rank == 0)
	{
		par.len.clear();
		par.disp.clear();

		int nPerP = static_cast<int>(n / np);

		for (int s = 0; s < np - 1; ++s)
		{
			par.len.push_back(nPerP);
			par.disp.push_back(s*nPerP);
		}
		
		par.len.push_back(static_cast<int>(n) - nPerP * (np - 1));
		par.disp.push_back(nPerP * (np - 1));
	}

	MPI_Scatter(par.len.data(), 1, MPI_INT, &par.myLen, 1, MPI_INT, 0, comm);
	MPI_Scatter(par.disp.data(), 1, MPI_INT, &par.myDisp, 1, MPI_INT, 0, comm);
	
	if (bcastAll)
	{
	    if (rank != 0)
	    { 
		par.len.resize(np);
		par.disp.resize(np);
	    }
	    MPI_Bcast(par.len.data(),  np, MPI_INT, 0, comm);
	    MPI_Bcast(par.disp.data(), np, MPI_INT, 0, comm);
	}
	
	return par;
	
}//SplitMPI(...)

void Parallel::CreateMpiType()
{
	{
		int          len[1] = { 3 };
		MPI_Aint     pos[1] = { 0 };
		MPI_Datatype typd[1] = { MPI_DOUBLE };
		MPI_Datatype typi[1] = { MPI_INT };

		MPI_Type_create_struct(1, len, pos, typd, &MPI_V3D);
		MPI_Type_create_struct(1, len, pos, typi, &MPI_I3D);
		MPI_Type_commit(&MPI_V3D);
		MPI_Type_commit(&MPI_I3D);
	}

	{
		int          len[1] = { 2 };
		MPI_Aint     pos[1] = { 0 };
		MPI_Datatype typd[1] = { MPI_DOUBLE };
		MPI_Datatype typi[1] = { MPI_INT };

		MPI_Type_create_struct(1, len, pos, typd, &MPI_V2D);
		MPI_Type_create_struct(1, len, pos, typi, &MPI_I2D);
		MPI_Type_commit(&MPI_V2D);
		MPI_Type_commit(&MPI_I2D);
	}

	{
		MPI_Datatype MPI_PAIRIIshort;
		std::pair<int, int> testPair({ 1, 1 });
		int          len[2] = { 1, 1 };
		MPI_Aint     pos[2] = { 0, reinterpret_cast<long long>(&testPair.second) - reinterpret_cast<long long>(&testPair.first) };
		MPI_Datatype typ[2] = { MPI_INT, MPI_INT };
		MPI_Type_create_struct(2, len, pos, typ, &MPI_PAIRIIshort);
		MPI_Type_commit(&MPI_PAIRIIshort);

		//На всякий случай вычисляем extent и удлинняем тип
		MPI_Type_create_resized(MPI_PAIRIIshort, 0, sizeof(std::pair<int, int>), &MPI_PAIRII);
		MPI_Type_commit(&MPI_PAIRII);
	}
}
