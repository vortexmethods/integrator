// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией явных специализаций класса Database
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 марта 2022 г.
\version 0.2
*/

#include "Database.h"

#include <fstream>
#include <iostream>

template <>
void Database<3>::readNodeTopoFromFile(const std::string& fileName)
{
	int NP_, NTr_;

	std::ifstream file;
	file.open(fileName, std::ios::in);
	if (file.is_open())
	{
		file >> NP_ >> NTr_;
		double x, y, z;
		int A, B, C, Number, Trigger;

		for (int i = 0; i < NP_; ++i)
		{
			file >> x >> x >> y >> z;
			v3D P{ x, y, z };
			node.push_back(P);
		};

		while (!file.eof())
		{
			file >> Number >> Trigger;
			if (Trigger == 203)
			{
				file >> A >> B >> C;
				topo.push_back({ A - 1, B - 1, C - 1 });
			}
			else
			{
				--NTr_;
				file >> Number >> Number;
			};
		};
	}
	else
		std::cout << "Error! File isn't opened.";

	file.close();
}



template <>
void Database<3>::fillNbh()
{
	std::vector<std::set<int>> db;
	db.resize(node.size());
	std::set<int> memb;

#pragma omp parallel for default(none) shared(db) private(memb)
	for (int nd = 0; nd < node.size(); ++nd)
	{
		memb.clear();
		for (int trg = 0; trg < topo.size(); ++trg)
		{
			if (topo[trg].member(nd) != -1)
				memb.insert(trg);
		}
		db[nd] = memb;
	}

	numvector <std::vector<std::set<int>>, 3> contacti;
	std::set<int> cont;

	for (int trg = 0; trg < topo.size(); ++trg)
	{
		for (int vert = 0; vert < 3; ++vert)
		{
			int nd = topo[trg][vert];
			cont = db[nd];
			cont.erase(trg);
			contacti[vert].push_back(cont);
		}
	}

	std::vector<numvector<int, 3>> sosedi;
	sosedi.resize(topo.size());
	for (int trg = 0; trg < topo.size(); ++trg)
	{
		for (int vert = 0; vert < 3; ++vert)
		{
			int vert2 = vert < 2 ? vert + 1 : 0;

			std::set<int>& cnt = contacti[vert][trg];
			for (auto othertrg = cnt.begin(); othertrg != cnt.end(); ++othertrg)
			{
				if (topo[*othertrg].member(topo[trg][vert2]) != -1)
				{
					sosedi[trg][vert] = *othertrg;
					break;
				}
			}
		}
	}

	for (int trg = 0; trg < topo.size(); ++trg)
		for (int vert = 0; vert < 3; ++vert)
			for (auto& it : contacti[vert][trg])
				nbh[{trg, it}] = contact;


	for (int trg = 0; trg < topo.size(); ++trg)
		for (int vert = 0; vert < 3; ++vert)
			nbh[{trg, sosedi[trg][vert]}] = sosed;

	/*
		std::string name = "G1";
		std::ofstream fsosed(name + "-sosedi.dat");
		fsosed << topo.size() << std::endl;
		for (int trg = 0; trg < topo.size(); ++trg)
		{
			fsosed << trg + 1 << " " << sosedi[trg][0] + 1 << " " << sosedi[trg][1] + 1 << " " << sosedi[trg][2] + 1 << std::endl;
		}
		fsosed.close();


		std::ofstream fcont(name + "-contacti.dat");
		fcont << topo.size() << std::endl;
		for (int trg = 0; trg < topo.size(); ++trg)
		{
			fcont << trg + 1 << std::endl;
			for (int vert = 0; vert < 3; ++vert)
			{
				contacti[vert][trg].erase(sosedi[trg][0]);
				contacti[vert][trg].erase(sosedi[trg][1]);
				contacti[vert][trg].erase(sosedi[trg][2]);

				fcont << contacti[vert][trg].size();
				for (auto it = contacti[vert][trg].begin(); it != contacti[vert][trg].end(); ++it)
					fcont << " " << (*it) + 1;
				fcont << std::endl;
			}
		}
		fcont.close();
	*/
}

template <>
void Database<3>::calcNrm()
{
	nrm.clear();
	nrm.resize(0);

	for (auto& trg : topo)
	{
		v3D v1 = node[trg[1]] - node[trg[0]];
		v3D v2 = node[trg[2]] - node[trg[0]];
		nrm.push_back((v1^v2).unit());
	}
}

template<> 
void Database<3>::calcCnt()
{
	cnt.clear();
	cnt.resize(0);

	cnt.reserve(topo.size());
	for (const auto& trg : topo)
		cnt.push_back((node[trg[0]] + node[trg[1]] + node[trg[2]])*0.3333333333333333);
}