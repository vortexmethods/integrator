// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией явных специализаций класса Database
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#include "Database.h"

#include <fstream>
#include <iostream>

template <>
void Database<3>::readNodeTopoFromFile(const std::string& fileName, double scale)
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
			node.push_back(scale * P);
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
void Database<2>::readNodeTopoFromFile(const std::string& fileName, double scale)
{
	int NP_, NTr_;

	std::ifstream file;
	file.open(fileName, std::ios::in);
	if (file.is_open())
	{
		file >> NP_ >> NTr_;
		double x, y, z;
		int A, B, Number, Trigger;

		for (int i = 0; i < NP_; ++i)
		{
			file >> x >> x >> y >> z;
			v3D P{ x, y, z };
			node.push_back(scale * v2D({P[0], P[1]}));
		};

		while (!file.eof())
		{
			file >> Number >> Trigger;
			if (Trigger == 102)
			{
				file >> A >> B;
				topo.push_back({ A - 1, B - 1 });
			}
			/*else
			{
				--NTr_;
				file >> Number >> Number;
			};*/
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
	sosedi.resize(topo.size(), {-1, -1, -1});
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

	//std::cout << "-----------1-------------" << std::endl;
	//for (auto& n : nbh)
	//	std::cout << n.first.first << ", " << n.first.second << ", " << (int)n.second << std::endl;

	for (int trg = 0; trg < topo.size(); ++trg)
		for (int vert = 0; vert < 3; ++vert)
			for (auto& it : contacti[vert][trg])
				nbh[{trg, it}] = contact;

	//std::cout << "-----------2-------------" << std::endl;
	//for (auto& n : nbh)
	//	std::cout << n.first.first << ", " << n.first.second << ", " << (int)n.second << std::endl;

	
	for (int trg = 0; trg < topo.size(); ++trg)
		for (int vert = 0; vert < 3; ++vert)
			if (sosedi[trg][vert] != -1)
				nbh[{trg, sosedi[trg][vert]}] = sosed;

	//std::cout << "-----------3-------------" << std::endl;
	//for (auto& n : nbh)
	//	std::cout << n.first.first << ", " << n.first.second << ", " << (int)n.second << std::endl;


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
void Database<2>::fillNbh()
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

	for (int i = 0; i < topo.size(); ++i)
	{
		std::set<int> ssd;
		for (int nd = 0; nd < node.size(); ++nd)
		{		
			auto it = db[nd].find(i);
			if (it != db[nd].end())
			{
				ssd.insert(db[nd].begin(), db[nd].end());
			}
			ssd.erase(i);
		}

		//std::cout << "i = " << i << ", ssd = ";
		//for (auto q : ssd)
		//	std::cout << q << " ";
		//std::cout << std::endl;

		for (auto q : ssd)
			nbh[{i, q}] = sosed;
	}



	
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

template<>
void Database<3>::calcMeasure()
{
	measure.clear();
	measure.resize(0);

	measure.reserve(topo.size());
	for (const auto& trg : topo)
		measure.push_back( 0.5* ((node[trg[2]] - node[trg[0]]) ^ (node[trg[1]] - node[trg[0]])).length() );
}


template <>
void Database<2>::calcNrm()
{
	nrm.clear();
	nrm.resize(0);

	for (auto& trg : topo)
	{
		v2D v = node[trg[1]] - node[trg[0]];
		nrm.push_back(-v.kcross().unit());
	}
}

template<>
void Database<2>::calcCnt()
{
	cnt.clear();
	cnt.resize(0);

	cnt.reserve(topo.size());
	for (const auto& trg : topo)
		cnt.push_back((node[trg[0]] + node[trg[1]]) * 0.5);
}

template<> 
void Database<2>::calcMeasure()
{
	measure.clear();
	measure.resize(0);

	measure.reserve(topo.size());
	for (const auto& trg : topo)
		measure.push_back((node[trg[1]] - node[trg[0]]).length());
}