// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией явных специализаций класса Gausspoints
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#include "Gausspoints.h"


template<>
void Gausspoints<2>::refine(std::vector<simplex>& lst) const
{
	std::vector<simplex> newLst;
	newLst.reserve(2 * lst.size());
	for (const auto& nds : lst)
	{
		v2D add = 0.5 * (nds.first[0] + nds.first[1]);
		newLst.push_back({ {nds.first[0], add}, 0.5 * nds.second });
		newLst.push_back({ {add, nds.first[1]}, 0.5 * nds.second });
	}
	lst.swap(newLst);
}

template<>
void Gausspoints<3>::refine(std::vector<simplex>& lst) const
{
	std::vector<simplex> newLst;
	newLst.reserve(4 * lst.size());
	for (const auto& nds : lst)
	{
		numvector<v3D, 3> add = { 0.5 * (nds.first[1] + nds.first[2]), 0.5 * (nds.first[2] + nds.first[0]), 0.5 * (nds.first[0] + nds.first[1]) };
		newLst.push_back({ {add[2], nds.first[1], add[0]}, 0.25 * nds.second });
		newLst.push_back({ {add[0], nds.first[2], add[1]}, 0.25 * nds.second });
		newLst.push_back({ {add[1], nds.first[0], add[2]}, 0.25 * nds.second });
		newLst.push_back({ {add[0], add[1],       add[2]}, 0.25 * nds.second });
	}
	lst.swap(newLst);
}

