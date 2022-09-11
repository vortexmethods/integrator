// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса ModelAnalyze
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#pragma once

#include <algorithm>
#include <utility>

#include "Database.h"

/*!
\brief Класс, определяющий характеристики трехмерной модели

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/

#define M_PI 3.14159265358979323846

class ModelAnalyze
{
public:
	/// Константная ссылка на базу данных геометрических параметров
	const Database<3>& db;

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	ModelAnalyze(const Database<3>& db_) : db(db_) {};

	/// Деструктор
	virtual ~ModelAnalyze(){};

	/// \brief Количество ячеек в модели
	///
	size_t NumberOfCells() const { return db.topo.size();};

	/// \brief Наибольшее отношение площадей ячеек
	///
	/// \param[in] neib если true --- поиск наибольшего отношения площадей соседних ячеек
	/// \return наибольшее отношение площадей ячеек
	double MaxAreaNeib(bool neib = false)
	{
		std::vector<double> area_ratio;
		size_t msize = db.measure.size();
		area_ratio.reserve(msize * msize);


		for (size_t i = 0; i < msize - 1; ++i)
			for (size_t j = 0; j < msize - 1; ++j)
			{
				if (i != j)
				{
					if (neib)
					{
						if (db.ifSosed({ i,j }) || db.ifContact({ i,j }))
						{
							area_ratio.push_back(db.measure[i] / db.measure[j]);
						}

						continue;
					}

					area_ratio.push_back(db.measure[i] / db.measure[j]);
				}
			}

		return(*(std::max_element(area_ratio.begin(), area_ratio.end())));
	}

	/// \brief Наибольшее отношение длин сторон ячейки
	///
	/// \return наибольшее отношение длин сторон ячейки
	double MaxSide()
	{
		std::vector<double> area_ratio;
		area_ratio.reserve(6 * db.topo.size());

		numvector<int, 3> trg;
		double st1, st2, st3, r1, r2, r3;

		for (size_t i = 0; i < db.topo.size(); ++i)
		{
			trg = db.topo[i];
			
			st1 = (db.node[trg[0]] - db.node[trg[1]]).length();
			st2 = (db.node[trg[0]] - db.node[trg[2]]).length();
			st3 = (db.node[trg[1]] - db.node[trg[2]]).length();

			r1 = st1 / st2;
			r2 = st1 / st3;
			r3 = st2 / st3;

			area_ratio.push_back(r1);
			area_ratio.push_back(1.0 / r1);
			area_ratio.push_back(r2);
			area_ratio.push_back(1.0 / r2);
			area_ratio.push_back(r3);
			area_ratio.push_back(1.0 / r3);
		}

		return(*(std::max_element(area_ratio.begin(), area_ratio.end())));
	}

	/// \brief Наименьший/наибольший угол в ячейке
	///
	/// \return пара наименьший и наибольший угол в ячейке
	std::pair< double, double > MinMaxAngle()
	{
		std::vector<double> area_ratio;
		/*area_ratio.resize(db.topo.size());*/
		
		numvector<int, 3> trg;
		double st1, st2, st3;

		for (size_t i = 0; i < db.topo.size(); ++i)
		{
			trg = db.topo[i];

			st1 = angle(db.node[trg[1]] - db.node[trg[0]], db.node[trg[2]] - db.node[trg[0]]);
			st2 = angle(db.node[trg[2]] - db.node[trg[1]], db.node[trg[0]] - db.node[trg[1]]);
			st3 = angle(db.node[trg[0]] - db.node[trg[2]], db.node[trg[1]] - db.node[trg[2]]);
						
			area_ratio.push_back(st1);
			//area_ratio.push_back(1.0 / r1);
			area_ratio.push_back(st2);
			//area_ratio.push_back(1.0 / r2);
			area_ratio.push_back(st3);
			//area_ratio.push_back(1.0 / r3);
		}

		return{(*(std::min_element(area_ratio.begin(), area_ratio.end()))), (*(std::max_element(area_ratio.begin(), area_ratio.end())))};
	}


	

};//class ModelAnalyze

