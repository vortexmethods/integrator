// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса CompTest (пример, наследник ComputerScalar)
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 01 марта 2022 г.
\version 0.1
*/

#pragma once
#include "ComputerScalar.h"

/*!
\brief Пример класса -- вычислителя скалярнозначных интегралов для пространственного случая
\n Наследован от Computer<3>
\n Для примера вычисялет расстояние между центрами панелей

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.1
\date 01 марта 2022 г.
*/

class CompTest :
	public ComputerScalar<3>
{
public:

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	CompTest(const Database<3>& db_);
	
	/// Деструктор
	~CompTest();

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \param[in] i индекс контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return скалярный результат --- расстояние между центрами панелей
	virtual double scalarEvaluate(int i, int j) override;	
};

