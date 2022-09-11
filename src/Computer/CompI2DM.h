// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса CompI2DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#pragma once
#include "Computer.h"

/*!
\brief Класс -- вычислитель для 2D случая однократного интеграла от функции Грина
\n Наследован от Computer<double, 2>

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/

class CompI2DM :
	public Computer<double, 2>
{
public:

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	CompI2DM(const Database<2>& db_, const Parallel& par_);

	/// Деструктор
	~CompI2DM();

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \param[in] i индекс контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return скалярный результат --- расстояние между центрами панелей
	virtual inline double evaluate(int i, int j) override;
};

