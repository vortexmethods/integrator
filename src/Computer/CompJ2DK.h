// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса CompJ2DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#pragma once
#include "Computer.h"

/*!
\brief Класс -- вычислитель для 2D случая однократного интеграла от градиента функции Грина
\n Наследован от Computer<v2D, 3>

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.3
\date 02 апреля 2022 г.
*/

class CompJ2DK :
	public Computer<v2D, 2>
{
public:

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	CompJ2DK(const Database<2>& db_, const Parallel& par_);

	/// Деструктор
	~CompJ2DK();

	double Theta(const v2D& v1, const v2D& v2);

	v2D omega(const v2D& v, const v2D& ti, const v2D& tj);

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \param[in] i индекс контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return скалярный результат --- расстояние между центрами панелей
	virtual inline v2D evaluate(int i, int j) override;
};

