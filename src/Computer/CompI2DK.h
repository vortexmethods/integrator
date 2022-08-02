// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса CompI2DK
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 августа 2022 г.
\version 0.4
*/

#pragma once
#include "Computer.h"

class CompI2DK :
	public Computer<double, 2>
{
public:

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	CompI2DK(const Database<2>& db_, const Parallel& par_);

	/// Деструктор
	~CompI2DK();

	double Theta(const v2D& v1, const v2D& v2);

	double psi(const v2D& v, const v2D& ti, const v2D& tj);

	//double psi(const v2D& v1, int i, int j);

	//double psi(const v2D& v1, const v2D& t1, const v2D& t2);

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \param[in] i индекс контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return скалярный результат --- расстояние между центрами панелей
	virtual inline double evaluate(int i, int j) override;

};


