// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием класса CompJ3DK
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#pragma once
#include "Computer.h"

/*!
\brief Класс -- вычислитель для 3D случая двукратного интеграла от градиента функции Грина
\n Наследован от Computer<v3D, 3>

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/

class CompJ3DK :
	public Computer<v3D, 3>
{
public:

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	/// \param[in] gp_ константный указатель на класс, обеспечивающий интегрирование по гауссовым точкам
	CompJ3DK(const Database<3>& db_, const Parallel& par_, const Gausspoints<3>* const gp_);

	/// Деструктор
	~CompJ3DK();

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \param[in] i индекс контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return скалярный результат --- расстояние между центрами панелей
	virtual inline v3D evaluate(int i, int j) override;
	
	//////////////////////

	p13D ThetaPsi(const v3D& pt, int j);
	v3D J3D(const v3D& pt, int j);

	//////////////////////
	
	i2D RenumerationSosed(int i, int j);
	i2D RenumerationContact(int i, int j);
		
	p13D SingSosed(const v3D& M, const i2D& jj);
	p13D SingContact(const v3D& M, const i2D& ii, const i2D& jj);
	
	v3D JSingSosed(const v3D& M, const i2D& jj);
	v3D JRegSosed(const v3D& M, const i2D& jj);
	p13D RegContact(const v3D& M, const i2D& ii, const i2D& jj);
	
	v3D IntJSingSosed(const i2D& ii, const i2D& jj);
	v3D IntJRegSosed(int i, const i2D& jj, size_t refineLevel = 0);
	std::pair<v3D, int> IntJRegSosedEpsRel(int i, const i2D& jj);
	
	p13D IntSingContact(const i2D& ii, const i2D& jj);
	p13D IntRegContact(const i2D& ii, const i2D& jj, size_t refineLevel = 0);
	std::pair<p13D, int> IntRegContactEpsRel(const i2D& ii, const i2D& jj);

	//std::vector<double> anglesSosed(int i, int j);	 // И.К.
		
	// И.К. Нижние функции не нужны
	/*
	v3D JSingContact(const v3D& M, const i2D& ii, const i2D& jj);
	v3D JRegContact(const v3D& M, const i2D& ii, const i2D& jj);
	
	v3D IntJRegContact(const i2D& ii, const i2D& jj, size_t refineLevel = 0);	
	v3D IntJSingContact(const i2D& ii, const i2D& jj);
	//*/

	std::vector<int> refines;
};

