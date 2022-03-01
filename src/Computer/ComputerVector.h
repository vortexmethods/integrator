// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием абстрактного шаблонного класса ComputerVector (наследник Computer)
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 01 марта 2022 г.
\version 0.1
*/

#pragma once
#include "Computer.h"

/*!
\brief Абстрактный шаблонный класс -- вычислитель векторнозначных интегралов
\n Наследован от Computer<dim>

\tparam dim размерность задачи

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.1
\date 01 марта 2022 г.
*/

template <int dim>
class ComputerVector :
	public Computer<dim>
{
public:
	using Computer<dim>::vectorResult;
	using Computer<dim>::task;
	
	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	ComputerVector(const Database<dim>& db_) : Computer(db_) {};
	
	/// Деструктор
	virtual ~ComputerVector() {};

	/// Перегрузка функции выполнения всего объема расчетов
	virtual void run() override
	{
		vectorResult.clear();
		vectorResult.resize(task.size());
		for (int tsk = 0; tsk < task.size(); ++tsk)
			vectorResult[tsk] = vectorEvaluate(task[tsk].first, task[tsk].second);
	}//run()

	/// \brief Перегрузка функции для выполнения одного скалярнозначного вычисления
	///
	/// \warning Кидает исключение при попытке обращения
	virtual double scalarEvaluate(int i, int j) override
	{
		throw (123);
	}//scalarEvaluate(...)

};//ComputerVector

