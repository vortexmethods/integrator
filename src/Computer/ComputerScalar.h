// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием абстрактного шаблонного класса ComputerScalar (наследник Computer)
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 01 марта 2022 г.
\version 0.1
*/

#pragma once
#include "Computer.h"

/*!
\brief Абстрактный шаблонный класс -- вычислитель скалярнозначных интегралов
\n Наследован от Computer<dim>

\tparam dim размерность задачи

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.1
\date 01 марта 2022 г.
*/

template <int dim>
class ComputerScalar :
	public Computer<dim>
{	
public:
	using Computer<dim>::scalarResult;
	using Computer<dim>::task;

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	ComputerScalar(const Database<dim>& db_) : Computer<dim>(db_) {};
	
	/// Деструктор
	virtual ~ComputerScalar() {};

	/// Перегрузка функции выполнения всего объема расчетов
	virtual void run() override
	{
		scalarResult.clear();
		scalarResult.resize(task.size());


#pragma omp parallel for schedule(dynamic, 1000)
		for (int tsk = 0; tsk < task.size(); ++tsk)
			scalarResult[tsk] = this->scalarEvaluate(task[tsk].first, task[tsk].second);
	}//run()

	/// \brief Перегрузка функции для выполнения одного векторнозначного вычисления
	///
	/// \warning Кидает исключение при попытке обращения
	virtual numvector<double, dim> vectorEvaluate(int i, int j) override
	{
		throw (123);
	}//vectorEvaluate(...)

};//ComputerScalar

