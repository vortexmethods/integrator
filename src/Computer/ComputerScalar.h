// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием абстрактного шаблонного класса ComputerScalar (наследник Computer)
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 марта 2022 г.
\version 0.2
*/

#pragma once
#include "Computer.h"

#include <iostream>
#include "omp.h"

/*!
\brief Абстрактный шаблонный класс -- вычислитель скалярнозначных интегралов
\n Наследован от Computer<dim>

\tparam dim размерность задачи

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.2
\date 11 марта 2022 г.
*/

template <int dim>
class ComputerScalar :
	public Computer<dim>
{	
public:
	using Computer<dim>::scalarResult;
	using Computer<dim>::task;
	using Computer<dim>::par;

	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	ComputerScalar(const Database<dim>& db_, const Parallel& par_) : Computer<dim>(db_, par_) {};
	
	/// Деструктор
	virtual ~ComputerScalar() {};

	/// Перегрузка функции выполнения всего объема расчетов	
	virtual void run(bool split) override
	{				
		scalarResult.clear();
		scalarResult.resize(task.size());

		if (!split)
		{		
			auto parall = par.SplitMPI(task.size(), true);
			//std::cout << "id = " << this->par.myidWork << ", len = " << parprop.myLen << ", disp = " << parprop.myDisp << std::endl;
						
			std::vector<std::pair<int, int>> locTask;
			parall.ScattervVector(task, locTask);
			
			std::vector<double> locScalarResult;
			locScalarResult.resize(parall.myLen);
					   
#pragma omp parallel for schedule(dynamic, OMP_SCHEDULE_BLOCK_SIZE)
			for (int tsk = 0; tsk < locTask.size(); ++tsk)
				locScalarResult[tsk] = this->scalarEvaluate(locTask[tsk].first, locTask[tsk].second);

			parall.GathervVector(locScalarResult, scalarResult);			
		}//if !split

		if (split)
		{
			std::vector<numvector<int, 3>> taskFar, taskSosed, taskContact;
			taskFar.reserve(task.size());
			taskSosed.reserve(task.size());
			taskContact.reserve(task.size());

			double tt1 = omp_get_wtime();

			for (int tsk = 0; tsk < task.size(); ++tsk)
			{
				auto it = db.nbh.find(task[tsk]);
				if (it != db.nbh.end())

					switch (it->second)
					{
					case sosed:
						taskSosed.push_back({ task[tsk].first, task[tsk].second, tsk });
						break;

					case contact:
						taskContact.push_back({ task[tsk].first, task[tsk].second, tsk });
						break;
					}

				else
					taskFar.push_back({ task[tsk].first, task[tsk].second, tsk });
				
			}

			double tt2 = omp_get_wtime();

			std::cout << "dt = " << tt2 - tt1 << std::endl;

#pragma omp parallel for schedule(dynamic, OMP_SCHEDULE_BLOCK_SIZE)
			for (int tsk = 0; tsk < taskFar.size(); ++tsk)
				scalarResult[taskFar[tsk][2]] = this->scalarEvaluate(taskFar[tsk][0], taskFar[tsk][1]);

#pragma omp parallel for schedule(dynamic, OMP_SCHEDULE_BLOCK_SIZE)
			for (int tsk = 0; tsk < taskSosed.size(); ++tsk)
				scalarResult[taskSosed[tsk][2]] = this->scalarEvaluate(taskSosed[tsk][0], taskSosed[tsk][1]);

#pragma omp parallel for schedule(dynamic, OMP_SCHEDULE_BLOCK_SIZE)
			for (int tsk = 0; tsk < taskContact.size(); ++tsk)
				scalarResult[taskContact[tsk][2]] = this->scalarEvaluate(taskContact[tsk][0], taskContact[tsk][1]);

		}//if split

	}//run()

	/// \brief Перегрузка функции для выполнения одного векторнозначного вычисления
	///
	/// \warning Кидает исключение при попытке обращения
	virtual numvector<double, dim> vectorEvaluate(int i, int j) override
	{
		throw (123);
	}//vectorEvaluate(...)

};//ComputerScalar

