// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл с описанием абстрактного шаблонного класса Computer
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#pragma once

#include "Database.h"
#include "Parallel.h"
#include "GaussPoints.h"

/*!
\brief Абстрактный шаблонный класс -- вычислитель интегралов

\tparam T тип результата
\tparam dim размерность задачи

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/


#define OMP_SCHEDULE_BLOCK_SIZE 10000

#define M_PI 3.14159265358979323846

template <typename T, int dim>
class Computer
{
public:
	/// Константная ссылка на базу данных геометрических параметров
	const Database<dim>& db;

	/// Константная ссылка на класс, управляющий распараллеливанием по MPI
	const Parallel& par;
	
	/// Константный указатель на интегратор по гауссовым точкам (может быть пустым)
	const Gausspoints<dim>* const gp;

	/// Вектор пар, определяющий необходимые вычисления
	std::vector<std::pair<int, int>> task;

	/// Вектор для хранения результатов вычислений
	std::vector<T> result;
		
	/// \brief Конструктор
	/// 	
	/// \param[in] db_ константная ссылка на базу данных геометрических параметров
	/// \param[in] par_ константная ссылка на класс, управляющий распараллеливанием по MPI
	/// \param[in] gp_ константный указатель на класс, обеспечивающий интегрирование по гауссовым точкам
	Computer(const Database<dim>& db_, const Parallel& par_, const Gausspoints<dim>* const gp_ = nullptr) : db(db_), par(par_), gp(gp_) {};
	
	/// Деструктор
	virtual ~Computer() {};

	/// \brief Виртуальная функция выполнения всего объема расчетов
	///
	/// \param[in] split признак разделения задач на 3 подсписка: для дальних пар ячеек, для ячеек с общим ребром и для ячеек с общей вершиной
	void run(bool split);

	/// \brief Виртуальная функция для выполнения одного вычисления
	///
	/// \param[in] i индекс контрольной точки или контрольной панели в базе данных
	/// \param[in] j индекс влияющей панели в базе данных
	/// \return результат интегрирования
	virtual inline T evaluate(int i, int j) = 0;	
}; //class Computer


template<typename T, int dim>
void Computer<T,dim>::run(bool split)
{
	result.clear();
	result.resize(task.size());

	if (!split)
	{
		auto parall = par.SplitMPI(task.size(), true);
		//std::cout << "id = " << this->par.myidWork << ", len = " << parprop.myLen << ", disp = " << parprop.myDisp << std::endl;

		std::vector<std::pair<int, int>> locTask;
		parall.ScattervVector(task, locTask);

		std::vector<T> locResult;
		locResult.resize(parall.myLen);

		MPI_Barrier(MPI_COMM_WORLD);
		double timempi = -omp_get_wtime();

		std::cout << "thread_max = " << omp_get_max_threads() << std::endl;
#pragma omp parallel for schedule(dynamic, OMP_SCHEDULE_BLOCK_SIZE)
		for (int tsk = 0; tsk < locTask.size(); ++tsk)
			locResult[tsk] = this->evaluate(locTask[tsk].first, locTask[tsk].second);
				
		parall.GathervVector(locResult, result, false);
		MPI_Barrier(MPI_COMM_WORLD);
		timempi += omp_get_wtime();
		if (parall.rank == 0)
		    std::cout << "time = " << timempi << std::endl;

	}//if !split

	if (split)
	{
		std::vector<numvector<int, 3>> taskFar, taskSosed, taskContact;
		taskFar.reserve(task.size());
		taskSosed.reserve(task.size());
		taskContact.reserve(task.size());

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

		double tfar = -omp_get_wtime(); 
#pragma omp parallel for schedule(dynamic, 50000) //OMP_SCHEDULE_BLOCK_SIZE)
		for (int tsk = 0; tsk < taskFar.size(); ++tsk)
			result[taskFar[tsk][2]] = this->evaluate(taskFar[tsk][0], taskFar[tsk][1]);
		tfar += omp_get_wtime();

		double tsosed = -omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 100)//OMP_SCHEDULE_BLOCK_SIZE)
		for (int tsk = 0; tsk < taskSosed.size(); ++tsk)
			result[taskSosed[tsk][2]] = this->evaluate(taskSosed[tsk][0], taskSosed[tsk][1]);
		tsosed += omp_get_wtime();

		double tcontact = -omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 300)//OMP_SCHEDULE_BLOCK_SIZE)
		for (int tsk = 0; tsk < taskContact.size(); ++tsk)
			result[taskContact[tsk][2]] = this->evaluate(taskContact[tsk][0], taskContact[tsk][1]);
		tcontact += omp_get_wtime();

		std::cout << "time_far     = " << tfar / taskFar.size() << " " << taskFar.size() << std::endl ;
		std::cout << "time_sosed   = " << tsosed / taskSosed.size() << " " << taskSosed.size() << std::endl;
		std::cout << "time_contact = " << tcontact / taskContact.size() << " " << taskContact.size() << std::endl;
		std::cout << "Total_time   = " << tfar + tsosed + tcontact << std::endl;


	}//if split
}