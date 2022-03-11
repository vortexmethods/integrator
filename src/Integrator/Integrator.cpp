// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Основной файл c функцией main
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 марта 2022 г.
\version 0.2
*/

#include <fstream>

#include "CompTest.h"
#include "CompI3DM.h"

/*!
\mainpage Вычисление интегралов от фундаментального решения уравнения Лапласа и его градиента
Данный программный модуль реализует алгоритмы вычислния однократных и повторных интегралов от 
фундаментального решения уравнения Лапласа и его градиента по панелям
(в двумерных задачах -- прямым, в трехмерных задачах -- треугольным).
Подобные интегралы возникают при решении граничных интегральных уравнений при кусочно-постоянном 
представлении решения на панелях

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 марта 2022 г.
\version 0.2
*/


int main(int argc, char** argv)
{
	Parallel par;
	
	// Инициализация подсистемы MPI
	MPI_Init(&argc, &argv);
	// Получить размер коммуникатора MPI_COMM_WORLD
	// (общее число процессов в рамках задачи)
	MPI_Comm_size(MPI_COMM_WORLD, &par.np);
	// Получить номер текущего процесса в рамках 
	// коммуникатора MPI_COMM_WORLD
	MPI_Comm_rank(MPI_COMM_WORLD, &par.rank);	
	
	Parallel::CreateMpiType();

	Database<3> db3;
	int nNode, nTopo;

	if (par.rank == 0)
	{
		std::cout << "Hello World!\n";	
		db3.readNodeTopoFromFile("G1.dat");
		nNode = (int)db3.node.size();
		nTopo = (int)db3.topo.size();
	}

	//Рассылка базы данных
	par.SplitMPI(db3.node.size()).BcastVector(db3.node);
	par.SplitMPI(db3.topo.size()).BcastVector(db3.topo);

	db3.fillNbh();
	db3.calcNrm();
	db3.calcCnt();
	db3.point = db3.cnt;		

	CompI3DM cmp(db3, par);
	if (par.rank == 0)
	{	
		for (int i = 0; i < 1 * db3.topo.size(); ++i)
			for (int j = 0; j < 1 * db3.topo.size(); ++j)
				cmp.task.push_back({ i % (int)db3.topo.size(), j % (int)db3.topo.size() });
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double t1 = omp_get_wtime();
	
	cmp.run(false);
	
	MPI_Barrier(MPI_COMM_WORLD);
	double t2 = omp_get_wtime();
	
	std::cout << t2 - t1 << " sec." << std::endl;

	/*
	if (par.rank == 0)
	{
		std::ofstream of("result_false.txt");
		for (const auto& res : cmp.scalarResult)
			of << res << '\n';
		of.close();
	}
	//*/

	MPI_Finalize();
	
}
