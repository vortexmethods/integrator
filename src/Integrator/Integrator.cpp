// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Основной файл c функцией main
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 01 марта 2022 г.
\version 0.1
*/

#include <iostream>
#include <fstream>
#include "omp.h"

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

\date 01 марта 2022 г.
\version 0.1
*/


int main(int argc, char** argv)
{
    std::cout << "Hello World!\n"; 

	Database<3> db3;
	db3.readNodeTopoFromFile("G1.dat");
	db3.fillNbh();
	db3.calcNrm();
	db3.calcCnt();

	db3.point = db3.cnt;

	
	CompI3DM cmp(db3);
	for (int i = 0; i < 100*db3.topo.size(); ++i)
		for (int j = 0; j < 100*db3.topo.size(); ++j)
			cmp.task.push_back({ i % (int)db3.topo.size(), j % (int)db3.topo.size() });
	
	double t1 = omp_get_wtime();
	cmp.run();
	double t2 = omp_get_wtime();
	std::cout << t2 - t1 << " sec." << std::endl;

	/*
	std::ofstream of("result.txt");
	for (const auto& res : cmp.scalarResult)
		of << res << '\n';
	of.close();
	*/

	
}
