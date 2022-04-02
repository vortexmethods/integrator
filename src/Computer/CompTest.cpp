// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией примера класса ComputerScalar
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#include "CompTest.h"

CompTest::CompTest(const Database<3>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompTest::~CompTest() {};

double CompTest::evaluate(int i, int j)
{
	v3D pi = (db.node[db.topo[i][0]] + db.node[db.topo[i][1]] + db.node[db.topo[i][2]]) * 0.33333333333333333;
	v3D pj = (db.node[db.topo[j][0]] + db.node[db.topo[j][1]] + db.node[db.topo[j][2]]) * 0.33333333333333333;
	return (pi - pj).length();
}
