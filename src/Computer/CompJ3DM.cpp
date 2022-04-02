// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompJ3DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#include "CompJ3DM.h"

#define _USE_MATH_DEFINES
#include <math.h>

CompJ3DM::CompJ3DM(const Database<3>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompJ3DM::~CompJ3DM() {};

//J_3D(M_i, K_j) - интеграл от градиента ньютоновского потенциала 

v3D CompJ3DM::evaluate(int i, int j)
{
	//Точка наблюдения
	const v3D& pt = db.point[i];

	//Векторы из вершин панели в точку наблюдения
	const v3D& va = pt - db.node[db.topo[j][0]];
	const v3D& vb = pt - db.node[db.topo[j][1]];
	const v3D& vc = pt - db.node[db.topo[j][2]];

	v3D ova(va), ovb(vb), ovc(vc);

	double lva = ova.normalize();
	double lvb = ovb.normalize();
	double lvc = ovc.normalize();

	//Орты векторов, направленных вдоль сторон треугольной панели K_j, лежащих против соответствующих вершин
	v3D taua = db.node[db.topo[j][2]] - db.node[db.topo[j][1]];
	v3D taub = db.node[db.topo[j][0]] - db.node[db.topo[j][2]];
	v3D tauc = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];

	taua.normalize();
	taub.normalize();
	tauc.normalize();

	v3D Psi = log((lva * (1.0 + (ova & tauc))) / (lvb * (1.0 + (ovb & tauc)))) * tauc +
		log((lvb * (1.0 + (ovb & taua))) / (lvc * (1.0 + (ovc & taua)))) * taua +
		log((lvc * (1.0 + (ovc & taub))) / (lva * (1.0 + (ova & taub)))) * taub;

	//Величина, равная по модулю телесному углу, под которым видна панель из точки наблюдения
	double Theta = 2.0 * atan2((ova ^ ovb) & ovc, 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova));

	return (0.25 / M_PI) * ((Theta * db.nrm[j])  + (Psi ^ db.nrm[j])) ;
}

// & - скалярное произведение
// ^ - вектороне произведение 

