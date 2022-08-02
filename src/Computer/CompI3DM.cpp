// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompI3DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 августа 2022 г.
\version 0.4
*/

#include "CompI3DM.h"

#define _USE_MATH_DEFINES
#include <math.h>

CompI3DM::CompI3DM(const Database<3>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompI3DM::~CompI3DM() {};

////CompI3DM::CompI3DM(const Database<3>& db_) : ComputerScalar(db_) {};
////CompI3DM::~CompI3DM() {};

//I_3D(M_i, K_j) - интеграл от ньютоновского потенциала

double CompI3DM::evaluate(int i, int j)
{
	//Точка наблюдения
	const v3D& pt = db.point[i];

	//Векторы из вершин панели в точку наблюдения
	const v3D& va = pt - db.node[db.topo[j][0]];
	const v3D& vb = pt - db.node[db.topo[j][1]];
	const v3D& vc = pt - db.node[db.topo[j][2]];

	v3D ova(va), ovb(vb), ovc(vc);

	//
	double lva = ova.normalize();
	double lvb = ovb.normalize();
	double lvc = ovc.normalize();
	
	//Орты векторов, направленных вдоль сторон треугольной панели K_j, лежащих против соответствующих вершин
	v3D taua = db.node[db.topo[j][2]] - db.node[db.topo[j][1]];
	v3D taub = db.node[db.topo[j][0]] - db.node[db.topo[j][2]];
	v3D tauc = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];

	//Длины сторон треугольника, лежащие против соответствующих вершин
	double La = taua.normalize();
	double Lb = taub.normalize();
	double Lc = tauc.normalize();

	//phi^p_q - плоский угол при основании боковой грани получающегося тетраэдра 
	//(между вектором v_p и стороной L_q)
	double fca = ova & tauc;
	double fcb = -(ovb & tauc);
	double fab = ovb & taua;
	double fac = -(ovc & taua);
	double fbc = ovc & taub;
	double fba = -(ova & taub);

	v3D Phi = (va ^ vb) * (log((lva * (1.0 + fca)) / (lvb * (1.0 - fcb))) / Lc) 
		    + (vb ^ vc) * (log((lvb * (1.0 + fab)) / (lvc * (1.0 - fac))) / La) 
		    + (vc ^ va) * (log((lvc * (1.0 + fbc)) / (lva * (1.0 - fba))) / Lb);

	//Величина, равная по модулю телесному углу, под которым видна панель из точки наблюдения
	double Theta = 2.0 * atan2((ova ^ ovb) & ovc, 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova));
	
	return 0.25/M_PI * ((Phi - va * Theta) & db.nrm[j]);
}
