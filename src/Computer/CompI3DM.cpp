// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией примера класса ComputerScalar
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 01 марта 2022 г.
\version 0.1
*/

#include "CompI3DM.h"

#define _USE_MATH_DEFINES
#include <math.h>

CompI3DM::CompI3DM(const Database<3>& db_) : ComputerScalar(db_) {};
CompI3DM::~CompI3DM() {};

//double Theta(const v3D& va, const v3D& vb, const v3D& vc)
//{
//	const v3D& ova = va.unit();
//	const v3D& ovb = vb.unit();
//	const v3D& ovc = vc.unit();
//
//	return 2.0 * atan2((ova ^ ovb) & ovc, 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova));
//
//}

double CompI3DM::scalarEvaluate(int i, int j)
{
	const v3D& pt = db.point[i];

	const v3D& va = pt - db.node[db.topo[j][0]];
	const v3D& vb = pt - db.node[db.topo[j][1]];
	const v3D& vc = pt - db.node[db.topo[j][2]];

	v3D ova(va), ovb(vb), ovc(vc);

	double lva = ova.normalize();
	double lvb = ovb.normalize();
	double lvc = ovc.normalize();
		
	v3D taua = db.node[db.topo[j][2]] - db.node[db.topo[j][1]];
	v3D taub = db.node[db.topo[j][0]] - db.node[db.topo[j][2]];
	v3D tauc = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];

	double La = taua.normalize();
	double Lb = taub.normalize();
	double Lc = tauc.normalize();

	
	double fca = ova & tauc;
	double fcb = -(ovb & tauc);
	double fab = ovb & taua;
	double fac = -(ovc & taua);
	double fbc = ovc & taub;
	double fba = -(ova & taub);

	v3D Phi = (va ^ vb) * (log((lva * (1.0 + fca)) / (lvb * (1.0 - fcb))) / Lc) 
		    + (vb ^ vc) * (log((lvb * (1.0 + fab)) / (lvc * (1.0 - fac))) / La) 
		    + (vc ^ va) * (log((lvc * (1.0 + fbc)) / (lva * (1.0 - fba))) / Lb);

	double Theta = 2.0 * atan2((ova ^ ovb) & ovc, 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova));
	
	return 0.25/M_PI * ((Phi - va * Theta) & db.nrm[j]);
}
