// ������ Integrator 
// (c) �.�. ��������, �.�. ����������, �.�. ����������, 2022

/*!
\file
\brief ���� c ����������� ������ CompJ2DM
\author �������� ���� ����������
\author ���������� ���� ��������������
\author ���������� ����� ���������

\date 02 ������ 2022 �.
\version 0.3
*/

#include "CompJ2DK.h"

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

CompJ2DK::CompJ2DK(const Database<2>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompJ2DK::~CompJ2DK() {};

const double eps_zero = 1e-6;

//J_2D(K_i, K_j) - ��������� �������� �� ��������� ���������������� ���������� 

double CompJ2DK::Theta(const v2D& v1, const v2D& v2)
{
	return atan2((v1 ^ v2), (v1 & v2));
}

v2D CompJ2DK::omega(const v2D& v, const v2D& ti, const v2D& tj)
{
	return (v & ti) * tj + (v & tj) * ti - (ti & tj) * v;
}

v2D CompJ2DK::evaluate(int i, int j)
{
	//����� ����������
	//const v2D& pt = db.point[i];

	const v2D& pt1 = db.node[db.topo[i][0]];
	const v2D& pt2 = db.node[db.topo[i][1]];

	//������� �� ������ ������ � ����� ����������
	const v2D& va = pt1 - db.node[db.topo[j][0]];
	const v2D& vb = pt1 - db.node[db.topo[j][1]];

	const v2D& wa = pt2 - db.node[db.topo[j][0]];
	const v2D& wb = pt2 - db.node[db.topo[j][1]];

	//
	v2D ova(va), ovb(vb), owa(wa), owb(wb);

	//������ �������� v_a � v_b ��������������
	double lva = ova.normalize();
	double lvb = ovb.normalize();
	double lwa = owa.normalize();
	double lwb = owb.normalize();

	//tau_j - ��������� ������������ ������ ������ 
	v2D tauj = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];

	v2D taui = db.node[db.topo[i][1]] - db.node[db.topo[i][0]];

	//std::cout << "i = " << i << " " << "j = " << j << " " << tau[0] << " " << tau[1] << std::endl;

	//std::cout << "----------------------------------------------------------------------------------------" << std::endl;

	//L_j - ����� ������������� ������
	double Lj = tauj.normalize();

	double Li = taui.normalize();

	//std::cout << "i = " << i << " " << "j = " << j << " " << tau[0] << " " << tau[1] << std::endl;

	//h - ������������ ���������� ������
	//const double mh = 1.0;

	v2D br1 = Li * Theta(wb, va) * tauj;
	v2D br2 = Theta(va, wa) * omega(wa, taui, tauj);
	v2D br3 = Theta(wb, vb) * omega(vb, taui, tauj);

	if ((lva / Li < eps_zero) && (lwb / Li < eps_zero))
		return { 0.0, 0.0 };
	
	if (lvb/Li < eps_zero)
		return 0.5 / M_PI * (-(br1+br2).kcross() + (Li * log(lva / lwb) * tauj + log(lwa / lva) * omega(wa, taui, tauj) ));
	
	if (lwa / Li < eps_zero)
		return 0.5 / M_PI * (-(br1+br3).kcross() + (Li * log(lva / lwb) * tauj + log(lvb / lwb) * omega(vb, taui, tauj)));
	
	return 0.5 / M_PI * ( -(br1+br2+br3).kcross() + (Li * log(lva / lwb) * tauj + log(lwa / lva) * omega(wa, taui, tauj) +
						  log(lvb / lwb) * omega(vb, taui, tauj)) );
}

