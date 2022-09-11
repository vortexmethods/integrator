// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompJ2DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompJ2DM.h"

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

const double eps_zero = 1e-6;

CompJ2DM::CompJ2DM(const Database<2>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompJ2DM::~CompJ2DM() {};

//J_2D(M_i, K_j) - интеграл от градиента логарифмического потенциала 

v2D CompJ2DM::evaluate(int i, int j)
{
	//Точка наблюдения
	const v2D& pt = db.point[i];

	//Векторы из вершин панели в точку наблюдения
	const v2D& va = pt - db.node[db.topo[j][0]];
	const v2D& vb = pt - db.node[db.topo[j][1]];

	//
	v2D ova(va), ovb(vb);

	//Модули векторов v_a и v_b соответственно
	double lva = ova.normalize();
	double lvb = ovb.normalize();

	//tau_j - единичный направляющий вектор панели 
	v2D tau = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];

	//std::cout << "i = " << i << " " << "j = " << j << " " << tau[0] << " " << tau[1] << std::endl;

	//L_j - длина прямолинейной панели
	double L = tau.normalize();

	//tau.normalize();

	//std::cout << "---------------------------------------------------------" << std::endl;
	//std::cout << "i = " << i << " " << "j = " << j << " " << tau[0] << " " << tau[1] << std::endl;

	v2D tau_k = { tau[1], -tau[0] };

	//Величина, равная по модулю углу, под которым видна панель из точки наблюдения
	double Theta = atan2((ova ^ ovb), (ova & ovb));

	//НОВОЕ
	//Проверка частного случая

	auto sqr = [](double x) {return x * x; };

	double nrmva = sqr(va[0] * va[0] + va[1] * va[1]);

	//std::cout << " nrmva = " << nrmva << std::endl;

	double nrmvb = sqr(vb[0] * vb[0] + vb[1] * vb[1]);

	//std::cout << " nrmvb = " << nrmvb << std::endl;

	//std::cout << " lva = " << lva << " " << "lvb = " << lvb << std::endl;

	//if ((nrmva > eps_zero) && (nrmvb > eps_zero) && (fabs((nrmva + nrmvb) / L - 1) < eps_zero))
	//ВЫДАЕТ НУЛИ!
	if (/*(lva/L > eps_zero) && (lvb/L > eps_zero) &&*/ (fabs((lva + lvb) / L - 1.0) < eps_zero))
	{
		//std::cout << "i = " << i << " " << " , " << " " << "j = " << j << std::endl;
		//std::cout << "log(lva / lvb) = " << log(lva / lvb) << std::endl;
		//std::cout << "(log(lva / lvb) * tau) = " << (log(lva / lvb) * tau) << std::endl;
		//std::cout << "res = " << 0.5 / M_PI * (log(lva / lvb) * tau) << std::endl;

		return 0.5 / M_PI * (log(lva / lvb) * tau);
	}
	else
	{
		//std::cout << "i = " << i << " " << " , " << " " << "j = " << j << std::endl;
		//std::cout << "res = " << 0.5 / M_PI * (log(lva / lvb) * tau - Theta * tau_k) << std::endl;

		return 0.5 / M_PI * (log(lva / lvb) * tau - Theta * tau_k);
	}


}
