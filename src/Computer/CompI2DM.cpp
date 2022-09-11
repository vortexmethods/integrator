// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompI2DM
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompI2DM.h"

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

const double eps_zero = 1e-6;

CompI2DM::CompI2DM(const Database<2>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompI2DM::~CompI2DM() {};

//I_2D(M_i, K_j) - интеграл от логарифмического потенциала 

double CompI2DM::evaluate(int i, int j)
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
	
	//std::cout << "----------------------------------------------------------------------------------------" << std::endl;

	//L_j - длина прямолинейной панели
	double L = tau.normalize();

	//std::cout << "i = " << i << " " << "j = " << j << " " << tau[0] << " " << tau[1] << std::endl;

	//h - произвольный постоянный вектор
	const double mh = 1.0;

	//Величина, равная по модулю углу, под которым видна панель из точки наблюдения
	double Theta = atan2((ova ^ ovb), (ova & ovb));

	//НОВОЕ!
	//Проверка частного случая 
	if ((va.length()/L < eps_zero) || (vb.length()/L < eps_zero))
	{	
		return 0.5 / M_PI * ( 1.0 - log(L / mh) ) * L;
			
	}
	else
		return 0.5 / M_PI * ( ( (log(lvb/mh) * vb - log(lva/mh) * va) & tau )  + L + Theta * ( va ^ tau ) );


	//return 0.5/M_PI * ( ( (log(lvb/mh) * vb - log(lva/mh) * va) & tau )  + L + Theta * ( va ^ tau ) );
}
