// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompI2DK
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompI2DK.h"

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

const double eps_zero = 1e-6;

CompI2DK::CompI2DK(const Database<2>& db_, const Parallel& par_) : Computer(db_, par_) {};
CompI2DK::~CompI2DK() {};  


//Величина, равная по модулю углу, под которым видна панель из точки наблюдения
double CompI2DK::Theta(const v2D& v1, const v2D& v2)
{
	return atan2((v1 ^ v2), (v1 & v2));
}

double CompI2DK::psi(const v2D& v, const v2D& ti, const v2D& tj)
{
	return ((v & ti) * (v & tj)) - ((v ^ ti) * (v ^ tj));
}

double CompI2DK::evaluate(int i, int j)
{
	const v2D& pt1 = db.node[db.topo[i][0]];
	const v2D& pt2 = db.node[db.topo[i][1]];

	//Векторы из вершин панели в точку наблюдения
	const v2D& va = pt1 - db.node[db.topo[j][0]];
	const v2D& vb = pt1 - db.node[db.topo[j][1]];

	const v2D& wa = pt2 - db.node[db.topo[j][0]];
	const v2D& wb = pt2 - db.node[db.topo[j][1]];

	//
	v2D ova(va), ovb(vb), owa(wa), owb(wb);

	//Модули векторов v_a и v_b соответственно
	double lva = ova.normalize();
	double lvb = ovb.normalize();
	double lwa = owa.normalize();
	double lwb = owb.normalize();

	//tau_j - единичный направляющий вектор панели 
	v2D tauj = db.node[db.topo[j][1]] - db.node[db.topo[j][0]];
	v2D taui = db.node[db.topo[i][1]] - db.node[db.topo[i][0]];

	//L_j - длина прямолинейной панели
	double Lj = tauj.normalize();
	double Li = taui.normalize();

	double LiEps_zero = Li * eps_zero;

	//h - произвольный постоянный вектор
	//const double mh = 1.0;



	//if ((lva < LiEps_zero) && (lwb < LiEps_zero))
	if (i == j)
		return 0.25 / M_PI * (3.0 - 2.0 * log(Lj /*/ mh*/)) * Lj * Lj;

	v2D c = (0.5 * (db.node[db.topo[i][1]] + db.node[db.topo[i][0]])) - (0.5 * (db.node[db.topo[j][1]] + db.node[db.topo[j][0]]));
	v2D c_k = -c.kcross();
		
	v2D tauj_k = -tauj.kcross();

	double br1, br21;
		
	double br22 = (((c & taui) * tauj) + ((c & tauj) * taui)) & c_k;
	double br23 = 0.25 * (Lj * Lj - Li * Li) * (taui & tauj_k);
		
	double br3;
	double Delta = 0.0;

	//if (db.ifSosed({ i, j }) ) 
	bool zerowa = (lwa < LiEps_zero);
	bool zerovb = (lvb < LiEps_zero);

	if (zerowa || zerovb)
	{
		if (zerowa)
		{
			double thVaVb = Theta(va, vb);			
			double thVbWb = Theta(vb, wb);

			br1 = ((Lj * thVbWb) * taui - (Li * thVaVb) * tauj) & c_k;
			br21 = -0.5 * (thVbWb - thVaVb);
			br3 = psi(va, taui, tauj) * log(lva /*/mh*/) - psi(vb, taui, tauj) * log(lvb /*/mh*/) + \
				  psi(wb, taui, tauj) * log(lwb /*/mh*/);
		}
		else if (zerovb)
		{
			double thWaWb = Theta(wa, wb);
			double thVaWa = Theta(va, wa);
			
			br1 = ((Lj * thVaWa) * taui - (Li * thWaWb) * tauj) & c_k;
			br21 = 0.5 * (thVaWa - thWaWb);
			br3 = psi(va, taui, tauj) * log(lva /*/mh*/) - \
				  psi(wa, taui, tauj) * log(lwa /*/mh*/) + psi(wb, taui, tauj) * log(lwb /*/mh*/);
		}
		
		//br1 = (Lj * ((zerowa ? 0.0 : Theta(va, wa)) + (zerovb ? 0.0 : Theta(vb, wb))) * taui \
			 - Li * ((zerovb ? 0.0 : Theta(va, vb)) + (zerowa ? 0.0 : Theta(wa, wb))) * tauj) & c_k;

		//br21 = 0.5 * ((zerowa ? 0.0 : Theta(va, wa) - Theta(wa, wb)) - (zerovb ? 0.0 : Theta(vb, wb) - Theta(va, vb)));
		
		//br3 = psi(va, taui, tauj) * log(lva /*/mh*/) - (zerovb ? 0.0 : psi(vb, taui, tauj) * log(lvb /*/mh*/)) - \
					 (zerowa ? 0.0 : psi(wa, taui, tauj) * log(lwa /*/mh*/)) + psi(wb, taui, tauj) * log(lwb /*/mh*/);
		
		double thij = Theta(taui, tauj);	
		Delta = -0.0625 / M_PI * (Li * Li + Lj * Lj) * thij * sin(thij);
	}
	else
	{
		double thVaVb = Theta(va, vb);
		double thWaWb = Theta(wa, wb);
		double thVaWa = Theta(va, wa);
		double thVbWb = Theta(vb, wb);

		br1 = ((Lj * (thVaWa + thVbWb)) * taui - (Li * (thVaVb + thWaWb)) * tauj) & c_k;

		br21 = 0.5 * (thVaWa - thWaWb - thVbWb + thVaVb);

		br3 = psi(va, taui, tauj) * log(lva /*/mh*/) - psi(vb, taui, tauj) * log(lvb /*/mh*/) - \
		      psi(wa, taui, tauj) * log(lwa /*/mh*/) + psi(wb, taui, tauj) * log(lwb /*/mh*/);
	}
	
	return 0.25 / M_PI * (3.0 * Li * Lj + br1 + br21 * (br22 + br23) + br3) + Delta;
}

// & - скалярное произведение
// ^ - вектороне произведение 