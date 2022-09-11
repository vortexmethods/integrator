// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompI3DK_Duffy
\author Иванова Юлия Витальевна
\author Марчевский Илья Константинович
\author Хорошева Анна Александровна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompI3DK-Duffy.h"

#include <iomanip> 


#define AUTOSPLIT

const int nrefine = 3;
const double epsRel = 1e-7; 

CompI3DK_Duffy::CompI3DK_Duffy(const Database<3>& db_, const Parallel& par_, const Gausspoints<3>* const gp_, const GausspointsCube<2>* const gp2_, const GausspointsCube<3>* const gp3_)
	:  Computer(db_, par_, gp_), gp2(gp2_), gp3(gp3_) {};
CompI3DK_Duffy::~CompI3DK_Duffy() {};

//J_3D(K_i, K_j) - повторный интеграл от ньютоновского потенциала (по алгоритму Тейлора-Даффи)

//*********************************************************************************************************************


i2D CompI3DK_Duffy::RenumerationSosed(int i, int j)
{
	int outerVertexI = db.topo[i].difference(db.topo[j])[0]; 
	int outerVertexJ = db.topo[j].difference(db.topo[i])[0];

	int shiftI = (int)db.topo[i].member(outerVertexI);
	int shiftJ = (int)db.topo[j].member(outerVertexJ);

	return { shiftI, shiftJ };
}

i2D CompI3DK_Duffy::RenumerationContact(int i, int j)
{
	int comVert = db.topo[i].intersection(db.topo[j])[0];
	int shiftI = (int)db.topo[i].member(comVert);
	int shiftJ = (int)db.topo[j].member(comVert);

	return { shiftI, shiftJ };	
}


v3D CompI3DK_Duffy::DistContact(double xi1, double xi2, double eta1, double eta2, const i3D& nodes_i, const i3D& nodes_j)
{
	const v3D& P11 = db.node[nodes_i[0]];
	const v3D& P12 = db.node[nodes_i[1]];
	const v3D& P13 = db.node[nodes_i[2]];
	const v3D& P22 = db.node[nodes_j[1]];
	const v3D& P23 = db.node[nodes_j[2]];
	
	return (eta1 - xi1) * P11 + (xi1 - xi2) * P12 + xi2 * P13 - (eta1 - eta2) * P22 - eta2 * P23;
}

v3D CompI3DK_Duffy::DistSosed(double xi2, double u1, double u2, const i3D& nodes_i, const i3D& nodes_j)
{
	const v3D& P11 = db.node[nodes_i[0]];
	const v3D& P12 = db.node[nodes_i[1]];
	const v3D& P13 = db.node[nodes_i[2]];
	const v3D& P22 = db.node[nodes_j[1]];
	const v3D& P23 = db.node[nodes_j[2]];

	return u1 * (P11-P12) + u2*(P12-P23)+xi2*(P13-P23);
}

double CompI3DK_Duffy::IntDuffyContact(const i2D& ii, const i2D& jj, size_t& refineLevel)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);


	std::vector<std::function<double(v3D)>> fu = {
		[&](const v3D& r) {return r[1] / DistContact(1.0, r[0], r[1], r[1] * r[2], nodes_i, nodes_j).length(); } ,
		[&](const v3D& r) {return r[1] / DistContact(r[1], r[1] * r[2], 1.0, r[0], nodes_i, nodes_j).length(); } 
	};
	
	//int refineLevel = 3;
	
	double ints = 0.0;
	int refs = 0;

#ifdef AUTOSPLIT
	for (int i = 0; i < fu.size(); ++i)
	{
		auto [ti, tr] = gp3->integrateEpsRel<double>(fu[i], epsRel);
		ints += ti;
		refs = std::max(refs, tr);
	}
	ints *= 0.3333333333333333 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
	refineLevel = refs;
#else
	for (int i = 0; i < fu.size(); ++i)
	{
		double ti = gp2->integrate<double>(fu[i], refineLevel);
		ints += ti;
	}
	ints *= 0.3333333333333333 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
#endif
	return ints;

}

double CompI3DK_Duffy::IntDuffySosed(const i2D& ii, const i2D& jj, size_t& refineLevel)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]+1);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]+1);

	std::vector<std::function<double(v2D)>> fu = {
		[&](const v2D& r) {return r[0] / DistSosed(1.0 - r[0] + r[0] * r[1], -r[0], -r[0] * r[1], nodes_i, nodes_j).length(); },
		[&](const v2D& r) {return r[0] / DistSosed(1.0-r[0], r[0], r[0]*r[1], nodes_i, nodes_j).length(); },
	    [&](const v2D& r) {return r[0] / DistSosed(1.0-r[0], -r[1] * r[0], r[0]*(1.0-r[1]), nodes_i, nodes_j).length(); },
		[&](const v2D& r) {return r[0] / DistSosed(1.0-r[1]*r[0], r[1] * r[0], -r[0]*(1.0-r[1]), nodes_i, nodes_j).length(); },
		[&](const v2D& r) {return r[0] / DistSosed(1.0, -r[0] * r[1], -r[0], nodes_i, nodes_j).length(); },
		[&](const v2D& r) {return r[0] / DistSosed(1.0 - r[0], r[0] * r[1], r[0], nodes_i, nodes_j).length(); } 
	};

	//int refineLevel = 3;
	double ints = 0.0;
	int refs = 0;

#ifdef AUTOSPLIT
	for (int i = 0; i < fu.size(); ++i)
	{
		auto [ti, tr] = gp2->integrateEpsRel<double>(fu[i], epsRel);
		ints += ti;
		refs = std::max(refs, tr);
	}
	ints *= 0.1666666666666667 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
	refineLevel = refs;
#else
	for (int i = 0; i < fu.size(); ++i)
	{
		double ti = gp2->integrate<double>(fu[i], refineLevel);
		ints += ti;
	}
	ints *= 0.1666666666666667 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
#endif
	return ints;
}

/////////////////////////////////////////////////////////////////////////////////////////////


//*********************************************************************************************************************

double CompI3DK_Duffy::evaluate(int i, int j)
{
	size_t locRefine = nrefine;

	if (i==j)
	//if (&(db.topo[i]) == &(db.topo[j]))	
	{
		return 0.0;
	}
	if (db.ifSosed({ i, j }))
	{
		i2D shifts = RenumerationSosed(i, j);
		auto Integral = IntDuffySosed({ i, shifts[0] }, { j, shifts[1] }, locRefine);
		//if (locRefine == maxRefine)
		//	std::cout << "Sosed, (i,j) = " << i2D{ i,j } << std::endl;
		//std::cout << "locRefine = " << locRefine << std::endl;
		return Integral;
	}

	if (db.ifContact({ i, j }))
	{
		i2D shifts = RenumerationContact(i, j);	
		auto Integral = IntDuffyContact({ i, shifts[0] }, { j, shifts[1] }, locRefine);
		//if (locRefine == maxRefine)
		//	std::cout << "Contact, (i,j) = " << i2D{ i,j } << std::endl;
		return Integral;
	}


	// Дальние панели
	auto onePt = [i, j, this](v3D pt)
	{
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

		return 0.25 / M_PI * ((Phi - va * Theta) & db.nrm[j]);
	};


	const int constnrefine = 1;
	const double epsRel = 1e-5;
#define AUTOSPLIT

#ifdef AUTOSPLIT		
	auto [res, nrefine] = gp->integrateEpsRel<double>(onePt, i, epsRel);
	//refines.push_back(nrefine);
#else
	double res = gp->integrate<double>(onePt, i, constnrefine);
	//refines.push_back(constnrefine);
#endif

	return res;





}



// & - скалярное произведение
// ^ - вектороне произведение 

