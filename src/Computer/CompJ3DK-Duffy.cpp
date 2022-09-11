// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompJ3DK
\author Иванова Юлия Витальевна
\author Марчевский Илья Константинович
\author Хорошева Анна Александровна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompJ3DK-Duffy.h"

#include <iomanip> 


#define AUTOSPLIT

const int nrefine = 3;
const double epsRel = 1e-7; 

const double eps_zero = 2e-6;
const double eps_zero2 = eps_zero * eps_zero;

const double eps_psi_theta = eps_zero;
const double eps_psi_theta2 = eps_psi_theta * eps_psi_theta;

CompJ3DK_Duffy::CompJ3DK_Duffy(const Database<3>& db_, const Parallel& par_, const Gausspoints<3>* const gp_, const GausspointsCube<2>* const gp2_, const GausspointsCube<3>* const gp3_)
	:  Computer(db_, par_, gp_), gp2(gp2_), gp3(gp3_) {};
CompJ3DK_Duffy::~CompJ3DK_Duffy() {};

//J_3D(K_i, K_j) - повторный интеграл от ньютоновского потенциала (по алгоритму Тейлора-Даффи)

//*********************************************************************************************************************


i2D CompJ3DK_Duffy::RenumerationSosed(int i, int j)
{
	int outerVertexI = db.topo[i].difference(db.topo[j])[0]; 
	int outerVertexJ = db.topo[j].difference(db.topo[i])[0];

	int shiftI = (int)db.topo[i].member(outerVertexI);
	int shiftJ = (int)db.topo[j].member(outerVertexJ);

	return { shiftI, shiftJ };
}

i2D CompJ3DK_Duffy::RenumerationContact(int i, int j)
{
	int comVert = db.topo[i].intersection(db.topo[j])[0];
	int shiftI = (int)db.topo[i].member(comVert);
	int shiftJ = (int)db.topo[j].member(comVert);

	return { shiftI, shiftJ };	
}


v3D CompJ3DK_Duffy::DistContact(double xi1, double xi2, double eta1, double eta2, const i3D& nodes_i, const i3D& nodes_j)
{
	const v3D& P11 = db.node[nodes_i[0]];
	const v3D& P12 = db.node[nodes_i[1]];
	const v3D& P13 = db.node[nodes_i[2]];
	const v3D& P22 = db.node[nodes_j[1]];
	const v3D& P23 = db.node[nodes_j[2]];
	
	return (eta1 - xi1) * P11 + (xi1 - xi2) * P12 + xi2 * P13 - (eta1 - eta2) * P22 - eta2 * P23;
}

v3D CompJ3DK_Duffy::DistSosed(double xi2, double u1, double u2, const i3D& nodes_i, const i3D& nodes_j)
{
	const v3D& P11 = db.node[nodes_i[0]];
	const v3D& P12 = db.node[nodes_i[1]];
	const v3D& P13 = db.node[nodes_i[2]];
	const v3D& P23 = db.node[nodes_j[2]];

	return u1 * (P11 - P12) + u2 * (P12 - P23) + xi2 * (P13 - P23);
}

template<typename T>
T cubpow(T x)
{
	return x * x * x;
}

v3D CompJ3DK_Duffy::IntDuffyContact(const i2D& ii, const i2D& jj, size_t& refineLevel)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);


	std::vector<std::function<v3D(v3D)>> fu = {
		[&](const v3D& r) {
			auto DC1 = DistContact(1.0, r[0], r[1], r[1] * r[2], nodes_i, nodes_j);
			return (r[1] / cubpow(DC1.length()) ) * DC1; } ,
		[&](const v3D& r) {
			auto DC2 = DistContact(r[1], r[1] * r[2], 1.0, r[0], nodes_i, nodes_j);
			return ( r[1] / cubpow(DC2.length()) ) * DC2; }
	};
	
	//int refineLevel = 3;
	
	v3D ints = { 0.0,0.0,0.0 };
	int refs = 0;

#ifdef AUTOSPLIT
	for (int i = 0; i < fu.size(); ++i)
	{
		auto [ti, tr] = gp3->integrateEpsRel<v3D>(fu[i], epsRel);
		ints += ti;
		refs = std::max(refs, tr);
	}
	ints *= 0.5 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
	refineLevel = refs;
#else
	for (int i = 0; i < fu.size(); ++i)
	{
		double ti = gp2->integrate<v3D>(fu[i], refineLevel);
		ints += ti;
	}
	ints *= 0.5 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
#endif
	return ints;

}

v3D CompJ3DK_Duffy::IntDuffySosed(const i2D& ii, const i2D& jj, size_t& refineLevel)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]+1);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]+1);

	std::vector<std::function<v3D(v2D)>> fu = {
		
		[&](const v2D& r) {
			auto DS1 = DistSosed(1.0 - r[0] + r[0] * r[1], -r[0], -r[0] * r[1], nodes_i, nodes_j); 
			return (r[0] / cubpow(DS1.length())) * DS1; },
		[&](const v2D& r) {
			auto DS2 = DistSosed(1.0 - r[0], r[0], r[0] * r[1], nodes_i, nodes_j);
			return (r[0] / cubpow(DS2.length()))*DS2; },
	    [&](const v2D& r) {
			auto DS3 = DistSosed(1.0 - r[0], -r[1] * r[0], r[0] * (1.0 - r[1]), nodes_i, nodes_j);
			return (r[0] / cubpow(DS3.length()))*DS3; },
		[&](const v2D& r) {
			auto DS4 = DistSosed(1.0 - r[1] * r[0], r[1] * r[0], -r[0] * (1.0 - r[1]), nodes_i, nodes_j);
			return (r[0] / cubpow(DS4.length()))*DS4; },
		[&](const v2D& r) {
			auto DS5 = DistSosed(1.0, -r[0] * r[1], -r[0], nodes_i, nodes_j);
			return (r[0] / cubpow( DS5.length()))*DS5; },
		[&](const v2D& r) {
			auto DS6 = DistSosed(1.0 - r[0], r[0] * r[1], r[0], nodes_i, nodes_j);
			return (r[0] / cubpow(DS6.length()))*DS6; }
	};

	//int refineLevel = 3;
	v3D ints = { 0.0,0.0,0.0 };
	int refs = 0;

#ifdef AUTOSPLIT
	for (int i = 0; i < fu.size(); ++i)
	{
		auto [ti, tr] = gp2->integrateEpsRel<v3D>(fu[i], epsRel);
		ints += ti;
		refs = std::max(refs, tr);
	}
	ints *= 0.5/ M_PI * db.measure[ii[0]] * db.measure[jj[0]];
	refineLevel = refs;
#else
	for (int i = 0; i < fu.size(); ++i)
	{
		double ti = gp2->integrate<double>(fu[i], refineLevel);
		ints += ti;
	}
	ints *= 0.5 / M_PI * db.measure[ii[0]] * db.measure[jj[0]];
#endif
	return ints;
}

/////////////////////////////////////////////////////////////////////////////////////////////

p13D CompJ3DK_Duffy::ThetaPsi(const v3D& pt, int j)
{
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

	double rac = (ova & tauc), rbc = (ovb & tauc), rba = (ovb & taua), rca = (ovc & taua), rcb = (ovc & taub), rab = (ova & taub);

	//auto w1 = ova ^ ovb;
	//auto w2 = (ova ^ ovb) & ovc;
	//auto w3 = 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova);

	double term1, term2, term3;

	//Проверка частных случаев
	//if ((fabs(rac + 1.0) < 0.5*eps_psi_theta2) && (fabs(rbc + 1.0) < 0.5*eps_psi_theta2))
	if (fabs(rbc + 1.0) < 0.5 * eps_psi_theta2)
		term1 = log(lvb / lva);
	else
		term1 = log((lva * (1.0 + (ova & tauc))) / (lvb * (1.0 + (ovb & tauc))));

	//if ((fabs(rba + 1.0) < 0.5 * eps_psi_theta2) && (fabs(rca + 1.0) < 0.5 * eps_psi_theta2))
	if (fabs(rca + 1.0) < 0.5 * eps_psi_theta2)
		term2 = log(lvc / lvb);
	else
		term2 = log((lvb * (1.0 + (ovb & taua))) / (lvc * (1.0 + (ovc & taua))));

	//if ((fabs(rcb + 1.0) < 0.5 * eps_psi_theta2) && (fabs(rab + 1.0) < 0.5 * eps_psi_theta2))
	if (fabs(rab + 1.0) < 0.5 * eps_psi_theta2)
		term3 = log(lva / lvc);
	else
		term3 = log((lvc * (1.0 + (ovc & taub))) / (lva * (1.0 + (ova & taub))));

	const v3D& Psi = term1 * tauc + term2 * taua + term3 * taub;

	//Величина, равная по модулю телесному углу, под которым видна панель из точки наблюдения
	double Theta = 2.0 * atan2((ova ^ ovb) & ovc, 1.0 + (ova & ovb) + (ovb & ovc) + (ovc & ova));

	return { Theta, Psi };
}

v3D CompJ3DK_Duffy::J3D(const v3D& pt, int j)
{
	auto [Theta, Psi] = ThetaPsi(pt, j);
	return (0.25 / M_PI) * ((Theta * db.nrm[j]) + (Psi ^ db.nrm[j]));
}



//*********************************************************************************************************************

v3D CompJ3DK_Duffy::evaluate(int i, int j)
{
	size_t locRefine = nrefine;

	if (i==j)
	//if (&(db.topo[i]) == &(db.topo[j]))	
	{
		return { 0.0,0.0,0.0 };
	}
	if (db.ifSosed({ i, j }))
	{
		i2D shifts = RenumerationSosed(i, j);
		auto Integral = IntDuffySosed({ i, shifts[0] }, { j, shifts[1] }, locRefine);
		if (locRefine == maxRefine)
			std::cout << "Sosed, (i,j) = " << i2D{ i,j } << std::endl;
		//std::cout << "locRefine = " << locRefine << std::endl;
		return Integral;
	}

	if (db.ifContact({ i, j }))
	{
		i2D shifts = RenumerationContact(i, j);	
		auto Integral = IntDuffyContact({ i, shifts[0] }, { j, shifts[1] }, locRefine);
		if (locRefine == maxRefine)
			std::cout << "Contact, (i,j) = " << i2D{ i,j } << std::endl;
		return Integral;
	}

	auto func = [&](const v3D& r) { return J3D(r, j); };


#ifdef AUTOSPLIT		
	auto [res, nrefine] = gp->integrateEpsRel<v3D>(func, i, epsRel);
	//refines.push_back(nrefine);
#else
	v3D res = gp->integrate<v3D>(func, i, nrefine);
	//refines.push_back(nrefine);
#endif

	//std::cout << "res = " << res << std::endl;

	return res;

}



// & - скалярное произведение
// ^ - вектороне произведение 

