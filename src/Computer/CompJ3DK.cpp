// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c реализацией класса CompJ3DK
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#include "CompJ3DK.h"

#include <iomanip> 


#define AUTOSPLIT

const double eps_zero = 1e-6;
const double eps_zero2 = eps_zero * eps_zero;

const double eps_psi_theta = eps_zero;
const double eps_psi_theta2 = eps_psi_theta * eps_psi_theta;


const int nrefine = 0;
const double epsRel = 1e-5; 

CompJ3DK::CompJ3DK(const Database<3>& db_, const Parallel& par_, const Gausspoints<3>* const gp_) : Computer(db_, par_, gp_) {};
CompJ3DK::~CompJ3DK() {};

//J_3D(K_i, K_j) - повторный интеграл от градиента ньютоновского потенциала

//*********************************************************************************************************************


p13D CompJ3DK::ThetaPsi(const v3D& pt, int j)
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

v3D CompJ3DK::J3D(const v3D& pt, int j)
{
	auto [Theta, Psi] = ThetaPsi(pt, j);
	return (0.25 / M_PI) * ((Theta * db.nrm[j]) + (Psi ^ db.nrm[j]));
}


i2D CompJ3DK::RenumerationSosed(int i, int j)
{
	int outerVertexI = db.topo[i].difference(db.topo[j])[0]; 
	int outerVertexJ = db.topo[j].difference(db.topo[i])[0];

	int shiftI = (int)db.topo[i].member(outerVertexI);
	int shiftJ = (int)db.topo[j].member(outerVertexJ);

	return { shiftI, shiftJ };
}

i2D CompJ3DK::RenumerationContact(int i, int j)
{
	int comVert = db.topo[i].intersection(db.topo[j])[0];
	int shiftI = (int)db.topo[i].member(comVert);
	int shiftJ = (int)db.topo[j].member(comVert);

	return { shiftI, shiftJ };	
}


p13D CompJ3DK::SingSosed(const v3D& M, const i2D& jj)
{
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);

	//Векторы из вершин панели в точку наблюдения
	const v3D& va = M - db.node[nodes_j[1]];
	const v3D& vb = M - db.node[nodes_j[2]];

	v3D ova(va), ovb(vb);

	double lva = ova.normalize();
	double lvb = ovb.normalize();

	//Орты векторов, направленных вдоль сторон треугольной панели K_j, лежащих против соответствующих вершин
	v3D taua = (db.node[nodes_j[0]] - db.node[nodes_j[2]]).unit();
	v3D taub = (db.node[nodes_j[1]] - db.node[nodes_j[0]]).unit();

	v3D tauc = (db.node[nodes_j[2]] - db.node[nodes_j[1]]); //Здесь вектор еще не нормирован
	double iLc = 1.0 / tauc.normalize();                    //После этой строки tauc нормирован

	return { 2.0 * (atan2((ova ^ taub) & tauc, (taub - tauc) & (taub + ova)) - atan2((ovb ^ taua) & tauc, (taua - tauc) & (taua - ovb))),
			 log((lvb * (tauc & (tauc - ovb))) / (lva * (tauc & (tauc - ova)))) * tauc -
			 log(lva * taub & (taub + ova) * iLc) * taub -
			 log(lvb * taua & (taua - ovb) * iLc) * taua };
}



p13D CompJ3DK::SingContact(const v3D& M, const i2D& ii, const i2D& jj)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);

	//Векторы из вершин панели в точку наблюдения
	v3D ovc = M - db.node[nodes_j[0]];
	double lvc = ovc.normalize();

	//Орты векторов, направленных вдоль сторон треугольной панели K_j, лежащих против соответствующих вершин
	v3D taua = (db.node[nodes_j[0]] - db.node[nodes_j[2]]).unit();
	v3D taub = (db.node[nodes_j[1]] - db.node[nodes_j[0]]).unit();

	v3D e = (db.nrm[ii[0]] ^ db.nrm[jj[0]]);
	if (e.length2() < eps_zero2)
		e = taub;
	else
		e.normalize();

	const v3D& nrmi = db.nrm[ii[0]];
	const v3D& nrmj = db.nrm[jj[0]];

	double deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
	double deltaB = atan2((e ^ taub) & nrmj, e & taub);

	if (((M_PI - fabs(deltaA)) < eps_zero) || ((M_PI - fabs(deltaB)) < eps_zero))
	{
		e *= -1.0;
		deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
		deltaB = atan2((e ^ taub) & nrmj, e & taub);
	}

	if ((deltaA * deltaB < eps_zero) && (fabs(deltaA - deltaB) > M_PI))
	{
		e *= -1.0;
		deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
		deltaB = atan2((e ^ taub) & nrmj, e & taub);
	}

	double ari = sqrt(db.measure[ii[0]]);

	return { 2.0*(atan2((ovc ^ taua) & e, (e - ovc) & (e - taua)) + atan2( (ovc ^ taub) & e, (e - ovc) & (e + taub) )),
			-(taua * log((lvc * (1 + (taua & ovc))) / ari) + taub * log((lvc * (1 - (taub & ovc))) / ari)) };
}

v3D CompJ3DK::JSingSosed(const v3D& M, const i2D& jj)
{
	auto [thSingSosed, psSingSosed] = SingSosed(M, jj);

	const v3D nrmj = db.nrm[jj[0]];
	return (0.25 / M_PI) * ((thSingSosed * nrmj) + (psSingSosed ^ nrmj));
}

v3D CompJ3DK::JRegSosed(const v3D& M, const i2D& jj)
{
	return J3D(M, jj[0]) - JSingSosed(M, jj);
}

p13D CompJ3DK::RegContact(const v3D& M, const i2D& ii, const i2D& jj)
{
	auto [Theta, Psi] = ThetaPsi(M, jj[0]);
	auto [SingTheta, SingPsi] = SingContact(M, ii, jj);
	return { Theta - SingTheta, Psi - SingPsi };
}

v3D CompJ3DK::IntJSingSosed(const i2D& ii, const i2D& jj)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);

	const v3D& taua = (db.node[nodes_j[0]] - db.node[nodes_j[2]]).unit();
	const v3D& taub = (db.node[nodes_j[1]] - db.node[nodes_j[0]]).unit();
	const v3D& tauc = (db.node[nodes_j[2]] - db.node[nodes_j[1]]).unit();

	//std::cout << db.node[nodes_i[0]] << " " << db.node[nodes_i[1]] << " " << db.node[nodes_i[2]] << std::endl;

	double alpha = angle(db.node[nodes_i[0]] - db.node[nodes_i[2]], db.node[nodes_i[1]] - db.node[nodes_i[2]]);
	double beta  = angle(db.node[nodes_i[2]] - db.node[nodes_i[1]], db.node[nodes_i[0]] - db.node[nodes_i[1]]);
	double gamma = angle(db.node[nodes_j[2]] - db.node[nodes_j[1]], db.node[nodes_j[0]] - db.node[nodes_j[1]]);
	double delta = angle(db.node[nodes_j[1]] - db.node[nodes_j[2]], db.node[nodes_j[0]] - db.node[nodes_j[2]]);

	double nu = M_PI - alpha - beta;

	const v3D& nrmi = db.nrm[ii[0]];
	const v3D& nrmj = db.nrm[jj[0]];

	double xi = atan2((nrmi ^ nrmj) & tauc, nrmi & nrmj);

	double sinXi = sin(xi), cosXi = cos(xi);

	// И.К. Добавил предвычисление синусов-косинусов
	double sinAlpha = sin(alpha), cosAlpha = cos(alpha);
	double sinBeta  = sin(beta),  cosBeta  = cos(beta);

	// И.К. Вычислим их здесь, чтобы однократно
	double sinNu = sinAlpha * cosBeta + cosAlpha * sinBeta, cosNu = sinAlpha * sinBeta - cosAlpha * cosBeta;

	double sinGamma = sin(gamma), cosGamma = cos(gamma);
	double sinDelta = sin(delta), cosDelta = cos(delta);


	double cosSigma = -(cosAlpha * cosDelta + cosXi * sinAlpha * sinDelta);
	double cosMu = -(cosBeta * cosGamma + cosXi * sinBeta * sinGamma);
	double cosLambda = -(cosAlpha * cosGamma - cosXi * sinAlpha * sinGamma);
	double cosTheta = -(cosBeta * cosDelta - cosXi * sinBeta * sinDelta);

	double q_alpha_beta = sinNu * log(tan(alpha / 2) * tan(nu / 2)) / sin(beta) + sinNu * log(tan(beta / 2) * tan(nu / 2)) / sin(alpha) \
		                + log(tan(alpha / 2) * tan(beta / 2));

	auto phi = [sinXi](double sinAlpha, double cosAlpha, double sinGamma, double cosGamma, double cosLambda)
	{
		//std::cout << sinXi * sinAlpha * sinGamma / (1.0 - cosAlpha + cosGamma + cosLambda) << std::endl;
		return 2.0*atan2(sinXi * sinAlpha * sinGamma, 1.0 - cosAlpha + cosGamma + cosLambda);
	};// []phi(...)

	auto q_Theta_Psi = [phi, sinNu, cosNu, cosXi, sinXi](
		                     double sinAlpha, double cosAlpha, \
		                     double sinBeta,  double cosBeta, \
							 double sinGamma, double cosGamma, \
		                     double cosMu, \
		                     double cosLambda)
	{
		double phi1 = phi(sinAlpha,  cosAlpha, sinGamma,  cosGamma, cosLambda); // xi, alpha, gamma, lambda
		double phi2 = phi(sinAlpha, -cosAlpha, sinGamma, -cosGamma, cosLambda); //xi, pi - alpha, pi - gamma, lambda

		double iSinAlphaSin2mu = 1.0 / (sinAlpha * (1.0 - cosMu * cosMu));

		return v2D{
			phi1 + sinGamma * sinNu * iSinAlphaSin2mu * (\
			(cosBeta * sinGamma - cosXi * sinBeta * cosGamma) * phi2 + \
			sinXi * sinBeta * (0.5 * (1.0 + cosMu) * log((1.0 + cosBeta) / (1.0 - cosNu)) + \
								0.5 * (1.0 - cosMu) * log((1.0 - cosBeta) / (1.0 + cosNu)) + \
								log((1.0 + cosLambda) / (1.0 - cosGamma)))),

			1.5 - iSinAlphaSin2mu * ( \
				sinBeta * (cosNu + cosMu * cosLambda) * log(1.0 + cosLambda) + \
				sinNu * (cosBeta + cosMu * cosGamma) * log(1.0 - cosGamma) + \
				sinBeta * (1.0 - cosMu) * (cosNu - cosLambda) * log(sinBeta / sinNu) + \
				sinNu * sinBeta * (sinBeta * cosGamma - cosXi * sinGamma * cosBeta) * log((1.0 - cosNu) / (1.0 + cosBeta)) + \
				phi2 * sinXi * sinGamma * sinNu * sinBeta ) 
		};
	}; // []q_Theta_Psi(...)

	auto q_Theta_Psi_zero = [phi, sinNu, cosNu, cosXi, sinXi](
		double sinAlpha, double cosAlpha, \
		double sinBeta, double cosBeta, \
		double sinGamma, double cosGamma, \
		double cosMu, \
		double cosLambda)
	{
		//std::cout << "Exception" << std::endl;

		return v2D{
				0.0,
				//1.5 - (cosNu * sinBeta * log(1.0 + cosNu) + sinNu * cosBeta * log(1.0 - cosBeta) + sinBeta * sinNu * (sinBeta / (1.0-cosBeta) - sinNu / (1.0+cosNu)) + sinBeta * cosNu * log(sinBeta / sinNu)) / sinAlpha
				1.5 - (cosNu * sinBeta * log(1.0 + cosNu) + sinNu * cosBeta * log(1.0 - cosBeta) + sinAlpha - sinBeta + sinNu + sinBeta * cosNu * log(sinBeta / sinNu)) / sinAlpha
		};
	}; // []q_Theta_Psi_zero(...)

	v2D qTP_A, qTP_B;

	if ((fabs(xi) < eps_zero) && (fabs(beta - gamma) < eps_zero))
		qTP_A = q_Theta_Psi_zero(sinAlpha, cosAlpha, sinBeta, cosBeta, sinGamma, cosGamma, cosMu, cosLambda);
	else
		qTP_A = q_Theta_Psi(sinAlpha, cosAlpha, sinBeta, cosBeta, sinGamma, cosGamma, cosMu, cosLambda);
	
	if ((fabs(xi) < eps_zero) && (fabs(alpha - delta) < eps_zero))
		qTP_B = q_Theta_Psi_zero(sinBeta, cosBeta, sinAlpha, cosAlpha, sinDelta, cosDelta, cosSigma, cosTheta);
	else
		qTP_B = q_Theta_Psi(sinBeta,  cosBeta,  sinAlpha, cosAlpha, sinDelta, cosDelta, cosSigma, cosTheta );

	double IntegrateThSing = db.measure[ii[0]] * (qTP_A[0] + qTP_B[0]);
	v3D IntegratePsiSing = db.measure[ii[0]] * (qTP_A[1] * taub + qTP_B[1] * taua - q_alpha_beta * tauc);

	//cout << "IntSing = " << IntegrateThSing << " " << IntegratePsiSing << endl;

	return 0.25 / M_PI * (IntegrateThSing * nrmj + (IntegratePsiSing ^ nrmj));
}


v3D CompJ3DK::IntJRegSosed(int i, const i2D& jj, size_t refineLevel)
{
	auto fo = [this, jj](const v3D& r) { return JRegSosed(r, jj); };
	return gp->integrate<v3D>(fo, i, refineLevel);
}

std::pair<v3D, int> CompJ3DK::IntJRegSosedEpsRel(int i, const i2D& jj)
{
	auto fo = [this, jj](const v3D& r) { return JRegSosed(r, jj); };
	return gp->integrateEpsRel<v3D>(fo, i, epsRel);
}


p13D CompJ3DK::IntSingContact(const i2D& ii, const i2D& jj)
{
	const i3D& nodes_i = db.topo[ii[0]].rotateLeft(ii[1]);
	const i3D& nodes_j = db.topo[jj[0]].rotateLeft(jj[1]);

	//Орты векторов, направленных вдоль сторон треугольной панели K_j, лежащих против соответствующих вершин
	const v3D& taua = (db.node[nodes_j[0]] - db.node[nodes_j[2]]).unit();
	const v3D& taub = (db.node[nodes_j[1]] - db.node[nodes_j[0]]).unit();
		
	v3D e = (db.nrm[ii[0]] ^ db.nrm[jj[0]]);
	if (e.length2() < eps_zero2)
	{
		e = taub;
		if ((db.nrm[ii[0]] & db.nrm[jj[0]]) < 0)
			std::cout << "Orientation!" << std::endl;
	}
	else
		e.normalize();
	
	const v3D& nrmi = db.nrm[ii[0]];
	const v3D& nrmj = db.nrm[jj[0]];

	double deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
	double deltaB = atan2( (e ^ taub) & nrmj,  e & taub);

	if ( ((M_PI - fabs(deltaA)) < eps_zero) || ((M_PI - fabs(deltaB)) < eps_zero) )
	{
		e *= -1.0;
		deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
		deltaB = atan2(( e ^ taub) & nrmj,  e & taub);
	}

	if ((deltaA * deltaB < 0) && (fabs(deltaA - deltaB) > M_PI))
	{
		e *= -1.0;
		deltaA = atan2((-e ^ taua) & nrmj, -e & taua);
		deltaB = atan2( (e ^ taub) & nrmj,  e & taub);
	}

	double xi = atan2((nrmi ^ nrmj) & e, nrmi & nrmj);

	// И.К. Предвычисление всех синусов и косинусов
	double sinDeltaA = sin(deltaA), cosDeltaA = cos(deltaA);
	double sinDeltaB = sin(deltaB), cosDeltaB = cos(deltaB);

	double sinXi = sin(xi), cosXi = cos(xi);

	//std::cout << "e  = " << e << std::endl;
	//std::cout << "xi = " << xi << std::endl;

	v3D s = db.node[nodes_i[2]] - db.node[nodes_i[1]];

	double nu = angle(db.node[nodes_i[0]] - db.node[nodes_i[1]], db.node[nodes_i[2]] - db.node[nodes_i[1]]);
	double mu = angle(db.node[nodes_i[1]] - db.node[nodes_i[2]], db.node[nodes_i[0]] - db.node[nodes_i[2]]);
	double kappa = angle(db.node[nodes_i[1]] - db.node[nodes_i[0]], db.node[nodes_i[2]] - db.node[nodes_i[0]]);

	double sinMu = sin(mu), cosMu = cos(mu);
	double sinNu = sin(nu), cosNu = cos(nu);
	double logSinNu = log(sinNu);
	double logSinMu = log(sinMu);

	double sinKappa = sin(kappa);

	double psi = atan2((e ^ s) & nrmi, e & s);
	double sinPsi = sin(psi), cosPsi = cos(psi);

	double sinNuPsi = sin(nu + psi), cosNuPsi = cos(nu + psi);
	double sinMuPsi = sin(mu - psi), cosMuPsi = cos(mu - psi);

	double cosSigmaA = sinDeltaA * sinPsi * cosXi + cosDeltaA * cosPsi;
	double cosChiA = sinDeltaA * cosPsi * cosXi - cosDeltaA * sinPsi;
	double cosEtaA = cosDeltaA * sinPsi * cosXi - sinDeltaA * cosPsi;
	double cosThetaA = sinDeltaA * sinNuPsi * cosXi + cosDeltaA * cosNuPsi;
	double cosLambdaA = sinDeltaA * sinMuPsi * cosXi - cosDeltaA * cosMuPsi;

	double cosSigmaB = sinDeltaB * sinPsi * cosXi + cosDeltaB * cosPsi;
	double cosChiB = sinDeltaB * cosPsi * cosXi - cosDeltaB * sinPsi;
	double cosEtaB = cosDeltaB * sinPsi * cosXi - sinDeltaB * cosPsi;
	double cosThetaB = sinDeltaB * sinNuPsi * cosXi + cosDeltaB * cosNuPsi;
	double cosLambdaB = sinDeltaB * sinMuPsi * cosXi - cosDeltaB * cosMuPsi;

	auto sqr = [](double x) {return x * x; };
	auto sign = [](double x) {if (!x) return 0.0; return (x>0 ? 1.0 : -1.0); };
	auto arg = [](double x) {return (x > 0 ? 0.0 : M_PI); };

	auto q_ThetaPsiCont = [ii, jj, cosXi, sinXi, mu, sinMu, cosMu, logSinMu, psi, sinPsi, cosPsi, nu, sinNu, cosNu, logSinNu, kappa, sinKappa, sinNuPsi, sinMuPsi, sqr, sign, arg](double delta, double cosLambda, double cosTheta, double cosEta, double cosSigma, double cosChi)
	{
		double logOneCosTheta = log(1.0 + cosTheta);
		double logOneCosLambda = log(1.0 + cosLambda);

		double Lambda1 = logOneCosLambda - logOneCosTheta + logSinNu - logSinMu;
		double Lambda2 = log(tan(nu / 2) * tan(mu / 2));

		double sinDelta2 = sin(delta / 2), cosDelta2 = cos(delta / 2);
		double sinDelta = 2.0 * sinDelta2 * cosDelta2, cosDelta = cosDelta2 * cosDelta2 - sinDelta2 * sinDelta2;
				
		//double Amu = atan2(sinDelta2 * cos((mu - psi) / 2) * sinXi, sinDelta2 * cos((mu - psi) / 2) * cosXi + cosDelta2 * sin((mu - psi) / 2));
		double Amu = atan2(sinDelta2 / cosDelta2 * cos((mu - psi) / 2) * sinXi, sinDelta2 / cosDelta2 * cos((mu - psi) / 2) * cosXi + sin((mu - psi) / 2));
		//double Amuv = atan2(sinXi * sinDelta2/cosDelta2, cosXi * sinDelta2/cosDelta2 + sin((mu - psi) / 2)/cos((mu - psi) / 2));
		
		//double Anu = atan2(sinDelta2 * sin((nu + psi) / 2) * sinXi, sinDelta2 * sin((nu + psi) / 2) * cosXi + cosDelta2 * cos((nu + psi) / 2));
		double Anu = atan2(sinDelta2 / cosDelta2 * sin((nu + psi) / 2) * sinXi, sinDelta2 / cosDelta2 * sin((nu + psi) / 2) * cosXi + cos((nu + psi) / 2));
		//double Anuv = atan2(sinXi * sinDelta2/cosDelta2,  cosXi * sinDelta2/cosDelta2 + cos((nu + psi) / 2) / sin((nu + psi) / 2));
		
		double W = atan2(sinDelta * sin(kappa / 2) * sinXi, \
			cos(kappa / 2) + sin(delta) * cos((mu-psi)/2.0 - (nu+psi)/2.0) *cosXi + cos(delta) * sin((mu - psi) / 2.0 - (nu + psi) / 2.0));
				
		double D = 1.0 / (sqr(sin(delta - psi)) + sinDelta * sinPsi * (1.0 - cosXi) * (cos(delta-psi) + cosSigma));

		//double G = cosDelta * (cosPsi * cosChi - sinPsi * cosSigma) + sinPsi * cosPsi * (sqr(cosSigma) + sqr(cosChi)) + sinDelta * cosEta / sinPsi;
		double G = cosPsi * (sinDelta * cosSigma * cosXi + cosDelta * cosChi - sqr(sinDelta) / sinPsi);
		

		double gent = 2.0 * (Anu * sinMu * sinNuPsi - Amu * sinNu * sinMuPsi - \
			D * sinMu * sinNu * sinDelta * (W * cosEta + 0.5 * sinPsi * sinXi * (Lambda1 - Lambda2 * cosSigma))) / (sinPsi * sinKappa);

		double gens =
			0.5 * (3.0 - log(2.0)) + \
			(sinMu*sinNu/sinKappa)*((logOneCosLambda-logOneCosTheta)*cosPsi/sinPsi + D*(Lambda1 * sinDelta * cosEta / sinPsi + \
				Lambda2 * cosChi - 2.0 * W * sinDelta * sinXi - G * (logSinNu - logSinMu))) - \
			(cosMu*sinNu*(logOneCosLambda-logSinMu) + sinMu*cosNu*(logOneCosTheta-logSinNu))/sinKappa - \
			0.5*(logSinMu+logSinNu-log(sinKappa));

		double mulPsi   = sign(sinXi * sinPsi);
		double mulDelta = sign(sinXi * sinDelta);


		if ((fabs(sinXi) < eps_zero) && (1.0 - fabs(cosSigma) < 0.5 * eps_zero2) && (fabs(sin(psi)) > eps_zero))
		{
			double ara = arg(sin(0.5 * (nu + psi)) * cosSigma);
			double arb = arg(cos(0.5 * (mu - psi)) * cosSigma);

			return v2D{ //tests 1,2,3,4			
				2.0 * mulPsi * cosSigma * cosXi / (sinKappa * sinPsi) * (sinMu * sin(nu + psi) * ara - sinNu * sin(mu - psi) * arb) ,
				0.5 * (1.0 - log(2.0)) - 0.5 * log((1.0 - cosSigma*cosMu) * (1.0 + cosSigma*cosNu) / sinKappa) + \
					cosSigma * (sinMu - sinNu + 0.5 * sin(mu - nu) * Lambda2) / sinKappa
			}; 
		}
						
		if ((fabs(sinXi) < eps_zero) && (fabs(sin(psi)) > eps_zero))
		{ // tests 1,2,3,4,5,6			
			double are = arg(sin(0.5 * (nu + mu)) + sin(0.5 * (mu - nu) - psi + delta * cosXi));
			double arc = arg(1.0 / tan(0.5 * (nu + psi)) + tan(0.5 * delta) * cosXi);
			double ard = arg(      tan(0.5 * (mu - psi)) + tan(0.5 * delta) * cosXi);		

			return v2D{ 
				2.0 * mulDelta * (sin(delta) * sinMu * sinNu / cosChi * cosXi * are + \
				                                     sinMu * sin(nu + psi) * arc - sinNu * sin(mu - psi) * ard ) / (sinKappa * sinPsi),
				gens 
			};
		}	

		if (fabs(sin(psi)) < eps_zero)
		{
			if (fabs(delta) > eps_zero)
				return v2D{ // tests 5,6,7,8
					2.0 * (sinMu * sinNu * (W * cosDelta * cosXi - 0.5 * (Lambda1 - Lambda2 * cosSigma) * sinXi) / sinDelta + \
						(Anu + mulDelta*arg(sin(0.5*(nu+psi)))) * sinMu * cosNu + \
					       (Amu + mulDelta*arg(cos(0.5*(mu-psi)))) * cosMu * sinNu ) / sinKappa,
					0.5 * (3.0 - log(2.0)) - 0.5 * (log((1.0 + cosLambda) * (1.0 + cosTheta) / sinKappa) - sin(mu - nu)  / sinKappa * Lambda1) - \
						sinMu * sinNu * ((Lambda1 * cosDelta - Lambda2 * cos(psi)) * cosXi + 2.0 * W * sinXi) / (sinKappa * sinDelta)
				};				
			else
				return v2D{ // tests 5,6,7,8
					2.0*arg(cosPsi),
					0.5*(1.0-log(2.0))-0.5*log((1.0-cosPsi*cosMu)*(1.0+cosPsi*cosNu)/sinKappa) + (0.5*sin(mu-nu)*Lambda1 + cosPsi*(sinMu-sinNu))/sinKappa
			};				
		}

		return v2D{
			gent, gens
		};
	};
	
	const v2D& qTP_A = q_ThetaPsiCont(deltaA, cosLambdaA, cosThetaA, cosEtaA, cosSigmaA, cosChiA);
	const v2D& qTP_B = q_ThetaPsiCont(deltaB, cosLambdaB, cosThetaB, cosEtaB, cosSigmaB, cosChiB);
		
	double IntThetaSing = db.measure[ii[0]] * (qTP_A[0] - qTP_B[0]);
	v3D IntPsiSing = db.measure[ii[0]] * (qTP_A[1] * taua + qTP_B[1] * taub);

	return { (fabs(xi) < 1e-5) ? 0.0 : IntThetaSing, IntPsiSing};
}

p13D CompJ3DK::IntRegContact(const i2D& ii, const i2D& jj, size_t refineLevel)
{
	auto fo = [this, ii, jj](const v3D& r) { return RegContact(r, ii, jj); };
	return gp->integrate<p13D>(fo, ii[0], refineLevel);
}

std::pair<p13D, int> CompJ3DK::IntRegContactEpsRel(const i2D& ii, const i2D& jj)
{
	auto fo = [this, ii, jj](const v3D& r) { return RegContact(r, ii, jj); };
	return gp->integrateEpsRel<p13D>(fo, ii[0], epsRel);
}

/////////////////////////////////////////////////////////////////////////////////////////////

//*********************************************************************************************************************

v3D CompJ3DK::evaluate(int i, int j)
{
	if (i==j)
	//if (&(db.topo[i]) == &(db.topo[j]))	
	{
		return { 0.0, 0.0, 0.0 };
	}

	if (db.ifSosed({ i, j }))
	{
		//std::cout << "YES!" << std::endl;

		//std::cout << "topo_i = " << db.topo[i] << std::endl;
		//std::cout << "topo_j = " << db.topo[j] << std::endl;

		i2D shifts = RenumerationSosed(i, j);

		//std::cout << "topo_i_rot = " << db.topo[i].rotateLeft(shifts[0]) << std::endl;
		//std::cout << "topo_j_rot = " << db.topo[j].rotateLeft(shifts[1]) << std::endl;

		//auto funcReg = [&](const v3D& r) { return JRegSosed(r, { j, shifts[1] }); };
		//auto IntegralReg = gp->integrate<v3D>(funcReg, i);

#ifdef AUTOSPLIT		
		auto [IntegralReg, nrefine] = IntJRegSosedEpsRel(i, { j, shifts[1] });
		//refines.push_back(nrefine);
#else
		v3D IntegralReg = IntJRegSosed(i, { j, shifts[1] }, nrefine);
		//refines.push_back(nrefine);
#endif


		v3D IntegralSing = IntJSingSosed({ i, shifts[0] }, { j, shifts[1] });

		v3D res = IntegralReg + IntegralSing;

		//cout << "AA = " << db.node[db.topo[i].rotateLeft(shifts[0])[2]] << endl;
		//cout << "BB = " << db.node[db.topo[i].rotateLeft(shifts[0])[1]] << endl;
		//cout << "MM = " << db.node[db.topo[i].rotateLeft(shifts[0])[0]] << endl;
		//cout << "CC = " << db.node[db.topo[j].rotateLeft(shifts[1])[0]] << endl;

		//auto [thSingSosed, psSingSosed] = SingSosed({ 0.1, -0.2, 0.0 }, { j, shifts[1] });
	    //cout << "sing = " << thSingSosed << " " << psSingSosed << endl;


		//std::cout << "res = " << res << std::endl;

		return IntegralReg + IntegralSing;
	}

	
	if (db.ifContact({ i, j }))
	{
		//std::cout << "CONTACT!" << std::endl;

		//std::cout << "i = " << i << " , j = " << j << std::endl;

		//std::cout << "topo_i = " << db.topo[i] << std::endl;
		//std::cout << "topo_j = " << db.topo[j] << std::endl;

		i2D shifts = RenumerationContact(i, j);

		//std::cout << "topo_i_rot = " << db.topo[i].rotateLeft(shifts[0]) << std::endl;
		//std::cout << "topo_j_rot = " << db.topo[j].rotateLeft(shifts[1]) << std::endl;
		//std::cout << std::endl;
		//std::cout << std::endl;

#ifdef AUTOSPLIT		
		auto[IntegralReg, nrefine] = IntRegContactEpsRel({ i, shifts[1] }, { j, shifts[1] });
		//refines.push_back(nrefine);
#else
		p13D IntegralReg = IntRegContact({ i, shifts[1] }, { j, shifts[1] }, nrefine);
		//refines.push_back(nrefine);
#endif	
		
		auto IntegralSing = IntSingContact({ i, shifts[0] }, { j, shifts[1] });

		double intTheta = IntegralReg.first + IntegralSing.first;

		

		int p = 0;		
		double refTheta = 2.0 * M_PI * db.measure[i];
		if (intTheta > refTheta)
			p = -((int)trunc((intTheta - refTheta) / (2.0 * refTheta)) + 1);
		else if (intTheta < -refTheta)
			p = ((int)trunc((-refTheta - intTheta) / (2.0 * refTheta)) + 1);
		
		//std::cout << "intTheta = " << std::setprecision(10) << intTheta + 2.0 * p * refTheta << std::endl;

		//if (p) 
		//	std::cout << "p = " << p << std::endl;

		v3D res = 0.25/M_PI *( (IntegralReg.first + IntegralSing.first + 2.0 * p * refTheta) * db.nrm[j]
			+ ((IntegralReg.second + IntegralSing.second) ^ db.nrm[j]));

		//cout << "AA = " << db.node[db.topo[i].rotateLeft(shifts[0])[2]] << endl;
		//cout << "BB = " << db.node[db.topo[i].rotateLeft(shifts[0])[1]] << endl;
		//cout << "MM = " << db.node[db.topo[i].rotateLeft(shifts[0])[0]] << endl;
		//cout << "CC = " << db.node[db.topo[j].rotateLeft(shifts[1])[0]] << endl;

		//auto [thSingSosed, psSingSosed] = SingSosed({ 0.1, -0.2, 0.0 }, { j, shifts[1] });
	    //cout << "sing = " << thSingSosed << " " << psSingSosed << endl;

		//std::cout << "res = " << res << std::endl;	
		
		return res;

	}
	

	auto func = [&](const v3D& r) { return J3D(r, j); };
		

#ifdef AUTOSPLIT		
	auto[res, nrefine] = gp->integrateEpsRel<v3D>(func, i, epsRel);
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

