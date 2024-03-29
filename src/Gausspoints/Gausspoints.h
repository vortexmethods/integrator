// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c описанием шаблонного класса Gausspoints
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 11 сентября 2022 г.
\version 0.5
*/

#pragma once

#include <functional> 
#include <algorithm>
#include <iostream>
#include <typeinfo>

#include "Database.h"

const int maxRefine = 7;

#define DIVZERO 2e-6

/*!
\brief Структура -- положения и веса гауссовых точек
\n Формулы для интегрирования по треугольнику взяты из статьи
\n Cowper, G. R. (1973). Gaussian quadrature formulas for triangles. 
\n International Journal for Numerical Methods in Engineering, 7(3), 405–408. doi:10.1002/nme.1620070316

\tparam dim размерность задачи (dim=3 - интегрирование по треугольнику, dim = 2 - интегрирование по отрезку)

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/
template<size_t dim>
struct gpPosWeightOrd
{
	/// L-координаты гауссовых точек (последняя не задается, получается дополнением до единицы)
	const std::vector<numvector<double, dim-1>> pos;

	/// Веса гауссовых точек
	const std::vector<double> weight;

	/// Порядок формулы
	const int ord;
};

/// Квадратурная формула с 1-й гауссовой точкой (2-й порядок точности)
const gpPosWeightOrd<2> gp2D1{
{
	{ 0.5 },
},
{
	1.0
},
	2
};


/// Квадратурная формула с 2-мя гауссовыми точками (4-й порядок точности)
const gpPosWeightOrd<2> gp2D2{
{
	{0.211324865405187}, 
	{0.788675134594813},
},
{
	0.500000000000000, 0.500000000000000
},
	4
};

/// Квадратурная формула с 4-мя гауссовыми точками (8-й порядок точности)
const gpPosWeightOrd<2> gp2D4{
{
	{0.069431844202974}, 
	{0.33000947820757},
	{0.66999052179243},
	{0.93056815579703}
},
{
	0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
},
	8
};

/// Квадратурная формула с 6-мя гауссовыми точками (12-й порядок точности)
const gpPosWeightOrd<2> gp2D6{
{
	{0.033765242898424},
	{0.16939530676687},
	{0.38069040695840},
	{0.61930959304160},
	{0.83060469323313},
	{0.96623475710158}
},
{
	0.085662246189585, 0.180380786524069, 0.233956967286346, \
0.233956967286346, 0.180380786524069, 0.085662246189585
},
	12
};

/// Квадратурная формула с 10-мя гауссовыми точками (20-й порядок точности)
const gpPosWeightOrd<2> gp2D10{
{
	{0.013046735741414},
	{0.06746831665551},
	{0.16029521585049},
	{0.28330230293538},
	{0.42556283050918},
	{0.5744371694908},
	{0.7166976970646},
	{0.8397047841495},
	{0.9325316833445}, 
	{0.9869532642586}
},
{
	0.0333356721543441, 0.074725674575290, 0.109543181257991, \
0.134633359654998, 0.147762112357376, 0.147762112357376, \
0.134633359654998, 0.109543181257991, 0.074725674575290, \
0.0333356721543441
},
	20
};

/// Квадратурная формула с 8-мя гауссовыми точками (16-й порядок точности)
const gpPosWeightOrd<2> gp2D8{
{
	{0.019855071751232},
	{0.10166676129319},
	{0.23723379504184},
	{0.40828267875218},
	{0.59171732124782}, 
	{0.7627662049582}, 
	{0.8983332387068},
	{0.9801449282488}
},
{
	0.0506142681451881, 0.111190517226687, 0.156853322938944, \
	0.181341891689181, 0.181341891689181, 0.156853322938944, \
	0.111190517226687, 0.0506142681451881
},
	16
};

/// Квадратурная формула с 12-мя гауссовыми точками (24-й порядок точности)
const gpPosWeightOrd<2> gp2D12{
{
	{0.009219682876640},
	{0.04794137181476},
	{0.11504866290285},
	{0.20634102285669},
	{0.3160842505009},
	{0.4373832957443},
	{0.5626167042557},
	{0.6839157494991},
	{0.7936589771433},
	{0.8849513370972},
	{0.9520586281852},
	{0.9907803171234}
},
{
	0.0235876681932559, 0.0534696629976592, 0.080039164271673, \
	0.101583713361533, 0.116746268269177, 0.124573522906701, \
	0.124573522906701, 0.116746268269177, 0.101583713361533, \
	0.080039164271673, 0.0534696629976592, 0.0235876681932559
},
	24
};

/// Квадратурная формула с 1-й гауссовой точкой (2-й порядок точности)
const gpPosWeightOrd<3> gp3D1{
{
	{ 0.3333333333333333, 0.3333333333333333 },	
},
{
	1.0	
},
	2
};

/// Квадратурная формула с 3-мя гауссовыми точками (2-й порядок точности)
const gpPosWeightOrd<3> gp3D3{
{
	{ 0.1666666666666667, 0.1666666666666667 },
	{ 0.6666666666666667, 0.1666666666666667 },
	{ 0.1666666666666667, 0.6666666666666667 },
},
{
	0.3333333333333333,
	0.3333333333333333,
	0.3333333333333333
},
	2
};

/// Квадратурная формула с 4-мя гауссовыми точками (3-й порядок точности)
const gpPosWeightOrd<3> gp3D4{
	//Первые две L-координаты
{
	{ 0.3333333333333333, 0.3333333333333333 },
	{ 0.6, 0.2 },
	{ 0.2, 0.6 },
	{ 0.2, 0.2 }
},
//Веса 
{
	-0.5625,
	0.5208333333333333,
	0.5208333333333333,
	0.5208333333333333
},
	3
};

/// Квадратурная формула с 6-ю гауссовыми точками (4-й порядок точности)
const gpPosWeightOrd<3> gp3D6{
{
	{ 0.816847572980459, 0.091576213509771 },
	{ 0.091576213509771, 0.816847572980459 },
	{ 0.091576213509771, 0.091576213509771 },
	{ 0.445948490915965, 0.108103018168070 },
	{ 0.108103018168070, 0.445948490915965 },
	{ 0.445948490915965, 0.445948490915965 }	
},
{
	0.109951743655322,
	0.109951743655322,
	0.109951743655322,
	0.223381589678011,
	0.223381589678011,
	0.223381589678011
},
	4
};

/// Квадратурная формула с 7-ю гауссовыми точками (5-й порядок точности)
const gpPosWeightOrd<3> gp3D7{
{
	{ 0.3333333333333333, 0.3333333333333333 },
	{ 0.101286507323456, 0.797426985353087 },
	{ 0.797426985353087, 0.101286507323456 },
	{ 0.101286507323456, 0.101286507323456 },
	{ 0.470142064105115, 0.059715871789770 },
	{ 0.059715871789770, 0.470142064105115 },
	{ 0.470142064105115, 0.470142064105115 }
},
{
	0.225,
	0.125939180544827,
	0.125939180544827,
	0.125939180544827,
	0.132394152788506,
	0.132394152788506,
	0.132394152788506
},
	5
};

/// Квадратурная формула с 9-ю гауссовыми точками (5-й порядок точности)
const gpPosWeightOrd<3> gp3D9{
{
	{ 0.437525248383384, 0.124949503233232 },
	{ 0.124949503233232, 0.437525248383384 },
	{ 0.437525248383384, 0.437525248383384 },
	{ 0.797112651860071, 0.165409927389841 },
	{ 0.165409927389841, 0.797112651860071 },
	{ 0.797112651860071, 0.037477420750088 },
	{ 0.037477420750088, 0.797112651860071 },
	{ 0.165409927389841, 0.037477420750088 },
	{ 0.037477420750088, 0.165409927389841 }
},
{
	0.205950504760887,
	0.205950504760887,
	0.205950504760887,
	0.063691414286223,
	0.063691414286223,
	0.063691414286223,
	0.063691414286223,
	0.063691414286223,
	0.063691414286223,	
},
	5
};

/// Квадратурная формула с 12-ю гауссовыми точками (6-й порядок точности)
const gpPosWeightOrd<3> gp3D12{
{
	{ 0.873821971016996, 0.063089014491502 },
	{ 0.063089014491502, 0.873821971016996 },
	{ 0.063089014491502, 0.063089014491502 },
	{ 0.501426509658179, 0.249286745170910 },
	{ 0.249286745170910, 0.501426509658179 },
	{ 0.249286745170910, 0.249286745170910 },
	{ 0.636502499121399, 0.310352451033785 },
	{ 0.310352451033785, 0.636502499121399 },
	{ 0.636502499121399, 0.053145049844816 },
	{ 0.053145049844816, 0.636502499121399 },
	{ 0.310352451033785, 0.053145049844816 }, 
	{ 0.053145049844816, 0.310352451033785 }	
},
{
	0.050844906370207,
	0.050844906370207,
	0.050844906370207,
	0.116786275726379,
	0.116786275726379,
	0.116786275726379,
	0.082851075618374,
	0.082851075618374,
	0.082851075618374,
	0.082851075618374,
	0.082851075618374,
	0.082851075618374	
},
	6
};

/// Квадратурная формула с 13-ю гауссовыми точками (7-й порядок точности)
const gpPosWeightOrd<3> gp3D13{
{
	{ 0.333333333333333, 0.333333333333333 },
	{ 0.479308067841923, 0.260345966079038 },
	{ 0.260345966079038, 0.479308067841923 },
	{ 0.260345966079038, 0.260345966079038 },
	{ 0.869739794195598, 0.065130102902216 },
	{ 0.065130102902216, 0.869739794195598 },
	{ 0.065130102902216, 0.065130102902216 },
	{ 0.638444188569809, 0.312865496004875 },
	{ 0.312865496004875, 0.638444188569809 },
	{ 0.638444188569809, 0.048690315425316 },
	{ 0.048690315425316, 0.638444188569809 },
	{ 0.312865496004875, 0.048690315425316 },
	{ 0.048690315425316, 0.312865496004875 }
}, 
{
   -0.149570044467670,
	0.175615257433204,
	0.175615257433204,
	0.175615257433204,
	0.053347235608839,
	0.053347235608839,
	0.053347235608839,
	0.077113760890257,
	0.077113760890257,
	0.077113760890257,
	0.077113760890257,
	0.077113760890257,
	0.077113760890257
},
	7
};
	

/*!
\brief Класс -- интегратор по гауссовым точкам по симплексу

\tparam dim размерность задачи (dim=3 - интегрирование по треугольнику, dim = 2 - интегрирование по отрезку)

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.5
\date 11 сентября 2022 г.
*/

template<size_t dim>
class Gausspoints
{
private:
	/// Количество гауссовых точек
	size_t ngp;

	/// L-координаты, вычисленные полностью (dim-мерный вектор для каждой гауссовой точки)
	std::vector<numvector<double, dim>> Lcoord;
	
	/// Константная ссылка на базу данных геометрии тела
	const Database<dim>& db;

	/// Константная ссылка на положения и веса гауссовых точек
	const gpPosWeightOrd<dim>& pw;

public:
	/// \brief Конструктор интегратора по гауссовым точкам
	///
	/// \param[in] db_ константная ссылка на базу данных геометрии тела
	/// \param[in] pw_ константная ссылка на структуру с положениями и весами гауссовых точек
	Gausspoints(const Database<dim>& db_, const gpPosWeightOrd<dim>& pw_);
	
	/// Пустой деструктор
	~Gausspoints() {};
		
	/// \brief Шаблонная функция вычисления интеграла по гауссовым точкам
	///
	/// \tparam T тип результата интегрируемой функции
	/// 
	/// \param[in] fo интегрируемая функция
	/// \param[in] pnl номер панели, по которой производится интегрирование
	/// \param[in] refineLevel количество дополнительных итераций разбиения панели
	template<typename T>
	T integrate(std::function<T(const numvector<double, dim>&)> fo, size_t pnl, size_t refineLevel) const;
	
	template<typename T>
	std::pair<T, int> integrateEpsRel(std::function<T(const numvector<double, dim>&)> fo, size_t pnl, double epsRel) const;

	using simplex = std::pair<numvector<numvector<double, dim>, dim>, double>;


	void refine(std::vector<simplex>& lst) const;
	
}; //Gausspoints

template<> 
void Gausspoints<2>::refine(std::vector<simplex>& lst) const;

template<>
void Gausspoints<3>::refine(std::vector<simplex>& lst) const;




template<size_t dim>
Gausspoints<dim>::Gausspoints(const Database<dim>& db_, const gpPosWeightOrd<dim>& pw_)
	: ngp(pw_.weight.size()), db(db_), pw(pw_)
{
	Lcoord.resize(ngp);

	for (size_t i = 0; i < ngp; ++i)
	{
		for (size_t s = 0; s < dim - 1; ++s)
			Lcoord[i][s] = pw.pos[i][s];

		double lastLcoord = 1.0;
		for (size_t s = 0; s < dim - 1; ++s)
			lastLcoord -= Lcoord[i][s];
		Lcoord[i][dim - 1] = lastLcoord;
	}
};//Gausspoints(...)


template<size_t dim> template<typename T>
T Gausspoints<dim>::integrate(std::function<T(const numvector<double, dim>&)> fo, size_t pnl, size_t refineLevel) const
{	
	const auto& topo = db.topo[pnl];

	simplex nds;
	for (size_t q = 0; q < dim; ++q)
		nds.first[q] = db.node[topo[q]];
	nds.second = db.measure[pnl];

	std::vector trgs { nds };
	for (int r = 0; r < refineLevel; ++r)
		refine(trgs);

	numvector<double, dim> pt;
	T res(0.0);

	for (auto& tr : trgs)
	{
		for (size_t i = 0; i < ngp; ++i)
		{
			pt.toZero();
			for (size_t s = 0; s < dim; ++s)
				pt += tr.first[s] * Lcoord[i][s];

			res += (tr.second * pw.weight[i]) * fo(pt);
			//std::cout << "pt = " << pt << ",  f(pt) = " << fo(pt) << std::endl;
		}//for i
	}
	return res;
}//integrate(...)

inline double newdiv(double x, double y)
{
	return ((fabs(x)< DIVZERO) && (fabs(y)< DIVZERO)) ? 0.0 : x/y;
}

inline v3D newdiv(const v3D& x, const v3D& y)
{
	return v3D{ 
		((fabs(x[0]) < DIVZERO) && (fabs(y[0]) < DIVZERO)) ? 0.0 : x[0] / y[0],
		((fabs(x[1]) < DIVZERO) && (fabs(y[1]) < DIVZERO)) ? 0.0 : x[1] / y[1],
		((fabs(x[2]) < DIVZERO) && (fabs(y[2]) < DIVZERO)) ? 0.0 : x[2] / y[2] };
}

inline p13D newdiv(const p13D& x, const p13D& y)
{
	return p13D{
		((fabs(x.first) < 1e-14) && (fabs(y.first) < 1e-14)) ? 0.0 : x.first / y.first,

	   {((fabs(x.second[0]) < DIVZERO) && (fabs(y.second[0]) < DIVZERO)) ? 0.0 : x.second[0] / y.second[0],
		((fabs(x.second[1]) < DIVZERO) && (fabs(y.second[1]) < DIVZERO)) ? 0.0 : x.second[1] / y.second[1],
		((fabs(x.second[2]) < DIVZERO) && (fabs(y.second[2]) < DIVZERO)) ? 0.0 : x.second[2] / y.second[2]} };
}



inline double fabs(const std::pair<double, double>& x)
{
	return fabs(x.first) + fabs(x.second);
};

inline double fabs(const p13D& x)
{
	return fabs(x.first) + fabs(x.second);
};

inline std::ostream& operator<<(std::ostream& str, const std::pair<double, double>& pr)
{
	str << "{ " << pr.first << ", " << pr.second << " }";
	return str;
}

template<size_t dim> template<typename T>
std::pair<T, int> Gausspoints<dim>::integrateEpsRel(std::function<T(const numvector<double, dim>&)> fo, size_t pnl, double epsRel) const
{
	int refine = -1;
	double deltah2;
	T Ih, Ih2;
	int p = pw.ord;
	do
	{		
		++refine;
		Ih = integrate(fo, pnl, refine);
		Ih2 = integrate(fo, pnl, refine + 1);

		auto num = (Ih - Ih2);
		auto den = (Ih2*(1 << (p)) - Ih);
		auto div = newdiv(num, den);
		deltah2 = fabs(div);

		//std::cout << "num = " << num << ", den = " << den << ", div = " << div << std::endl;
		//std::cout << "Ih = " << Ih << ", Ih2 = " << Ih2 << ", deltah2 = " << deltah2 << std::endl;

		/*
		auto num = fabs(Ih - Ih2);
		auto den = fabs(Ih2*(1 << (p + 1)) - Ih);

		//std::cout << typeid(num).name() << std::endl;
		//std::cout << typeid(den).name() << std::endl;
		auto div = newdiv(num, den);
		deltah2 = fabs(div);
		std::cout << "num = " << num << ", den = " << den << ", div = " << div << std::endl;
		std::cout << "Ih = " << Ih << ", Ih2 = " << Ih2 << ", deltah2 = " << deltah2 << std::endl;
		*/
	} while ((deltah2 > epsRel) && (refine < maxRefine));
	
	if (refine == maxRefine)
	{
		std::cout << "Doesn't converge!" << std::endl;
	}

	return { Ih2, refine };
}






/*!
\brief Класс -- интегратор по гауссовым точкам по квадрату

\tparam dim размерность задачи (dim=3 - интегрирование по треугольнику, dim = 2 - интегрирование по отрезку)

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.3
\date 02 апреля 2022 г.
*/



template<size_t dim>
class GausspointsCube
{
private:
	/// Константная ссылка на положения и веса гауссовых точек
	const gpPosWeightOrd<2>& pw;

public:
	/// \brief Конструктор интегратора по гауссовым точкам
	///	
	/// \param[in] pw_ константная ссылка на структуру с положениями и весами гауссовых точек
	GausspointsCube(const gpPosWeightOrd<2>& pw_) : pw(pw_) {};

	/// Пустой деструктор
	~GausspointsCube() {};

	using cube = std::pair<numvector<double, dim>, double>;

	/// \brief Шаблонная функция вычисления интеграла по гауссовым точкам
	///
	/// \tparam T тип результата интегрируемой функции
	/// 
	/// \param[in] fo интегрируемая функция
	/// \param[in] refineLevel количество дополнительных итераций разбиения панели
	template<typename T>
	T integrate(std::function<T(const numvector<double, dim>&)> fo, size_t refineLevel) const;
	
	template<typename T>
	std::pair<T, int> integrateEpsRel(std::function<T(const numvector<double, dim>&)> fo, double epsRel) const;
	
	void refineCube(std::vector<cube>& lst) const;

}; //GausspointsCube

template<>
void GausspointsCube<2>::refineCube(std::vector<cube>& lst) const;

template<>
void GausspointsCube<3>::refineCube(std::vector<cube>& lst) const;


template<size_t dim> template<typename T>
T GausspointsCube<dim>::integrate(std::function<T(const numvector<double, dim>&)> fo, size_t refineLevel) const
{


	numvector<double, dim> zeroVec;
	zeroVec.toZero();

	std::vector<std::pair<numvector<double, dim>, double>> cubs{ { zeroVec, 1.0 } };
	for (int r = 0; r < refineLevel; ++r)
		refineCube(cubs);


	numvector<T, dim> resijk;
	numvector<double, dim> posijk;
	numvector<double, dim> weightijk;

	T result(0);

	double h;

	if (dim == 2)
	{
		for (auto& cb : cubs)
		{
			h = cb.second;
			const numvector<double, dim>& a = cb.first;
			

			resijk[0] = 0.0;
			for (size_t i = 0; i < pw.pos.size(); ++i)
			{
				posijk[0] = a[0] + h * pw.pos[i][0];
				weightijk[0] = pw.weight[i];

				resijk[1] = 0.0;
				for (size_t j = 0; j < pw.pos.size(); ++j)
				{
					posijk[1] = a[1] + h * pw.pos[j][0];
					weightijk[1] = pw.weight[j];

					resijk[1] += weightijk[1] * fo(posijk);
				}//for j

				resijk[0] += weightijk[0] * resijk[1];
			}//for i

			result += resijk[0];
		}

		result *= h * h;
	}

	if (dim == 3)
	{
		for (auto& cb : cubs)
		{
			h = cb.second;
			const numvector<double, dim>& a = cb.first;
			
			resijk[0] = 0.0;
			for (size_t i = 0; i < pw.pos.size(); ++i)
			{
				posijk[0] = a[0] + h * pw.pos[i][0];
				weightijk[0] = pw.weight[i];

				resijk[1] = 0.0;
				for (size_t j = 0; j < pw.pos.size(); ++j)
				{
					posijk[1] = a[1] + h * pw.pos[j][0];
					weightijk[1] = pw.weight[j];

					resijk[2] = 0.0;
					for (size_t k = 0; k < pw.pos.size(); ++k)
					{
						posijk[2] = a[2] + h * pw.pos[k][0];
						weightijk[2] = pw.weight[k];

						resijk[2] += weightijk[2] * fo(posijk);
					}//for k

					resijk[1] += weightijk[1] * resijk[2];
				}//for j

				resijk[0] += weightijk[0] * resijk[1];
			}//for i

			result += resijk[0];
		}
		result *= h * h * h;
	}

	return result;
}//integrate(...)


template<size_t dim> template<typename T>
std::pair<T, int> GausspointsCube<dim>::integrateEpsRel(std::function<T(const numvector<double, dim>&)> fo, double epsRel) const
{
	int refine = -1;
	double deltah2;
	T Ih, Ih2;
	int p = pw.ord;
	do
	{
		++refine;
		Ih = integrate(fo, refine);
		Ih2 = integrate(fo,  refine + 1);

		auto num = (Ih - Ih2);
		auto den = (Ih2 * (1 << (p)) - Ih);
		auto div = newdiv(num, den);
		deltah2 = fabs(div);

		//std::cout << "refine = " << refine << std::endl;
		//std::cout << "num = " << num << ", den = " << den << ", div = " << div << std::endl;
		//std::cout << "Ih = " << Ih << ", Ih2 = " << Ih2 << ", deltah2 = " << deltah2 << std::endl << std::endl;


	} while ((deltah2 > epsRel) && (refine < maxRefine));

	if (refine == maxRefine)
	{
		std::cout << "Doesn't converge!" << std::endl;
	}

	return { Ih2, refine };
}


