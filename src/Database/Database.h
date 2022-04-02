// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Файл c описанием шаблонного класса Database
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#pragma once

#include <map>

#include "numvector.h"

/// Признак соседних панелей
enum nbh_t {
	/// соседняя через грань (3D) или через граничную точку (2D)
	sosed,  //(0)
	/// соседняя через вершину (3D)
	contact //(1)
	/*, self, far*/
}; //enum nbh_t

/*!
\brief Шаблонный класс, определяющий базу данных геометрических объектов --- панелей и точек наблюдения

\tparam dim размерность задачи

\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна
\version 0.3
\date 02 апреля 2022 г.
*/
template <int dim>
class Database
{
public:
	/// Координаты узлов сетки
	std::vector<numvector<double, dim>> node;

	/// Топология сетки
	std::vector<numvector<int, dim>> topo;

	/// Нормали к ячейкам сетки
	std::vector<numvector<double, dim>> nrm;

	/// Центры ячеек сетки
	std::vector<numvector<double, dim>> cnt;
		
	///  Меры ячеек (площади или длины)
	std::vector<double> measure;

	/// Координаты точек вычисления однократных интегралов
	std::vector<numvector<double, dim>> point;

	/// Признаки соседства
	std::map<std::pair<int, int>, nbh_t> nbh;

	bool ifSosed(std::pair<int, int> pr) const;
	bool ifContact(std::pair<int, int> pr) const;

	/// Конструктор (пустой)
	Database() {};

	/// Деструктор (пустой)
	~Database() {};

	/// \brief Считывание базы данных из файла
	/// 
	/// \param[in] fileName имя файла в формате Salome
	/// \warning Базовая реализация отсутствует, см. явные специализации для конкретных dim
	void readNodeTopoFromFile(const std::string& fileName, double scale = 1.0);
		
	/// Заполнение сведений о соседних панелях
	/// \warning Базовая реализация отсутствует, см. явные специализации для конкретных dim
	void fillNbh();

	/// Вычисление нормалей
	/// \warning Базовая реализация отсутствует, см. явные специализации для конкретных dim
	void calcNrm();

	/// Вычисление центров панелей
	/// \warning Базовая реализация отсутствует, см. явные специализации для конкретных dim
	void calcCnt();

	/// Вычисление мер
	/// \warning Базовая реализация отсутствует, см. явные специализации для конкретных dim
	void calcMeasure();

};//class Database

//Специализации функций-членов шаблонного класса для dim=3

template<> void Database<3>::readNodeTopoFromFile(const std::string& fileName, double scale);
template<> void Database<2>::readNodeTopoFromFile(const std::string& fileName, double scale);

template<> void Database<3>::fillNbh();
template<> void Database<2>::fillNbh();

template<> void Database<3>::calcNrm();
template<> void Database<2>::calcNrm();

template<> void Database<3>::calcCnt();
template<> void Database<2>::calcCnt();

template<> void Database<3>::calcMeasure();
template<> void Database<2>::calcMeasure();

template<int dim>
bool Database<dim>::ifSosed(std::pair<int, int> pr) const
{
	auto it = nbh.find(pr);
	if ((it != nbh.end()) && (it->second == sosed))
		return true;
	return false;
}

template<int dim>
bool Database<dim>::ifContact(std::pair<int, int> pr) const
{
	auto it = nbh.find(pr);
	if ((it != nbh.end()) && (it->second == contact))
		return true;
	return false;
}
