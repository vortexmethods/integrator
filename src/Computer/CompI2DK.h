// ������ Integrator 
// (c) �.�. ��������, �.�. ����������, �.�. ����������, 2022

/*!
\file
\brief ���� � ��������� ������ CompI2DK
\author �������� ���� ����������
\author ���������� ���� ��������������
\author ���������� ����� ���������

\date 02 ������ 2022 �.
\version 0.3
*/

#pragma once
#include "Computer.h"

class CompI2DK :
	public Computer<double, 2>
{
public:

	/// \brief �����������
	/// 	
	/// \param[in] db_ ����������� ������ �� ���� ������ �������������� ����������
	/// \param[in] par_ ����������� ������ �� �����, ����������� ������������������ �� MPI
	CompI2DK(const Database<2>& db_, const Parallel& par_);

	/// ����������
	~CompI2DK();

	double Theta(const v2D& v1, const v2D& v2);

	double psi(const v2D& v, const v2D& ti, const v2D& tj);

	//double psi(const v2D& v1, int i, int j);

	//double psi(const v2D& v1, const v2D& t1, const v2D& t2);

	/// \brief ���������� ������� ��� ���������� ������ ���������������� ����������
	///
	/// \param[in] i ������ ����������� ������ � ���� ������
	/// \param[in] j ������ �������� ������ � ���� ������
	/// \return ��������� ��������� --- ���������� ����� �������� �������
	virtual inline double evaluate(int i, int j) override;

};


