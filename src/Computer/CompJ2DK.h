// ������ Integrator 
// (c) �.�. ��������, �.�. ����������, �.�. ����������, 2022

/*!
\file
\brief ���� � ��������� ������ CompJ2DM
\author �������� ���� ����������
\author ���������� ���� ��������������
\author ���������� ����� ���������

\date 02 ������ 2022 �.
\version 0.3
*/

#pragma once
#include "Computer.h"

/*!
\brief ����� -- ����������� ��� 2D ������ ������������ ��������� �� ��������� ������� �����
\n ���������� �� Computer<v2D, 3>

\author �������� ���� ����������
\author ���������� ���� ��������������
\author ���������� ����� ���������
\version 0.3
\date 02 ������ 2022 �.
*/

class CompJ2DK :
	public Computer<v2D, 2>
{
public:

	/// \brief �����������
	/// 	
	/// \param[in] db_ ����������� ������ �� ���� ������ �������������� ����������
	/// \param[in] par_ ����������� ������ �� �����, ����������� ������������������ �� MPI
	CompJ2DK(const Database<2>& db_, const Parallel& par_);

	/// ����������
	~CompJ2DK();

	double Theta(const v2D& v1, const v2D& v2);

	v2D omega(const v2D& v, const v2D& ti, const v2D& tj);

	/// \brief ���������� ������� ��� ���������� ������ ���������������� ����������
	///
	/// \param[in] i ������ ����������� ������ � ���� ������
	/// \param[in] j ������ �������� ������ � ���� ������
	/// \return ��������� ��������� --- ���������� ����� �������� �������
	virtual inline v2D evaluate(int i, int j) override;
};

