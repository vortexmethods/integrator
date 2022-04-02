/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.10   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2021/05/17     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: parallel.h                                                       |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса Parallel и структуры parProp
\author Марчевский Илья Константинович
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef PARALLEL_H
#define PARALLEL_H

#include "mpi.h"
#include <assert.h>

#include "numvector.h"

//namespace VMlib
//{

template <typename T>
/*[[nodiscard]]*/ constexpr MPI_Datatype mpi_get_type() noexcept;


	/*!
	\brief Стрктура, содержащая параметры исполнения задачи в параллельном MPI-режиме
	\author Марчевский Илья Константинович
	\version 1.10
	\date 17 мая 2021 г.
	*/
struct parProp
{
    /// Коммуникатор для решения конкретной задачи
    MPI_Comm comm;

    /// Локальный номер процессора, решающего конкретную задачу
    int rank;

    /// Число процессоров, решающих конкретную задачу
    int np;

    /// Список из чисел витков циклов, предназначенных для каждого процессора
    std::vector<int> len;

    /// Список, определяющий номер витка цикла, с которого должен начинать работу данный процессор
    std::vector<int> disp;

    /// Число витков, предназначенное текущему процессору
    int myLen;

    /// Индекс первого витка из числа витков, предназначенных текущему процессору
    int myDisp;

    /// Общее число витков, разделенное между всеми процессорами
    int totalLen;


    /// \brief Обертка функции MPI_Bcast(...) для рассылки (синхронизации) вектора
    /// 
    /// \tparam T тип данных компонент вектора
    /// \param[in,out] vec ссылка на синхронизируемый вектор 
    /// \param[in] resize признак выполнения resize для вектора (по умолчанию true)   
    template<typename T>
    void BcastVector(std::vector<T>& vec, bool resize = true)
    {
        if ((rank > 0) && (resize))
            vec.resize(totalLen);
        MPI_Bcast(vec.data(), totalLen, mpi_get_type<T>(), 0, comm);
    }//BcastVector(...)

    /// \brief Обертка функции MPI_Scatterv(...) для распределения вектора по процессорам
    /// 
    /// \tparam T тип данных компонент вектора
    /// \param[in] globVec константная ссылка на рассылаемый вектор  
    /// \param[out] locVec ссылка на заполняемые (локальные) векторы
    /// \param[in] locResize признак выполнения resize для локального вектора (по умолчанию true)
    template<typename T>
    void ScattervVector(std::vector<T>& globVec, std::vector<T>& locVec, bool locResize = true)
    {
        if (locResize)
            locVec.resize(myLen);
               
        MPI_Scatterv(const_cast<T*>(globVec.data()), len.data(), disp.data(), mpi_get_type<T>(), locVec.data(), myLen, mpi_get_type<T>(), 0, comm);
    };

    /// \brief Обертка функции MPI_Gatherv(...) для сборки вектора с процессоров
    /// 
    /// \tparam T тип данных компонент вектора
    /// \param[out] locVec константная ссылка на локальные векторы
    /// \param[in] globVec ссылка на глобальный (собираемый) вектор 
    /// \param[in] globResize признак выполнения resize для глобального вектора (по умолчанию true)
    template<typename T>
    void GathervVector(std::vector<T>& locVec, std::vector<T>& globVec, bool globResize = true)
    {
        if (globResize)
            globVec.resize(totalLen);

        MPI_Gatherv(const_cast<T*>(locVec.data()), myLen, mpi_get_type<T>(), globVec.data(), len.data(), disp.data(), mpi_get_type<T>(), 0, comm);
    };
};


	/*!
	\brief Класс, опеделяющий параметры исполнения задачи в параллельном MPI-режиме

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Parallel
	{
	public:
		/// Конструктор по умолчанию
		Parallel() : comm(MPI_COMM_WORLD), rank(0), np(1) {};
		
		/// Коммуникатор для решения конкретной задачи
		MPI_Comm comm;

		/// Локальный номер процессора, решающего конкретную задачу
		int rank;

		/// Число процессоров, решающих конкретную задачу
		int np;

		/// \brief Распределение задач по процессорам
		///
		/// \param[in] n число распределяемых витков цикла
		/// \param[in] bcastAll признак рассылки всей информации всем процессорам (по умолчанию false)
		/// \return структуру типа parProp, заполненную для текущего процессора
		parProp SplitMPIone(size_t n, bool bcastAll = false) const;

		/// \brief Распределение задач по процессорам
		///
		/// \param[in] n число распределяемых витков цикла
		/// \param[in] bcastAll признак рассылки всей информации всем процессорам (по умолчанию false)
		/// \return структуру типа parProp, заполненную для текущего процессора
		parProp SplitMPI(size_t n, bool bcastAll = false) const;

        /// MPI-описатель типа v3D (numvector<double, 3>)
        static MPI_Datatype MPI_V3D;

        /// MPI-описатель типа i3D (numvector<int, 3>)
        static MPI_Datatype MPI_I3D;

        /// MPI-описатель типа v2D (numvector<double, 2>)
        static MPI_Datatype MPI_V2D;

        /// MPI-описатель типа i2D (numvector<int, 2>)
        static MPI_Datatype MPI_I2D;

        /// MPI-описатель типа std::pair<int, int>
        static MPI_Datatype MPI_PAIRII;

        /// Формирование MPI-описателей пользовательских типов
        static void CreateMpiType();

	};


    /// \brief Шаблонная функция автоподбора MPI-описателя типа
/// 
/// \warning Требует добавления пользовательских типов MPI
    template <typename T>
    /*[[nodiscard]]*/ constexpr MPI_Datatype mpi_get_type() noexcept
    {
        MPI_Datatype mpi_type = MPI_DATATYPE_NULL;

        if constexpr (std::is_same_v<T, char>) {
            mpi_type = MPI_CHAR;
        }
        else if constexpr (std::is_same_v<T, signed char>) {
            mpi_type = MPI_SIGNED_CHAR;
        }
        else if constexpr (std::is_same_v<T, unsigned char>) {
            mpi_type = MPI_UNSIGNED_CHAR;
        }
        else if constexpr (std::is_same_v<T, wchar_t>) {
            mpi_type = MPI_WCHAR;
        }
        else if constexpr (std::is_same_v<T, signed short>) {
            mpi_type = MPI_SHORT;
        }
        else if constexpr (std::is_same_v<T, unsigned short>) {
            mpi_type = MPI_UNSIGNED_SHORT;
        }
        else if constexpr (std::is_same_v<T, signed int>) {
            mpi_type = MPI_INT;
        }
        else if constexpr (std::is_same_v<T, unsigned int>) {
            mpi_type = MPI_UNSIGNED;
        }
        else if constexpr (std::is_same_v<T, signed long int>) {
            mpi_type = MPI_LONG;
        }
        else if constexpr (std::is_same_v<T, unsigned long int>) {
            mpi_type = MPI_UNSIGNED_LONG;
        }
        else if constexpr (std::is_same_v<T, signed long long int>) {
            mpi_type = MPI_LONG_LONG;
        }
        else if constexpr (std::is_same_v<T, unsigned long long int>) {
            mpi_type = MPI_UNSIGNED_LONG_LONG;
        }
        else if constexpr (std::is_same_v<T, float>) {
            mpi_type = MPI_FLOAT;
        }
        else if constexpr (std::is_same_v<T, double>) {
            mpi_type = MPI_DOUBLE;
        }
        else if constexpr (std::is_same_v<T, long double>) {
            mpi_type = MPI_LONG_DOUBLE;
        }
        else if constexpr (std::is_same_v<T, std::int8_t>) {
            mpi_type = MPI_INT8_T;
        }
        else if constexpr (std::is_same_v<T, std::int16_t>) {
            mpi_type = MPI_INT16_T;
        }
        else if constexpr (std::is_same_v<T, std::int32_t>) {
            mpi_type = MPI_INT32_T;
        }
        else if constexpr (std::is_same_v<T, std::int64_t>) {
            mpi_type = MPI_INT64_T;
        }
        else if constexpr (std::is_same_v<T, std::uint8_t>) {
            mpi_type = MPI_UINT8_T;
        }
        else if constexpr (std::is_same_v<T, std::uint16_t>) {
            mpi_type = MPI_UINT16_T;
        }
        else if constexpr (std::is_same_v<T, std::uint32_t>) {
            mpi_type = MPI_UINT32_T;
        }
        else if constexpr (std::is_same_v<T, std::uint64_t>) {
            mpi_type = MPI_UINT64_T;
        }
        else if constexpr (std::is_same_v<T, bool>) {
            mpi_type = MPI_C_BOOL;
        }
        //    else if constexpr (std::is_same_v<T, std::complex<float>>) {
        //        mpi_type = MPI_C_COMPLEX;
        //    }
        //    else if constexpr (std::is_same_v<T, std::complex<double>>) {
        //        mpi_type = MPI_C_DOUBLE_COMPLEX;
        //    }
        //    else if constexpr (std::is_same_v<T, std::complex<long double>>) {
        //        mpi_type = MPI_C_LONG_DOUBLE_COMPLEX;
        //    }
        else if constexpr (std::is_same_v<T, v3D>) {
            mpi_type = Parallel::MPI_V3D;
        }
        else if constexpr (std::is_same_v<T, i3D>) {
            mpi_type = Parallel::MPI_I3D;
        }
        else if constexpr (std::is_same_v<T, v2D>) {
            mpi_type = Parallel::MPI_V2D;
        }
        else if constexpr (std::is_same_v<T, i2D>) {
            mpi_type = Parallel::MPI_I2D;
        }
        else if constexpr (std::is_same_v<T, std::pair<int, int>>) {
            mpi_type = Parallel::MPI_PAIRII;
        }

        assert(mpi_type != MPI_DATATYPE_NULL);

        return mpi_type;
    }



//}//namespace VMlib

#endif
