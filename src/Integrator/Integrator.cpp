// Проект Integrator 
// (c) А.И. Гумирова, И.К. Марчевский, С.Р. Серафимова, 2022

/*!
\file
\brief Основной файл c функцией main
\author Гумирова Алия Ильдусовна
\author Марчевский Илья Константинович
\author Серафимова София Романовна

\date 02 апреля 2022 г.
\version 0.3
*/

#include <omp.h>
#include <fstream>

#include "CompTest.h"
#include "CompI3DM.h"
#include "CompJ3DM.h"
#include "CompI2DM.h"
#include "CompJ2DM.h"
#include "CompJ3DK.h"
#include "CompI2DK.h"
#include "CompJ2DK.h"

#include <algorithm>

int main(int argc, char** argv)
{
	Parallel par;

	// Инициализация подсистемы MPI
	MPI_Init(&argc, &argv);
	// Получить размер коммуникатора MPI_COMM_WORLD
	// (общее число процессов в рамках задачи)
	MPI_Comm_size(MPI_COMM_WORLD, &par.np);
	// Получить номер текущего процесса в рамках 
	// коммуникатора MPI_COMM_WORLD
	MPI_Comm_rank(MPI_COMM_WORLD, &par.rank);

	//Имя узла
	char proc_name[255];
	int namelen;
	MPI_Get_processor_name(proc_name, &namelen);

	Parallel::CreateMpiType();

	

	
	int nNode, nTopo;

	//std::cout << "Hello World!\n";
	//std::cout << "proc = " << proc_name << "\n";
	//std::cout << "omp = " << omp_get_max_threads() << "\n";

	/*
	Database<3> db3;
	if (par.rank == 0)
	{		
		//db3.readNodeTopoFromFile("G1.dat");
		//db3.readNodeTopoFromFile("13bad");
		db3.readNodeTopoFromFile("Girja");
		//db3.readNodeTopoFromFile("G1contactR.dat");
		//db3.readNodeTopoFromFile("Krylo01.dat");
		//db3.readNodeTopoFromFile("cubehole.dat", 20.0);
		//db3.readNodeTopoFromFile("Fish.dat", 100.0);
		//db3.readNodeTopoFromFile("s5m.dat", 0.0005); ///// <--- (359, 360)
		//db3.readNodeTopoFromFile("s5m2.dat", 0.0005);
		//db3.readNodeTopoFromFile("MeshScreen.dat", 0.1);
		//db3.readNodeTopoFromFile("0012e2", 1.0);
		//db3.readNodeTopoFromFile("1x1x1_extrafine", 1.0);
		
		//db3.readNodeTopoFromFile("Case1.dat", 1.0);
		//db3.readNodeTopoFromFile("Case-6-4.dat", 1.0);

		nNode = (int)db3.node.size();
		nTopo = (int)db3.topo.size();
	}

	//Рассылка базы данных
	par.SplitMPI(db3.node.size()).BcastVector(db3.node);
	par.SplitMPI(db3.topo.size()).BcastVector(db3.topo);

	db3.fillNbh();
	db3.calcNrm();
	db3.calcCnt();
	db3.calcMeasure();

	db3.point = db3.cnt;

	//CompJ3DK cmp(db3, par);
	//if (par.rank == 0)
	//{
	//	for (int i = 0; i < 1 * db3.topo.size(); ++i)
	//		for (int j = 0; j < 1 * db3.topo.size(); ++j)
	//			cmp.task.push_back({ i % (int)db3.topo.size(), j % (int)db3.topo.size() });
	//}
	
	Gausspoints<3> gaussianQuadratures(db3, gp3D13);  //доступны ngp = 1 3 4 6 7 9 12 13 гауссовых точек

	std::cout << "db_size = " << db3.topo.size() << std::endl;

	CompJ3DK cmp(db3, par, &gaussianQuadratures);
	
	//std::cout << db3.topo[5826] << std::endl;
	
	if (par.rank == 0)
	{
		//for (int i = 0; i < db3.topo.size(); ++i)
		int stp = 2000;
		int q = 0;
		
		std::cout << "range : [ " << q * stp << "..." << std::min((q + 1) * stp, (int)db3.topo.size()) << " )" << std::endl;
		//int i = 1724;
		for (int i = q* stp; i < std::min((q+1)*stp, (int)db3.topo.size()); ++i)
		for (int j = 0; j < db3.topo.size(); ++j)
		{
			//int i = 60;
			//int j = 43;
			
			cmp.task.push_back({ i % (int)db3.topo.size(), j % (int)db3.topo.size() });
			cmp.task.push_back({ j % (int)db3.topo.size(), i % (int)db3.topo.size() });
		}
	}
	*/

	

	


	Database<2> db2;
	if (par.rank == 0)
	{
		db2.readNodeTopoFromFile("G2-1600.dat");
		nNode = (int)db2.node.size();
		nTopo = (int)db2.topo.size();
	}
	
	
	////Рассылка базы данных
	par.SplitMPI(db2.node.size()).BcastVector(db2.node);
	par.SplitMPI(db2.topo.size()).BcastVector(db2.topo);

	db2.fillNbh();
	db2.calcNrm();
	db2.calcCnt();
	db2.calcMeasure();
	db2.point = db2.cnt;

	//Тест
	//Gausspoints<2> gaussianQuadratures(db2, gp2D1);  //доступны ngp = 1 гауссовых точек
	//auto q = gaussianQuadratures.integrate<double, 10>([&](const v2D& pt) { return (pt-db2.node[0]).length2(); }, 0);
	//std::cout << q << std::endl;
	
	
	CompI2DK cmp(db2, par);
	if (par.rank == 0)
	{
		for (int i = 0; i < 1*db2.topo.size(); ++i)
			for (int j = 0; j < 1*db2.topo.size(); ++j)		
				cmp.task.push_back({ i % (int)db2.topo.size(), j % (int)db2.topo.size() });
	}
	

	MPI_Barrier(MPI_COMM_WORLD);
	double t1 = omp_get_wtime();

	cmp.run(false);

	//std::cout << cmp.result[0] << std::endl;

	/* //3D
	size_t numpairs = cmp.result.size() / 2;
	for (size_t i = 0; i < numpairs; ++i)
	{
		v3D A = cmp.result[2 * i + 0];
		v3D B = cmp.result[2 * i + 1];
		double diff = (A + B).length() / std::max(A.length(), B.length());
		//if (diff > 1e-4)
		if ((cmp.task[2 * i].first!= cmp.task[2 * i].second) && (diff > 1e-3 || std::isnan(diff)))
		{

			if (db3.ifSosed({ cmp.task[2 * i].first, cmp.task[2 * i].second }))
			{
				std::cout << "i,j = {" << cmp.task[2 * i].first << ", " << cmp.task[2 * i].second << "}, " << diff << ",  " << A << " " << B << std::endl;
				std::cout << "SOSED!" << std::endl;

			}
			else if (db3.ifContact({ cmp.task[2 * i].first, cmp.task[2 * i].second }))
			{
				std::cout << "i,j = {" << cmp.task[2 * i].first << ", " << cmp.task[2 * i].second << "}, " << diff << ",  " << A << " " << B << std::endl;
				if (std::max(A.length(), B.length()) > 1e-5)
					std::cout << "CONTACT!" << std::endl;
				else
					std::cout << "CONTACT, but value = 0!" << std::endl;
			}
			else
			{
				std::cout << "i,j = {" << cmp.task[2 * i].first << ", " << cmp.task[2 * i].second << "}, " << diff << ",  " << A << " " << B << std::endl;
				std::cout << "FAR!" << std::endl;
			}
		}
		//if
		//	std::cout << "ContactOK" << std::endl;
	}
	*/



	MPI_Barrier(MPI_COMM_WORLD);
	double t2 = omp_get_wtime();

	std::cout << t2 - t1 << " sec." << std::endl;

	//std::cout << cmp.result[0] << " " << cmp.result[1] << std::endl;



	//*
	if (par.rank == 0)
	{
		std::ofstream of("result.txt");
		of.precision(17);
		for (const auto& res : cmp.result)
			of << res << '\n';
			//of << res[0] << " " << res[1] << '\n';
		of.close();
	}
	//*/

	//std::cout << cmp.scalarResult[1 * db2.topo.size() + 25] << std::endl;

	MPI_Finalize();

}

