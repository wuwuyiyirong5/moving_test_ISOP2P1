/**
 * @file   ISOP2P1.cpp
 * @author Heyu Wang <scshw@cslin107.csunix.comp.leeds.ac.uk>
 * @date   Tue Nov  4 14:33:17 2014
 * 
 * @brief ISO P2P1 有限元的网格, 空间, 矩阵信息集中封装类的实现.
 * 
 * 
 */

#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::initialize()
{
	/// 读取配置文件.
	config("config");
	/// 构建网格.
	buildMesh();
	/// 构建混合有限元 ISOP1P2 空间. 
	buildFEMSpace();
	/// 构建稀疏矩阵模板.
	buildMatrixStruct();
	/// 构建矩阵.
	buildMatrix();
	if (Stokes_init == 1 || NS_init == 1)
		solveStokes();
};

void ISOP2P1::run()
{
	initialize();
	std::cout << "Initialize mesh ... " << std::endl;
	if (isMoving == 1)
	{
		double scale_step = 0.2;
		scale = scale_step;
		do {
			v_h[0].reinit(fem_space_v);
			v_h[1].reinit(fem_space_v);
			p_h.reinit(fem_space_p);
			buildMatrix();
			// stepForwardEuler();
			solveStokes();
			v_h[0].scale(scale);
			v_h[1].scale(scale);
			movingMesh();
			std::cout << "\r\tscale = " << scale << std::endl;
			scale += scale_step;
		} while (scale <= 1.0);
	}
	outputSolution();
	v_h[0].reinit(fem_space_v);
	v_h[1].reinit(fem_space_v);
	p_h.reinit(fem_space_p);
	solveStokes();
	outputTecplotP("P0");
	outputTecplot("initial_value0");
	getchar();

	int steps = 0;
	bool isOutput = false;

	while (t < t1)
	{
		/// 准备每步输出一个 tecplot 数据, 注意别把硬盘写爆了!
		std::stringstream ss;
		ss << "NS_Euler";

		if (scheme == 1)
		{
			if(isMoving == 1)
				buildMatrix();
			stepForwardEuler();
			time_step();
		}
		else if (scheme == 2)
			stepForwardLinearizedEuler();
		else
			break;
		t += dt;
		/// 网格移动.
		if (isMoving == 1)
			movingMesh();
		/// 输出.
		if (isOutput)
		{
			steps++;
			ss << steps;
			outputTecplot(ss.str());
			std::cout << "Data outputed!" << std::endl;
			isOutput = false;
		}
		double t_stop = int((t / 0.01) + 0.5) * 0.01;
		if (t > t_stop - dt && t < t_stop + dt)
		{
			t = t_stop;
			isOutput = true;
		}
		std::cout << "t = " << t << std::endl;
	}

}
#undef DIM
