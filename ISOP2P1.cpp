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
		solveStokes();
	    v_h[0].scale(scale);
	    v_h[1].scale(scale);
	    movingMesh();
	    std::cout << "\r\tscale = " << scale << std::endl;
	    scale += scale_step;
	} while (scale <= 1.0);
    }
    outputSolution();
    outputTecplotP("P0");
    outputTecplot("initial_value0");
    getchar();
    do {
	buildMatrix();
	stepForwardEuler();
	time_step();
	t += dt;
	if (isMoving == 1)
	{
	    scale = 1.0;
	    movingMesh();
	    outputSolution();
	}
	std::cout << "t  = " << t << std::endl;
    } while (t < t1);
    outputTecplotP("P1");
    outputTecplot("V1");

};

#undef DIM
