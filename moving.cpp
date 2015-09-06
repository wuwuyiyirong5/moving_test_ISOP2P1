#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"

#define DIM 2

void ISOP2P1::syncMesh()
{
    int n_point = this->n_point();

    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();

    for (int i = 0; i < mesh_p.n_geometry(0); ++i)
    {
    	(*mesh_p.h_geometry<0>(i))[0] = this->point(i)[0];
    	(*mesh_p.h_geometry<0>(i))[1] = this->point(i)[1];
    }

    /// 更新p网格的中点. 学习这个过程, 似乎可以去掉单元对应? 
    for (int j = 0; j < mesh_p.n_geometry(1); ++j)
    {
    	GeometryBM &bnd = mesh_p.geometry(1, j);
    	(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[0]
    	    = 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[0] + 
    		     (*mesh_p.h_geometry<0>(bnd.vertex(1)))[0]);
    	(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[1]
    	    = 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[1] + 
    		     (*mesh_p.h_geometry<0>(bnd.vertex(1)))[1]);
    }

    irregular_mesh_p->semiregularize();
    irregular_mesh_p->regularize(false);
    irregular_mesh_v->semiregularize();
    irregular_mesh_v->regularize(false);
};

void ISOP2P1::getMonitor()
{
	/// 限制最大迭代步数.
	maxStep() = max_step;
	/// 同步网格.
	syncMesh();
	// // /// 输出一下.
	// outputTecplot("NS_Euler");
	// getchar();

	FEMFunction<double, DIM> _u_h(fem_space_p);
	FEMFunction<double, DIM> _v_h(fem_space_p);
	Operator::L2Interpolate(v_h[0], _u_h);
	Operator::L2Interpolate(v_h[1], _v_h);
    
	FEMSpace<double, DIM>::ElementIterator the_element = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element = fem_space_p.endElement();
	for (int i = 0; the_element != end_element; ++the_element) 
	{
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(2);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
		std::vector<double> u_h_value = _u_h.value(q_point, *the_element);
		std::vector<double> v_h_value = _v_h.value(q_point, *the_element);
		std::vector<std::vector<double> > u_h_gradient = _u_h.gradient(q_point, *the_element);
		std::vector<std::vector<double> > v_h_gradient = _v_h.gradient(q_point, *the_element);
		float d = 0, area = 0, norm = 0;
		for (int l = 0; l < n_quadrature_point; ++l) 
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			area += Jxw;
			norm += u_h_value[l] * u_h_value[l] + v_h_value[l] * v_h_value[l];
			d += Jxw * fabs(v_h_gradient[l][0] - u_h_gradient[l][1]);
		}
		norm = 1.0 / (eps + sqrt(norm));
		monitor(i++) = d * norm / area;
	}
	std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
		  << "\tmin monitor=" << *std::min_element(monitor().begin(), monitor().end())
		  << std::endl;
	double max_monitor = *std::max_element(monitor().begin(), monitor().end());
	smoothMonitor(2);
	for (int i = 0; i < n_geometry(2); ++i)
		monitor(i) = 1.0 / sqrt(1.0 + alpha * monitor(i));
};

void ISOP2P1::updateSolution()
{
	/// 更新插值点.
	fem_space_p.updateDofInterpPoint();
	fem_space_v.updateDofInterpPoint();
	
	/// 备份数值解.
	FEMFunction<double, DIM> _u_h(v_h[0]);
	FEMFunction<double, DIM> _v_h(v_h[1]);
	FEMFunction<double, DIM> _p_h(p_h);
	const double& msl = moveStepLength();
	/// 因为有限元空间插值点变化, 重新构造矩阵.
	buildMatrixStruct();
	buildMatrix();
	/// 因为网格移动量为小量, 因此时间步长可以相对取的大些.
	double _dt = 1.0;
    
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;
	int n_total_dof_v = 2 * n_dof_v;
	
	/// 一步Euler.
	for (int m = 3; m > 0; --m)
	{
		/// 系数矩阵直接使用 Stokes 矩阵结构.
		SparseMatrix<double> mat_moving;
		mat_moving.reinit(sp_stokes);

		/// (0, 0) 
		for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_vxvx[i]) = mat_v_mass.global_entry(i); 
		/// (1, 1) 这两个对角块仅有质量块.
		for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_vyvy[i]) = mat_v_mass.global_entry(i);
		/// (0, 2) 这个不是方阵. 在矩阵结构定义的时候已经直接排除了对角元优
		/// 先.
		for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_pvx[i]) = (1.0 / m) * mat_pvx_divT.global_entry(i);

		/// (1, 2)
		for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_pvy[i]) = (1.0 / m) * mat_pvy_divT.global_entry(i);

		/// (2, 0)
		for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_vxp[i]) = (1.0 / m) * mat_vxp_div.global_entry(i);
	
		/// (2, 1) 这四块直接复制散度矩阵. 
		for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
			mat_moving.global_entry(index_vyp[i]) = (1.0 / m) * mat_vyp_div.global_entry(i);

		/// 问题右端项.
		rhs.reinit(n_total_dof);

		FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
		FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
		/// 遍历速度单元, 拼装相关系数矩阵和右端项.
		for (the_element_v = fem_space_v.beginElement(); 
		     the_element_v != end_element_v; ++the_element_v) 
		{
			/// 当前单元信息.
			double volume = the_element_v->templateElement().volume();
			/// 积分精度至少为2, u 和 p 都是 1 次, 梯度和散度 u 都是常数.
			const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(2);
			std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
 			/// 速度单元信息.
			std::vector<std::vector<double> > basis_value_v = the_element_v->basis_function_value(q_point);
			std::vector<double> vx_value = _u_h.value(q_point, *the_element_v);
			std::vector<double> vy_value = _v_h.value(q_point, *the_element_v);
			std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
			std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
			const std::vector<int>& element_dof_v = the_element_v->dof();
			int n_element_dof_v = the_element_v->n_dof();
			/// 速度积分点上的移动方向, 注意是速度单元的积分点, 在压力单元上的移动方向.
			/// 会不会是在这里出问题: 移动的是P网格, 但是对应4个速度单元上的积分点也位于大的压力单元上, 
			/// 可以看作是压力单元上的点, 所以下面一行程序, 想想应该是没有问题.
			std::vector<std::vector<double> > move_vector = moveDirection(q_point, index_v2p[the_element_v->index()]);
			for (int l = 0; l < n_quadrature_point; ++l)
			{
				double Jxw = quad_info.weight(l) * jacobian[l] * volume;
				for (int i = 0; i < n_element_dof_v; ++i)
				{
					double rhs_cont = (vx_value[l] + (1.0 / m) * msl * innerProduct(move_vector[l], vx_gradient[l])) * basis_value_v[i][l];
					rhs_cont *= Jxw;
					rhs(element_dof_v[i]) += rhs_cont;

					rhs_cont = (vy_value[l] + (1.0 / m) * msl * innerProduct(move_vector[l], vy_gradient[l])) * basis_value_v[i][l];
					rhs_cont *= Jxw;
					rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
				}
			}
		}
		// /// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
		// AccuracyVx accuracy_vx(viscosity, t);
		// AccuracyVy accuracy_vy(viscosity, t);

		/// 构建系数矩阵和右端项.
		Vector<double> x(n_total_dof);
		// AccuracyVx real_vx(viscosity, t);
		// AccuracyVy real_vy(viscosity, t);
		// Operator::L2Project(real_vx, v_h[0], Operator::LOCAL_LEAST_SQUARE, 3);
		// Operator::L2Project(real_vy, v_h[1], Operator::LOCAL_LEAST_SQUARE, 3);
		
		// /// 边界处理.
		const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
		const unsigned int * colnum = sp_stokes.get_column_numbers();

		/// 遍历全部维度的速度节点.
		for (unsigned int i = 0; i < n_total_dof_v; ++i)
		{
			/// 边界标志.
			int bm = -1;
			/// 判断一下是 x 方向还是 y 方向. 分别读取标志.
			if (i < n_dof_v)
				bm = fem_space_v.dofInfo(i).boundary_mark;
			else
				bm = fem_space_v.dofInfo(i - n_dof_v).boundary_mark;

			if (bm == 0)
				continue;
			/// 对 Dirichelet 边界根据边界分别赋值. 注意同时还要区别 x 和
			if (bm == 1 || bm == 2 || bm == 4)
				x(i) = 0.0;
			if (bm == 3)
				if (i < n_dof_v)
				{
					Regularized regularize;
					x(i) = regularize.value(fem_space_v.dofInfo(i).interp_point);
				}
				else
					x(i) = 0.0;
			// if (bm == 1 || bm == 2 || bm == 3 || bm == 4)
			// 	if (i < n_dof_v)
			// 		x(i) = real_vx.value(fem_space_v.dofInfo(i).interp_point);
			// 	else
			// 		x(i) = real_vy.value(fem_space_v.dofInfo(i - n_dof_v).interp_point);
			/// 右端项这样改, 如果该行和列其余元素均为零, 则在迭代中确
			/// 保该数值解和边界一致.
			if (bm == 1 || bm == 2 || bm == 3 || bm == 4)
			{
				rhs(i) = mat_moving.diag_element(i) * x(i);
				/// 遍历 i 行.
				for (unsigned int j = rowstart[i] + 1; j < rowstart[i + 1]; ++j)
				{
					/// 第 j 个元素消成零(不是第 j 列!). 注意避开了对角元.
					mat_moving.global_entry(j) -= mat_moving.global_entry(j);
					/// 第 j 个元素是第 k 列.
					unsigned int k = colnum[j];
					/// 看看第 k 行的 i 列是否为零元.
					const unsigned int *p = std::find(&colnum[rowstart[k] + 1],
									  &colnum[rowstart[k + 1]],
									  i);
					/// 如果是非零元. 则需要将这一项移动到右端项. 因为第 i 个未知量已知.
					if (p != &colnum[rowstart[k + 1]])
					{
						/// 计算 k 行 i 列的存储位置.
						unsigned int l = p - &colnum[rowstart[0]];
						/// 移动到右端项. 等价于 r(k) = r(k) - x(i) * A(k, i).
						rhs(k) -= mat_moving.global_entry(l) * x(i);
						/// 移完此项自然是零.
						mat_moving.global_entry(l) -= mat_moving.global_entry(l);
					}
				}
			}
		}
		std::cout << "boundary values for updateSolution OK!" << std::endl;

		clock_t t_cost = clock();
		/// 预处理矩阵.
		/// 不完全LU分解.     
		dealii::SparseILU <double> preconditioner;
		preconditioner.initialize(mat_moving);	
		/// 矩阵求解. 
		dealii::SolverControl solver_control (4000000, l_tol, check);

		SolverMinRes<Vector<double> > minres (solver_control);
		SolverGMRES<Vector<double> >::AdditionalData para(5000, false, true);
		SolverGMRES<Vector<double> > gmres (solver_control, para);
		

		minres.solve (mat_moving, x, rhs, PreconditionIdentity());
		// gmres.solve(mat_moving, x, rhs, preconditioner);
		t_cost = clock() - t_cost;	
		std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;
		for (int i = 0; i < n_dof_v; ++i)
		{
			v_h[0](i) = x(i);
			v_h[1](i) = x(i + n_dof_v);
		}
		for (int i = 0; i < n_dof_p; ++i)
			p_h(i) = x(i + 2 * n_dof_v);
	
		// /// debug
		// std::ofstream mat_deb;
		// // rowstart = sp_pvx.get_rowstart_indices();
		// // colnum = sp_pvx.get_column_numbers();
		// mat_deb.open("mat.m", std::ofstream::out);
		// mat_deb.setf(std::ios::fixed);
		// mat_deb.precision(20);
	    
		// for (int i = 0; i < n_total_dof; ++i)
		// {
		// 	for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
		// 	{
		// 		mat_deb << "A(" << i + 1 << ", " << colnum[j] + 1 << ")=" 
		// 			<< mat_moving.global_entry(j) << ";" << std::endl;
		// 	}
		// 	mat_deb << "x(" << i + 1<< ") = " << x(i) << ";" << std::endl;
		// 	mat_deb << "rhs(" << i + 1 << ") = " << rhs(i) << ";" << std::endl;
		// }
		// // for(int i = 0; i < n_dof_p; ++i)
		// // 	mat_deb << "x(" << i + 1<< ") = " << p_h(i) << ";" << std::endl;
		// mat_deb.close();
		// std::cout << "mat output" << std::endl;
		Vector<double> res(n_total_dof);
		mat_moving.vmult(res, x);
		res *= -1;
		res += rhs;
		std::cout << "res_l2norm =" << res.l2_norm() << std::endl;
                /// debug
	}
}; 

void ISOP2P1::outputSolution()
{
	outputPhysicalMesh("E");
	p_h.writeOpenDXData("p_h.dx");
};


void ISOP2P1::movingMesh()
{
	moveMesh();
	syncMesh();
	fem_space_p.updateDofInterpPoint();
	fem_space_v.updateDofInterpPoint();
};

