#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::solveStokes()
{
	buildStokesSys();
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// 构建系数矩阵和右端项.
	/// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
	Vector<double> x(n_total_dof);
	rhs.reinit(n_total_dof);
	/// 边界条件一起处理了.
	boundaryValueStokes(x);

	/// 矩阵求解. 
	dealii::SolverControl solver_control (400000, l_tol, check);
	SolverMinRes<Vector<double> > minres (solver_control);
	StokesPreconditioner preconditioner;
	/// 预处理矩阵.
	SparseMatrix<double> matrix_vxvx(sp_vxvx);
	SparseMatrix<double> matrix_vyvy(sp_vyvy);
	/// 这里从 Stokes 取是因为加了边界条件.
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix_vxvx.global_entry(i) = matrix.global_entry(index_vxvx[i]);
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix_vyvy.global_entry(i) = matrix.global_entry(index_vyvy[i]);

	preconditioner.initialize(matrix_vxvx, matrix_vyvy, mat_p_mass);

	clock_t t_cost = clock();
	minres.solve (matrix, x, rhs, preconditioner);
	// minres.solve (matrix, x, rhs, PreconditionIdentity());
	t_cost = clock() - t_cost;

	std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;

	/// 将整体数值解分割成速度和压力.
	for (int i = 0; i < n_dof_v; ++i)
	{
		v_h[0](i) = x(i);
		v_h[1](i) = x(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
		p_h(i) =  x(i + 2 * n_dof_v);

	/// 记录计算结果参数.
	std::ofstream output;
	output.open("record", std::ofstream::out | std::ofstream::app);
	output.setf(std::ios::fixed);
	output.precision(20);
	output << "nu = " << n_dof_v << std::endl;
	output << "np = " << n_dof_p << std::endl;
	output << "steps = " << solver_control.last_step() << std::endl;
	output << "linear tol = " << l_tol << std::endl;
	if (error_check == true)
	{
		// RealVx real_Vx;
		// RealVy real_Vy;    
		// double mean_p_h= Functional::meanValue(p_h, 3); 
		// RealP real_P(mean_p_h);
		AccuracyVx accuracy_vx(viscosity, t + dt);
		AccuracyVy accuracy_vy(viscosity, t + dt);


		FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
		FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();
		/// 误差.
		double H1_err = 0.0;
		double L2_err = 0.0;

		/// 遍历速度单元, 拼装相关系数矩阵和右端项.
		for (; the_element_v != end_element_v; ++the_element_v) 
		{
			/// 当前单元信息.
			double volume = the_element_v->templateElement().volume();
			/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
			/// 装时积分精度不用超过 1 次. (验证一下!)
			const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(1);
			std::vector<double> jacobian 
				= the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<Point<DIM> > q_point 
				= the_element_v->local_to_global(quad_info.quadraturePoint());
			/// 速度信息.
			std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
			std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
			std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
			std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
			/// 实际拼装.
			for (int l = 0; l < n_quadrature_point; ++l)
			{
				double Jxw = quad_info.weight(l) * jacobian[l] * volume;
				double dx_value = vx_value[l] - accuracy_vx.value(q_point[l]);
				double dy_value = vy_value[l] - accuracy_vy.value(q_point[l]);
				L2_err += Jxw * (dx_value * dx_value + dy_value * dy_value);
				std::vector<double> real_vx_gradient = accuracy_vx.gradient(q_point[l]);
				std::vector<double> real_vy_gradient = accuracy_vy.gradient(q_point[l]);

				for (int i = 0; i < DIM; ++i)
				{
					dx_value = vx_gradient[l][i] - real_vx_gradient[i];
					dy_value = vy_gradient[l][i] - real_vy_gradient[i];
					H1_err += Jxw * (dx_value * dx_value + dy_value * dy_value); 
				}
			}
		}
		H1_err = sqrt(H1_err);
		L2_err = sqrt(L2_err);

		double error;

		error = Functional::L2Error(v_h[0], accuracy_vx, 2);
		std::cout << "|| u - u_h ||_L2 = " << error << std::endl;
		// output << "|| u - u_h ||_L2 = " << error << std::endl;

		error = Functional::H1SemiError(v_h[0], accuracy_vx, 1);
		std::cout << "|| u - u_h ||_H1 = " << error << std::endl;
		// error = Functional::L2Error(v_h[0], real_Vx, 2);
		// std::cout << "|| u - u_h ||_L2 = " << error << std::endl;
		// // output << "|| u - u_h ||_L2 = " << error << std::endl;

		// error = Functional::H1SemiError(v_h[0], real_Vy, 1);
		// std::cout << "|| u - u_h ||_H1 = " << error << std::endl;
		// output << "|| u - u_h ||_H1 = " << error << std::endl;

		// error = Functional::H1SemiError(p_h, real_P, 1);
		// std::cout << "|| p - p_h ||_H1 = " << error << std::endl;
		// // output << "ph_H1err(1) = " << error << std::endl;

		// error = Functional::H1SemiError(p_h_refine, real_P, 1);
		// std::cout << "|| p - p_h_refine ||_H1 = " << error << std::endl;
		// output << "ph_refine_H1err(1) = " << error << std::endl;

		// error = Functional::L2Error(p_h, real_P, 2);
		// std::cout << "|| p - p_h ||_L^2 = " << error << std::endl;
		// // output << "ph_L2err(1) = " << error << std::endl;

		// error = Functional::L2Error(p_h_refine, real_P, 2);
		// std::cout << "|| p - p_h_refine ||_L^2 = " << error << std::endl;
		// output << "ph_refine_L2err(1) = " << error << std::endl;

		std::cout << "uh_L2err() = " << L2_err << std::endl;
		std::cout << "uh_H1err() = " << H1_err << std::endl;	
		output << "uh_L2err(1) = " << L2_err << std::endl;
		output << "uh_H1err(1) = " << H1_err << std::endl;	
	}
	output.close();
};

void ISOP2P1::solveNS(int method)
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = n_dof_v * 2 + n_dof_p;

	/// 开始迭代.
	double error_N = 1.0;
	int iteration_times = 0;
	while (error_N > n_tol)
	{
		std::ofstream output;
		output.open("record", std::ofstream::out);
		output.setf(std::ios::fixed);
		output.precision(20);
		std::cout.setf(std::ios::scientific);
		std::cout.precision(20);

		/// Newton 迭代或 Picard 迭代.
		/// 先更新和速度场有关的矩阵块. 
		updateNonlinearMatrix();
		/// 构建迭代矩阵.
		if (method == 1)
			buildNewtonSys4NS();
		else if (method == 2)
			buildPicardSys4NS();
		else if (method == 3)
			if (iteration_times < 2)
				buildPicardSys4NS();
			else
				buildNewtonSys4NS();
		else
		{
			std::cout << "Newton: 1, Picard: 2, Hybrid: 3." << std::endl;
			exit(1);
		}

		/// 建立右端项.
	
		rhs.reinit(n_total_dof);

		FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
		FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();
		FEMSpace<double,2>::ElementIterator the_element_p = fem_space_p.beginElement();
		FEMSpace<double,2>::ElementIterator end_element_p = fem_space_p.endElement();
		/// 遍历速度单元, 拼装相关系数矩阵和右端项.
		for (the_element_v = fem_space_v.beginElement(); 
		     the_element_v != end_element_v; ++the_element_v) 
		{
			/// 当前单元信息.
			double volume = the_element_v->templateElement().volume();
			/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
			/// 装时积分精度不用超过 1 次. (验证一下!)
			const QuadratureInfo<2>& quad_info = the_element_v->findQuadratureInfo(4);
			std::vector<double> jacobian 
				= the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<Point<2> > q_point 
				= the_element_v->local_to_global(quad_info.quadraturePoint());
			/// 速度单元信息.
			std::vector<std::vector<std::vector<double> > > basis_gradient_v 
				= the_element_v->basis_function_gradient(q_point);
			std::vector<std::vector<double> >  basis_value_v 
				= the_element_v->basis_function_value(q_point);
			const std::vector<int>& element_dof_v = the_element_v->dof();
			std::vector<double> fx_value = source_v[0].value(q_point, *the_element_v);
			std::vector<double> fy_value = source_v[1].value(q_point, *the_element_v);
			int n_element_dof_v = the_element_v->n_dof();
			std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
			std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
			std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
			std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
			/// 压力单元信息.
			Element<double, 2> &p_element = fem_space_p.element(index_v2p[the_element_v->index()]);
			const std::vector<int>& element_dof_p = p_element.dof();
			std::vector<std::vector<std::vector<double> > > basis_gradient_p 
				= p_element.basis_function_gradient(q_point);
			std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
			std::vector<double> p_value = p_h.value(q_point, p_element);
			int n_element_dof_p = p_element.n_dof();
			/// 实际拼装.
			for (int l = 0; l < n_quadrature_point; ++l)
			{
				double Jxw = quad_info.weight(l) * jacobian[l] * volume;
				for (int i = 0; i < n_element_dof_v; ++i)
				{
					double rhs_cont = fx_value[l] * basis_value_v[i][l];
					rhs_cont -= (vx_value[l] * vx_gradient[l][0] + 
						     vy_value[l] * vx_gradient[l][1]) * basis_value_v[i][l];
					rhs_cont -= viscosity * innerProduct(basis_gradient_v[i][l], vx_gradient[l]);
					rhs_cont += p_value[l] * basis_gradient_v[i][l][0];
					rhs_cont *= Jxw;
					rhs(element_dof_v[i]) += rhs_cont;

					rhs_cont = fy_value[l] * basis_value_v[i][l];
					rhs_cont -= (vx_value[l] * vy_gradient[l][0] + 
						     vy_value[l] * vy_gradient[l][1]) * basis_value_v[i][l];
					rhs_cont -= viscosity * innerProduct(basis_gradient_v[i][l], vy_gradient[l]);
					rhs_cont += p_value[l] * basis_gradient_v[i][l][1];
					rhs_cont *= Jxw;
					rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
				}
			}
		}

		/// 遍历压力单元. 拼装矩阵和右端项.
		for (the_element_p = fem_space_p.beginElement(); 
		     the_element_p != end_element_p; ++the_element_p) 
		{
			const std::vector<int>& element_dof_p = the_element_p->dof();
			int n_element_dof_p = the_element_p->n_dof();
			for (int i = 0; i < n_element_dof_p; ++i)
			{
				int idx_p = the_element_p->index();
				int n_chi = index_p2v[idx_p].size();
				for (int k = 0; k < n_chi; k++)
				{
					/// 速度单元信息.
					Element<double, 2> &v_element = fem_space_v.element(index_p2v[idx_p][k]);
					/// 几何信息.
					double volume = v_element.templateElement().volume();
					const QuadratureInfo<2>& quad_info = v_element.findQuadratureInfo(4);
					std::vector<double> jacobian 
						= v_element.local_to_global_jacobian(quad_info.quadraturePoint());
					int n_quadrature_point = quad_info.n_quadraturePoint();
					std::vector<Point<2> > q_point 
						= v_element.local_to_global(quad_info.quadraturePoint());
					std::vector<std::vector<double> > vx_gradient 
						= v_h[0].gradient(q_point, v_element);
					std::vector<std::vector<double> > vy_gradient 
						= v_h[1].gradient(q_point, v_element);
					std::vector<double> vx_value = v_h[0].value(q_point, v_element);
					std::vector<double> vy_value = v_h[1].value(q_point, v_element);
					/// 压力单元信息.
					std::vector<std::vector<double> >  basis_value_p 
						= the_element_p->basis_function_value(q_point);
					/// 具体拼装.
					for (int l = 0; l < n_quadrature_point; ++l)
					{
						double Jxw = quad_info.weight(l) * jacobian[l] * volume;

						/// 右端项还是零. 源项和 Neumann 条件.
						double rhs_cont = Jxw * basis_value_p[i][l] 
							* (vx_gradient[l][0] + vy_gradient[l][1]);
						rhs(2 * n_dof_v + element_dof_p[i]) += rhs_cont;
					}
				}
			}

		}

	
		/// 初始化未知量.
		Vector<double> x(n_total_dof);

		/// 边界条件处理.
		boundaryValueNS(x);


		/// debug
		// const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
		// const unsigned int * colnum = sp_stokes.get_column_numbers();
		// std::ofstream mat_deb;
		// mat_deb.open("mat_p.m", std::ofstream::out | std::ofstream::app);
		// mat_deb.setf(std::ios::fixed);
		// mat_deb.precision(20);
	    
		// mat_deb << "A = sparse(" << n_total_dof << ", " << n_total_dof << ");" << std::endl;
		// for (int i = 0; i < n_total_dof; ++i)
		// {
		// 	for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
		// 	{
		// 	    mat_deb << "A(" << i + 1 << ", " << colnum[j] + 1 << ")=" 
		// 		    << matrix.global_entry(j) << ";" << std::endl;
		// 	}
		// }
		// mat_deb.close();
		// std::cout << "mat output" << std::endl;
		// getchar();
		/// debug


		std::cout << "nonlinear res:" << std::endl;
		double revx = 0.0;
		for (int i = 0; i < n_dof_v; ++i)
			revx += rhs(i) * rhs(i);
		std::cout << "vx re: " << sqrt(revx) << std::endl;
		output << "vx re: " << sqrt(revx) << std::endl;
		double revy = 0.0;
		for (int i = 0; i < n_dof_v; ++i)
			revy += rhs(i + n_dof_v) * rhs(i + n_dof_v);
		std::cout << "vy re: " << sqrt(revy) << std::endl;
		output << "vy re: " << sqrt(revy) << std::endl;
		double rep = 0.0;
		for (int i = 0; i < n_dof_p; ++i)
			rep += rhs(i + 2 * n_dof_v) * rhs(i + 2 * n_dof_v);
		std::cout << "p re: " << sqrt(rep) << std::endl;
		output << "p re: " << sqrt(rep) << std::endl;
	
		double re = revx + revy +rep;
		std::cout << "total re: " << sqrt(re) << std::endl;
		output << "total re: " << sqrt(re) << std::endl;
		std::cout << "pause ..." << std::endl;
		// getchar();
		if (sqrt(re) < n_tol)
		{
			std::cout << "Covergence with residual: " << sqrt(re) 
				  << " in step " << iteration_times << std::endl;
			break;
		}
		std::cout << "Building precondition matrix ..." << std::endl;

		SparseMatrix<double> mat_Axx(sp_vxvx);
		SparseMatrix<double> mat_Ayy(sp_vyvy);
		SparseMatrix<double> mat_Wxy(sp_vyvx);
		SparseMatrix<double> mat_Wyx(sp_vxvy);
		SparseMatrix<double> mat_BTx(sp_pvx);
		SparseMatrix<double> mat_BTy(sp_pvy);

                //	updatePCDMatrix();

		for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
			mat_Axx.global_entry(i) = matrix.global_entry(index_vxvx[i]);
		for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
			mat_Ayy.global_entry(i) = matrix.global_entry(index_vyvy[i]);
		for (int i = 0; i < sp_vyvx.n_nonzero_elements(); ++i)
			mat_Wxy.global_entry(i) = matrix.global_entry(index_vyvx[i]);
		for (int i = 0; i < sp_vxvy.n_nonzero_elements(); ++i)
			mat_Wyx.global_entry(i) = matrix.global_entry(index_vxvy[i]);
		for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
			mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
		for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
			mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
		std::cout << "Precondition matrix builded!" << std::endl;

		// NSPreconditioner navierstokes_preconditioner;
		// navierstokes_preconditioner.initialize(mat_Axx, 
		// 				       mat_Ayy,
		// 				       mat_BTx,
		// 				       mat_BTy,
		// 				       mat_p_mass);
		std::cout << "Precondition initialized!" << std::endl;

		/// 矩阵求解. 
		dealii::SolverControl solver_control (4000000, l_tol, check);

		SolverGMRES<Vector<double> >::AdditionalData para(2000, false, true);
		SolverGMRES<Vector<double> > gmres (solver_control, para);
		std::cout << "Begin to solve linear system ..." << std::endl;
		// gmres.solve (matrix, x, rhs, navierstokes_preconditioner);
		gmres.solve (matrix, x, rhs, PreconditionIdentity());

		/// 调试块: 直接观测真实残量.
//	Vector<double> tmp(n_total_dof);
//	matrix.vmult(tmp, x);
//	tmp -= rhs;
//	std::cout << "linear residual: " << tmp.l2_norm() << std::endl;
//	getchar();

		FEMFunction<double, DIM> res_vx(fem_space_v);
		FEMFunction<double, DIM> res_vy(fem_space_v);
		FEMFunction<double, DIM> res_p(fem_space_p);

		/// 更新数值解.
		for (int i = 0; i < n_dof_v; ++i)
		{
			v_h[0](i) += x(i);
			res_vx(i) = x(i);
			v_h[1](i) += x(i + n_dof_v);
			res_vy(i) = x(i+ n_dof_v);
		}
		for (int i = 0; i < n_dof_p; ++i)
		{
			p_h(i) += x(i + 2 * n_dof_v);
			res_p(i) = x(i + 2 * n_dof_v);
		}

		double r_vx = Functional::L2Norm(res_vx, 1);
		double r_vy = Functional::L2Norm(res_vy, 1);
		double r_p = Functional::L2Norm(res_p, 1);
		/// 这个其实是更新...
		error_N = r_vx + r_vy + r_p;

		std::cout.setf(std::ios::fixed);
		std::cout.precision(20);
		std::cout << "updated vx: " << r_vx << std::endl;
		std::cout << "updated vy: " << r_vy << std::endl;
		std::cout << "updated p: " << r_p << std::endl;
		std::cout << "total updated: " << error_N << std::endl;
		std::cout << "step " << iteration_times << ", total updated: " << error_N 
			  << ", GMRES stpes: " << solver_control.last_step() << std::endl;

		output << "updated vx: " << r_vx << std::endl;
		output << "updated vy: " << r_vy << std::endl;
		output << "updated p: " << r_p << std::endl;
		output << "total updated: " << error_N << std::endl;
		output << "step " << iteration_times  
		       << ", GMRES stpes: " << solver_control.last_step() << std::endl;
		iteration_times++;
		if (iteration_times > 10)
		{
			std::cout << "Disconvergence at step 10." << std::endl;
			break;
		}
		output.close();
	}
};

void ISOP2P1::buildStokesSys()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;
    
	matrix.reinit(sp_stokes);

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i);
	/// (0, 2)
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvx[i]) = mat_pvx_divT.global_entry(i);

	/// (1, 2)
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvy[i]) = mat_pvy_divT.global_entry(i);

	/// (2, 0)
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxp[i]) =  mat_vxp_div.global_entry(i);
	
	/// (2, 1)
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyp[i]) =  mat_vyp_div.global_entry(i);

};


void ISOP2P1::buildNewtonSys4NS()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i) 
			+ mat_v_convection.global_entry(i)
			+ mat_v_Jacobi_xx.global_entry(i);
	/// (0, 1)
	for (int i = 0; i < sp_vyvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvx[i]) = mat_v_Jacobi_xy.global_entry(i);
	/// (1, 0)
	for (int i = 0; i < sp_vxvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvy[i]) = mat_v_Jacobi_yx.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i)
			+ mat_v_Jacobi_yy.global_entry(i);
};


void ISOP2P1::buildPicardSys4NS()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i) 
			+ mat_v_convection.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i);
};

#undef DIM
