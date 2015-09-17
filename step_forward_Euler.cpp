#include "ISOP2P1.h"
#include "functions.h"
#include "preconditioner.h"


#define DIM 2

void ISOP2P1::stepForwardEuler()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// 系数矩阵直接使用 Stokes 矩阵结构.
	matrix.reinit(sp_stokes);
	std::cout << "dt * viscosity = " << dt * viscosity << std::endl;
	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = dt * viscosity * mat_v_stiff.global_entry(i)
		+ mat_v_mass.global_entry(i);
	/// (1, 1) 这两个对角块对应扩散算子和质量算子, dt 直接乘上, 以避免
	/// 矩阵系数过大.
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = dt * viscosity * mat_v_stiff.global_entry(i)
		+ mat_v_mass.global_entry(i);
	/// (0, 2) 这个不是方阵. 在矩阵结构定义的时候已经直接排除了对角元优
	/// 先.
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvx[i]) = dt * mat_pvx_divT.global_entry(i);

	/// (1, 2)
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvy[i]) = dt * mat_pvy_divT.global_entry(i);

	/// (2, 0)
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxp[i]) = dt * mat_vxp_div.global_entry(i);

	/// (2, 1) 这四块直接复制散度矩阵.
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyp[i]) = dt * mat_vyp_div.global_entry(i);

	/// 问题右端项.
	rhs.reinit(n_total_dof);

	FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
	FEMSpace<double, DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_p = fem_space_p.endElement();

	/// 遍历速度单元, 拼装相关系数矩阵和右端项.
	for (the_element_v = fem_space_v.beginElement();
			the_element_v != end_element_v; ++the_element_v)
	{
		/// 当前单元信息.
		double volume = the_element_v->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
		/// 装时积分精度不用超过 1 次. (验证一下!)
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
		/// 速度单元信息.
		std::vector<std::vector<double> > basis_value_v = the_element_v->basis_function_value(q_point);
		std::vector<std::vector<std::vector<double> > > basis_gradient_v = the_element_v->basis_function_gradient(q_point);
		std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
		std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
		std::vector<double> fx_value = source_v[0].value(q_point, *the_element_v);
		std::vector<double> fy_value = source_v[1].value(q_point, *the_element_v);
		std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
		std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
		const std::vector<int>& element_dof_v = the_element_v->dof();
		int n_element_dof_v = the_element_v->n_dof();
		Element<double, DIM> &p_element = fem_space_p.element(index_v2p[the_element_v->index()]);
		const std::vector<int>& element_dof_p = p_element.dof();
		std::vector<std::vector<std::vector<double> > > basis_gradient_p = p_element.basis_function_gradient(q_point);
		std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
		int n_element_dof_p = p_element.n_dof();
		std::vector<double> p_value = p_h.value(q_point, p_element);
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			for (int i = 0; i < n_element_dof_v; ++i)
			{
				double rhs_cont = fx_value[l] * basis_value_v[i][l] + vx_value[l] * basis_value_v[i][l];
				rhs_cont -= dt * (vx_value[l] * vx_gradient[l][0] +
						vy_value[l] * vx_gradient[l][1]) * basis_value_v[i][l];
				rhs_cont *= Jxw;
				rhs(element_dof_v[i]) += rhs_cont;

				rhs_cont = fy_value[l] * basis_value_v[i][l] + vy_value[l] * basis_value_v[i][l];
				rhs_cont -= dt * (vx_value[l] * vy_gradient[l][0] +
						vy_value[l] * vy_gradient[l][1]) * basis_value_v[i][l];
				rhs_cont *= Jxw;
				rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
			}
		}
	}

	/// 构建系数矩阵和右端项.
	/// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
	Vector<double> x(n_total_dof);

	// for (int i = 0; i < n_dof_v; ++i)
	// {
	// 	x(i) = v_h[0](i);
	// 	x(i + n_dof_v) = v_h[1](i);
	// }

	// for (int i = 0; i < n_dof_p; ++i)
	// 	x(i + 2 * n_dof_v) = p_h(i);

	/// 边界条件一起处理了. 这里需要传递 x 因为 x 是临时的. 这里似乎应
	/// 该把 v_h, p_h 和 x 统一起来, 避免冗余错误.
	boundaryValueStokes(x);

	clock_t t_cost = clock();
	/// 矩阵求解. 
	dealii::SolverControl solver_control (4000000, l_tol, check);
	SolverMinRes<Vector<double> > minres (solver_control);
	
	minres.solve (matrix, x, rhs, PreconditionIdentity());
	
	for (int i = 0; i < n_dof_v; ++i)
	{
	    v_h[0](i) = x(i);
	    v_h[1](i) = x(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
	    p_h(i) = x(i + 2 * n_dof_v);

	// /// debug
	// const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
	// const unsigned int * colnum = sp_stokes.get_column_numbers();
	// std::ofstream mat_deb;
	// mat_deb.open("mat.m", std::ofstream::out);
	// mat_deb.setf(std::ios::fixed);
	// mat_deb.precision(20);
	    
	// mat_deb << "A = sparse(" << n_total_dof << ", " << n_total_dof << ");" << std::endl;
	// for (int i = 0; i < n_total_dof; ++i)
	// {
	//     for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
	//     {
	// 	mat_deb << "A(" << i + 1 << ", " << colnum[j] + 1 << ")=" 
	// 		<< matrix.global_entry(j) << ";" << std::endl;
	//     }
	//     mat_deb << "x(" << i + 1<< ") = " << x(i) << ";" << std::endl;
	//     mat_deb << "rhs(" << i + 1 << ") = " << rhs(i) << ";" << std::endl;
	// }
	// mat_deb.close();
	// std::cout << "mat output" << std::endl;
	Vector<double> res(n_total_dof);
	matrix.vmult(res, x);
	res *= -1;
	res += rhs;
	std::cout << "res_l2norm =" << res.l2_norm() << std::endl;
	/// debug

	/// 记录计算结果参数.
	std::ofstream output;
	output.open("record", std::ofstream::out | std::ofstream::app);
	output.setf(std::ios::fixed);
	output.precision(20);
	output << "nu = " << n_dof_v << std::endl;
	output << "np = " << n_dof_p << std::endl;

	if (error_check == true)
	{
	    // RealVx real_Vx;
	    // RealVy real_Vy;    
	    // double mean_p_h= Functional::meanValue(p_h, 3); 
	    // RealP real_P(mean_p_h);
		PoiseuilleVx accuracy_vx(-1.0, 1.0);
		PoiseuilleVy accuracy_vy;
		PoiseuilleP poiseuille_p(0.0, viscosity);

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
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(3);
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
	    error = Functional::L2Error(v_h[0], accuracy_vx, 3);
	    std::cout << "|| u - u_h ||_L2 = " << error << std::endl;

	    error = Functional::H1SemiError(v_h[0], accuracy_vx, 3);
	    std::cout << "|| u - u_h ||_H1 = " << error << std::endl;
	    error = Functional::L2Error(p_h, poiseuille_p, 3);
	    std::cout << "|| p - p_h ||_L2 = " << error << std::endl;

	    std::cout << "uh_L2err() = " << L2_err << std::endl;
	    std::cout << "uh_H1err() = " << H1_err << std::endl;	
	    output << "uh_L2err(1) = " << L2_err << std::endl;
	    output << "uh_H1err(1) = " << H1_err << std::endl;	
	}
	output.close();
	double mean_p_h = Functional::meanValue(p_h, 3);
	FEMFunction<double, DIM> _p_h(p_h);
	mean_p_h = -mean_p_h;
	_p_h.add(mean_p_h);
	std::cout << "|| p - mean_p_h ||_L2 = " << _p_h.l2_norm() << std::endl;
	
};

#undef DIM
