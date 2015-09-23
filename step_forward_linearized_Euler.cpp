#include "ISOP2P1.h"
#include "preconditioner.h"
#define DIM 2

void ISOP2P1::stepForwardLinearizedEuler()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = DIM * n_dof_v + n_dof_p;

	matrix.reinit(sp_stokes);
	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ (1.0 / dt) * mat_v_mass.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ (1.0 / dt) * mat_v_mass.global_entry(i);
	// /// (0, 2)
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvx[i]) = mat_pvx_divT.global_entry(i);

	// /// (1, 2)
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvy[i]) = mat_pvy_divT.global_entry(i);

	// /// (2, 0)
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxp[i]) = mat_vxp_div.global_entry(i);
	
	// /// (2, 1)
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyp[i]) = mat_vyp_div.global_entry(i);

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
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(4);
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
				for (int j = 0; j < n_element_dof_v; ++j)
				{
					double cont = Jxw * (vx_value[l] * basis_gradient_v[j][l][0] +
							     vy_value[l] * basis_gradient_v[j][l][1]) * basis_value_v[i][l];
					matrix.add(element_dof_v[i], element_dof_v[j], cont);
					matrix.add(n_dof_v + element_dof_v[i], n_dof_v + element_dof_v[j], cont);
				}
				/// 右端项. 这里可施加源项和 Neumann 条件.
				double rhs_cont = fx_value[l] * basis_value_v[i][l] + (1.0 / dt) * vx_value[l] * basis_value_v[i][l];
				rhs_cont *= Jxw;
				rhs(element_dof_v[i]) += rhs_cont;
				rhs_cont = fy_value[l] * basis_value_v[i][l] + (1.0 / dt ) * vy_value[l] * basis_value_v[i][l];
				rhs_cont *= Jxw;
				rhs(n_dof_v + element_dof_v[i]) += rhs_cont;

			}
		}
	}
	/// 构建系数矩阵和右端项.
	/// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
	Vector<double> x(n_total_dof);

	for (int i = 0; i < n_dof_v; ++i)
	{
		x(i) = v_h[0](i);
		x(i + n_dof_v) = v_h[1](i);
	}

	for (int i = 0; i < n_dof_p; ++i)
		x(i + 2 * n_dof_v) = p_h(i);

	boundaryValueStokes(x);

	/// 矩阵求解. 
	SparseMatrix<double> mat_Axx(sp_vxvx);
	SparseMatrix<double> mat_Ayy(sp_vyvy);
	SparseMatrix<double> mat_BTx(sp_pvx);
	SparseMatrix<double> mat_BTy(sp_pvy);

	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		mat_Axx.global_entry(i) = matrix.global_entry(index_vxvx[i]); 
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		mat_Ayy.global_entry(i) = matrix.global_entry(index_vyvy[i]);
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
	std::cout << "Precondition matrix builded!" << std::endl;
	
	// NSPreconditioner preconditioner_ns;
	// preconditioner_ns.initialize(mat_Axx, 
	// 				 mat_Ayy,
	// 				 mat_BTx,
	// 				 mat_BTy,
	// 				 mat_p_mass);

	clock_t t_cost = clock();
	dealii::SolverControl solver_control (4000000, l_tol, check);

	SolverGMRES<Vector<double> >::AdditionalData para(2000, false, true);
	SolverGMRES<Vector<double> > gmres (solver_control, para);
	gmres.solve (matrix, x, rhs, PreconditionIdentity());
	t_cost = clock() - t_cost;

	std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;

	/// 将整体数值解分割成速度和压力.
	for (int i = 0; i < n_dof_v; ++i)
	{
		v_h[0](i) = x(i);
		v_h[1](i) = x(i + n_dof_v);
	}

	for (int i = 0; i < n_dof_p; ++i)
		p_h(i) = x(i + 2 * n_dof_v);
};

#undef DIM
