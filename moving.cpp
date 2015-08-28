#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"

#define DIM 2

void ISOP2P1::syncMesh()
{
    int n_point = this->n_point();
    const double & msl = moveStepLength();
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
    /// 输出一下.
    outputTecplotP("P0");
    outputTecplot("NS_Euler");

    FEMFunction<double, DIM> _u_h(fem_space_p);
    FEMFunction<double, DIM> _v_h(fem_space_p);
    Operator::L2Interpolate(v_h[0], _u_h);
    Operator::L2Interpolate(v_h[1], _v_h);
    
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space_p.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space_p.endElement();
    for (int i = 0; the_element != end_element; ++the_element, ++i) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(1);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	std::vector<std::vector<double> > u_h_gradient = _u_h.gradient(q_point, *the_element);
	std::vector<std::vector<double> > v_h_gradient = _v_h.gradient(q_point, *the_element);
	float d = 0, area = 0;
	for (int l = 0; l < n_quadrature_point; ++l) 
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    area += Jxw;
	    d += Jxw * (innerProduct(u_h_gradient[l], u_h_gradient[l]) + innerProduct(v_h_gradient[l], v_h_gradient[l]));
	}
	monitor(i) = d / area;
    }
    std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
	      << "\tmin monitor=" << *std::min_element(monitor().begin(), monitor().end())
	      << std::endl;
    double max_monitor = *std::max_element(monitor().begin(), monitor().end());
    smoothMonitor(5);
    for (int i = 0; i < n_geometry(2); ++i)
	monitor(i) = 1.0 / sqrt(1.0 + alpha * monitor(i));

};

void ISOP2P1::updateSolution()
{
    fem_space_p.updateDofInterpPoint();
    fem_space_v.updateDofInterpPoint();

    FEMFunction<double, DIM> _u_h(v_h[0]);
    FEMFunction<double, DIM> _v_h(v_h[1]);
    FEMFunction<double, DIM> _p_h(p_h);
    const double& msl = moveStepLength();
    
    /// 重新构造矩阵.
    buildMatrixStruct();
    buildMatrix();
    /// 虚拟时间.
    double _dt = 1.0;
    
    int n_dof_v = fem_space_v.n_dof();
    int n_dof_p = fem_space_p.n_dof();
    int n_total_dof = 2 * n_dof_v + n_dof_p;
    int n_total_dof_v = 2 * n_dof_v;

    for (int m = 1; m > 0; --m)
    {
	/// 系数矩阵直接使用 Stokes 矩阵结构.
	matrix.reinit(sp_stokes);

	/// (0, 0) 
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_vxvx[i]) = mat_v_mass.global_entry(i); 
	/// (1, 1) 这两个对角块对应扩散算子和质量算子, dt 直接乘上, 以避免
	/// 矩阵系数过大.
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_vyvy[i]) = mat_v_mass.global_entry(i);
	/// (0, 2) 这个不是方阵. 在矩阵结构定义的时候已经直接排除了对角元优
	/// 先.
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_pvx[i]) = (_dt / m) * mat_pvx_divT.global_entry(i);

	/// (1, 2)
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_pvy[i]) = (_dt / m) * mat_pvy_divT.global_entry(i);

	/// (2, 0)
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_vxp[i]) = (_dt / m) * mat_vxp_div.global_entry(i);
	
	/// (2, 1) 这四块直接复制散度矩阵. 
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
	    matrix.global_entry(index_vyp[i]) = (_dt / m) * mat_vyp_div.global_entry(i);

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
	    /// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数.
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
	    /// 速度积分点上的移动方向.
	    std::vector<std::vector<double> > move_vector = moveDirection(q_point, index_v2p[the_element_v->index()]);
	    for (int l = 0; l < n_quadrature_point; ++l)
	    {
		double Jxw = quad_info.weight(l) * jacobian[l] * volume;
		for (int i = 0; i < n_element_dof_v; ++i)
		{
		    double rhs_cont = vx_value[l] * basis_value_v[i][l];
		    rhs_cont += (_dt / m) * msl * innerProduct(move_vector[l], vx_gradient[l]) * basis_value_v[i][l];
		    rhs_cont *= Jxw;
		    rhs(element_dof_v[i]) += rhs_cont;

		    rhs_cont = vy_value[l] * basis_value_v[i][l];
		    rhs_cont += (_dt / m) * msl * innerProduct(move_vector[l], vy_gradient[l]) * basis_value_v[i][l];
		    rhs_cont *= Jxw;
		    rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
		}
	    }
	}
	// /// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
	// AccuracyVx accuracy_vx(viscosity, t);
	// AccuracyVy accuracy_vy(viscosity, t);

	// Operator::L2Project(accuracy_vx, v_h[0], Operator::LOCAL_LEAST_SQUARE, 3);
	// Operator::L2Project(accuracy_vy, v_h[1], Operator::LOCAL_LEAST_SQUARE, 3);
	/// 构建系数矩阵和右端项.
	Vector<double> x(n_total_dof);
	for (int i = 0; i < n_dof_v; ++i)
	{
	    x(i) = v_h[0](i);
	    x(i + n_dof_v) = v_h[1](i);
	}
	
	for (int i = 0; i < n_dof_p; ++i)
	    x(i + 2 * n_dof_v) = p_h(i);
	/// 边界处理.
	boundaryValueStokes(x);

	clock_t t_cost = clock();
	/// 矩阵求解. 
	dealii::SolverControl solver_control (400000, l_tol, check);
	SolverMinRes<Vector<double> > minres (solver_control);
	
	minres.solve (matrix, x, rhs, PreconditionIdentity());

	t_cost = clock() - t_cost;	
	std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;
	for (int i = 0; i < n_dof_v; ++i)
	{
	    v_h[0](i) = x(i);
	    v_h[1](i) = x(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
	    p_h(i) = x(i + 2 * n_dof_v);
	
	/// debug
	const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
	const unsigned int * colnum = sp_stokes.get_column_numbers();
	std::ofstream mat_deb;
	mat_deb.open("mat.m", std::ofstream::out);
	mat_deb.setf(std::ios::fixed);
	mat_deb.precision(20);
	    
	mat_deb << "A = sparse(" << n_total_dof << ", " << n_total_dof << ");" << std::endl;
	for (int i = 0; i < n_total_dof; ++i)
	{
	    for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
	    {
		mat_deb << "A(" << i + 1 << ", " << colnum[j] + 1 << ")=" 
			<< matrix.global_entry(j) << ";" << std::endl;
	    }
	    mat_deb << "x(" << i + 1<< ") = " << x(i) << ";" << std::endl;
	    mat_deb << "rhs(" << i + 1 << ") = " << rhs(i) << ";" << std::endl;
	}
	mat_deb.close();
	std::cout << "mat output" << std::endl;
	Vector<double> res(n_total_dof);
	matrix.vmult(res, x);
	res *= -1;
	res += rhs;
	std::cout << "res_l2norm =" << res.l2_norm() << std::endl;
	/// debug
    }
    
    /// 为了调试, 输出一下数值解.
    outputSolution();
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

