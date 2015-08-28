#include "ISOP2P1.h"

#define DIM 2

void ISOP2P1::time_step()
{
    int n_dof_v = fem_space_v.n_dof(); 
    int n_dof_p = fem_space_p.n_dof(); 

    double min_v_norm = 10;

    FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
    FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();

    for (; the_element_v != end_element_v; ++the_element_v)
    {
	/// 当前单元信息.
	double volume = the_element_v->templateElement().volume();
	/// 算个面积, 积分精度 0 够了.
	const QuadratureInfo<2>& quad_info = the_element_v->findQuadratureInfo(0);
	std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<2> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
	/// 速度单元信息.
	std::vector<std::vector<double> > basis_value_v = the_element_v->basis_function_value(q_point);
	std::vector<std::vector<std::vector<double> > > basis_gradient_v = the_element_v->basis_function_gradient(q_point);
	std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
	std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
	
	// /// 获取单元最小边长.
	// const GeometryBM& geometry = the_element_v->geometry();
	// for (int i = 0; i < geometry.n_boundary(); ++i) 
	// {
	//     int j = geometry.boundary(i);
	//     const GeometryBM& side = ele_mesh.geometry(1, j);
	//     const Point<DIM>& p0 = ele_mesh.point(ele_mesh.geometry(0, side.boundary(0)).vertex(0));
	//     const Point<DIM>& p1 = ele_mesh.point(ele_mesh.geometry(0, side.boundary(1)).vertex(0));
	//     double side_length = fabs(distance(p0, p1));

	//     if (min_side_length > side_length)
	// 	min_side_length = side_length;
	// }
	const std::vector<int>& element_dof_v = the_element_v->dof();
	int n_element_dof_v = the_element_v->n_dof();
	for (int i =0; i < n_element_dof_v; ++i)
	{
	    /// local_v_norm中有单元面积开平方,其实这样获得的时间步长dt是与网格的长度同阶的,
	    /// 因此不需要再乘上最小的网格长度.
	    double local_v_norm 
	    	= sqrt(2.0 * jacobian[0] * volume 
	    	       / (vx_value[0] * vx_value[0] +  vy_value[0] * vy_value[0] + eps));
	    if (min_v_norm > local_v_norm)
		min_v_norm = local_v_norm;
	}
    }
    /// 时间步长满足CFL条件.
    if (t < 1)
    	dt = 0.1 * min_v_norm;
    else
	dt = 0.1 * min_v_norm;
    std::cout << "dt * max_v / min_h < 0.5 " << std::endl; 
    std::cout<< "=> dt = 0.5 * min_h / max_v -> dt = " << dt << std::endl;
};

#undef DIM
