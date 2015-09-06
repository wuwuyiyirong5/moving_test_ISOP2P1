#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::outputTecplot(const std::string &prefix)
{
    if (output_vorticity == true)
	computVorticity();
    if (output_divergence == true)
	computDivergence();

    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    int n_node = mesh_v.n_geometry(0);
    int n_ele = mesh_v.n_geometry(2);

    FEMFunction <double, DIM> p_h_refine(fem_space_v);

    Operator::L2Interpolate(p_h, p_h_refine);
    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(4);
    result << prefix << ".dat";
    std::ofstream tecplot(result.str().c_str()); 
    tecplot.setf(std::ios::fixed);
    tecplot.precision(20);
    tecplot << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\"";
    if (output_vorticity == true) 
	tecplot << ", \"W\"";
    if (output_divergence == true)
	    tecplot << ", \"D\"";
    tecplot << std::endl;
    tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
    tecplot << "ZONETYPE=FETRIANGLE" << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_v.point(i)[0] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_v.point(i)[1] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << p_h_refine(i) << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << v_h[0](i) << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << v_h[1](i) << "\n";
    tecplot << std::endl;
    if (output_vorticity == true)
    {
	for (int i = 0; i < n_node; ++i)
		tecplot << vot(i) << "\n";
	tecplot << std::endl;
    }
    if (output_divergence == true)
	    for (int i = 0; i < n_node; ++i)
		    tecplot << div(i) << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_ele; ++i)
    {
	std::vector<int> &vtx =  fem_space_v.element(i).geometry().vertex();
        tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
    }
    tecplot.close();
};

void ISOP2P1::outputTecplotP(const std::string &prefix)
{
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
    int n_node = mesh_p.n_geometry(0);
    int n_ele = mesh_p.n_geometry(2);

    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(4);
    result << prefix << ".dat";
    std::ofstream tecplot(result.str().c_str()); 
    tecplot.setf(std::ios::fixed);
    tecplot.precision(20);
    tecplot << "VARIABLES = \"X\", \"Y\", \"P\"";
    tecplot << std::endl;
    tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
    tecplot << "ZONETYPE=FETRIANGLE" << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_p.point(i)[0] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_p.point(i)[1] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << p_h(i) << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_ele; ++i)
    {
	std::vector<int> &vtx =  fem_space_p.element(i).geometry().vertex();
        tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
    }
    tecplot.close();
};

void ISOP2P1::computVorticity()
{
    int n_dof_v = fem_space_v.n_dof();
    vot.reinit(fem_space_v);
 
    /// 准备一个遍历全部单元的迭代器. 包括 v 和 p .
    FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
    Vector<double> rhs_vot(n_dof_v);

    /// 遍历速度单元, 将压力值插值到速度空间.
    for (the_element_v = fem_space_v.beginElement(); 
    	 the_element_v != end_element_v; ++the_element_v) 
    {
    	/// 当前单元信息.
    	double volume = the_element_v->templateElement().volume();
    	/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
    	/// 装时积分精度不用超过 1 次. (验证一下!)
    	const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(2);
    	std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
    	int n_quadrature_point = quad_info.n_quadraturePoint();
    	std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
    	const std::vector<int>& element_dof_v = the_element_v->dof();
    	int n_element_dof_v = the_element_v->n_dof();
    	std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
    	std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);

    	for (int l = 0; l < n_quadrature_point; ++l)
    	{
    	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
    	    for (int i = 0; i < n_element_dof_v; ++i)
    	    {
    		double cont = (vy_gradient[l][0] - vx_gradient[l][1]) * Jxw; 
    		rhs_vot(element_dof_v[i]) += cont;
    	    }
    	}
    }
    AMGSolver solver(mat_v_mass);
    solver.solve(vot, rhs_vot, 1.0e-08, 200);	
};

void ISOP2P1::computDivergence()
{
    int n_dof_v = fem_space_v.n_dof();
    div.reinit(fem_space_v);
 
    /// 准备一个遍历全部单元的迭代器. 包括 v 和 p .
    FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
    Vector<double> rhs_div(n_dof_v);

    /// 遍历速度单元, 将压力值插值到速度空间.
    for (the_element_v = fem_space_v.beginElement(); 
    	 the_element_v != end_element_v; ++the_element_v) 
    {
    	/// 当前单元信息.
    	double volume = the_element_v->templateElement().volume();
    	/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
    	/// 装时积分精度不用超过 1 次. (验证一下!)
    	const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(2);
    	std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
    	int n_quadrature_point = quad_info.n_quadraturePoint();
    	std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
    	const std::vector<int>& element_dof_v = the_element_v->dof();
    	int n_element_dof_v = the_element_v->n_dof();
    	std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
    	std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);

    	for (int l = 0; l < n_quadrature_point; ++l)
    	{
    	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
    	    for (int i = 0; i < n_element_dof_v; ++i)
    	    {
    		double cont = (vx_gradient[l][0] + vy_gradient[l][1]) * Jxw; 
    		rhs_div(element_dof_v[i]) += cont;
    	    }
    	}
    }
    AMGSolver solver(mat_v_mass);
    solver.solve(div, rhs_div, 1.0e-08, 200);
    std::cout << "divergence norm : " << div.l2_norm() << std::endl;
};

#undef DIM
