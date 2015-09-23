#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::buildFEMSpace()
{
    /// 组装有限元空间, 同时给宏单元网格也组装了一个空间. 暂时没有用.

    /// 参考几何体一致, 为三角形剖分.
    template_geometry.readData("triangle.tmp_geo");
    coord_transform.readData("triangle.crd_trs");
    template_dof.reinit(template_geometry);
    
    /// 均为分片线性三角插值. 即 P1 空间.
    template_dof.readData("triangle.1.tmp_dof"); 
    basis_function.reinit(template_dof);
    basis_function.readData("triangle.1.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry, 
    				 template_dof,
    				 coord_transform, 
    				 basis_function);
    /// 本地引用正则化网格. 正则化网格只能本地引用. 它实际存储在非正则
    /// 化网格中.
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    
    /// 压力和速度, 根据各自的网格, 建立各自的空间.
    fem_space_v.reinit(mesh_v, template_element);
    fem_space_p.reinit(mesh_p, template_element);

    /// 建立各自的有限元空间.
    unsigned int n_ele = mesh_v.n_geometry(DIM);
    unsigned int n_vtx_v = mesh_v.n_geometry(0);
    unsigned int n_macro = mesh_p.n_geometry(2);
    unsigned int n_vtx_p = mesh_p.n_geometry(0);
    fem_space_v.element().resize(n_ele);
    fem_space_p.element().resize(n_macro);
    for (int i = 0; i < n_ele; ++i)
    	fem_space_v.element(i).reinit(fem_space_v, i, 0);
    for (int i = 0; i < n_macro; ++i)
    	fem_space_p.element(i).reinit(fem_space_p, i, 0);
    fem_space_v.buildElement();
    fem_space_v.buildDof();
    fem_space_v.buildDofBoundaryMark();
    fem_space_p.buildElement();
    fem_space_p.buildDof();
    fem_space_p.buildDofBoundaryMark();

    /// 建立数值解, 源项向量
    v_h.resize(DIM);
    for (int i = 0; i < DIM; ++i)
	    v_h[i].reinit(fem_space_v);
    p_h.reinit(fem_space_p);
    source_v.resize(DIM);
    for (int i = 0; i < DIM; ++i)
	    source_v[i].reinit(fem_space_v);
    source_p.reinit(fem_space_p);

    /// 这里施加外力项. 
    ForceX force_x(body_force, angle);
    ForceY force_y(body_force, angle);
    Operator::L2Interpolate(force_x, source_v[0]);
    Operator::L2Interpolate(force_y, source_v[1]);
};

#undef DIM
