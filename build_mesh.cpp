#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::buildMesh()
{
    std::cout << "The mesh tree is " << mesh_file << std::endl;
    /// 读入网格.
    h_tree.readEasyMesh(mesh_file);

    irregular_mesh_p = new IrregularMesh<DIM>;
    /// 产生宏单元网格.
    irregular_mesh_p->reinit(h_tree);
    irregular_mesh_p->semiregularize();
    irregular_mesh_p->regularize(false);
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();

    /// 初始化网格.
    std::cout << "Initialize mesh ... " << std::endl;
    readDomain(mesh_file);

    irregular_mesh_v = new IrregularMesh<DIM>(*irregular_mesh_p);
    /// 产生单元网格.
    irregular_mesh_v->globalRefine(1);
    irregular_mesh_v->semiregularize();
    irregular_mesh_v->regularize(false);
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    
    unsigned int n_ele_v = mesh_v.n_geometry(2);
    unsigned int n_ele_p = mesh_p.n_geometry(2);
    index_v2p.resize(n_ele_v);
    index_p2v.resize(n_ele_p);
    IrregularMesh<DIM, DIM>::RootIterator root_iterator = irregular_mesh_v->beginRootElement();
    for (; root_iterator != irregular_mesh_v->endRootElement(); ++root_iterator)
    {
	int ele_p_idx = root_iterator->index;
	int n_chi = root_iterator->n_child;
	index_p2v[ele_p_idx].resize(n_chi);
        for (unsigned int i = 0; i < root_iterator->n_child; ++i)
    	{
	    HElement<DIM, DIM> *chi = root_iterator->child[i]; 
	    int ele_v_idx = chi->index;
	    index_p2v[ele_p_idx][i] = ele_v_idx;
	    index_v2p[ele_v_idx] = ele_p_idx;
	}
    }
};

#undef DIM
