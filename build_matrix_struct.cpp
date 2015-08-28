#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::buildMatrixStruct()
{
    /// 构建系数矩阵模板.
    /// 速度空间自由度.
    int n_dof_v = fem_space_v.n_dof();
    /// 压力空间自由度.
    int n_dof_p = fem_space_p.n_dof();
    /// 总自由度.
    int n_total_dof = 2 * n_dof_v + n_dof_p;
    /// 准备统计系数矩阵的每一行有多少个非零元.
    std::vector<unsigned int> n_non_zero_per_row(n_total_dof);
    /// 设置各块大小.
    int row_vx = n_dof_v;	/**< (0, 0) */
    int col_vx = n_dof_v;
    std::vector<unsigned int> n_non_zero_per_row_vxvx(n_dof_v);
    int col_vy = n_dof_v;       /**< (0, 1) */
    std::vector<unsigned int> n_non_zero_per_row_vyvx(n_dof_v);
    int col_p = n_dof_p;        /**< (0, 2) */
    std::vector<unsigned int> n_non_zero_per_row_pvx(n_dof_v);
    int row_vy = n_dof_v;	/**< (1, 0) */
    std::vector<unsigned int> n_non_zero_per_row_vxvy(n_dof_v); 
    /**< (1, 1) */
    std::vector<unsigned int> n_non_zero_per_row_vyvy(n_dof_v); 
    /**< (1, 2) */
    std::vector<unsigned int> n_non_zero_per_row_pvy(n_dof_v);
    int row_p = n_dof_p;	/**< (2, 0) */
    std::vector<unsigned int> n_non_zero_per_row_vxp(n_dof_p); 
    /**< (2, 1) */
    std::vector<unsigned int> n_non_zero_per_row_vyp(n_dof_p); 
    /**< (2, 2) */
    std::vector<unsigned int> n_non_zero_per_row_penalty(n_dof_p);
    /// 压力空间质量矩阵块. 用于预处理. 
    std::vector<unsigned int> n_non_zero_per_row_mass_p(n_dof_p);

    /// 准备一个遍历全部单元的迭代器. 包括 v 和 p .
    FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
    FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();
    FEMSpace<double,2>::ElementIterator the_element_p = fem_space_p.beginElement();
    FEMSpace<double,2>::ElementIterator end_element_p = fem_space_p.endElement();

    /// 第一次循环遍历全部单元, 只是为了统计每一行的非零元个数.
    for (; the_element_v != end_element_v; ++the_element_v) 
    {
	/// 在每一速度单元中, 不但需要和自己的局部自由度拼装, 也需要和
	/// 包含自己的压力单元中的全部自由度拼装. 用这个方式简化基础拼
	/// 装. 现在的拼装就是线性三角单元的标准拼装, 唯一的区别就是需
	/// 要有一个压力和速度单元匹配, 这个在生成有限元空间已经做了.
	/// 实际上的单元应该是压力单元, 里面包含四个速度单元也可以看作
	/// 独立的六个自由度, 这样自由度总数和 Taylor-Hood 元是一样的. 
    	const std::vector<int>& element_dof = the_element_v->dof();
    	/// 由索引找到对应的压力空间单元.
    	Element<double, DIM> &p_element = fem_space_p.element(index_v2p[the_element_v->index()]);
    	const std::vector<int>& element_dof_p = p_element.dof();
    	int n_element_dof_v = the_element_v->n_dof();
    	int n_element_dof_p = p_element.n_dof();
	/// j 是检测函数指标. 行指标. 在整个计算中, 默认行是检测函数而
	/// 列是基函数.
    	for (int j = 0; j < n_element_dof_v; ++j)
    	{
	    /// k 是基函数指标. 列指标.
    	    for (int k = 0; k < n_element_dof_v; ++k)
    	    {
    		/// element_dof[j] 行两个.
    		n_non_zero_per_row[element_dof[j]] += 2;
    		/// 一个在 (0, 0) 块.
    		n_non_zero_per_row_vxvx[element_dof[j]]++;
    		/// 一个在 (0, 1) 块.
    		n_non_zero_per_row_vyvx[element_dof[j]]++;
    		/// n_dof_v + element_dof[j] 行两个.
    		n_non_zero_per_row[element_dof[j] + n_dof_v] += 2;
    		/// 这个在 (1, 0) 块.
    		n_non_zero_per_row_vxvy[element_dof[j]]++;
    		/// 这个在 (1, 1) 块.
    		n_non_zero_per_row_vyvy[element_dof[j]]++;
    	    }
	    /// 和 j 行有关系的列是套住 j 所在单元的那个宏单元的全部自由度.
    	    for (int k = 0; k < n_element_dof_p; ++k)
    	    { 
    		/// element_dof[j] 行一个.
    		n_non_zero_per_row[element_dof[j]]++;
    		/// 这个在 (0, 2) 块.
    		n_non_zero_per_row_pvx[element_dof[j]]++;
    		/// n_dof_v + element_dof[j] 行一个. 
    		n_non_zero_per_row[element_dof[j] + n_dof_v]++;
    		/// 这个在 (1, 2) 块.
    		n_non_zero_per_row_pvy[element_dof[j]]++;
    	    }
    	}
    }

    for (; the_element_p != end_element_p; ++the_element_p)
    {
    	const std::vector<int>& element_dof = the_element_p->dof();
    	int n_element_dof_p = the_element_p->n_dof();
    	for (int j = 0; j < n_element_dof_p; ++j)  
    	{    
	    int idx_p = the_element_p->index();
	    int n_chi = index_p2v[idx_p].size();
	    for (int i = 0; i < n_chi; ++i)
	    {
		/// 由索引找到对应的速度空间单元.
		Element<double, DIM> &v_element = fem_space_v.element(index_p2v[idx_p][i]);
		const std::vector<int>& element_dof_v = v_element.dof();
		int n_element_dof_v = v_element.n_dof();
		for (int k = 0; k < n_element_dof_v; ++k)
		{
		    /// 2 * n_dof_v + element_dof[j] 行 2 个. 
		    n_non_zero_per_row[element_dof[j] + 2 * n_dof_v] += 2;
		    /// 一个在 (2, 0) 块.
		    n_non_zero_per_row_vxp[element_dof[j]]++;
		    /// 另一个在 (2, 1) 块.
		    n_non_zero_per_row_vyp[element_dof[j]]++;
		}
    	    }
    	}
    }
    /// 用于预处理的压力空间质量矩阵.
    for (the_element_p = fem_space_p.beginElement(); 
	 the_element_p != end_element_p; ++the_element_p)
    {
    	const std::vector<int>& element_dof = the_element_p->dof();
    	int n_element_dof_p = the_element_p->n_dof();
    	for (int j = 0; j < n_element_dof_p; ++j)  
	    for (int k = 0; k < n_element_dof_p; ++k)
		n_non_zero_per_row_mass_p[element_dof[j]]++;
    }

    /// 给右下角对角元留带宽. 因为对角元实际是 0 .
    for (int i = 0; i < n_dof_p; ++i)
    {
	n_non_zero_per_row[i + 2 * n_dof_v]++;
	n_non_zero_per_row_penalty[i]++;
    }

    /// 指定矩阵模板带宽.
    sp_stokes.reinit(n_total_dof, n_total_dof, n_non_zero_per_row, true);
    /// 应用各块带宽.
    sp_vxvx.reinit(row_vx, col_vx, n_non_zero_per_row_vxvx, true);	/**< (0, 0) */
    sp_vyvx.reinit(row_vx, col_vy, n_non_zero_per_row_vyvx, false);	/**< (0, 1) */
    sp_pvx.reinit(row_vx, col_p, n_non_zero_per_row_pvx, false);	/**< (0, 2) */
    sp_vxvy.reinit(row_vy, col_vx, n_non_zero_per_row_vxvy, false);	/**< (1, 0) */
    sp_vyvy.reinit(row_vy, col_vy, n_non_zero_per_row_vyvy, true);	/**< (1, 1) */
    sp_pvy.reinit(row_vy, col_p, n_non_zero_per_row_pvy, false);	/**< (1, 2) */
    sp_vxp.reinit(row_p, col_vx, n_non_zero_per_row_vxp, false);	/**< (2, 0) */
    sp_vyp.reinit(row_p, col_vy, n_non_zero_per_row_vyp, false);	/**< (2, 1) */
    sp_penalty.reinit(row_p, col_p, n_non_zero_per_row_penalty, true);/**< (2, 2) */

    sp_mass_p.reinit(row_p, col_p, n_non_zero_per_row_mass_p);

    /// 第二次遍历, 指定每个非零元的坐标.
    for (the_element_v = fem_space_v.beginElement(); 
    	 the_element_v != end_element_v; ++the_element_v) 
    {
    	const std::vector<int>& element_dof_v = the_element_v->dof();
    	/// 由索引找到对应的压力空间单元.
    	Element<double, DIM> &p_element = fem_space_p.element(index_v2p[the_element_v->index()]);
    	const std::vector<int>& element_dof_p = p_element.dof();
    	int n_element_dof_v = the_element_v->n_dof();
    	int n_element_dof_p = p_element.n_dof();
	/// j 是检测函数指标. 行指标.
    	for (int j = 0; j < n_element_dof_v; ++j)
    	{
    	    for (int k = 0; k < n_element_dof_v; ++k)
    	    {
		/// (0, 0)
    		sp_stokes.add(element_dof_v[j], element_dof_v[k]);
    		sp_vxvx.add(element_dof_v[j], element_dof_v[k]);
		/// (0, 1)
    		sp_stokes.add(element_dof_v[j], element_dof_v[k] + n_dof_v);
    		sp_vyvx.add(element_dof_v[j], element_dof_v[k]);
		/// (1, 0)
    		sp_stokes.add(element_dof_v[j] + n_dof_v, element_dof_v[k]);
    		sp_vxvy.add(element_dof_v[j], element_dof_v[k]);
		/// (1, 1)
    		sp_stokes.add(element_dof_v[j] + n_dof_v, element_dof_v[k] + n_dof_v);
    		sp_vyvy.add(element_dof_v[j], element_dof_v[k]);
    	    }
    	    for (int k = 0; k < n_element_dof_p; ++k)
    	    { 
    		sp_stokes.add(element_dof_v[j], element_dof_p[k] + 2 * n_dof_v);
    		sp_pvx.add(element_dof_v[j], element_dof_p[k]);
    		sp_stokes.add(element_dof_v[j] + n_dof_v, element_dof_p[k] + 2 * n_dof_v);
    		sp_pvy.add(element_dof_v[j], element_dof_p[k]);
    	    }
    	}
    }

    for (the_element_p = fem_space_p.beginElement(); 
    	 the_element_p != end_element_p; ++the_element_p) 
    {
    	const std::vector<int>& element_dof_p = the_element_p->dof();
    	int n_element_dof_p = the_element_p->n_dof();
    	for (int j = 0; j < n_element_dof_p; ++j)
    	{
	    int idx_p = the_element_p->index();
	    int n_chi = index_p2v[idx_p].size();
	    for (int i = 0; i < n_chi; ++i)
	    {
		/// 由索引找到对应的速度空间单元.
		Element<double, DIM> &v_element = fem_space_v.element(index_p2v[idx_p][i]);
		const std::vector<int>& element_dof_v = v_element.dof();
		int n_element_dof_v = v_element.n_dof();
		for (int k = 0; k < n_element_dof_v; ++k)
		{
		    sp_stokes.add(element_dof_p[j] + 2 * n_dof_v, element_dof_v[k]);
		    sp_vxp.add(element_dof_p[j], element_dof_v[k]);
		    sp_stokes.add(element_dof_p[j] + 2 * n_dof_v, element_dof_v[k] + n_dof_v);
		    sp_vyp.add(element_dof_p[j], element_dof_v[k]);
		}
	    }
    	}
    }
    /// 构建压力空间质量矩阵结构, 用于预处理.
    for (the_element_p = fem_space_p.beginElement(); 
    	 the_element_p != end_element_p; ++the_element_p)
    {
    	const std::vector<int>& element_dof_p = the_element_p->dof();
    	int n_element_dof_p = the_element_p->n_dof();
    	for (int j = 0; j < n_element_dof_p; ++j)  
    	    for (int k = 0; k < n_element_dof_p; ++k)
    		sp_mass_p.add(element_dof_p[j], element_dof_p[k]);	    
    }


    /// 给右下角对角元留位置.
    for (int i = 0; i < n_dof_p; ++i)
    {
	sp_stokes.add(i + 2 * n_dof_v, i + 2 * n_dof_v);
	sp_penalty.add(i, i);
    }

    /// 矩阵模板压缩. 创建矩阵.
    sp_stokes.compress();
    /// 压缩各块模板.
    sp_vxvx.compress();	/**< (0, 0) */
    sp_vyvx.compress();	/**< (0, 1) */
    sp_pvx.compress();  /**< (0, 2) */
    sp_vxvy.compress();	/**< (1, 0) */
    sp_vyvy.compress();	/**< (1, 1) */
    sp_pvy.compress();	/**< (1, 2) */
    sp_vxp.compress();	/**< (2, 0) */
    sp_vyp.compress();	/**< (2, 1) */
    sp_penalty.compress();  /**< (2, 2) */

    sp_mass_p.compress();

    index_vxvx.resize(sp_vxvx.n_nonzero_elements());
    index_vyvx.resize(sp_vyvx.n_nonzero_elements());
    index_pvx.resize(sp_pvx.n_nonzero_elements());
    index_vxvy.resize(sp_vxvy.n_nonzero_elements());
    index_vyvy.resize(sp_vyvy.n_nonzero_elements());
    index_pvy.resize(sp_pvy.n_nonzero_elements());
    index_vxp.resize(sp_vxp.n_nonzero_elements());
    index_vyp.resize(sp_vyp.n_nonzero_elements());
    index_penalty.resize(sp_penalty.n_nonzero_elements());
    
    std::vector<int>::iterator index_vxvx_iterator = index_vxvx.begin();
    std::vector<int>::iterator index_vyvx_iterator = index_vyvx.begin();
    std::vector<int>::iterator index_pvx_iterator = index_pvx.begin();
    std::vector<int>::iterator index_vxvy_iterator = index_vxvy.begin();
    std::vector<int>::iterator index_vyvy_iterator = index_vyvy.begin();
    std::vector<int>::iterator index_pvy_iterator = index_pvy.begin();
    std::vector<int>::iterator index_vxp_iterator = index_vxp.begin();
    std::vector<int>::iterator index_vyp_iterator = index_vyp.begin();
    std::vector<int>::iterator index_penalty_iterator = index_penalty.begin();

    /// 准备整体矩阵索引.
    const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
    const unsigned int * column = sp_stokes.get_column_numbers();
    /// 准备各块列索引.
    const unsigned int * column_vxvx = sp_vxvx.get_column_numbers();
    const unsigned int * column_vyvx = sp_vyvx.get_column_numbers();
    const unsigned int * column_pvx = sp_pvx.get_column_numbers();
    const unsigned int * column_vxvy = sp_vxvy.get_column_numbers();
    const unsigned int * column_vyvy = sp_vyvy.get_column_numbers();
    const unsigned int * column_pvy = sp_pvy.get_column_numbers();
    const unsigned int * column_vxp = sp_vxp.get_column_numbers();
    const unsigned int * column_vyp = sp_vyp.get_column_numbers();
    const unsigned int * column_penalty = sp_penalty.get_column_numbers();

    for (int i = 0; i < n_total_dof; ++i)
    {
    	int l_i1 = row_vx;
    	int l_i2 = row_vx + row_vy;
    	int l_i3 = row_vx + row_vy + row_p;
    	int l_j1 = col_vx;
    	int l_j2 = col_vx + col_vy;
    	int l_j3 = col_vx + col_vy + col_p;
    	for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
    	{
    	    if (i >= 0 && i < l_i1)
    	    {
    		int l_i = i;
    		if (column[j] >= 0 && column[j] < l_j1)
		{
    		    if (column[j] !=  *column_vxvx)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vxvx_iterator = j;
    			index_vxvx_iterator++;
    			column_vxvx++;
    		    }
		}
    		else if (column[j] >= l_j1 && column[j] < l_j2)
		{
    		    if (column[j] - l_j1 !=  *column_vyvx)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vyvx_iterator = j;
    			index_vyvx_iterator++;
    			column_vyvx++;
    		    }
		}
    		else if (column[j] >= l_j2 && column[j] < l_j3)
		{
    		    if (column[j] - l_j2 !=  *column_pvx)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_pvx_iterator = j;
    			index_pvx_iterator++;
    			column_pvx++;
    		    }
		}
    	    }
    	    else if (i >= l_i1 && i < l_i2)
    	    {
		if (column[j] == i)
		{
		    if (column[j] - l_j1 != *column_vyvy)
		    {
			std::cout << "Matrix structure error!" << std::endl;
			exit(-1);
		    }
		    else
		    {
			*index_vyvy_iterator = j;
			index_vyvy_iterator++;
			column_vyvy++;
		    }

		}
    		else if (column[j] >= 0 && column[j] < l_j1)
		{
    		    if (column[j] !=  *column_vxvy)
    		    {
    			std::cout << "Matrix structure error!!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vxvy_iterator = j;
    			index_vxvy_iterator++;
    			column_vxvy++;
    		    }
		} 
    		else if (column[j] >= l_j1 && column[j] < l_j2)
    		    if (column[j] - l_j1 !=  *column_vyvy)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vyvy_iterator = j;
    			index_vyvy_iterator++;
    			column_vyvy++;
    		    }
    		else if (column[j] >= l_j2 && column[j] < l_j3)
    		    if (column[j] - l_j2 !=  *column_pvy)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_pvy_iterator = j;
    			index_pvy_iterator++;
    			column_pvy++;
    		    }
    	    }
    	    else if (i >= l_i2 && i < l_i3)
    	    {
    		if (column[j] >= 0 && column[j] < l_j1)
    		    if (column[j] !=  *column_vxp)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vxp_iterator = j;
    			index_vxp_iterator++;
    			column_vxp++;
    		    }
    		else if (column[j] >= l_j1 && column[j] < l_j2)
    		    if (column[j] - l_j1 !=  *column_vyp)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_vyp_iterator = j;
    			index_vyp_iterator++;
    			column_vyp++;
    		    }
    		else if (column[j] >= l_j2 && column[j] < l_j3)
    		    if (column[j] - l_j2 !=  *column_penalty)
    		    {
    			std::cout << "Matrix structure error!" << std::endl;
    			exit(-1);
    		    }
    		    else
    		    {
    			*index_penalty_iterator = j;
    			index_penalty_iterator++;
    			column_penalty++;
    		    }		   
    	    }
    	}
    }
    std::cout << "Matrix struct partion OK!" << std::endl;
};

#undef DIM
