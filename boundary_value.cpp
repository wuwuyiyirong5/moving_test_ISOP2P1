#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::boundaryValueStokes(Vector<double> &x)
{
	/// 各空间自由度.
	unsigned int n_dof_v = fem_space_v.n_dof();
	unsigned int n_dof_p = fem_space_p.n_dof();
	unsigned int n_total_dof_v = 2 * n_dof_v;
	const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
	const unsigned int * colnum = sp_stokes.get_column_numbers();
	std::cout << "n_dof_v: " << n_dof_v << ", n_dof_p: " << n_dof_p << std::endl;
	std::cout << "n_A: " << sp_stokes.n_rows() << ", m_A: " << sp_stokes.n_cols() << std::endl;

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

		/// 方腔流边界条件.
		if (bm == 1 || bm == 2 || bm == 4)
			x(i) = 0.0;
		else if (bm == 3)
			if (i < n_dof_v)
				x(i) = 1.0;
			else
				x(i) = 0.0;

		/// 右端项这样改, 如果该行和列其余元素均为零, 则在迭代中确
		/// 保该数值解和边界一致.
		if (bm == 1 || bm == 2 || bm == 3 || bm == 4)
		{
			rhs(i) = matrix.diag_element(i) * x(i);
			/// 遍历 i 行.
			for (unsigned int j = rowstart[i] + 1;
					j < rowstart[i + 1]; ++j)
			{
				/// 第 j 个元素消成零(不是第 j 列!). 注意避开了对角元.
				matrix.global_entry(j) -= matrix.global_entry(j);
				/// 第 j 个元素是第 k 列.
				unsigned int k = colnum[j];
				/// 看看第 k 行的 i 列是否easymesh 为零元.
				const unsigned int *p = std::find(&colnum[rowstart[k] + 1],
						&colnum[rowstart[k + 1]],
						i);
				/// 如果是非零元. 则需要将这一项移动到右端项. 因为第 i 个未知量已知.
				if (p != &colnum[rowstart[k + 1]])
				{
					/// 计算 k 行 i 列的存储位置.
					unsigned int l = p - &colnum[rowstart[0]];
					/// 移动到右端项. 等价于 r(k) = r(k) - x(i) * A(k, i).
					rhs(k) -= matrix.global_entry(l)
					* x(i);
					/// 移完此项自然是零.
					matrix.global_entry(l) -= matrix.global_entry(l);
				}
			}
		}
	}
	std::cout << "boundary values for Stokes OK!" << std::endl;
};

void ISOP2P1::boundaryValueNS(Vector<double> &x)
{
	/// 各空间自由度.
	unsigned int n_dof_v = fem_space_v.n_dof();
	unsigned int n_dof_p = fem_space_p.n_dof();
	unsigned int n_total_dof_v = 2 * n_dof_v;
	const std::size_t * rowstart = sp_stokes.get_rowstart_indices();
	const unsigned int * colnum = sp_stokes.get_column_numbers();

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
		// /// 对全部 Dirichlet 边界.
		if (bm  == 2 || bm  == 3 || bm == 5 || bm == 1 || bm == 4 || bm == 11)
			// if (bm < 10 && bm > 0 && bm != 6)
			// if (bm == 1 || bm == 2 || bm == 3 || bm == 4)
		{
			/// 数值解对应点按成驱动速度.
			x(i) = 0.0;
			/// 右端项这样改, 如果该行和列其余元素均为零, 则在迭代中确
			/// 保该数值解和边界一致.
			rhs(i) = matrix.diag_element(i) * x(i);
			/// 遍历 i 行.
			for (unsigned int j = rowstart[i] + 1;
					j < rowstart[i + 1]; ++j)
			{
				/// 第 j 个元素消成零(不是第 j 列!). 注意避开了对角元.
				matrix.global_entry(j) -= matrix.global_entry(j);
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
					rhs(k) -= matrix.global_entry(l)
					* x(i);
					/// 移完此项自然是零.
					matrix.global_entry(l) -= matrix.global_entry(l);
				}
			}
		}
	}
	std::cout << "boundary apply to NS OK!" << std::endl;
};


#undef DIM
