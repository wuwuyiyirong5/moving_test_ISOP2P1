/**
 * @file   ISOP2P1.h
 * @author Heyu Wang <scshw@cslin107.csunix.comp.leeds.ac.uk>
 * @date   Tue Nov  4 14:19:55 2014
 * 
 * @brief 将 ISO P2P1 有限元的网格, 空间, 矩阵信息集中封装, 提供给
 * Stokes 和 Navier-Stokes 求解器使用.
 * 
 * 
 */

#ifndef __CRAZYFISH__ISOP2P1__
#define __CRAZYFISH__ISOP2P1__

#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Operator.h>
#include <AFEPack/MovingMesh2D.h>


#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_gmres.h>
#include <lac/solver_minres.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_mic.h>

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#define DIM 2

#define PI (atan(1.0) * 4.0)


/**
 * 主类. 针对无结构三角形网格. 带有两层网格结构. 细网格是粗网格三角形四
 * 等分形成. 原始一个粗网格单元构成的四个细网格单元形成一个宏单元片. 细
 * 网格将用于离散速度, 粗网格是压力. 
 * 
 */
class ISOP2P1 : public MovingMesh2D
{
private:
	/// 网格和几何信息.
	HGeometryTree<DIM> h_tree; /**< 网格树. */
	IrregularMesh<DIM> *irregular_mesh_p;  /**< P 空间不规则网格. (链表) */
	IrregularMesh<DIM> *irregular_mesh_v;  /**< V 空间不规则网格. (链表) */
	/* MovingMesh2D::Domain domain; /\**< 计算区域信息. *\/ */
	std::vector<Point<DIM> > logical_node_mesh_v; /**< 逻辑网格节点座标. */
	std::vector<Point<DIM> > displacement4mesh_p; /**< 记录压力移动位移, 为了更新速度网格. */
	std::vector<Point<DIM> > move_direction4mesh_v; /**< 速度网格移动方向. */

	TemplateGeometry<DIM> template_geometry; /**< 参考几何信息. */
	CoordTransform<DIM, DIM> coord_transform; /**< 座标变换信息. */

	/// 有限元空间构成. 这里本质上只有一种有限元, P1元.
	TemplateDOF<DIM> template_dof; /**< 参考自由度. */
	BasisFunctionAdmin<double, DIM, DIM> basis_function; /**< 基函数. */
	std::vector<TemplateElement<double, DIM, DIM> > template_element; /**< 参考单元. */

	FEMSpace<double, DIM> fem_space_v; /**< 速度空间. */
	std::vector<FEMFunction<double, DIM> > v_h;	/**< 速度数值解. (向量型) */
	std::vector<FEMFunction<double, DIM> > source_v; /**< 速度对应源项. (向量型) */

	FEMSpace<double, DIM> fem_space_p; /**< 压力空间. */
	FEMFunction<double, DIM> p_h; /**< 压力数值解. (标量型) */
	FEMFunction<double, DIM> source_p; /**< 压力对应源项. (标量型) */

	FEMFunction<double, DIM> wxx; /**< 外场速度, 用于调试. */
	FEMFunction<double, DIM> wxy; /**< 外场速度, 用于调试. */

	FEMFunction<double, DIM> vot; /**< 涡量, 用于输出. */
	FEMFunction<double, DIM> div; /**< 散度, 用于输出. */

	/// 速度和压力空间单元索引.
	std::vector<std::vector<int> > index_p2v;
	std::vector<int> index_v2p;

	/// 矩阵模板.
	SparsityPattern sp_stokes; /**< 离散 Stokes 方程系数矩阵模板. */

	/// 将整体的系数矩阵分裂成 3 x 3 的矩阵块, 建立分块操作索引.
	SparsityPattern sp_vxvx;	/**< (0, 0) */
	std::vector<int> index_vxvx; 
	SparsityPattern sp_vyvx;	/**< (0, 1) */
	std::vector<int> index_vyvx; 
	SparsityPattern sp_vxvy;	/**< (1, 0)*/
	std::vector<int> index_vxvy; 
	SparsityPattern sp_vyvy;	/**< (1, 1) */
	std::vector<int> index_vyvy; 
	SparsityPattern sp_pvx;	/**< (0, 2) */
	std::vector<int> index_pvx; 
	SparsityPattern sp_pvy;	/**< (1, 2) */
	std::vector<int> index_pvy; 
	SparsityPattern sp_vxp;	/**< (2, 0) */
	std::vector<int> index_vxp; 
	SparsityPattern sp_vyp;	/**< (2, 1) */
	std::vector<int> index_vyp; 
	SparsityPattern sp_penalty;	/**< (2, 2) */
	std::vector<int> index_penalty; 

	SparsityPattern sp_mass_p;	/**< 预处理压力质量矩阵模板. */
	SparsityPattern sp_div_p;	/**< 预处理压力散度矩阵模板. */

	SparseMatrix<double> mat_v_stiff; /**< 速度空间 x 方向刚度矩阵. */
	SparseMatrix<double> mat_v_mass; /**< 速度空间 x 方向质量矩阵. */
	SparseMatrix<double> mat_vxp_div; /**< 混合空间 x 方向散度矩阵. */
	SparseMatrix<double> mat_vyp_div; /**< 混合空间 y 方向散度矩阵. */
	SparseMatrix<double> mat_vzp_div; /**< 混合空间 z 方向散度矩阵. */
	SparseMatrix<double> mat_pvx_divT; /**< 混合空间 x 方向散度矩阵转置. */
	SparseMatrix<double> mat_pvy_divT; /**< 混合空间 y 方向散度矩阵转置. */
	SparseMatrix<double> mat_pvz_divT; /**< 混合空间 z 方向散度矩阵转置. */

	SparseMatrix<double> mat_p_stiff; /**< 压力空间刚度矩阵. */
	SparseMatrix<double> mat_p_mass; /**< 压力空间质量矩阵. */

	SparseMatrix<double> mat_v_convection; /**< 速度空间对流矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_xx; /**< 速度空间 Jacobi xx 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_xy; /**< 速度空间 Jacobi xy 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_xz; /**< 速度空间 Jacobi xz 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_yx; /**< 速度空间 Jacobi yx 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_yy; /**< 速度空间 Jacobi yy 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_yz; /**< 速度空间 Jacobi yz 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_zx; /**< 速度空间 Jacobi zx 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_zy; /**< 速度空间 Jacobi zy 矩阵块. */
	SparseMatrix<double> mat_v_Jacobi_zz; /**< 速度空间 Jacobi zz 矩阵块. */

	SparseMatrix<double> mat_pcd; /**< 预处理, 对流扩散矩阵在压力空间的投影. */

	SparseMatrix<double> matrix; /**< 问题的总系数矩阵. */

	Vector<double> rhs;		/**< 离散方程右端项. */

	double eps;			/**< 机器浮点精度. */
	std::string mesh_file;	/**< 网格树文件名, easymesh 格式. */
	double l_tol;			/**< 线性求解容许残量. */
	double l_Euler_tol;     /**< 显示欧拉方法容许残量. */
	double n_tol;			/**< 非线性求解容许残量. */
	bool record;		/**< 是否记录结果. */
	bool check;			/**< 是否显示迭代过程. */
	double error_check;		/**< 是否分析后验误差. */
	double viscosity;		/**< 方程粘性系数. */

	double t0;			/**< 起始时间. */
	double t1;			/**< 终止时间. */
	double t;			/**< 当前时间. */
	double dt;			/**< 下一时间步. */

	int n_method;		/**< 非线性迭代选择, 1: Newton; 2: Picard; 3: Hybrid. */
	bool time_step_control;	/**< 是否要自适应时间步长. */
	bool Stokes_init;		/**< 是否要先运行稳态 Stokes 做初值. */
	bool NS_init;		/**< 是否要先运行稳态 Navier-Stokes 做初值. */
	bool isMoving;             /**< 网格是否进行移动. */
	double scale;           /**< 初始网格参数. */
	int max_step;           /**< 最大移动次数. */
	double alpha;		/**< 计算monitor的参数. */
	double beta;            /**< 计算monitor的参数. */
	int scheme;			/**< 时间发展格式. */
	bool output_vorticity;	/**< 是否要输出涡量. */
	bool output_divergence;	/**< 是否要输出散度. */
	double body_force;          /**< 外力. */
	double angle;		/**< 倾角. x 轴正方向为 0 度.*/

	/**
	 * 关于时间发展格式: 
	 *
	 * = 1: Back Euler, 显式处理非线性项, 由于稳定性
	 * 问题, time_step_control 必须为真. 
	 *
	 * = 2: Back Euler, 线性化非线性项, 无条件稳定, 时间步长不宜太小,
	 * RE 不宜太大 (< 100). 
	 *
	 * 如想添加其他格式可以在这个框架下增加.
	 */

public:

	class Matrix : public StiffMatrix<2,double>
	{
	private:
		double dt, a;
	public:
		Matrix(FEMSpace<double,2>& sp, const double& _dt, const double& _a) :
			dt(_dt), a(_a),
			StiffMatrix<2,double>(sp) {};
		virtual ~Matrix() {};
	public:
		virtual void getElementMatrix(const Element<double,2>& e0,
				const Element<double,2>& e1,
				const ActiveElementPairIterator<2>::State state);
	};



	/** 
	 * 
	 * 缺省构造. 
	 */
	ISOP2P1()
		: scheme(1),
		  NS_init(false),
		  Stokes_init(false),
		  time_step_control(false),
		  angle(0),
		  l_Euler_tol(1.0e-12),
		  n_tol(1.0e-12),
		  viscosity(1.0),
		  dt(1.0e-3),
		  n_method(1),
		  body_force(0.0),
		  l_tol(1.0e-12),
		  t1(1.0),
		  t0(0.0),
		  t(0.0)
	{
		irregular_mesh_v = NULL;
		irregular_mesh_p = NULL;
		eps = std::numeric_limits<double>::epsilon();
	};

	/** 
	 * 缺省析构.
	 * 
	 */
	~ISOP2P1()
	{
		if (irregular_mesh_p != NULL)
			delete irregular_mesh_p;
		if (irregular_mesh_v != NULL)
			delete irregular_mesh_v;
	};
     
	/** 
	 * 主流程.
	 * 
	 */
	void initialize();

	/** 
	 * 构建计算网格和宏单元.
	 * 
	 */
	void buildMesh();

	/** 
	 * 构建有限元空间.
	 * 
	 */
	void buildFEMSpace();
	
	/** 
	 * 当树结构改变的时候,需要重新构造有限元空间,现在在调试.
	 * 重新构造网格也在这里操作了,懒得再写重新构造网格的函数了.
	 */
	void rebuildFEMSpace();

	/** 
	 * 构建稀疏矩阵, 块矩阵结构.
	 * 
	 */
	void buildMatrixStruct();

	/** 
	 * 构建线性矩阵和右端项.
	 * 
	 */
	void buildMatrix();

	/** 
	 * 构建用于预处理的速度在压力空间投射构成的对流扩散矩阵.
	 * 
	 */
	void updatePCDMatrix();

	/** 
	 * 更新非线性残量 ( 右端项 ) 和系数矩阵块. 
	 * 
	 * @param mode 0: 非线性残量; 1: Linearized back Euler 右端项.
	 */
	void updateNonlinearMatrix();

	/** 
	 * 构建 Stokes 矩阵.
	 * 
	 */
	void buildStokesSys();

	/** 
	 * 构建 Newton 法求解 NS 问题的迭代矩阵.
	 * 
	 */
	void buildNewtonSys4NS();

	/** 
	 * 构建 Picard 迭代求解 NS 问题的迭代矩阵.
	 * 
	 */
	void buildPicardSys4NS();

	/** 
	 * 求解 Stokes .
	 * 
	 */
	void solveStokes();

	/** 
	 * 求解 Navier-Stokes 方程的非线性求解器.
	 * 
	 * @param method 迭代方法, 1: Newton; 2: Picard; 3: Hybrid.
	 */
	void solveNS(int method);

	/** 
	 * 求解一步 linearized back Euler 问题.
	 * 
	 */
	void solveLinearizedBackEuler();

	/** 
	 * 求解一个问题的主流程.
	 * 
	 */
	void run();

	/** 
	 * 读入配置文件.
	 * 
	 * @param _config_file 配置文件名.
	 */
	void config(std::string _config_file);
	
	/** 
	 * debug config函数.
	 * 
	 */
	void config_debug();
	/** 
	 * Stokes 问题的边界条件处理.
	 * 
	 * @param x 线性方程组未知量.
	 */
	void boundaryValueStokes(Vector<double> &x);

	/** 
	 * updateSolution()中边界条件处理.
	 * 
	 */
	void boundaryValueUpdateSolution(Vector<double> &x);
	/** 
	 * Stokes 问题的边界条件处理.
	 * 
	 * @param x 线性方程组未知量.
	 */
	void boundaryValueNS(Vector<double> &x);

	/** 
	 * Linearized back Euler 的系统构建.
	 * 
	 */
	void buildLinearizedBackEulerScheme();

	/** 
	 * 发展一步线性化 Euler .
	 * 
	 */
	void stepForwardLinearizedEuler();

	/** 
	 * 发展一步显式 Euler (非线性项显式处理).
	 * 
	 */
	void stepForwardEuler();

	/** 
	 * 时间步长判定.
	 * 
	 */
	void time_step();
	
	/** 
	 * 计算散度.
	 * 
	 */
	void computDivergence();

	/** 
	 * 计算涡量.
	 * 
	 */
	void computVorticity();
    
	/** 
	 * 输出 tecplot 格式的数值解.
	 * 
	 * @param prefix 文件名前缀.
	 */
	void outputTecplot(const std::string &prefix);

	/** 
	 * 输出 tecplot 格式的 P 网格, 用于调试.
	 * 
	 * @param prefix 
	 */
	void outputTecplotP(const std::string &prefix);

	/** 
	 * 输出 tecplot 格式的调试数据.
	 * 
	 * @param prefix 文件名前缀.
	 */
	void outputTecplot_debug(const std::string &prefix);

	/** 
	 * 移动网格时的 monitor, 暂时取全 1. 因为是移动区域.
	 * 
	 */
	virtual void getMonitor(); 

	/** 
	 * 移动网格时的数值解更新.
	 * 
	 */
	virtual void updateSolution(); 

	/** 
	 * 暂时给一个,基类中的纯虚函数必须在派生类中实现.
	 * 
	 */
	virtual void outputSolution();

	/** 
	 * 同步全部网格. 用于移动网格.
	 * 
	 */
	void syncMesh();

	/** 
	 * 一步网格移动.
	 * 
	 */
	void movingMesh();

	void stepForward();

	
	/** 
	 * 调试测试程序.
	 * 
	 */
	void debug();

	void assemble();

	void buildJacobi();

	void solve();
};
#undef DIM
#endif

//
// end of file
//////////////////////////////////////////////////////////////////////////////
