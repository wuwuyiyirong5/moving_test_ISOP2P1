#include <AFEPack/AMGSolver.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <vector>

#define DIM 2

/**
 * 对应的预处理.
 *
 */
class StokesPreconditioner
{
private:
	const SparseMatrix<double> *Ax; /**< 预处理矩阵各分块. */
	const SparseMatrix<double> *Ay;
	const SparseMatrix<double> *Q;

public:
	StokesPreconditioner();

	~StokesPreconditioner();

	/**
	 * 预处理子初始化.
	 *
	 * @param _stiff_vx vx 空间的刚度矩阵.
	 * @param _stiff_vy vy 空间的刚度矩阵.
	 * @param _mass_p_diag p 空间的质量矩阵的对角元.
	 */
	void initialize (const SparseMatrix<double> &_stiff_vx,
			const SparseMatrix<double> &_stiff_vy,
			const SparseMatrix<double> &_mass_p_diag);

	/**
	 * 实际估值 dst = M^{-1}src.
	 *
	 * @param dst
	 * @param src
	 */
	void vmult (Vector<double> &dst,
			const Vector<double> &src) const;
};

class InverseMatrix : public Subscriptor
{
public:
	InverseMatrix (const SparseMatrix<double> &m, const AMGSolver &a);
	~InverseMatrix ();

	void vmult (Vector<double>       &dst,
			const Vector<double> &src) const;

private:
	const SmartPointer<const SparseMatrix<double> > matrix;
	const AMGSolver *amg;
};



class SchurComplement : public Subscriptor
{
public:
	SchurComplement(SparseMatrix<double> &_BTx,
			SparseMatrix<double> &_BTy,
			SparseMatrix<double> &_Bx,
			SparseMatrix<double> &_By,
			InverseMatrix &_AInvx,
			InverseMatrix &_AInvy);

	void vmult (Vector<double>       &dst,
			const Vector<double> &src) const;

private:
	const SmartPointer<const SparseMatrix<double> > BTx;
	const SmartPointer<const SparseMatrix<double> > BTy;
	const SmartPointer<const SparseMatrix<double> > Bx;
	const SmartPointer<const SparseMatrix<double> > By;
	const SmartPointer<const InverseMatrix> AInvx;
	const SmartPointer<const InverseMatrix> AInvy;

	mutable Vector<double> tmp11, tmp12, tmp21, tmp22;
};

#undef DIM
