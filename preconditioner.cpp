#include "preconditioner.h"
#define DIM 2

SchurComplement::SchurComplement(SparseMatrix<double> &_BTx,
				 SparseMatrix<double> &_BTy,
				 SparseMatrix<double> &_Bx,
				 SparseMatrix<double> &_By,
				 InverseMatrix &_Ax,
				 InverseMatrix &_Ay)
		:
	BTx(&_BTx),
	BTy(&_BTy),
	Bx(&_Bx),
	By(&_By),
	AInvx(&_Ax),
	AInvy(&_Ay),
	tmp11(BTx->m()),
	tmp12(BTy->m()),
	tmp21(BTx->m()),
	tmp22(BTy->m())
{
};


void SchurComplement::vmult (Vector<double>       &dst,
			     const Vector<double> &src) const
{
	BTx->vmult (tmp11, src);
	BTy->vmult (tmp12, src);
	AInvx->vmult(tmp21, tmp11);
	AInvy->vmult(tmp22, tmp12);
	Bx->vmult(dst, tmp21);
	By->vmult_add(dst,tmp22);
};

StokesPreconditioner::StokesPreconditioner()
{
	Ax = NULL;
	Ay = NULL;
	Q = NULL;
};

StokesPreconditioner::~StokesPreconditioner()
{};

void StokesPreconditioner::initialize (const SparseMatrix<double> &_stiff_vx, 
				       const SparseMatrix<double> &_stiff_vy,
				       const SparseMatrix<double> &_mass_p_diag)
{
	Ax = &_stiff_vx;
	Ay = &_stiff_vy;
	Q = &_mass_p_diag;
	AMGx.reinit(*Ax);
	AMGy.reinit(*Ax);
};

void StokesPreconditioner::vmult (Vector<double> &dst,
				  const Vector<double> &src) const
{
	int n_dof_v = Ax->n();
	int n_dof_p = Q->n();
	Vector<double> d0(n_dof_v);
	Vector<double> d1(n_dof_v);
	Vector<double> s0(n_dof_v);
	Vector<double> s1(n_dof_v);

	for (int i = 0; i < n_dof_v; ++i)
		s0(i) = src(i);
	for (int i = 0; i < n_dof_v; ++i)
		s1(i) = src(n_dof_v + i);
	for (int i = 0; i < n_dof_p; ++i)
		dst(2 * n_dof_v + i) = src(2 * n_dof_v + i) / (*Q).diag_element(i);

	AMGx.solve(d0, s0, 1e-8, 1, 1);

	AMGy.solve(d1, s1, 1e-8, 1, 1);

	for (int i = 0; i < n_dof_v; ++i)
		dst(i) = d0(i);
	for (int i = 0; i < n_dof_v; ++i)
		dst(n_dof_v + i) = d1(i);
};

InverseMatrix::InverseMatrix (const SparseMatrix<double> &m,  const AMGSolver &a)
 :matrix (&m), amg(&a)
{
};

InverseMatrix::~InverseMatrix()
{
};

void InverseMatrix::vmult(Vector<double> &dst,
			  const Vector<double> &src) const
{
	dst = 0;
	amg->solve(dst, src, 1e-8*src.l2_norm(), 0);
};

#undef DIM

