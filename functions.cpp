#include "functions.h"

#define DIM 2
#define EPS 1e-12

double BoundaryMapping::value(const double *p) const
{
    for (int i = 0; i < domain->n_vertex; ++i)
    {
	/// 第 i 个节点
	if (fabs(p[0] - domain->physical_domain_vertex[i][0]) < EPS && 
	    fabs(p[1] - domain->physical_domain_vertex[i][1]) < EPS)
	{
	    if (mode == 1)
		return domain->logical_domain_vertex[i][0];
	    else if (mode == 2) 
		return domain->logical_domain_vertex[i][1];
	    else
	    {
		std::cerr << "No boundary mapping mode " << mode << std::endl;
		abort();
	    }
	}
	/// 由第 i - i + 1 个节点构成的边. 当 i = n_vertex 时, 是 i -
	/// 0.
	
    }
    return 0;
};


std::vector<double> BoundaryMapping::gradient(const double *p) const
{
    std::vector<double> result(DIM);
    result[0] = 0.0;
    result[1] = 0.0;
    return result;
};

#undef DIM
#undef EPS
