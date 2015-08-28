#include <AFEPack/Functional.h>
#include <vector>

#include "ISOP2P1.h"

#define DIM 2

class BurgV : public Function<double>
{
private:
    double a;
    double t;
public:
    BurgV(double _a, double _t) : a(_a), t(_t) 
    {};
    ~BurgV()
    {};
public:
    double value(const double *p) const
    {
	return 1.0 / (1.0 + exp((p[0] + p[1] - t) / (2.0 * a)));
    };
    std::vector<double> gradient(const double *p) const
    {
	std::vector<double> result(2);
	return result;
    };
};

class RealVx : public Function<double>
{
public:
	RealVx()
	{};
	~RealVx()
	{};
public:
	double value(const double *p) const
	{
		return 20.0 * p[0] * p[1] * p[1] * p[1];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 20.0 * p[1] * p[1] * p[1];
		result[1] = 60.0 * p[0] * p[1] * p[1];
		return result;
	};
};

class RealVy : public Function<double>
{
public:
	RealVy()
	{};
	~RealVy()
	{};
public:
	double value(const double *p) const
	{
		return 5.0 * p[0] * p[0] * p[0] * p[0] - 5.0 * p[1] * p[1] * p[1] * p[1];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 20.0 * p[0] * p[0] * p[0];
		result[1] = -20.0 * p[1] * p[1] * p[1];
		return result;
	};
};

class AccuracyVx: public Function<double>
{
private:
	double nu;
	double t;
public:
AccuracyVx(double _nu, double _t) : nu(_nu), t(_t)
	{};
	~AccuracyVx()
	{};
public:
	double value(const double *p) const
	{
		return -cos(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = sin(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
		result[1] = -cos(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
		return result;
	};
};

class AccuracyVy: public Function<double>
{
private:
	double nu;
	double t;
public:
AccuracyVy(double _nu, double _t) : nu(_nu), t(_t)
	{};
	~AccuracyVy()
	{};
public:
	double value(const double *p) const
	{
		return sin(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = cos(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
		result[1] = -sin(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
		return result;
	};
};


class RealP : public Function<double>
{
private:
	double average;
public:
RealP(double _a) : average(_a)
	{};
	~RealP()
	{};
public:
	double value(const double *p) const
	{
		return 60.0 * p[0] * p[0] * p[1] - 20.0 * p[1] * p[1] * p[1] + average;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 120.0 * p[0] * p[1];
		result[1] = 60.0 * p[0] * p[0] - 60.0 * p[1] * p[1];
		return result;
	};
};

class PoiseuilleVx : public Function<double>
{
private:
	double y0;
	double y1;
	double c;
	double v;
public:
PoiseuilleVx(double _y0, double _y1) : y0(_y0), y1(_y1)
	{
		c = (y0 + y1) * 0.5;
		v = (y0 - y1) * 0.5;
	};
	~PoiseuilleVx()
	{};
public:
	double value(const double *p) const
	{
		return   (1.0 - (p[1] - c) * (p[1] - c)  / (v * v));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = -2.0 * p[1];
		return result;
	};
};

class PoiseuilleVy : public Function<double>
{
public:
	PoiseuilleVy()
	{};
	~PoiseuilleVy() 
	{};
public:
	double value(const double *p) const
	{
		return 0.0;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

class PoiseuilleP : public Function<double>
{
private:
	double average;
public:
PoiseuilleP(double _a) : average(_a)
	{};
	~PoiseuilleP() 
	{};
public:
	double value(const double *p) const
	{
		return -2.0 * p[0] + average;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = -2.0;
		result[1] = 0.0;
		return result;
	};
};

class Regularized : public Function<double>
{
public:
	Regularized() 
	{};
	~Regularized() 
	{};
public:
	double value(const double *p) const
	{
		return 1.0 - p[0] * p[0] * p[0] * p[0];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = -4.0 * p[0] * p[0] * p[0];
		result[1] = 0.0;
		return result;
	};
};

class ForceX : public Function<double>
{
private:
	double body_force;
	double angle;
public:
ForceX(double _bf, double _an) : body_force(_bf),
		angle(_an)
		{};
	~ForceX()
	{};
public:
	double value(const double *p) const
	{
		double deg = angle / 180 * PI;
		return (body_force * cos(deg));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

class ForceY : public Function<double>
{
private:
	double body_force;
	double angle;
public:
        ForceY(double _bf, double _an) : body_force(_bf),
	                               	angle(_an)
		{};
	~ForceY()
	{};
public:
	double value(const double *p) const
	{
		double deg = angle / 180 * PI;
		return (body_force * sin(deg));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

class BoundaryMapping : public Function<double>
{
private:
	const MovingMesh2D::Domain *domain;
	const Mesh<DIM> *mesh;
	int mode;
public:
	/** 
	 * 构建边界映射.
	 * 
	 * @param _domain 计算和逻辑区域信息.
	 * @param _mesh 网格信息.
	 * @param _mode 模式: 1. xy2xi; 2. xy2eta; 暂时只有这些.
	 * 
	 * @return 
	 */
        BoundaryMapping(MovingMesh2D::Domain &_domain, 
	                	Mesh<DIM> &_mesh,
		int _mode) 
	: domain(&_domain), mesh(&_mesh), mode(_mode)
	{}; 
	~BoundaryMapping()
	{};
public:
	double value(const double *p) const;
	std::vector<double> gradient(const double *p) const;
};


#undef DIM
