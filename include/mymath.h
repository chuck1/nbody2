#ifndef MYMATH_H
#define MYMATH_H

#include <exception>
#include <cmath>

#define M_PI (3.1415926535897)
#define M_G (6.67408E-11)

namespace mymath
{

	long		fact2(long x);

	double		K(double k);
	double		E(double k);

	double		volume_sphere_sphere_intersection(double d, double r0, double r1);

	double		sphere_radius(double volume);

	double		gravity_uniform_disc(double density, double p, double r);
}

#endif

