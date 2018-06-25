#include "stdafx.h"

#include "mymath.h"

long	mymath::fact2(long x)
{
	if (x < -1) throw std::exception();

	if (x < 1) return 1;

	long ret = 1;

	while (x > 1)
	{
		ret *= x;

		x -= 2;
	}

	return ret;
}

double	mymath::K(double k)
{
	// complete elliptic integral of the first kind

	double pi = 3.1415926535897;

	double ret = 0;

	for (int i = 0; i < 10; ++i)
	{
		double a = (double)fact2(2 * i - 1) / (double)fact2(2 * i);
		double b = pow(a, 2.0);
		double c = pow(k, 2.0 * i);

		ret += b * c;
	}

	ret *= pi / 2.0;

	return ret;
}
double	mymath::E(double k)
{
	// complete elliptic integral of the second kind

	double pi = 3.1415926535897;

	double ret = 0;

	for (int i = 1; i < 10; ++i)
	{
		ret += pow((double)fact2(2 * i - 1) / (double)fact2(2 * i), 2.0) * pow(k, 2.0 * i) / (double)(2 * i - 1);
	}

	ret = 1 - ret;

	ret *= pi / 2.0;

	return ret;
}

double	mymath::volume_sphere_sphere_intersection(double d, double r0, double r1)
{
	double p = r0 + r1 - d;
	double h0 = (r1 - r0 + d) * p / 2.0 / d;
	double h1 = (r0 - r1 + d) * p / 2.0 / d;
	double V0 = M_PI / 3.0 * pow(h0, 2.0) * (3.0 * r0 - h0);
	double V1 = M_PI / 3.0 * pow(h1, 2.0) * (3.0 * r1 - h1);
	return V0 + V1;
}


double	mymath::gravity_uniform_disc(double density, double p, double r)
{
	double k = 4.0 * r * p / (pow(r, 2.0) + pow(p, 2.0) + 2.0 * r * p);
	double K = mymath::K(k);
	double E = mymath::E(k);
	double a = (1.0 - 0.5 * pow(k, 2.0)) * K - E;
	return 4.0 * M_G * density * pow(p, 0.5) / k / pow(r, 0.5) * a;
}

double		mymath::sphere_radius(double volume)
{
	return pow(3.0 * volume / 4.0 / M_PI, 1.0 / 3.0);
}