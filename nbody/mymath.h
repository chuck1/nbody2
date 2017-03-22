
#include <exception>
#include <cmath>

#define M_PI (3.1415926535897)

namespace mymath
{

	long fact2(long x)
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

	double	K(double k)
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
	double	E(double k)
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

	double			volume_sphere_sphere_intersection(double d, double r0, double r1)
	{
		double p = r0 + r1 - d;
		double h0 = (r1 - r0 + d) * p / 2.0 / d;
		double h1 = (r0 - r1 + d) * p / 2.0 / d;
		double V0 = M_PI / 3.0 * pow(h0, 2.0) * (3.0 * r0 - h0);
		double V1 = M_PI / 3.0 * pow(h1, 2.0) * (3.0 * r1 - h1);
		return V0 + V1;
	}
}