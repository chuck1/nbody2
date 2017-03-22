
#ifndef KERNEL_H
#define KERNEL_H

struct Pair
{
	int i;
	int j;

	// distance
	double d;

	// relative speed
	double s;

	double dt;

	// radius force
	double F;
};
struct Header
{
	int bodies_size;
	double t;
	double dt;

	int count_pen;
};
struct Vec3
{
	double v[3];
};
struct Vec4
{
	double v[4];
};
struct Quat
{
	double v[4];
};
struct Mat4
{
	double v[4][4];
};
struct Mat3
{
	double v[3][3];
};
struct Mat43
{
	double v[4][3];
};

struct Body
{
	struct Vec3 pos;
	struct Vec3 vel;
	struct Vec3 acc;

	struct Quat q;
	struct Vec3 w;
	struct Vec3 a;

	struct Vec3 force;
	struct Vec3 torque;

	// debugging
	struct Vec4 q1;

	double radius;
	double mass;
	double density;
};

#endif


