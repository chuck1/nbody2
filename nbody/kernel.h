
#ifndef KERNEL_H
#define KERNEL_H

struct Pair
{
	int i;
	int j;

#ifdef KDEBUG
	// distance
	double d;
	// relative speed
	double s;
	// radial force
	double F;
#endif

	//double dt;
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

struct Body0
{
	struct Vec3 pos;
	
	struct Quat q;

	double radius;
};
struct Body1
{
	struct Vec3 vel;
	struct Vec3 acc;

	struct Vec3 w;
	struct Vec3 a;
	
	struct Vec3 force;
	struct Vec3 torque;

	struct Vec4 q1;

	double mass;
	double density;
};

#endif


