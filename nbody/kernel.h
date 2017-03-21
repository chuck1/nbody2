
struct Pair
{
	int i;
	int j;

	// distance
	double d;

	// relative speed
	double s;

	double dt;
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
struct Mat4
{
	double v[4][4];
};


struct Body
{
	struct Vec3 pos;
	struct Vec3 vel;
	struct Vec3 acc;

	struct Vec4 q0;
	struct Vec4 q1;
	struct Vec4 q2;

	double radius;
	double mass;
};