
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
	double x, y, z;
};

struct Body
{
	struct Vec3 pos;
	struct Vec3 vel;
	struct Vec3 acc;
	double radius;
	double mass;
};