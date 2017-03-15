
struct Pair
{
	int i;
	int j;
};
struct Header
{
	int bodies_size;
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