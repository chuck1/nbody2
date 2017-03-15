
#include "kernel.h"



__kernel void calc_acc(
	__global struct Header * header,
	__global struct Body * bodies,
	__global struct Pair * pairs,
	volatile __global unsigned int * counter)
{
	if (get_global_id(0) == 0) *counter = 0;
	barrier(CLK_GLOBAL_MEM_FENCE);

	int n = header->bodies_size;
	
	unsigned int len = n * (n - 1) / 2;

	double G = 6.67408E-11;

	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		int i = pairs[ind].i;
		int j = pairs[ind].j;

		// debugging
		//if (i > n) break;
		//if (j > n) break;

		__global struct Body * b1 = &bodies[i];
		__global struct Body * b2 = &bodies[j];

		double x = b1->pos.x - b2->pos.x;
		double y = b1->pos.y - b2->pos.y;
		double z = b1->pos.z - b2->pos.z;

		double rs = x*x + y*y + z*z;

		double r = sqrt(rs);

		double F = G / rs;

		double m0 = b1->mass;
		double m1 = b2->mass;

		double pen = b1.radius + b2.radius - r;

		if (pen > 0)
		{
			double d0 = pen * (m1 / (m0 + m1)) * 0.5;
			b1->pos.x += (x / rs) * d0;
			b1->pos.y += (y / rs) * d0;
			b1->pos.z += (z / rs) * d0;

			double d1 = pen * (m0 / (m0 + m1)) * 0.5;
			b2->pos.x += (x / rs) * d1;
			b2->pos.y += (y / rs) * d1;
			b2->pos.z += (z / rs) * d1;

			conserve momentum!
		}

		b1->acc.x -= x / r * F * b2->mass;
		b1->acc.y -= y / r * F * b2->mass;
		b1->acc.z -= z / r * F * b2->mass;

		b2->acc.x += x / r * F * b1->mass;
		b2->acc.y += y / r * F * b1->mass;
		b2->acc.z += z / r * F * b1->mass;
	}

}



__kernel void step_pos(
	__global struct Header * header,
	__global struct Body * bodies,
	__global struct Pair * pairs,
	volatile __global unsigned int * counter)
{
	if (get_global_id(0) == 0) *counter = 0;
	barrier(CLK_GLOBAL_MEM_FENCE);

	int len = header->bodies_size;

	double dt = 1.0E-1;

	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		int i = ind;

		__global struct Body * b = &bodies[i];

		b->vel.x += b->acc.x * dt;
		b->vel.y += b->acc.y * dt;
		b->vel.z += b->acc.z * dt;

		b->pos.x += b->vel.x * dt;
		b->pos.y += b->vel.y * dt;
		b->pos.z += b->vel.z * dt;

		b->acc.x = 0;
		b->acc.y = 0;
		b->acc.z = 0;
	}
}

