
#include "kernel.h"

#define TIME_STEP_FACTOR (0.001)

struct Vec3 vector_sub(
	__global struct Vec3 * a, 
	__global struct Vec3 * b)
{
	struct Vec3 ret;
	ret.x = a->x - b->x;
	ret.y = a->y - b->y;
	ret.z = a->z - b->z;
	return ret;
}

double vector_length(struct Vec3 * v)
{
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

double CDF_uniform(double a, double b, double x)
{
	if (x < a) return 0;
	if (x > b) return 1;
	return (x - a) / (b - a);
}

__kernel void k_min(
	__global const double * input,
	__global unsigned int * len,
	__global double * partial,
	__local double * loc)
{
	/*
	length of partial must be at least number of work groups
	final result is in partial[0]
	*/

	uint local_id = get_local_id(0);
	uint global_id = get_global_id(0);
	uint group_size = get_local_size(0);
	uint global_size = get_global_size(0);

	// Copy from global to local memory
	if (global_id < *len)
		loc[local_id] = input[global_id];
	else
		loc[local_id] = 0;

	// Loop for computing localSums : divide WorkGroup into 2 parts
	for (uint stride = group_size / 2; stride > 0; stride /= 2)
	{
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);

		// Add elements 2 by 2 between local_id and local_id + stride
		if (local_id < stride)
		{
			if (loc[local_id] == 0)
			{
				loc[local_id] = loc[local_id + stride];
			}
			else if (loc[local_id + stride] == 0)
			{
			}
			else 
			{
				loc[local_id] = min(loc[local_id], loc[local_id + stride]);
			}
		}
	}

	// Write result into partialSums[nWorkGroups]
	if (local_id == 0)
		partial[get_group_id(0)] = loc[0];
	
	barrier(CLK_GLOBAL_MEM_FENCE);

	if (global_id == 0)
	{
		for (int i = 1; i < get_num_groups(0); ++i)
		{
			if (partial[0] == 0)
			{
				partial[0] = partial[i];
			}
			else if (partial[i] == 0)
			{
			}
			else
			{
				partial[0] = min(partial[0], partial[i]);
			}
		}
	}
}

__kernel void store_dt(
	__global double * dt_partial,
	__global struct Header * header
	)
{
	if (get_global_id(0) == 0)
	{
		header->dt = dt_partial[0];
	}
}

__kernel void calc_acc(
	__global struct Header * header,
	__global struct Body * bodies,
	__global struct Pair * pairs,
	__global double * dt_input,
	volatile __global unsigned int * counter)
{
	if (get_global_id(0) == 0) {
		header->dt = 0;
		header->count_pen = 0;
		*counter = 0;
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

	int n = header->bodies_size;
	
	unsigned int len = n * (n - 1) / 2;

	double G = 6.67408E-11;

	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		__global struct Pair * p = &pairs[ind];

		int i = p->i;
		int j = p->j;

		// debugging
		//if (i > n) break;
		//if (j > n) break;

		__global struct Body * b1 = &bodies[i];
		__global struct Body * b2 = &bodies[j];

		struct Vec3 R = vector_sub(&b1->pos, &b2->pos);

		double r = vector_length(&R);
		
		struct Vec3 v = vector_sub(&b2->vel, &b1->vel);
		double s = vector_length(&v);

		// save
		p->d = r;
		p->s = s;

		// calc dt
		dt_input[ind] = TIME_STEP_FACTOR * r / s;
		
		double m0 = b1->mass;
		double m1 = b2->mass;

		double F = G * b1->mass * b2->mass / r / r;
		
		double a0 = G * b2->mass / r/r;
		double a1 = G * b1->mass / r/r;

		// positive penetration is no penetration

		double pen = p->d - b1->radius - b2->radius;

		double pen_frac = max(pen / b1->radius, pen / b2->radius);

		if (pen < 0)
		{
			atomic_inc(&header->count_pen);
			
			// damping force
			// must have body rotation in order to allow two touching bodies to spin together

			if (1)
			{
				// repulsive force shall be equal and opposite to gravitaional force

				double c = CDF_uniform(-0.1, 0.1, pen_frac);

				a0 = a0 * (2.0 * c - 1.0);
				a1 = a1 * (2.0 * c - 1.0);

				/*
				double d0 = pen * (m1 / (m0 + m1)) * 0.5;
				b1->pos.x += (x / rs) * d0;
				b1->pos.y += (y / rs) * d0;
				b1->pos.z += (z / rs) * d0;

				double d1 = pen * (m0 / (m0 + m1)) * 0.5;
				b2->pos.x += (x / rs) * d1;
				b2->pos.y += (y / rs) * d1;
				b2->pos.z += (z / rs) * d1;
				*/
				//conserve momentum!
			}
		}

		b1->acc.x -= R.x / r * a0;
		b1->acc.y -= R.y / r * a0;
		b1->acc.z -= R.z / r * a0;

		b2->acc.x += R.x / r * a1;
		b2->acc.y += R.y / r * a1;
		b2->acc.z += R.z / r * a1;
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

	double dt = header->dt;

	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		int i = ind;

		__global struct Body * b = &bodies[i];

		struct Vec3 v0 = b->vel;
		struct Vec3 a0 = b->acc;

		b->vel.x += b->acc.x * dt;
		b->vel.y += b->acc.y * dt;
		b->vel.z += b->acc.z * dt;

		if (0)
		{
			b->pos.x += b->vel.x * dt;
			b->pos.y += b->vel.y * dt;
			b->pos.z += b->vel.z * dt;
		}
		else{
			b->pos.x += v0.x * dt + a0.x / 2 * dt * dt;
			b->pos.y += v0.y * dt + a0.y / 2 * dt * dt;
			b->pos.z += v0.z * dt + a0.z / 2 * dt * dt;
		}

		b->acc.x = 0;
		b->acc.y = 0;
		b->acc.z = 0;
	}
}





