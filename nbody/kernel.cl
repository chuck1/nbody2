
#include "kernel.h"

#define TIME_STEP_FACTOR (0.0001)

struct Vec3 *	vector_cross_ppp(struct Vec3 * ret, struct Vec3 * a, struct Vec3 * b)
{
	ret->v[0] = a->v[1] * b->v[2] - a->v[2] * b->v[1];
	ret->v[1] = a->v[2] * b->v[0] - a->v[0] * b->v[2];
	ret->v[2] = a->v[0] * b->v[1] - a->v[1] * b->v[0];
	return ret;
}
struct Vec3 *	vector_cross_pgp(struct Vec3 * ret, __global struct Vec3 * a, struct Vec3 * b)
{
	ret->v[0] = a->v[1] * b->v[2] - a->v[2] * b->v[1];
	ret->v[1] = a->v[2] * b->v[0] - a->v[0] * b->v[2];
	ret->v[2] = a->v[0] * b->v[1] - a->v[1] * b->v[0];
	return ret;
}

struct Vec3 *	vector_sub_pgg(struct Vec3 * ret, __global struct Vec3 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] - b->v[i];
	return ret;
}
struct Vec3 *	vector_sub_pgp(struct Vec3 * ret, __global struct Vec3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] - b->v[i];
	return ret;
}
struct Vec3 *	vector_sub_ppg(struct Vec3 * ret, struct Vec3 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] - b->v[i];
	return ret;
}
struct Vec3 *	vector_sub_ppp(struct Vec3 * ret, struct Vec3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] - b->v[i];
	return ret;
}

void			vector_sub_self_gg(__global struct Vec3 * a,__global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) a->v[i] -= b->v[i];
}
void			vector_sub_self_gp(__global struct Vec3 * a,struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) a->v[i] -= b->v[i];
}

struct Vec3 *	vector_add_ppp(struct Vec3 * ret, struct Vec3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] + b->v[i];
	return ret;
}
struct Vec3 *	vector_add_ppg(struct Vec3 * ret, struct Vec3 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] + b->v[i];
	return ret;
}

void			vector_add_self_gg(__global struct Vec3 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) a->v[i] += b->v[i];
}
void			vector_add_self_gp(__global struct Vec3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i) a->v[i] += b->v[i];
}

__global struct Vec3 *	vector_mul_gg(__global struct Vec3 * ret, __global struct Vec3 * a, double b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] * b;
	return ret;
}
struct Vec3 *			vector_mul_pg(struct Vec3 * ret, __global struct Vec3 * a, double b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] * b;
	return ret;
}
struct Vec3 *			vector_mul_pp( struct Vec3 * ret, struct Vec3 * a, double b)
{
	for (int i = 0; i < 3; ++i) ret->v[i] = a->v[i] * b;
	return ret;
}

struct Vec4 *			vec4_mul_pp(struct Vec4 * ret, struct Vec4 * a, double b)
{
	for (int i = 0; i < 4; ++i) ret->v[i] = a->v[i] * b;
	return ret;
}
__global struct Vec4 *	vec4_mul_gp(__global struct Vec4 * ret, struct Vec4 * a, double b)
{
	for (int i = 0; i < 4; ++i) ret->v[i] = a->v[i] * b;
	return ret;
}

void			quat_add_self_gp(__global struct Quat * a, struct Quat * b)
{
	for (int i = 0; i < 4; ++i) a->v[i] += b->v[i];
}

struct Quat *	quat_mul_ppg(struct Quat * ret, struct Quat * a, __global struct Quat * b)
{
	ret->v[0] = a->v[0] * b->v[0] - a->v[1] * b->v[1] - a->v[2] * b->v[2] - a->v[3] * b->v[3];
	ret->v[1] = a->v[0] * b->v[1] + a->v[1] * b->v[0] + a->v[2] * b->v[3] - a->v[3] * b->v[2];
	ret->v[2] = a->v[0] * b->v[2] - a->v[1] * b->v[3] + a->v[2] * b->v[0] + a->v[3] * b->v[1];
	ret->v[3] = a->v[0] * b->v[3] + a->v[1] * b->v[2] - a->v[2] * b->v[1] + a->v[3] * b->v[0];
	return ret;
}

void			quat_mul_self_g(__global struct Quat * a, double b)
{
	for(int i = 0; i < 4; ++i) a->v[i] *= b;
}

struct Quat *	cast_vec3_to_quat_pp( struct Quat * ret, struct Vec3 * a)
{
	ret->v[0] = 0;
	ret->v[1] = a->v[0];
	ret->v[2] = a->v[1];
	ret->v[3] = a->v[2];
	return ret;
}

double			vector_length(struct Vec3 * v)
{
	double l = 0;
	for (int i = 0; i < 3; ++i) l += v->v[i] * v->v[i];
	return sqrt(l);
}

double			quat_length_g(__global struct Quat * v)
{
	double l = 0;
	for (int i = 0; i < 4; ++i) l += v->v[i] * v->v[i];
	return sqrt(l);
}
void			quat_normalize_g(__global struct Quat * v)
{
	double l = quat_length_g(v);
	if (l == 0) return;
	quat_mul_self_g(v, 1.0 / l);
}

struct Vec4 *	mat43_mul_ppg(struct Vec4 * ret, struct Mat43 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 4; ++i)
	{
		ret->v[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			ret->v[i] += a->v[i][j] * b->v[j];
		}
	}
	return ret;
}
struct Vec4 *	mat43_mul_ppp(struct Vec4 * ret, struct Mat43 * a, struct Vec3 * b)
{
	for (int i = 0; i < 4; ++i)
	{
		ret->v[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			ret->v[i] += a->v[i][j] * b->v[j];
		}
	}
	return ret;
}

struct Vec3 *	mat3_mul_pgp(struct Vec3 * ret, __global struct Mat3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i)
	{
		ret->v[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			ret->v[i] += a->v[i][j] * b->v[j];
		}
	}
	return ret;
}
struct Vec3 *	mat3_mul_ppg(struct Vec3 * ret, struct Mat3 * a, __global struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i)
	{
		ret->v[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			ret->v[i] += a->v[i][j] * b->v[j];
		}
	}
	return ret;
}
struct Vec3 *	mat3_mul_ppp(struct Vec3 * ret, struct Mat3 * a, struct Vec3 * b)
{
	for (int i = 0; i < 3; ++i)
	{
		ret->v[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			ret->v[i] += a->v[i][j] * b->v[j];
		}
	}
	return ret;
}

struct Mat43 *	construct_B_quat_pg(struct Mat43 * ret, __global struct Quat * q)
{
	ret->v[0][0] = -q->v[1];
	ret->v[0][1] = -q->v[2];
	ret->v[0][2] = -q->v[3];
	
	ret->v[1][0] = +q->v[0];
	ret->v[1][1] = +q->v[3];
	ret->v[1][2] = -q->v[2];
	
	ret->v[2][0] = -q->v[3];
	ret->v[2][1] = +q->v[0];
	ret->v[2][2] = +q->v[1];

	ret->v[3][0] = +q->v[2];
	ret->v[3][1] = -q->v[1];
	ret->v[3][2] = +q->v[0];

	return ret;
}
struct Mat43 *	construct_B_vec4_pg(struct Mat43 * ret, __global struct Vec4 * q)
{
	ret->v[0][0] = -q->v[1];
	ret->v[0][1] = -q->v[2];
	ret->v[0][2] = -q->v[3];

	ret->v[1][0] = +q->v[0];
	ret->v[1][1] = +q->v[3];
	ret->v[1][2] = -q->v[2];

	ret->v[2][0] = -q->v[3];
	ret->v[2][1] = +q->v[0];
	ret->v[2][2] = +q->v[1];

	ret->v[3][0] = +q->v[2];
	ret->v[3][1] = -q->v[1];
	ret->v[3][2] = +q->v[0];

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

double			CDF_uniform(double a, double b, double x)
{
	if (x < a) return 0;
	if (x > b) return 1;
	return (x - a) / (b - a);
}

void			mat3_diagonal(struct Mat3 *	a, double b)
{
	a->v[0][0] = b;
	a->v[0][1] = 0;
	a->v[0][2] = 0;
	a->v[1][0] = 0;
	a->v[1][1] = b;
	a->v[1][2] = 0;
	a->v[2][0] = 0;
	a->v[2][1] = 0;
	a->v[2][2] = b;
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

		// override

		//header->dt = 1.0E0;
	}
}

__kernel void calc_acc(
	__global struct Header * header,
	__global struct Body * bodies,
	__global struct Pair * pairs,
	__global double * dt_input,
	volatile __global unsigned int * counter)
{
	struct Vec3 temp0;

	if (get_global_id(0) == 0) {
		header->dt = 0;
		header->count_pen = 0;
		*counter = 0;
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

	int n = header->bodies_size;
	
	unsigned int len = n * (n - 1) / 2;

	double G = 6.67408E-11;

	double k_drag = 1.0E-2;
	
	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		__global struct Pair * p = &pairs[ind];

		int i = p->i;
		int j = p->j;

		__global struct Body * b0 = &bodies[i];
		__global struct Body * b1 = &bodies[j];

		struct Vec3 R;
		
		vector_sub_pgg(&R, &b1->pos, &b0->pos);

		double r = vector_length(&R);
		
		struct Vec3 v;
		
		vector_sub_pgg(&v, &b1->vel, &b0->vel);
		
		double s = vector_length(&v);

		// save
		p->d = r;
		p->s = s;

		// calc dt
		if (s == 0)
		{
			dt_input[ind] = 1.0;
		}
		else
		{
			dt_input[ind] = min(TIME_STEP_FACTOR * r / s, 100.0);
		}

		double m0 = b0->mass;
		double m1 = b1->mass;

		double F = G * m0 * m1 / r / r;
		
		// positive pen is penetration
		// negative pen is no penetration

		double pen = b0->radius + b1->radius - p->d;

		if (pen > 0)
		{
			atomic_inc(&header->count_pen);
			
			// volume and mass of interseting region
			double V = volume_sphere_sphere_intersection(p->d, b0->radius, b1->radius);
			double m = V * (b0->density + b1->density);

			// point at which drag force will act
			struct Vec3 o;
			vector_add_ppg(&o, vector_mul_pp(&temp0, &R, (b0->radius - pen / 2.0) / r), &b0->pos);

			// drag force
			struct Vec3 F_drag0, F_drag1;

			struct Vec3 R0, R1;

			vector_sub_ppg(&R0, &o, &b0->pos);
			vector_sub_ppg(&R1, &o, &b1->pos);

			struct Vec3 vo0, vo1, vo_rel;

			vector_add_ppg(&vo0, vector_cross_pgp(&temp0, &b0->w, &R0), &b0->vel);
			vector_add_ppg(&vo1, vector_cross_pgp(&temp0, &b1->w, &R1), &b1->vel);

			vector_sub_ppp(&vo_rel, &vo0, &vo1);

			

			vector_mul_pp(&F_drag0, &vo_rel, -k_drag * m);
			vector_mul_pp(&F_drag1, &vo_rel, +k_drag * m);

			// torque
			struct Vec3 T0, T1;
			
			vector_cross_ppp(&T0, vector_sub_ppg(&temp0, &o, &b0->pos), &F_drag0);
			vector_cross_ppp(&T1, vector_sub_ppg(&temp0, &o, &b1->pos), &F_drag1);

			 
			// repulsion

			double k = 0.01;
			
			if (1) {
				F -= k * m;
			}

			vector_add_self_gp(&b0->force, &F_drag0);
			vector_add_self_gp(&b1->force, &F_drag1);

			vector_add_self_gp(&b0->torque, &T0);
			vector_add_self_gp(&b1->torque, &T1);
		}

		// gravity and repulsion
		vector_add_self_gp(&b0->force, vector_mul_pp(&temp0, &R, F / r));
		vector_sub_self_gp(&b1->force, vector_mul_pp(&temp0, &R, F / r));

		p->F = F;
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

		struct Vec3 temp0, temp1, temp2;

		vector_mul_gg(&b->acc, &b->force, 1.0 / b->mass);

		vector_add_self_gp(&b->vel, vector_mul_pg(&temp0, &b->acc, dt));

		vector_add_self_gp(&b->pos, vector_add_ppp(&temp0, vector_mul_pp(&temp1, &v0, dt), vector_mul_pp(&temp2, &a0, dt * dt / 2.0)));

		b->force.v[0] = 0;
		b->force.v[1] = 0;
		b->force.v[2] = 0;

		// rotation

		struct Quat tempq0, tempq1;

		struct Mat43 B0;
		struct Mat43 B1;
		
		struct Vec4 tempv40, tempv41, tempv42;

		construct_B_quat_pg(&B0, &b->q);
		
		vec4_mul_gp(&b->q1, mat43_mul_ppg(&tempv40, &B0, &b->w), 0.5);

		construct_B_vec4_pg(&B1, &b->q1);

		struct Vec3 w1;

		//mat43_mul_ppg(&tempv41, &B1, &b->w);
		//mat43_mul_ppp(&tempv42, &B0, &w1);
		//vec4_mul_gp(&b->q2, vec4_add_ppp(&tempv40, &tempv41, &tempv42), 0.5);

		for (int i = 0; i < 4; ++i) b->q.v[i] += b->q1.v[i] * dt;

		quat_normalize_g(&b->q);

		// update rotational velocity

		__global struct Vec3 * w = &b->w;
		struct Vec3 alpha;
		struct Vec3 * a = &alpha;

		struct Mat3 I, I_inv;
		double I_scalar = 2.0 / 5.0 * b->mass * pow(b->radius, 2.0);
		mat3_diagonal(&I, I_scalar);
		mat3_diagonal(&I_inv, 1.0 / I_scalar);

		mat3_mul_ppp(a, &I_inv, vector_sub_pgp(&temp0, &b->torque, vector_cross_pgp(&temp1, w, mat3_mul_ppg(&temp2, &I, w))));

		vector_add_self_gp(w, vector_mul_pp(&temp0, a, dt));

		

		b->torque.v[0] = 0;
		b->torque.v[1] = 0;
		b->torque.v[2] = 0;
	}
}





