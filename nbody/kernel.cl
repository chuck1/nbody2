
#include "kernel.h"

#define TIME_STEP_FACTOR_PROXIMITY (0.001)
#define TIME_STEP_FACTOR_ACC (0.1)

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

double			vector_length_p(struct Vec3 * v)
{
	double l = 0;
	for (int i = 0; i < 3; ++i) l += v->v[i] * v->v[i];
	return sqrt(l);
}
double			vector_length_g(__global struct Vec3 * v)
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

__kernel void k_sum_uint(
	__global const unsigned int * input,
	__global unsigned int * plen,
	__global unsigned int * output,
	__local unsigned int * loc)
{
	/*
	length of partial must be at least number of work groups
	final result is in partial[0]
	*/

	uint local_id = get_local_id(0);
	uint global_id = get_global_id(0);
	uint local_size = get_local_size(0);
	uint global_size = get_global_size(0);

	unsigned int len = *plen;

	if (global_id > (len - 1)) return;

	// Copy from global to local memory
	loc[local_id] = input[global_id];

	// Loop for computing localSums : divide WorkGroup into 2 parts
	for (uint stride = local_size / 2; stride > 0; stride /= 2)
	{
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (local_id > (stride - 1)) break;

		if ((global_id + stride) > (len - 1)) continue;

		loc[local_id] = loc[local_id] + loc[local_id + stride];
	}

	//barrier(CLK_LOCAL_MEM_FENCE);

	// write result into partial[nWorkGroups]
	if (local_id == 0) output[get_group_id(0)] = loc[0];
}

__kernel void k_min_double(
	__global const double * input,
	__global unsigned int * plen,
	__global double * output,
	__local double * loc)
{
	/*
	length of partial must be at least number of work groups
	final result is in partial[0]
	*/

	uint local_id = get_local_id(0);
	uint global_id = get_global_id(0);
	uint local_size = get_local_size(0);
	uint global_size = get_global_size(0);

	unsigned int len = *plen;

	if (global_id > (len - 1)) return;

	// Copy from global to local memory
	loc[local_id] = input[global_id];

	// Loop for computing localSums : divide WorkGroup into 2 parts
	for (uint stride = local_size / 2; stride > 0; stride /= 2)
	{
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);

		if (local_id > (stride - 1)) break;

		if ((global_id + stride) > (len - 1)) continue;

		loc[local_id] = min(loc[local_id], loc[local_id + stride]);
	}

	//barrier(CLK_LOCAL_MEM_FENCE);

	// write result into partial[nWorkGroups]
	if (local_id == 0) output[get_group_id(0)] = loc[0];
}



__kernel void store_dt(
	__global double * dt_output,
	__global struct Header * header
	)
{
	if (get_global_id(0) == 0)
	{
		header->dt = dt_output[0];

		// override

		//header->dt = 1.0E0;
	}
}


void calc_acc_sub(
	__global struct Pair * p, 
	__global struct Header * header, 
	__global struct Body0 * bodies0,
	__global struct Body1 * bodies1,
	__global double * dt_input, 
	int ind)
{
	struct Vec3 temp0;

	double G = 6.67408E-11;

	double k_drag = 1.0E+0;
	double k_repulsion = 0.01;

	__global struct Body0 * b0_0 = &bodies0[p->i];
	__global struct Body1 * b1_0 = &bodies1[p->i];

	__global struct Body0 * b0_1 = &bodies0[p->j];
	__global struct Body1 * b1_1 = &bodies1[p->j];

	struct Vec3 R;

	vector_sub_pgg(&R, &b0_1->pos, &b0_0->pos);

	double r = vector_length_p(&R);

	struct Vec3 v;

	vector_sub_pgg(&v, &b1_1->vel, &b1_0->vel);

	double s = vector_length_p(&v);

	// save
#ifdef DEBUG
	p->d = r;
	p->s = s;
#endif

	// calc dt
	if (s == 0)
	{
		dt_input[ind] = 1.0;
	}
	else
	{
		dt_input[ind] = min(TIME_STEP_FACTOR_PROXIMITY * r / s, 100.0);
	}

	double m0 = b1_0->mass;
	double m1 = b1_1->mass;

	double F = G * m0 * m1 / r / r;

	// positive pen is penetration
	// negative pen is no penetration

	double pen = b0_0->radius + b0_1->radius - r;

	if (pen > 0)
	{
		atomic_inc(&header->count_pen);

		// volume and mass of interseting region
		double V = volume_sphere_sphere_intersection(r, b0_0->radius, b0_1->radius);
		double m = V * (b1_0->density + b1_1->density);

		// point at which drag force will act
		struct Vec3 o;
		vector_add_ppg(&o, vector_mul_pp(&temp0, &R, (b0_0->radius - pen / 2.0) / r), &b0_0->pos);

		// drag force
		struct Vec3 F_drag0, F_drag1;

		struct Vec3 R0, R1;

		vector_sub_ppg(&R0, &o, &b0_0->pos);
		vector_sub_ppg(&R1, &o, &b0_1->pos);

		struct Vec3 vo0, vo1, vo_rel;

		vector_add_ppg(&vo0, vector_cross_pgp(&temp0, &b1_0->w, &R0), &b1_0->vel);
		vector_add_ppg(&vo1, vector_cross_pgp(&temp0, &b1_1->w, &R1), &b1_1->vel);

		vector_sub_ppp(&vo_rel, &vo0, &vo1);



		vector_mul_pp(&F_drag0, &vo_rel, -k_drag * m);
		vector_mul_pp(&F_drag1, &vo_rel, +k_drag * m);

		// torque
		struct Vec3 T0, T1;

		vector_cross_ppp(&T0, &R0, &F_drag0);
		vector_cross_ppp(&T1, &R1, &F_drag1);


		// repulsion

		if (1) {
			F -= k_repulsion * m;
		}

		vector_add_self_gp(&b1_0->force, &F_drag0);
		vector_add_self_gp(&b1_1->force, &F_drag1);

		vector_add_self_gp(&b1_0->torque, &T0);
		vector_add_self_gp(&b1_1->torque, &T1);
	}

	// gravity and repulsion
	vector_add_self_gp(&b1_0->force, vector_mul_pp(&temp0, &R, F / r));
	vector_sub_self_gp(&b1_1->force, vector_mul_pp(&temp0, &R, F / r));

#ifdef KDEBUG
	p->F = F;
#endif
}

int get_pair(__global struct Pair * pairs, int p, int i, int j)
{
	for (int k = 0; k < p; ++k)
	{
		if ((pairs[k].i == i) && (pairs[k].j == j)) return k;
		if ((pairs[k].i == j) && (pairs[k].j == i)) return k;
	}
	return -1;
}

__kernel void calc_acc(
	__global struct Header * header,
	__global struct Body0 * bodies0,
	__global struct Body1 * bodies1,
	__global struct Pair * pairs,
	__global double * dt_input)
{
	/* MUST USE ONLY ONE WORK GROUP
	*/

	if (get_global_id(0) != get_local_id(0)) return;

	if (get_local_id(0) == 0) {
		header->dt = 0;
		header->count_pen = 0;
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

	int n = header->bodies_size;
	int p = n * (n - 1) / 2;

	/*unsigned int len = n * (n - 1) / 2;

	while (true)
	{
		unsigned int ind = atomic_inc(counter);
		if (ind > (len - 1)) break;

		calc_acc_sub(&pairs[ind], header, bodies, dt_input, ind);
	}*/

	/* NEW METHOD
	diagonals of the body-body table
	*/
	

	int local_id = get_local_id(0);
	int local_size = get_local_size(0);

	for (int i = 0; i < (n - 1); i += local_size)
	{
		for (int j = 1; j < n; ++j)
		{
			barrier(CLK_GLOBAL_MEM_FENCE);

			int n0 = min(local_size, n - j);

			if (local_id > (n0 - 1)) continue;

			int i0 = i + local_id;
			int j0 = j + local_id;

			int k = get_pair(pairs, p, i0, j0);

			calc_acc_sub(&pairs[k], header, bodies0, bodies1, dt_input, k);
		}
	}
}

__kernel void dt_calc(
	__global struct Header * header,
	__global struct Body0 * bodies0,
	__global struct Body1 * bodies1,
	__global double * dt_input)
{
	int n = header->bodies_size;
	int p = n * (n - 1) / 2;

	unsigned int ind = get_global_id(0);
	
	if (ind > (n - 1)) return;

	__global struct Body0 * b0 = &bodies0[ind];
	__global struct Body1 * b1 = &bodies1[ind];

	double f = vector_length_g(&b1->force);
	double a = f / b1->mass;
	double v = vector_length_g(&b1->vel);
	
	if (a == 0)
	{
		dt_input[p + ind] = 0.0;
	}
	else
	{
		dt_input[p + ind] = TIME_STEP_FACTOR_ACC * v / a;
			
		dt_input[p + ind] = max(dt_input[p + ind], 1.0E-5);
	}
}

__kernel void step_pos(
	__global struct Header * header,
	__global struct Body0 * bodies0,
	__global struct Body1 * bodies1,
	__global struct Pair * pairs)
{
	int len = header->bodies_size;

	double dt = header->dt;

	int ind = get_global_id(0);

	if (ind > (len - 1)) return;

	__global struct Body0 * b0 = &bodies0[ind];
	__global struct Body1 * b1 = &bodies1[ind];

	struct Vec3 v0 = b1->vel;
	struct Vec3 a0 = b1->acc;

	struct Vec3 temp0, temp1, temp2;

	vector_mul_gg(&b1->acc, &b1->force, 1.0 / b1->mass);

	vector_add_self_gp(&b1->vel, vector_mul_pp(&temp0, &a0, dt));

	vector_add_self_gp(&b0->pos, vector_add_ppp(&temp0, vector_mul_pp(&temp1, &v0, dt), vector_mul_pp(&temp2, &a0, dt * dt / 2.0)));

	b1->force.v[0] = 0;
	b1->force.v[1] = 0;
	b1->force.v[2] = 0;

	// rotation

	struct Quat tempq0, tempq1;

	struct Mat43 B0;
	struct Mat43 B1;
		
	struct Vec4 tempv40, tempv41, tempv42;

	construct_B_quat_pg(&B0, &b0->q);
		
	vec4_mul_gp(&b1->q1, mat43_mul_ppg(&tempv40, &B0, &b1->w), 0.5);

	construct_B_vec4_pg(&B1, &b1->q1);

	struct Vec3 w1;

	//mat43_mul_ppg(&tempv41, &B1, &b->w);
	//mat43_mul_ppp(&tempv42, &B0, &w1);
	//vec4_mul_gp(&b->q2, vec4_add_ppp(&tempv40, &tempv41, &tempv42), 0.5);

	for (int i = 0; i < 4; ++i) b0->q.v[i] += b1->q1.v[i] * dt;

	quat_normalize_g(&b0->q);

	// update rotational velocity

	__global struct Vec3 * w = &b1->w;
	struct Vec3 alpha;
	struct Vec3 * a = &alpha;

	struct Mat3 I, I_inv;
	double I_scalar = 2.0 / 5.0 * b1->mass * pow(b0->radius, 2.0);
	mat3_diagonal(&I, I_scalar);
	mat3_diagonal(&I_inv, 1.0 / I_scalar);

	mat3_mul_ppp(a, &I_inv, vector_sub_pgp(&temp0, &b1->torque, vector_cross_pgp(&temp1, w, mat3_mul_ppg(&temp2, &I, w))));

	vector_add_self_gp(w, vector_mul_pp(&temp0, a, dt));

		

	b1->torque.v[0] = 0;
	b1->torque.v[1] = 0;
	b1->torque.v[2] = 0;
	


	if (get_global_id(0) == 0)
	{
		header->t += header->dt;
	}
}





