// nbody.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <ctime>

#include <Windows.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>


#include "OCL.h"
#include "kernel.h"

double G = 6.67408E-11;
double pi = 3.1415926535897;

int m = 100000;
int n = 20;
int p = n * (n - 1) / 2;


class Frame
{
public:
	Header header;
	std::vector<Body> bodies;
	std::vector<Pair> pairs;
};

class History
{
public:
	void write()
	{
		std::ofstream myFile;
		
		myFile.open("P:/temp/nbody.bin", std::ios::out | std::ios::binary);

		unsigned int l = frames.size();

		myFile.write((const char *)&l, sizeof(unsigned int));
	}

	void push(Header & header, std::shared_ptr<OCL::MemObj> memobj_bodies, std::shared_ptr<OCL::MemObj> memobj_pairs)
	{
		frames.emplace_back();
		Frame & frame = frames.back();

		frame.header = header;
		frame.bodies.resize(n);
		frame.pairs.resize(p);
		memobj_bodies->EnqueueRead(&frame.bodies[0], n * sizeof(Body));
		memobj_pairs->EnqueueRead(&frame.pairs[0], p * sizeof(Pair));
	}

	std::vector<Frame> frames;
};

void generate_binary_system(
	std::vector<Body> & bodies,
	std::vector<Pair> & pairs)
{
	int n = bodies.size();

	Body & b0 = bodies[0];
	Body & b1 = bodies[1];

	double & x0 = b0.pos.x;
	double & x1 = b1.pos.x;

	double & m0 = b0.mass;
	double & m1 = b1.mass;

	m0 = 1.0E11;
	m1 = 1.0E11;

	double d = 200.0;

	x0 = d * m1 / (m0 + m1);
	x1 = x0 - d;

	double f = 0.2;

	b0.vel.y = f * sqrt(G * m1 * x0 / d / d);
	b1.vel.y = -f * sqrt(G * m1 * -x1 / d / d);
}

void generate_pairs(
	std::vector<Body> & bodies,
	std::vector<Pair> & pairs)
{
	int n = bodies.size();

	int k = 0;
	for (int i = 0; i < (n - 1); ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			pairs[k].i = i;
			pairs[k].j = j;
			++k;
		}
	}
}

void generate_bodies(
	std::vector<Body> & bodies,
	std::vector<Pair> & pairs)
{
	int n = bodies.size();

	for (int i = 0; i < n; ++i)
	{
		double scale = 500;
		double x = (float)rand() / (float)RAND_MAX - 0.5;
		double y = (float)rand() / (float)RAND_MAX - 0.5;
		x = x * scale;
		y = y * scale;

		double r = sqrt(x*x + y*y);

		Body & b = bodies[i];

		b.mass = 1E8;
		b.radius = pow(3.0 * b.mass / 4.0 / pi / 5000.0, 1.0/3.0);

		b.pos.x = x;
		b.pos.y = y;

	}


	double X = 0;
	double Y = 0;
	double Z = 0;
	double M = 0;
	for (int i = 0; i < n; ++i)
	{
		Body & b = bodies[i];
		X += b.mass * b.pos.x;
		Y += b.mass * b.pos.y;
		Z += b.mass * b.pos.z;
		M += b.mass;
	}
	X /= M;
	Y /= M;
	Z /= M;

	for (int i = 0; i < n; ++i)
	{
		Body & b = bodies[i];

		double x = b.pos.x - X;
		double y = b.pos.y - Y;

		double r = sqrt(x*x + y*y);

		if (r == 0)
		{
			b.vel.y = 0;
			b.vel.x = 0;
		}
		else
		{
			double s = 0.6 * sqrt(G * M / r);

			b.vel.y = s * x / r;
			b.vel.x = -s * y / r;
		}
	}

}

void draw_body(Body & b)
{
	glPushMatrix();
	{
		glTranslated(b.pos.x, b.pos.y, b.pos.z);

		double s = b.radius;
		glScaled(s,s,s);

		glColor3f(1.0f, 1.0f, 1.0f);

		glBegin(GL_QUADS);
		
		// top
		glNormal3f( 0.0f, 1.0f, 0.0f);
		glVertex3f(-1.0f, 1.0f, 1.0f);
		glVertex3f( 1.0f, 1.0f, 1.0f);
		glVertex3f( 1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);

		glNormal3f( 0.0f, -1.0f,  0.0f);
		glVertex3f(-1.0f, -1.0f,  1.0f);
		glVertex3f( 1.0f, -1.0f,  1.0f);
		glVertex3f( 1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);

		// front
		glNormal3f( 0.0f,  0.0f, 1.0f);
		glVertex3f( 1.0f, -1.0f, 1.0f);
		glVertex3f( 1.0f,  1.0f, 1.0f);
		glVertex3f(-1.0f,  1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);

		glNormal3f(0.0f, 0.0f, -1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);
		glVertex3f(1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);

		// right
		glNormal3f(1.0f,  0.0f,  0.0f);
		glVertex3f(1.0f,  1.0f, -1.0f);
		glVertex3f(1.0f,  1.0f,  1.0f);
		glVertex3f(1.0f, -1.0f,  1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);

		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);

		glEnd();
	}
	glPopMatrix();
}
void draw_frame(Frame & frame)
{
	for (int i = 0; i < frame.bodies.size(); ++i)
	{
		draw_body(frame.bodies[i]);
	}
}


class Plot
{
public:
	double t_factor;
	double t;
	unsigned int i;
	bool pause;
} plot;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_Q && action == GLFW_PRESS)
	{	plot.t_factor *= 10.0;
	}

	if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{
		plot.t_factor /= 10.0;
	}

	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
	{
		plot.t = 0;
		plot.i = 0;
	}

	if (key == GLFW_KEY_P && action == GLFW_PRESS)
	{
		plot.pause = true;
	}
}
void plotfunc(History & hist)
{
	plot.t_factor = 1.0E1;
	plot.t = 0;
	plot.i = 0;
	plot.pause = false;

	GLFWwindow* window;

	if (!glfwInit())
	{
		exit(EXIT_FAILURE);
	}

	//Create a window and create its OpenGL context
	window = glfwCreateWindow(200, 200, "Test Window", NULL, NULL);
	if (!window)
	{
		fprintf(stderr, "Failed to open GLFW window.\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwSetKeyCallback(window, key_callback);

	glfwMakeContextCurrent(window);

	gluPerspective(90, 1, 100, 10000);

	gluLookAt(
		0, 0, 400, 
		0, 0, 0, 
		0, 1, 0);

	double t0 = glfwGetTime();

	do
	{
		double t1 = glfwGetTime();
		double dt = t1 - t0;
		t0 = t1;

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		Frame & f = hist.frames[plot.i];

		glPushMatrix();
		{
			draw_frame(f);	
		}
		glPopMatrix();

		glfwSwapBuffers(window);
		glfwPollEvents();

		if (!plot.pause)
		{
			plot.t += dt * plot.t_factor;

			while (f.header.t < plot.t)
			{
				plot.i += 1;

				if (plot.i >= hist.frames.size()) {
					plot.i = 0;
					plot.t = 0;
				}

				f = hist.frames[plot.i];
			}
		}


		printf("t_sim=%8f i=%6i t=%16f dt=%16f tf=%8.2e\n", plot.t, plot.i, f.header.t, f.header.dt, plot.t_factor);

	} while (!glfwWindowShouldClose(window));

	
}



unsigned int next_power_of_two(unsigned int x)
{
	int ret = 2;
	while (ret < x)
	{
		ret *= 2;
	}
	return ret;
}

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned int file_size = m * (n * sizeof(Body) + p * sizeof(Pair));

	printf("file size = %u MB\n", file_size / 1024 / 1024);

	srand(time(NULL));

	double * dt_partial = new double[p];

	History hist;

	std::vector<Body> bodies(n);
	std::vector<Pair> pairs(n*(n-1)/2);

	//generate_binary_system(bodies, pairs);
	generate_bodies(bodies, pairs);
	generate_pairs(bodies, pairs);

	Header header;
	header.bodies_size = n;
	header.t = 0;
	header.dt = 0;

	auto ocl = std::make_shared<OCL::Manager>();
	ocl->init();

	char * files[2] = { "kernel.cl", "kernel.h" };

	auto program = ocl->create_program(files, 1, "");

	auto kernel_calc_acc = program->create_kernel("calc_acc");
	auto kernel_step = program->create_kernel("step_pos");
	auto kernel_dt_min = program->create_kernel("k_min");
	auto kernel_dt_store = program->create_kernel("store_dt");

	auto memobj_bodies = ocl->create_buffer(CL_MEM_READ_WRITE, bodies.size() * sizeof(Body));
	memobj_bodies->EnqueueWrite(&bodies[0], bodies.size() * sizeof(Body));

	auto memobj_pairs = ocl->create_buffer(CL_MEM_READ_WRITE, pairs.size() * sizeof(Pair));
	memobj_pairs->EnqueueWrite(&pairs[0], pairs.size() * sizeof(Pair));

	auto memobj_header = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(Header));
	memobj_header->EnqueueWrite(&header, sizeof(Header));

	unsigned int counter = 0;
	auto memobj_counter = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(unsigned int));
	memobj_counter->EnqueueWrite(&counter, sizeof(unsigned int));

	auto memobj_dt_input = ocl->create_buffer(CL_MEM_READ_WRITE, p * sizeof(double));

	auto memobj_dt_partial = ocl->create_buffer(CL_MEM_READ_WRITE, p * sizeof(double));

	auto memobj_dt_len = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(unsigned int));
	memobj_dt_len->EnqueueWrite(&p, sizeof(unsigned int));

	kernel_calc_acc->set_arg(memobj_header, 0);
	kernel_calc_acc->set_arg(memobj_bodies, 1);
	kernel_calc_acc->set_arg(memobj_pairs, 2);
	kernel_calc_acc->set_arg(memobj_dt_input, 3);
	kernel_calc_acc->set_arg(memobj_counter, 4);

	kernel_step->set_arg(memobj_header, 0);
	kernel_step->set_arg(memobj_bodies, 1);
	kernel_step->set_arg(memobj_pairs, 2);
	kernel_step->set_arg(memobj_counter, 3);

	/*
	__global const double * input,
		__global unsigned int * len,
		__global double * partial,
		__local double * loc)
		*/
	kernel_dt_min->set_arg(memobj_dt_input, 0);
	kernel_dt_min->set_arg(memobj_dt_len, 1);
	kernel_dt_min->set_arg(memobj_dt_partial, 2);
	kernel_dt_min->set_arg(3, 1024 * sizeof(double));

	kernel_dt_store->set_arg(memobj_dt_partial, 0);
	kernel_dt_store->set_arg(memobj_header, 1);

	std::cout << "kernel start" << std::endl;

	hist.push(header, memobj_bodies, memobj_pairs);

	

	for (int i = 0; i < m; ++i)
	{
		double t0 = header.t;

		kernel_calc_acc->enqueue_ND_range_kernel(next_power_of_two(p), 1);
		kernel_dt_min->enqueue_ND_range_kernel(next_power_of_two(p), 256);

		kernel_dt_store->enqueue_ND_range_kernel(1, 1);

		kernel_step->enqueue_ND_range_kernel(next_power_of_two(n), 1);


		memobj_header->EnqueueRead(&header, sizeof(Header));

		header.t = t0 + header.dt;

		if (i % (m / 100) == 0) printf("%8i t=%16f dt=%16f count_pen=%8i\n", i, header.t, header.dt, header.count_pen);

		

		// save
		hist.push(header, memobj_bodies, memobj_pairs);
	}

	hist.write();

	std::cout << "read" << std::endl;

	memobj_bodies->EnqueueRead(&bodies[0], n * sizeof(Body));

	for (int i = 0; i < n; ++i)
	{
		printf("%4i %16f %16f %16f\n     %16f %16f %16f\n     %16f %16f %16f\n", i, 
			bodies[i].pos.x, 
			bodies[i].pos.y,
			bodies[i].pos.z,
			bodies[i].vel.x,
			bodies[i].vel.y,
			bodies[i].vel.z,
			bodies[i].acc.x,
			bodies[i].acc.y,
			bodies[i].acc.z
			);
	}

	std::cout << "done" << std::endl;

	ocl->shutdown();

	plotfunc(hist);

	getchar();

	return 0;
}



