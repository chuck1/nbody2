// nbody.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>

#include <Windows.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>


#include "OCL.h"
#include "kernel.h"

double G = 6.67408E-11;
double pi = 3.1415926535897;

class Frame
{
public:
	std::vector<Body> bodies;
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

	b0.vel.y = sqrt(G * m1 * x0 / d / d);
	b1.vel.y = -sqrt(G * m1 * -x1 / d / d);
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

			printf("%4i %4i %4i\n", k, i, j);
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

		b.mass = 1E11;
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

		glScaled(10, 10, 10);

		glColor3f(1.0f, 1.0f, 1.0f);

		glBegin(GL_QUADS);
		
		// top
		glNormal3f( 0.0f, 1.0f, 0.0f);
		glVertex3f(-0.5f, 0.5f, 0.5f);
		glVertex3f( 0.5f, 0.5f, 0.5f);
		glVertex3f( 0.5f, 0.5f, -0.5f);
		glVertex3f(-0.5f, 0.5f, -0.5f);

		glNormal3f( 0.0f, -1.0f,  0.0f);
		glVertex3f(-0.5f, -0.5f,  0.5f);
		glVertex3f( 0.5f, -0.5f,  0.5f);
		glVertex3f( 0.5f, -0.5f, -0.5f);
		glVertex3f(-0.5f, -0.5f, -0.5f);

		// front
		glNormal3f( 0.0f,  0.0f, 1.0f);
		glVertex3f( 0.5f, -0.5f, 0.5f);
		glVertex3f( 0.5f,  0.5f, 0.5f);
		glVertex3f(-0.5f,  0.5f, 0.5f);
		glVertex3f(-0.5f, -0.5f, 0.5f);

		glNormal3f(0.0f, 0.0f, -1.0f);
		glVertex3f(0.5f, -0.5f, -0.5f);
		glVertex3f(0.5f, 0.5f, -0.5f);
		glVertex3f(-0.5f, 0.5f, -0.5f);
		glVertex3f(-0.5f, -0.5f, -0.5f);

		// right
		glNormal3f(1.0f,  0.0f,  0.0f);
		glVertex3f(0.5f,  0.5f, -0.5f);
		glVertex3f(0.5f,  0.5f,  0.5f);
		glVertex3f(0.5f, -0.5f,  0.5f);
		glVertex3f(0.5f, -0.5f, -0.5f);

		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-0.5f, 0.5f, -0.5f);
		glVertex3f(-0.5f, 0.5f, 0.5f);
		glVertex3f(-0.5f, -0.5f, 0.5f);
		glVertex3f(-0.5f, -0.5f, -0.5f);

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
void plot(History & hist)
{
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
	glfwMakeContextCurrent(window);

	//glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	int i = 0;

	gluPerspective(90, 1, 100, 10000);

	gluLookAt(0, 0, 400, 0, 0, 0, 0, 1, 0);

	do{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		

		glPushMatrix();
		{
			draw_frame(hist.frames[i]);	
		}
		glPopMatrix();

		

		glfwSwapBuffers(window);
		glfwPollEvents();

		printf("i=%i\n", i);

		i += 10;
		if (i >= hist.frames.size()) i = 0;

		

	} while (!glfwWindowShouldClose(window));

	
}

int _tmain(int argc, _TCHAR* argv[])
{
	int m = 10000;
	int n = 20;

	History hist;

	std::vector<Body> bodies(n);
	std::vector<Pair> pairs(n*(n-1)/2);

	//generate_binary_system(bodies, pairs);
	generate_bodies(bodies, pairs);
	generate_pairs(bodies, pairs);

	Header header;
	header.bodies_size = n;

	auto ocl = std::make_shared<OCL::Manager>();
	ocl->init();

	char * files[2] = { "kernel.cl", "kernel.h" };

	auto program = ocl->create_program(files, 1, "");

	auto kernel_calc_acc = program->create_kernel("calc_acc");
	auto kernel_step = program->create_kernel("step_pos");

	auto memobj_bodies = ocl->create_buffer(CL_MEM_READ_WRITE, bodies.size() * sizeof(Body));
	memobj_bodies->EnqueueWrite(&bodies[0], bodies.size() * sizeof(Body));

	auto memobj_pairs = ocl->create_buffer(CL_MEM_READ_WRITE, pairs.size() * sizeof(Pair));
	memobj_pairs->EnqueueWrite(&pairs[0], pairs.size() * sizeof(Pair));

	auto memobj_header = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(Header));
	memobj_header->EnqueueWrite(&header, sizeof(Header));

	unsigned int counter = 0;
	auto memobj_counter = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(unsigned int));
	memobj_counter->EnqueueWrite(&counter, sizeof(unsigned int));

	kernel_calc_acc->set_arg(memobj_header, 0);
	kernel_calc_acc->set_arg(memobj_bodies, 1);
	kernel_calc_acc->set_arg(memobj_pairs, 2);
	kernel_calc_acc->set_arg(memobj_counter, 3);

	kernel_step->set_arg(memobj_header, 0);
	kernel_step->set_arg(memobj_bodies, 1);
	kernel_step->set_arg(memobj_pairs, 2);
	kernel_step->set_arg(memobj_counter, 3);

	
	
	std::cout << "kernel start" << std::endl;


	for (int i = 0; i < m; ++i)
	{
		if (i % (m / 10) == 0) printf("%i\n", i);

		kernel_calc_acc->enqueue_ND_range_kernel(256, 1);
		kernel_step->enqueue_ND_range_kernel(256, 1);

		
		hist.frames.emplace_back();
		Frame & frame = hist.frames.back();
		frame.bodies.resize(n);
		memobj_bodies->EnqueueRead(&frame.bodies[0], n * sizeof(Body));
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

	plot(hist);

	getchar();

	return 0;
}



