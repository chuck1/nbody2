// nbody.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <algorithm>

#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include "OCL.h"
#include "kernel.h"
#include "RenderGeometry.h"

double G = 6.67408E-11;
double pi = 3.1415926535897;

double eye_z = 800.0;
double perspective_near = 1.0;
double perspective_far = 10000.0;

int m = 10000;
int n = 20;
int p = n * (n - 1) / 2;

Sphere sphere;
GLuint matrixID_mv;
GLuint matrixID_p;
glm::mat4 proj, view;

void gl_check_error()
{
	GLenum err(glGetError());
	if (err != GL_NO_ERROR) {
		printf("GL error %u %s\n", err, gluErrorString(err));
		getchar();
		exit(1);
	}
}


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

	b0.radius = pow(3.0 * b0.mass / 4.0 / pi / 5000.0, 1.0 / 3.0);
	b1.radius = pow(3.0 * b1.mass / 4.0 / pi / 5000.0, 1.0 / 3.0);

	printf("b0.radius = %f\n", b0.radius);
	printf("b1.radius = %f\n", b1.radius);

	double d = 1000.0;

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


void draw_cube()
{
	glBegin(GL_QUADS);

	// top
	glNormal3f(0.0f, 1.0f, 0.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);

	glNormal3f(0.0f, -1.0f, 0.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);

	// front
	glNormal3f(0.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);

	glNormal3f(0.0f, 0.0f, -1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);

	// right
	glNormal3f(1.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);

	glNormal3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);

	glEnd();
}
void draw_body(Body & b)
{
	/*glPushMatrix();
	{
		glTranslated(b.pos.x, b.pos.y, b.pos.z);

		double s = b.radius;
		glScaled(s,s,s);

		glColor3f(1.0f, 1.0f, 1.0f);

		draw_cube();
	}
	glPopMatrix();*/


	glm::mat4 model = glm::translate(glm::vec3(b.pos.x, b.pos.y, b.pos.z)) * glm::scale(glm::vec3(b.radius));

	glm::mat4 mv = view * model;

	glUniformMatrix4fv(matrixID_mv, 1, GL_FALSE, &mv[0][0]);
	gl_check_error();

	sphere.drawTriangles(); gl_check_error();


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

GLuint LoadShaders(const char * vertex_file_path, const char * fragment_file_path){

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if (VertexShaderStream.is_open()){
		std::string Line = "";
		while (getline(VertexShaderStream, Line))
			VertexShaderCode += "\n" + Line;
		VertexShaderStream.close();
	}
	else{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if (FragmentShaderStream.is_open()){
		std::string Line = "";
		while (getline(FragmentShaderStream, Line))
			FragmentShaderCode += "\n" + Line;
		FragmentShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;


	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer, NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}



	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer, NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}



	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0){
		std::vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}


	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	return ProgramID;
}

void print(glm::mat4 const & m)
{
	for (int i = 0; i < 4; ++i)
	{
		printf("    %8.2f%8.2f%8.2f%8.2f\n", m[i][0], m[i][1], m[i][2], m[i][3]);
	}
}
void print(glm::vec4 const & m)
{
	printf("    %8.2f%8.2f%8.2f%8.2f\n", m[0], m[1], m[2], m[3]);
}


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
void init_glew()
{
	//Initialize GLEW
	GLenum err = glewInit();

	//If GLEW hasn't initialized
	if (err != GLEW_OK)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		getchar();
		exit(-1);
	}
}
void plotfunc(History & hist)
{
	plot.t_factor = 1.0E1;
	plot.t = 0;
	plot.i = 0;
	plot.pause = false;

	if (!glfwInit())
	{
		exit(EXIT_FAILURE);
	}

	int width = 800;
	int height = 800;

	//Create a window and create its OpenGL context
	GLFWwindow * window = glfwCreateWindow(width, height, "Test Window", NULL, NULL);
	if (!window)
	{
		fprintf(stderr, "Failed to open GLFW window.\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwSetKeyCallback(window, key_callback);

	glfwMakeContextCurrent(window);

	init_glew();

	
	sphere.construct();
	sphere.setup();

	Rect rect;
	//rect.construct();
	//rect.setup();

	proj = glm::perspective<float>(
		glm::radians(90.0f),
		(float)width / (float)height,
		perspective_near,
		perspective_far);

	view = glm::lookAt(
		glm::vec3(0, 0, eye_z), // Camera is at (4,3,3), in World Space
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
		);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("SimpleVertexShader.vertexshader", "SimpleFragmentShader.fragmentshader");

	matrixID_mv = glGetUniformLocation(programID, "MV");
	matrixID_p = glGetUniformLocation(programID, "P");
	gl_check_error();

	printf("matrixID_mv = %u\n", matrixID_mv);
	printf("matrixID_p  = %u\n", matrixID_p);

	
	/*
	print(model); printf("\n");
	print(view); printf("\n");
	print(proj); printf("\n");
	print(mvp); printf("\n");*/

	glm::vec4 v(1, 1, 0, 1);

	//v = mvp * v;
	//print(v); printf("\n");

	double t0 = glfwGetTime();

	do
	{
		double t1 = glfwGetTime();
		double dt = t1 - t0;
		t0 = t1;

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glUseProgram(programID);

		glUniformMatrix4fv(matrixID_p, 1, GL_FALSE, &proj[0][0]);
		gl_check_error();

		Frame & f = hist.frames[plot.i];


		draw_frame(f);
		

		//sphere.drawTriangles();

		//glPushMatrix();
		{
			//glRotatef(t1 * pi, 1,0,0);

			//sphere.drawTriangles(); gl_check_error();
			
			//rect.drawTriangles();
			
			
			
			//draw_cube();
		}
		//glPopMatrix();



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
		kernel_dt_min->enqueue_ND_range_kernel(next_power_of_two(p), std::min<unsigned int>(next_power_of_two(p), 1024));

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



