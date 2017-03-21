// nbody.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <regex>

#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

#include "OCL.h"
#include "kernel.h"
#include "RenderGeometry.h"
#include "mymath.h"

double G = 6.67408E-11;
double pi = 3.1415926535897;

double eye_z = 800.0;
double perspective_near = 1.0;
double perspective_far = 10000.0;

int m = 0;
int n = 5;
int p = n * (n - 1) / 2;

int window_width = 200;
int window_height = 200;

Sphere sphere;
GLuint matrixID_mv;
GLuint matrixID_p;
glm::mat4 proj, view;

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
void print(Vec3 const & m)
{
	printf("    %8.2f%8.2f%8.2f\n", m.v[0], m.v[1], m.v[2]);
}

void gl_check_error()
{
	GLenum err(glGetError());
	if (err != GL_NO_ERROR) {
		printf("GL error %u %s\n", err, gluErrorString(err));
		getchar();
		exit(1);
	}
}

double gravity_uniform_disc(double density, double p, double r)
{
	double k = 4.0 * r * p / (pow(r, 2.0) + pow(p, 2.0) + 2.0 * r * p);
	double K = mymath::K(k);
	double E = mymath::E(k);
	double a = (1.0 - 0.5 * pow(k, 2.0)) * K - E;
	return 4.0 * G * density * pow(p, 0.5) / k / pow(r, 0.5) * a;
}

class Frame
{
public:
	void				load(std::string folder, unsigned int i)
	{
		char buffer[128];
		sprintf_s(buffer, "frame_%08u.bin", i);
		std::ifstream myFile(folder + "\\frames\\" + buffer, std::ios::in | std::ios::binary);

		load(myFile);

		myFile.close();
	}
	void				save(std::string folder, unsigned int i)
	{
		char buffer[128];
		sprintf_s(buffer, "frame_%08u.bin", i);

		std::string filename = folder + "\\frames\\" + buffer;

		std::ofstream myFile(filename, std::ios::out | std::ios::binary);

		save(myFile);

		myFile.close();
	}
	void				load(std::ifstream & myFile)
	{
		myFile.read((char*)&header, sizeof(Header));

		unsigned int s;

		// bodies

		myFile.read((char*)&s, sizeof(unsigned int));

		bodies.resize(s);

		for (int i = 0; i < s; ++i)
		{
			myFile.read((char*)&bodies[i], sizeof(Body));
		}

		// pairs

		myFile.read((char*)&s, sizeof(unsigned int));

		pairs.resize(s);

		for (int i = 0; i < s; ++i)
		{
			myFile.read((char*)&pairs[i], sizeof(Pair));
		}
	}
	void				save(std::ofstream & myFile)
	{
		myFile.write((char*)&header, sizeof(Header));

		unsigned int s;

		// bodies

		s = bodies.size();

		myFile.write((char*)&s, sizeof(unsigned int));

		for (int i = 0; i < s; ++i)
		{
			myFile.write((char*)&bodies[i], sizeof(Body));
		}

		// pairs

		s = pairs.size();

		myFile.write((char*)&s, sizeof(unsigned int));

		for (int i = 0; i < s; ++i)
		{
			myFile.write((char*)&pairs[i], sizeof(Pair));
		}
	}

	double				x_max()
	{
		double x = 0;

		for (int i = 0; i < bodies.size(); ++i)
		{
			x = std::max(x, bodies[i].pos.v[0]);
		}

		return x;
	}
	void				print()
	{
		for (int i = 0; i < n; ++i)
		{
			printf("%4i\n", i);
			::print(bodies[i].pos);
			::print(bodies[i].vel);
			::print(bodies[i].acc);
		}
	}

	Header				header;
	std::vector<Body>	bodies;
	std::vector<Pair>	pairs;
};

class History
{
public:
	History(std::string f) : folder(f)
	{}

	void write()
	{
		std::ofstream myFile;

		myFile.open(folder + "\\" + "frame_times.bin", std::ios::out | std::ios::binary);

		write(myFile);

		myFile.close();


		/*for (int i = 0; i < frames.size(); ++i)
		{
			frames[i].save(folder, i);
		}*/
	}
	void write(std::ofstream & myFile)
	{
		unsigned int l = frame_times.size();

		myFile.write((char *)&l, sizeof(unsigned int));

		myFile.write((char *)&frame_times[0], sizeof(double) * l);
	}

	void push(Header & header, std::shared_ptr<OCL::MemObj> memobj_bodies, std::shared_ptr<OCL::MemObj> memobj_pairs)
	{
		//frames.emplace_back();
		//Frame & frame = frames.back();
		Frame frame;

		frame.header = header;
		frame.bodies.resize(n);
		frame.pairs.resize(p);
		memobj_bodies->EnqueueRead(&frame.bodies[0], n * sizeof(Body));
		memobj_pairs->EnqueueRead(&frame.pairs[0], p * sizeof(Pair));

		frame_times.push_back(header.t);

		unsigned int i = frame_times.size() - 1;

		frame.save(folder, i);
	}


	void						load()
	{
		std::ifstream myFile(folder + "\\" + "frame_times.bin", std::ios::in | std::ios::binary);

		unsigned int s;

		myFile.read((char*)&s, sizeof(unsigned int));

		frame_times.resize(s);

		myFile.read((char *)&frame_times[0], sizeof(double) * s);

		myFile.close();
	}
	std::shared_ptr<Frame>		get_frame(unsigned int i)
	{
		auto it = frames.find(i);
		if (it != frames.end()) return it->second;

		return load_frame(i);
	}
	std::shared_ptr<Frame>		load_frame(unsigned int i)
	{
		auto f = std::make_shared<Frame>();

		f->load(folder, i);

		return f;
	}

private:
	std::map<unsigned int, std::shared_ptr<Frame>>	frames;

public:
	std::vector<double>		frame_times;

	std::string				folder;
};

void generate_binary_system(
	std::vector<Body> & bodies,
	std::vector<Pair> & pairs)
{
	int n = bodies.size();

	Body & b0 = bodies[0];
	Body & b1 = bodies[1];

	double & x0 = b0.pos.v[0];
	double & x1 = b1.pos.v[0];

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

	b0.vel.v[1] = f * sqrt(G * m1 * x0 / d / d);
	b1.vel.v[1] = -f * sqrt(G * m1 * -x1 / d / d);
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

void generate_disc(
	std::vector<Body> & bodies,
	std::vector<Pair> & pairs)
{
	int n = bodies.size();
	
	double x, y, r;

	double body_mass = 1E8;
	double body_radius = pow(3.0 * body_mass / 4.0 / pi / 5000.0, 1.0/3.0);

	double radius = sqrt(10.0 * (double)n) * body_radius;

	double density = body_mass * (double)n / (pi * pow(radius, 2.0));

	for (int i = 0; i < n; ++i)
	{
		while (true)
		{
			x = (float)rand() / (float)RAND_MAX * 2.0 - 1.0;
			y = (float)rand() / (float)RAND_MAX * 2.0 - 1.0;
			x = x * radius;
			y = y * radius;
			r = sqrt(x*x + y*y);

			if (r < radius) break;
		}

		Body & b = bodies[i];

		b.mass = body_mass;
		b.radius = body_radius;

		b.pos.v[0] = x;
		b.pos.v[1] = y;
	}

	double X = 0;
	double Y = 0;
	double Z = 0;
	double M = 0;
	for (int i = 0; i < n; ++i)
	{
		Body & b = bodies[i];
		X += b.mass * b.pos.v[0];
		Y += b.mass * b.pos.v[1];
		Z += b.mass * b.pos.v[2];
		M += b.mass;
	}
	X /= M;
	Y /= M;
	Z /= M;

	for (int i = 0; i < n; ++i)
	{
		Body & b = bodies[i];

		double x = b.pos.v[0] - X;
		double y = b.pos.v[1] - Y;

		double r = sqrt(x*x + y*y);

		if (r == 0)
		{
			b.vel.v[1] = 0;
			b.vel.v[0] = 0;
		}
		else
		{
			double k = 0.5;

			double acc1 = G * M / r / r;
			double acc2 = k * gravity_uniform_disc(density, radius, r);

			printf("grav acc %16e %16e\n", acc1, acc2);

			double s = sqrt(acc2 * r);

			b.vel.v[1] = s * x / r;
			b.vel.v[0] = -s * y / r;
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
		glTranslated(b.pos.v[0], b.pos.v[1], b.pos.v[2]);

		double s = b.radius;
		glScaled(s,s,s);

		glColor3f(1.0f, 1.0f, 1.0f);

		draw_cube();
	}
	glPopMatrix();*/

	glm::mat4 t = glm::translate(glm::vec3(b.pos.v[0], b.pos.v[1], b.pos.v[2]));
	
	glm::quat q(b.q0.v[0], b.q0.v[1], b.q0.v[2], b.q0.v[3]);

	glm::mat4 r = glm::mat4_cast(q);

	glm::mat4 s = glm::scale(glm::vec3(b.radius));

	glm::mat4 model = t * r * s;

	glm::mat4 mv = view * model;

	glUniformMatrix4fv(matrixID_mv, 1, GL_FALSE, &mv[0][0]);
	gl_check_error();

	sphere.drawTriangles(); gl_check_error();
}
void draw_frame(std::shared_ptr<Frame> & frame)
{
	for (int i = 0; i < frame->bodies.size(); ++i)
	{
		draw_body(frame->bodies[i]);
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
		plot.pause = true;

	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
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
	auto f = hist.get_frame(0);
	eye_z = 1.5 * f->x_max();
	
	plot.t_factor = 1.0E1;
	plot.t = 0;
	plot.i = 0;
	plot.pause = false;

	if (!glfwInit())
	{
		exit(EXIT_FAILURE);
	}

	

	//Create a window and create its OpenGL context
	GLFWwindow * window = glfwCreateWindow(window_width, window_height, "Test Window", NULL, NULL);
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
		(float)window_width / (float)window_height,
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

	glEnable(GL_DEPTH_TEST);

	do
	{
		double t1 = glfwGetTime();
		double dt = t1 - t0;
		t0 = t1;

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glUseProgram(programID);

		glUniformMatrix4fv(matrixID_p, 1, GL_FALSE, &proj[0][0]);
		gl_check_error();

		f = hist.get_frame(plot.i);


		draw_frame(f);
		

		glfwSwapBuffers(window);
		glfwPollEvents();

		if (!plot.pause)
		{
			plot.t += dt * plot.t_factor;

			while (hist.frame_times[plot.i] < plot.t)
			{
				plot.i += 1;

				if (plot.i >= hist.frame_times.size()) {
					plot.i = 0;
					plot.t = 0;
				}
			}

			f = hist.get_frame(plot.i);
		}


		printf("t_sim=%8f i=%6i t=%16f dt=%16f tf=%8.2e\n", plot.t, plot.i, f->header.t, f->header.dt, plot.t_factor);

	} while (!glfwWindowShouldClose(window));

	glfwDestroyWindow(window);
}

inline bool file_exists(const std::string& name)
{
	std::ifstream f(name.c_str());
	return f.good();
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

class OCL_manager
{
public:
	void	init(std::vector<Body> & bodies, std::vector<Pair> & pairs, Header & header)
	{
		ocl = std::make_shared<OCL::Manager>();
		ocl->init();

		char * files[2] = { "kernel.cl", "kernel.h" };

		auto program = ocl->create_program(files, 1, "");

		kernel_calc_acc = program->create_kernel("calc_acc");
		kernel_step = program->create_kernel("step_pos");
		kernel_dt_min = program->create_kernel("k_min");
		kernel_dt_store = program->create_kernel("store_dt");

		memobj_bodies = ocl->create_buffer(CL_MEM_READ_WRITE, bodies.size() * sizeof(Body));
		memobj_bodies->EnqueueWrite(&bodies[0], bodies.size() * sizeof(Body));

		memobj_pairs = ocl->create_buffer(CL_MEM_READ_WRITE, pairs.size() * sizeof(Pair));
		memobj_pairs->EnqueueWrite(&pairs[0], pairs.size() * sizeof(Pair));

		memobj_header = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(Header));
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

		kernel_dt_min->set_arg(memobj_dt_input, 0);
		kernel_dt_min->set_arg(memobj_dt_len, 1);
		kernel_dt_min->set_arg(memobj_dt_partial, 2);
		kernel_dt_min->set_arg(3, 1024 * sizeof(double));

		kernel_dt_store->set_arg(memobj_dt_partial, 0);
		kernel_dt_store->set_arg(memobj_header, 1);
	}
	void	shutdown()
	{
		ocl->shutdown();
	}

	std::shared_ptr<OCL::Manager>	ocl;

	std::shared_ptr<OCL::MemObj>	memobj_bodies;
	std::shared_ptr<OCL::MemObj>	memobj_pairs;
	std::shared_ptr<OCL::MemObj>	memobj_header;

	std::shared_ptr<OCL::Kernel>	kernel_calc_acc;
	std::shared_ptr<OCL::Kernel>	kernel_step;
	std::shared_ptr<OCL::Kernel>	kernel_dt_min;
	std::shared_ptr<OCL::Kernel>	kernel_dt_store;
};

int simulate()
{
	unsigned int file_size = m * (n * sizeof(Body)+p * sizeof(Pair));

	printf("file size = %u MB\n", file_size / 1024 / 1024);

	srand(time(NULL));

	double * dt_partial = new double[p];

	History hist("C:\\test");
	
	std::vector<Body> bodies;
	std::vector<Pair> pairs;

	Header header;
	header.bodies_size = n;
	header.t = 0;
	header.dt = 0;

	if (file_exists(hist.folder + "\\frame_times.bin"))
	{
		hist.load();

		unsigned int i = hist.frame_times.size() - 1;

		auto f = hist.get_frame(i);

		header = f->header;
		bodies = f->bodies;
		pairs = f->pairs;
	}
	else
	{
		bodies.resize(n);
		pairs.resize(p);

		//generate_binary_system(bodies, pairs);
		generate_disc(bodies, pairs);

		generate_pairs(bodies, pairs);
	}

	
	OCL_manager om;
	om.init(bodies, pairs, header);
	

	std::cout << "kernel start" << std::endl;

	// if this is a new simulation, save the zero frame
	if (hist.frame_times.empty())
	{
		hist.push(header, om.memobj_bodies, om.memobj_pairs);
	}
	else
	{
		header.t = hist.frame_times.back();
	}


	for (int i = 0; i < m; ++i)
	{
		double t0 = header.t;

		om.kernel_calc_acc->enqueue_ND_range_kernel(next_power_of_two(p), 1);
		om.kernel_dt_min->enqueue_ND_range_kernel(next_power_of_two(p), std::min<unsigned int>(next_power_of_two(p), 1024));

		om.kernel_dt_store->enqueue_ND_range_kernel(1, 1);

		om.kernel_step->enqueue_ND_range_kernel(next_power_of_two(n), 1);


		om.memobj_header->EnqueueRead(&header, sizeof(Header));

		header.t = t0 + header.dt;

		if (i % (m / 10) == 0) printf("%8i t=%16f dt=%16f count_pen=%8i\n", i, header.t, header.dt, header.count_pen);



		// save
		hist.push(header, om.memobj_bodies, om.memobj_pairs);
	}

	

	std::cout << "read" << std::endl;

	om.memobj_bodies->EnqueueRead(&bodies[0], n * sizeof(Body));

	

	hist.write();

	std::cout << "done" << std::endl;

	om.shutdown();

	getchar();

	return 0;
}

int render()
{
	History hist("C:\\test");
	hist.load();

	plotfunc(hist);

	return 0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	while (true){

		std::cout << "enter command" << std::endl;
		std::string s;
		std::getline(std::cin, s);

		std::smatch sm1;

		std::regex e1("s (\\d+)");

		if (std::regex_match(s, sm1, e1))
		{
			m = stoi(sm1[1]);

			simulate();
		}
		else if (s.compare("r") == 0)
		{
			render();
		}
		else{
			break;
		}
	}
	return 0;
}





