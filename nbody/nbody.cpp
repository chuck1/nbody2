// nbody.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Socket.h"

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
#include <signal.h>
#include <deque>

#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include "OCL.h"
#include "kernel.h"
#include "RenderGeometry.h"
#include "mymath.h"
#include "Frame.h"
#include "History.h"
#include "AppRender.h"
#include "AppSimulate.h"

int n_disc = 100;

//int n = 5;
//int p = n * (n - 1) / 2;

int frame_skip = 10;


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
void print(glm::quat const & m)
{
	printf("    %16.2e%16.2e%16.2e%16.2e\n", m[0], m[1], m[2], m[3]);
}
void print(Vec3 const & m)
{
	printf("    %16.2e%16.2e%16.2e\n", m.v[0], m.v[1], m.v[2]);
}
void print(Vec4 const & m)
{
	printf("    %16.2e%16.2e%16.2e%16.2e\n", m.v[0], m.v[1], m.v[2], m.v[3]);
}
void print(Quat const & m)
{
	printf("    %16.2e%16.2e%16.2e%16.2e\n", m.v[0], m.v[1], m.v[2], m.v[2]);
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

double body_surface_gravity(Body0 & b0, Body1 & b1)
{
	return M_G * b1.mass / pow(b0.radius, b0.radius);
}

void notification()
{

}

inline bool		file_exists(const std::string& name)
{
	std::ifstream f(name.c_str());
	return f.good();
}

std::deque<std::string> tokenize(const std::string& input)
{
	typedef boost::escaped_list_separator<char> separator_type;
	separator_type separator("\\",    // The escape characters.
		"= ",    // The separator characters.
		"\"\'"); // The quote characters.

	// Tokenize the intput.
	boost::tokenizer<separator_type> tokens(input, separator);

	// Copy non-empty tokens from the tokenizer into the result.
	std::deque<std::string> result;
	
	copy_if(tokens.begin(), tokens.end(), std::back_inserter(result), !boost::bind(&std::string::empty, _1));

	return result;
}

void print_usage(boost::program_options::options_description & desc)
{
	std::cout << "Usage: " << "<command> [options]" << std::endl;
	std::cout << desc << std::endl;
}

App * new_app_simulate(){ return new AppSimulate(); }
App * new_app_render(){ return new AppRender(); }

int _tmain(int argc, char* argv[])
{
	/*auto oclm = std::make_shared<OCL::Manager>();
	oclm->test();
	getchar();
	return 0;*/

	init_signal();

	while (true)
	{

		std::cout << "enter command" << std::endl;
		std::string s;
		std::getline(std::cin, s);

		auto tok = tokenize(s);

		if (tok.size() == 0) break;

		std::string command = tok.front();

		std::map<std::string, std::function<App*()>> app_map;
		app_map["s"] = &new_app_simulate;
		app_map["r"] = &new_app_render;

		auto it = app_map.find(command);

		if (it == app_map.end())
		{
			std::cout << "Invalid command" << std::endl;
			std::cout << "Commands" << std::endl;
			std::cout << "  s" << std::endl;
			std::cout << "  r" << std::endl;
			continue;
		}

		App * app = (it->second)();

		tok.pop_front();

		std::vector<std::string> vec(tok.begin(), tok.end());

		//app->run(vec);
		
		try
		{
			app->run(vec);
		}
		catch (OCL::Error & e)
		{

			std::cout << "OCL Error" << std::endl;
			break;
		}
	}

	printf("goodbye\n");
	getchar();
	return 0;
}


