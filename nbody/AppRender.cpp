#include "stdafx.h"

#include "mymath.h"
#include "AppRender.h"
#include "RenderGeometry.h"

#include <boost/program_options.hpp>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

void	callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
	AppRender::s_app_render->callback_scroll(window, xoffset, yoffset);
}
void	callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	AppRender::s_app_render->callback_key(window, key, scancode, action, mods);
}
void	callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
	AppRender::s_app_render->callback_cursor_position(window, xpos, ypos);
}
void	callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
	AppRender::s_app_render->callback_mouse_button(window, button, action, mods);
}

AppRender *	AppRender::s_app_render;

AppRender::AppRender()
{
	s_app_render = this;

	eye_z = 800.0;
	perspective_near = 1.0;
	perspective_far = 10000.0;

	window_width = 200;
	window_height = 200;

	q_view = glm::quat(1, 0, 0, 0);

	cursor_x = 0;
	cursor_y = 0;
	mouse_left_pressed = false;
}
void	AppRender::run(std::vector<std::string> vec)
{
	boost::program_options::options_description desc("render options");

	desc.add_options()
		("help", "produce help message");

	boost::program_options::variables_map vm;

	try{
		boost::program_options::store(boost::program_options::command_line_parser(vec).options(desc).run(), vm);
	}
	catch (std::exception & e)
	{
		std::cout << e.what() << std::endl;
		std::cout << desc << std::endl;
		return;
	}

	boost::program_options::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return;
	}

	render();
}
void	AppRender::render()
{
	hist.folder = "C:\\test";
	hist.load();

	plotfunc();
}
void	AppRender::plotfunc()
{
	auto f = hist.get_frame(0);
	//eye_z = 1.5 * f->x_max();

	t_factor = 1.0E1;
	t = 0;
	i = 0;
	pause = false;

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

	glfwSetKeyCallback(window, ::callback_key);
	glfwSetScrollCallback(window, ::callback_scroll);
	glfwSetCursorPosCallback(window, ::callback_cursor_position);
	glfwSetMouseButtonCallback(window, ::callback_mouse_button);

	glfwMakeContextCurrent(window);

	init_glew();

	sphere = new Sphere();
	sphere->construct();
	sphere->setup();

	proj = glm::perspective<float>(
		glm::radians(90.0f),
		(float)window_width / (float)window_height,
		perspective_near,
		perspective_far);



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
		view = glm::lookAt(
			glm::vec3(0, 0, eye_z), // Camera is at (4,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
			);

		glm::mat4 view_rot = glm::mat4_cast(q_view);

		view = view * view_rot;

		double t1 = glfwGetTime();
		double dt = t1 - t0;
		t0 = t1;

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glUseProgram(programID);

		glUniformMatrix4fv(matrixID_p, 1, GL_FALSE, &proj[0][0]);
		gl_check_error();

		f = hist.get_frame(i);


		draw_frame(f);


		glfwSwapBuffers(window);
		glfwPollEvents();

		if (!pause)
		{
			t += dt * t_factor;

			if (1) // constant time mode
			{

				vec_find(hist.frame_times, t, i);

				/*while (true)
				{
				if (plot.t_factor > 0)
				{
				if (hist.frame_times[plot.i] > plot.t) break;
				plot.i += 1;
				}
				else if (plot.t_factor < 0)
				{
				if (hist.frame_times[plot.i] < plot.t) break;
				plot.i -= 1;
				}



				if (plot.i >= hist.frame_times.size()) {
				plot.i = 0;
				plot.t = 0;
				}
				}*/
			}
			else //  constant frame rate mode
			{
				i += 1;

				if (i >= (int)hist.frame_times.size()) {
					i = 0;
				}
			}

			f = hist.get_frame(i);
		}

		//printf("i=%6i t=%16f\n", plot.i, f->header.t);
		//printf("t_sim=%8f i=%6i t=%16f dt=%16f tf=%8.2e\n", plot.t, plot.i, f->header.t, f->header.dt, plot.t_factor);

	} while (!glfwWindowShouldClose(window));

	glfwDestroyWindow(window);
}

void	AppRender::callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
	//printf("xoffset=%f yoffset=%f\n", xoffset, yoffset);

	eye_z *= (float)pow(1.1f, -yoffset);
}
void	AppRender::callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_Q && action == GLFW_PRESS)
	{
		t_factor *= 2.0;
		printf("time factor %f\n", AppRender::s_app_render->t_factor);
	}

	if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{
		t_factor /= 2.0;
		printf("time factor %f\n", AppRender::s_app_render->t_factor);
	}

	if (key == GLFW_KEY_R && action == GLFW_PRESS)
	{
		t_factor *= -1;
		printf("time factor %f\n", AppRender::s_app_render->t_factor);
	}

	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
	{
		t = 0;
		i = 0;
	}

	if (key == GLFW_KEY_P && action == GLFW_PRESS)
		pause = true;

	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
}
void	AppRender::callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
	float dx = (float)(xpos - cursor_x);
	float dy = (float)(ypos - cursor_y);
	cursor_x = xpos;
	cursor_y = ypos;

	if (mouse_left_pressed)
	{
		printf("cursor %16.2f %16.2f\n", dx, dy);

		glm::quat q = glm::angleAxis<float>(sqrt(dy*dy + dx*dx) / 100.0f * (float)M_PI, glm::normalize(glm::vec3(dy, dx, 0.0f)));

		q_view = q * q_view;
	}
}
void	AppRender::callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
	printf("mouse button %i %i %i\n", button, action, mods);

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		mouse_left_pressed = true;
		glfwGetCursorPos(window, &cursor_x, &cursor_y);
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
	{
		mouse_left_pressed = false;
	}
}

void	AppRender::draw_body(Body0 & b0)
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

	glm::mat4 t = glm::translate(glm::vec3(b0.pos.v[0], b0.pos.v[1], b0.pos.v[2]));

	glm::quat q(b0.q.v[0], b0.q.v[1], b0.q.v[2], b0.q.v[3]);

	glm::mat4 r = glm::mat4_cast(q);

	glm::mat4 s = glm::scale(glm::vec3(b0.radius));

	glm::mat4 model = t * r * s;

	glm::mat4 mv = view * model;

	glUniformMatrix4fv(matrixID_mv, 1, GL_FALSE, &mv[0][0]);
	gl_check_error();

	sphere->drawTriangles();
	gl_check_error();
}
void	AppRender::draw_frame(std::shared_ptr<Frame> & frame)
{
	for (unsigned int i = 0; i < frame->bodies0.size(); ++i)
	{
		draw_body(frame->bodies0[i]);
	}

#if 0
	// debug
	for (unsigned int i = 0; i < frame->pairs.size(); ++i)
	{
		printf("pair F = %16.2e\n", frame->pairs[i].F);
	}
#endif
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



void	init_glew()
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
GLuint	LoadShaders(const char * vertex_file_path, const char * fragment_file_path){

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
void	vec_find(std::vector<double> v, double & x, int & i)
{
	int d;
	int l = v.size();
	double s;
	if (x > v[i])
	{
		d = 1;
		s = 1;
	}
	else
	{
		d = -1;
		s = -1;
	}

	while (((x - v[i])*s) > 0)
	{
		i += d;

		if (i < 0) {
			i = l - 1;
			x = v[i];
			break;
		}

		if (i >(l - 1)) {
			i = 0;
			x = v[i];
			break;
		}
	}
}


