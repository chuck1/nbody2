#ifndef APP_RENDER_H
#define APP_RENDER_H

#include "App.h"
#include "History.h"
#include "decl.h"

#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

void	callback_scroll(GLFWwindow* window, double xoffset, double yoffset);
void	callback_key(GLFWwindow* window, int key, int scancode, int action, int mods);
void	callback_cursor_position(GLFWwindow* window, double xpos, double ypos);
void	callback_mouse_button(GLFWwindow* window, int button, int action, int mods);


void draw_cube();



void	init_glew();
GLuint	LoadShaders(const char * vertex_file_path, const char * fragment_file_path);
void	vec_find(std::vector<double> v, double & x, int & i);

class AppRender : public App
{
public:
	AppRender();
	
	virtual void	run(std::vector<std::string> vec);

	void			render();
	void			plotfunc();
	
	void			draw_frame(std::shared_ptr<Frame> & frame);
	void			draw_body(Body0 & b0);

	History			hist;

	float			eye_z;
	float			perspective_near;
	float			perspective_far;

	static AppRender *	s_app_render;

	void			callback_scroll(GLFWwindow* window, double xoffset, double yoffset);
	void			callback_key(GLFWwindow* window, int key, int scancode, int action, int mods);
	void			callback_cursor_position(GLFWwindow* window, double xpos, double ypos);
	void			callback_mouse_button(GLFWwindow* window, int button, int action, int mods);


	double	t_factor;
	double	t;
	int		i;
	bool	pause;


	int			window_width;
	int			window_height;

	Sphere *	sphere;
	GLuint		matrixID_mv;
	GLuint		matrixID_p;
	glm::mat4	proj, view;
	glm::quat	q_view;



	double cursor_x = 0;
	double cursor_y = 0;
	bool mouse_left_pressed = false;
};

#endif



