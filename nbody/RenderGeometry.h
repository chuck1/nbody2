
#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLFW/glfw3.h>

void gl_check_error();


class Geometry
{
public:
	void drawTriangles()
	{
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			4,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
			);
		gl_check_error();
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		if (1)
		{
			glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
			glVertexAttribPointer(
				1,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
				3,                  // size
				GL_FLOAT,           // type
				GL_FALSE,           // normalized?
				0,                  // stride
				(void*)0            // array buffer offset
				);
			gl_check_error();
			glBindBuffer(GL_ARRAY_BUFFER, 0);
		}

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

		glDrawElements(
			GL_TRIANGLES,      // mode
			indexArrayIndex,    // count
			GL_UNSIGNED_INT,   // type
			(void*)0           // element array buffer offset
			);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
	}
	void drawLines()
	{
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			4,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
			);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

		glDrawElements(
			GL_LINES,      // mode
			indexArrayIndex,    // count
			GL_UNSIGNED_INT,   // type
			(void*)0           // element array buffer offset
			);

		glDisableVertexAttribArray(0);
	}
	void setup()
	{
		glGenVertexArrays(1, &VertexArrayID);
		glBindVertexArray(VertexArrayID);

		// Generate 1 buffer, put the resulting identifier in vertexbuffer
		glGenBuffers(1, &vertexBuffer);
		glGenBuffers(1, &normalBuffer);
		glGenBuffers(1, &indexBuffer);

		// The following commands will talk about our 'vertexbuffer' buffer
		glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, verts * 4 * sizeof(GLfloat), p, GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
		glBufferData(GL_ARRAY_BUFFER, verts * 3 * sizeof(GLfloat), normal, GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexArrayIndex * sizeof(GLuint), indexArray, GL_STATIC_DRAW);
	}

	GLuint vertexBuffer;
	GLuint normalBuffer;
	GLuint indexBuffer;

	GLuint VertexArrayID;

	GLfloat * p;
	GLfloat * normal;
	int verts;

	GLuint * indexArray;
	int indexArrayIndex;
};

class Sphere : public Geometry
{
public:
	void construct()
	{
		double pi = 3.1415926535897;

		int m = 16;
		int n = 16;

		verts = (m - 1) * n;
		
		p = new GLfloat[verts * 4];
		normal = new GLfloat[verts * 3];

		for (int i = 0; i < (m - 1); ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				double a = (float(i + 1) / (float)m - 0.5) * pi;
				double b = (float)j / (float)n * 2.0 * pi;

				double x = cos(b) * cos(a);
				double y = sin(b) * cos(a);
				double z = sin(a);

				GLfloat * p0 = p + (i * n + j) * 4;
				GLfloat * n0 = normal + (i * n + j) * 3;

				p0[0] = x;
				p0[1] = y;
				p0[2] = z;
				p0[3] = 1;

				n0[0] = x;
				n0[1] = y;
				n0[2] = z;

				//printf("%8.2f%8.2f%8.2f\n", x, y, z);

			}
		}

		indexArrayIndex = (m - 2) * n * 6;

		indexArray = new GLuint[indexArrayIndex];

		int k = 0;
		for (int i = 0; i < (m - 2); ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				int j0 = j;
				int j1 = (j + 1) % n;

				indexArray[k++] = ((i + 0) * n + j0);
				indexArray[k++] = ((i + 0) * n + j1);
				indexArray[k++] = ((i + 1) * n + j0);

				indexArray[k++] = ((i + 0) * n + j1);
				indexArray[k++] = ((i + 1) * n + j0);
				indexArray[k++] = ((i + 1) * n + j1);
			}
		}

		
	}
};

class Rect : public Geometry
{
public:
	void construct()
	{
		static const GLfloat g_vertex_buffer_data[] = {
			-1.0f, -1.0f, 0.0f, 1.0f,
			1.0f, -1.0f, 0.0f, 1.0f,
			1.0f, 1.0f, 0.0f, 1.0f,
			-1.0f, 1.0f, 0.0f, 1.0f,
		};

		p = new GLfloat[16];

		memcpy(p, g_vertex_buffer_data, 16 * sizeof(GLfloat));

		verts = 4;

		indexArray = new GLuint[6];
		indexArray[0] = 0;
		indexArray[1] = 1;
		indexArray[2] = 3;
		indexArray[3] = 1;
		indexArray[4] = 2;
		indexArray[5] = 3;
		indexArrayIndex = 6;
	}
};
