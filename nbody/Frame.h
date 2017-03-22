

#include <string>
#include <fstream>
#include <vector>

#include "kernel.h"


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
		printf("bodies\n");
		for (int i = 0; i < bodies.size(); ++i)
		{
			printf("%4i\n", i);
			//::print(bodies[i].pos);
			//::print(bodies[i].vel);
			//::print(bodies[i].acc);
		}

		printf("pairs\n");
		for (int i = 0; i < pairs.size(); ++i)
		{
			printf("  d=%8f s=%8f\n", pairs[i].d, pairs[i].s);
		}
	}

	Header				header;
	std::vector<Body>	bodies;
	std::vector<Pair>	pairs;
};

