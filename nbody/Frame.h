#ifndef FRAME_H
#define FRAME_H

#include <algorithm>
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

		bodies0.resize(s);

		for (unsigned int i = 0; i < s; ++i)
		{
			myFile.read((char*)&bodies0[i], sizeof(Body0));
		}

		// pairs
#if 0
		myFile.read((char*)&s, sizeof(unsigned int));

		pairs.resize(s);

		for (unsigned int i = 0; i < s; ++i)
		{
			myFile.read((char*)&pairs[i], sizeof(Pair));
		}
#endif
	}
	void				save(std::ofstream & myFile)
	{
		myFile.write((char*)&header, sizeof(Header));

		unsigned int s;

		// bodies

		s = bodies0.size();

		myFile.write((char*)&s, sizeof(unsigned int));

		for (unsigned int i = 0; i < s; ++i)
		{
			myFile.write((char*)&bodies0[i], sizeof(Body0));
		}

		// pairs
#if 0
		s = pairs.size();

		myFile.write((char*)&s, sizeof(unsigned int));

		for (unsigned int i = 0; i < s; ++i)
		{
			myFile.write((char*)&pairs[i], sizeof(Pair));
		}
#endif
	}

	double				x_max()
	{
		double x = 0;

		for (unsigned int i = 0; i < bodies0.size(); ++i)
		{
			x = std::max(x, bodies0[i].pos.v[0]);
		}

		return x;
	}
	void				print()
	{
		printf("bodies\n");
		for (unsigned int i = 0; i < bodies0.size(); ++i)
		{
			printf("%4i\n", i);
			//::print(bodies[i].pos);
			//::print(bodies[i].vel);
			//::print(bodies[i].acc);
		}

#if 0
		printf("pairs\n");
		for (unsigned int i = 0; i < pairs.size(); ++i)
		{
			printf("  d=%8f s=%8f\n", pairs[i].d, pairs[i].s);
		}
#endif
	}

	Header				header;
	std::vector<Body0>	bodies0;
	//std::vector<Pair>	pairs;
};


#endif