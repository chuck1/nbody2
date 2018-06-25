#ifndef HISTORY_H
#define HISTORY_H

#include <string>
#include <fstream>
#include <map>

#include "OCL.h"
#include "Misc.h"
#include "Frame.h"

class History
{
public:
	History()
	{}

	void						write()
	{
		std::ofstream myFile;

		myFile.open(folder + "\\" + "frame_times.bin", std::ios::out | std::ios::binary);

		write(myFile);

		myFile.close();
	}
	void						write(std::ofstream & myFile)
	{
		// frame times

		unsigned int l = frame_times.size();

		myFile.write((char *)&l, sizeof(unsigned int));

		myFile.write((char *)&frame_times[0], sizeof(double)* l);

		// bodies1

		unsigned int s = bodies1.size();

		myFile.write((char*)&s, sizeof(unsigned int));

		myFile.write((char*)&bodies1[0], s * sizeof(Body1));
	}

	void						resize(int n)
	{
		frame.bodies0.resize(n);
	}
	void						push(
		std::shared_ptr<OCL::MemObj> memobj_header, 
		std::shared_ptr<OCL::MemObj> memobj_bodies0,
		std::shared_ptr<OCL::MemObj> memobj_pairs, 
		int n)
	{
		int p = n*(n - 1) / 2;

		//frame.header = header;
		
		// this is done in the resize function
		//frame.bodies0.resize(n);
		
		//frame.pairs.resize(p);

		//TimeMeasurement tm;
		//tm.start();
		memobj_bodies0->EnqueueRead(&frame.bodies0[0], n * sizeof(Body0));
		//tm.stop();
		//printf("read buffer %i bytes in %f seconds\n", n * sizeof(Body0), tm.dif);
		
		//memobj_pairs->EnqueueRead(&frame.pairs[0], p * sizeof(Pair));
		memobj_header->EnqueueRead(&frame.header, sizeof(Header));

		frame_times.push_back(frame.header.t);

		unsigned int i = frame_times.size() - 1;

		frame.save(folder, i);
	}
	void						load()
	{
		std::ifstream myFile(folder + "\\" + "frame_times.bin", std::ios::in | std::ios::binary);

		if (!myFile.good()) throw std::exception("error with file");

		load(myFile);

		myFile.close();

		printf("frame_times loaded. m = %u t = %f ... %f\n", frame_times.size(), frame_times.front(), frame_times.back());
	}
	void						load(std::ifstream & myFile)
	{// frame times

		unsigned int s;

		myFile.read((char*)&s, sizeof(unsigned int));

		frame_times.resize(s);

		myFile.read((char *)&frame_times[0], sizeof(double)* s);

		// bodies1

		myFile.read((char*)&s, sizeof(unsigned int));

		bodies1.resize(s);

		myFile.read((char*)&bodies1[0], s * sizeof(Body1));
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

	// buffer for reading from GPU and writing to file
	Frame						frame;

public:

	std::vector<double>		frame_times;

	// Body1 from last frame, needed to resume simulation
	std::vector<Body1>		bodies1;

	std::string				folder;
};




#endif