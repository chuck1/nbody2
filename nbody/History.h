
#include <string>
#include <fstream>



class History
{
public:
	History(std::string f) : folder(f)
	{}

	void						write()
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
	void						write(std::ofstream & myFile)
	{
		unsigned int l = frame_times.size();

		myFile.write((char *)&l, sizeof(unsigned int));

		myFile.write((char *)&frame_times[0], sizeof(double)* l);
	}

	void						push(Header & header, std::shared_ptr<OCL::MemObj> memobj_bodies, std::shared_ptr<OCL::MemObj> memobj_pairs, int n)
	{
		//frames.emplace_back();
		//Frame & frame = frames.back();
		Frame frame;

		int p = n*(n - 1) / 2;

		frame.header = header;
		frame.bodies.resize(n);
		frame.pairs.resize(p);
		memobj_bodies->EnqueueRead(&frame.bodies[0], n * sizeof(Body));
		memobj_pairs->EnqueueRead(&frame.pairs[0], p * sizeof(Pair));

		frame_times.push_back(header.t);

		unsigned int i = frame_times.size() - 1;

		frame.save(folder, i);

		//frame.print();
	}

	void						load()
	{
		std::ifstream myFile(folder + "\\" + "frame_times.bin", std::ios::in | std::ios::binary);

		unsigned int s;

		myFile.read((char*)&s, sizeof(unsigned int));

		frame_times.resize(s);

		myFile.read((char *)&frame_times[0], sizeof(double)* s);

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




