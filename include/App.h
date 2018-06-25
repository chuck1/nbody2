#ifndef APP_H
#define APP_H

#include <vector>

class App
{
public:
	virtual void run(std::vector<std::string>) = 0;
};

#endif