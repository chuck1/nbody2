#ifndef MY_SOCKET_H
#define MY_SOCKET_H

//
// blocking_tcp_echo_client.cpp
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <boost/asio.hpp>

//using boost::asio::ip::tcp;



int client(std::string host, std::string port, std::string message)
{
	const int max_length = 1024;

	try
	{
		boost::asio::io_service io_service;

		boost::asio::ip::tcp::resolver resolver(io_service);
		boost::asio::ip::tcp::resolver::query query(boost::asio::ip::tcp::v4(), host, port);
		boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve(query);

		boost::asio::ip::tcp::socket s(io_service);
		boost::asio::connect(s, iterator);

		std::cout << "connected" << std::endl;

		using namespace std; // For strlen.
		
		char request[max_length];
		
		strcpy_s(request, message.c_str());
		
		size_t request_length = strlen(request);
		
		boost::asio::write(s, boost::asio::buffer(request, request_length));

		std::cout << "message sent" << std::endl;

		/*char reply[max_length];
		
		size_t reply_length = boost::asio::read(s, boost::asio::buffer(reply, request_length));

		std::cout << "Reply is: ";
		std::cout.write(reply, reply_length);
		std::cout << "\n";*/
	}
	catch (std::exception& e)
	{
		std::cerr << "Exception: " << e.what() << "\n";
	}

	return 0;
}



#endif

