#
/*
 *    Copyright (C) 2016, 2017
 *    Jan van Katwijk (J.vanKatwijk@gmail.com)
 *    Lazy Chair Programming
 *
 *    This file is part of the  DAB-library
 *    DAB-library is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    DAB-library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with DAB-library; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *	Simple streaming server, for e.g. epg data and tpg data
 */

#ifndef	__TCP_SERVER__
#define	__TCP_SERVER__

#include	<stdint.h>
#include	<ringbuffer.h>
#include	<sys/types.h>
#include	<sys/socket.h>
#include	<netdb.h>
#include	<string>
#include	<thread>
#include	<unistd.h>
#include	<atomic>

class	tcpServer {
public:
		tcpServer	(std::string);
		~tcpServer	(void);
	void	sendData	(uint8_t *, int32_t);
	void	run		(std::string port);
private:
	std::thread	threadHandle;
	RingBuffer<uint8_t> 	*buffer;
	std::atomic<bool>	running;
	std::atomic<bool>	connected;
};
#endif

