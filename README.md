
DAB COMMAND LINE and DAB LIBRARY

======================================================================


=======================================================================

There is an obvious need - at least felt by me - to experiment with other (forms of) GUI(s) for a DAB handling program, using the same mechanism - preferably the same code - to handle the DAB data stream. That is why a choice was made to pack the full DAB handling as a library. 

The library provides entries for the functionality through some simple calls, while a few callback functions provide the communication back from the library to the gui.
The library interface is given in dab-api.h

To show the use of the library, several example programs are included:

	- The example programs example-1, example-2, example-3, example-4,
	  example-5 and example-6 are regular C (C++) programs. 
	  
	- simpleDab is a Qt based GUI program, linking to the library.
	  It shows the use of the library when handled from with a Qt GUI.

	- python-example contains an example program in python to use (an
	  extended form of) the library.
	  An additional file, dab-python.cpp, contains
	  the sources for binding C and Python for this linbrary.

Note that the basic idea of providing these example is to provide
examples of how to use the library and/or the sources of the library. 

THEY ARE NOT INTENDED TO BE FULL SOLUTIONS, BUT MERELY THERE TO GIVE YOU AN 
IDEA HOW TO USE THE LIBRARY OR ITS SOURCES.

Feel free to adapt each (all) of these programs to your own ideas.

=======================================================================

LIBRARY CONFIGURATION PARAMETER (IMPORTANT)

For e.g. the RPI 2, where the CPU has more cores, but the capacity
of each core is limited, the best option is to split the processing
of the frontend into two parts. 
A configuration parameter in the CMakeLists.txt file controls this

If __THREADED_DECODING is defined, either in the file
"library/includes/ofdm/ofdm-decoder.h" or with a line
add_definitions (-D__THREADED_DECODING) in the CMakeLists.txt file,

The "ofdm_decoding" (involving a lot of FFT's) will be done in a
separate thread. This is the way I run it on an RPI2.

The drawback is that synchronization might take longer.


========================================================================
 
=======================================================================

The C (C++) example programs
------------------------------------------------------------------------

As said, there are 6 versions of an example dab command line program,
they are written in C, communicate with the library functions through callbacks.
The programs are all more or less the same, with some differences in
the functionality provided.

For each of the programs, a CMakeLists.txt file exists with which a
Makefile can be generated using Cmake.

The standard way to create an executable is

	cd example-X   (replace X by the appropriate digit)
	mkdir build
        cd build
        cmake .. -DXXX=ON
        make
        sudo make install

where XXX is one of the supported input devices, i.e. SDRPLAY, AIRSPY,
RTLSDR, WAVFILES, or RTL_TCP. The name of the generated executable
is dab-xxx-y, where xxx is the device name, and y the number
of the example to which the executable belongs.

(A WAVFILE is an ".sdr" file, it is a regular ".wav" file, however, 
generated by the qt-dab program).

The executable will be installed in /usr/local/bin, so yo need to have
permissions (you can obviously also just run the generated program
from the build directory).

Invocation of the program, with some parameters specified, then is
     
	dab-sdrplay-x -M 1 -B "BAND III" -C 12C -P "Radio 4" -G 80 -A default

The example programs are different though:

	- example 1 is the example that dynamically links to the DAB library,
	  i.e. the DAB library should be pre-installed
	  example 1 does provide support for neither WAVFILES nor RTL_TCP.
	  Note that compiling the main program requires the availability
	  of both the portaudio library (version 19 rather than 18)
	  and libsamplerate. It goes without saying that the library
	  for supporting the selected device also should be available.

	- example 2 has the same functionality as example 1, the sources
          of the library are "compiled-in", however. Additionally it provides
	  support for sending the Tpeg data to a server (if so configured).
	  Note that compiling the main program requires the availability
	  of both the portaudio library (version 19 rather than 18)
	  and libsamplerate. It goes without saying that the library
	  for supporting the selected device also should be available.

	- example 3 has the same functionality as example 2, and here
	  the library sources are "compiled in" as well. However, the
	  PCM samples are being sent out to stdout.
	  One might use one of the available programs to make the sound
	  audible
	  dab-example-3 .... | aplay -r 48000 -f S16_LE -t raw -c 2

	- example 4 has the sample functionality as examples 2 and 3, and
	  here the library sources are "compiled in" as well. However,
	  no sound decoding takes place. The MP2 frames (in case of DAB)
	  or the AAC frames (in case of DAB+) are just emitted through stdout.
	  (Note that the AAC frames have 960 rather than 1024 samples)

	- example 5 is a small experimental extension to example 2,
	  It contains a simple "keyboard listener", that will react
	  on entering a stroke on the return key. It will cause the
	  "next" (audio) service to be selected.

	- example 6 is an experimental version where control is
	  through an IP port.

=======================================================================

A note on the callback functions

-----------------------------------------------------------------------

The library (whether separate or compiled in) sends its data to the
main program using callbacks. These callbacks, the specification of
which is given in the file dab-api.h, are implemented here as
simple C functions. WHAT MUST BE NOTED IS THAT THESE FUNCTIONS ARE
EXECUTED IN THE THREAD OF THE CALLER, and while the library is built around
more than a single thread, it is wise to add some locking when extending
the callback functions.

========================================================================

Creating the library
------------------------------------------------------------------------------

The library can be created by - if needed - adapting the
`CMakeLists.txt` file in the dab-library/library directory and running

	mkdir build 
	cd build 
	cmake .. 
	make 
	sudo make install
	
from within the dab-library directory.

IMPORTANT: YOU NEED C++11 SUPPORT FOR THIS

Note that contrary to earlier versions, the "device" is NOT part of the library,
the user program has to provide some functions to the library for getting samples.
The interface can be found in the file "device-handler.h". 

===============================================================================

Libraries (together with the "development" or ".h" files) needed for creating the library are

	libfaad
	libfftw3f
	libusb-1.0
	zlib


============================================================================


For "simpleDab" one uses qt-make, there is a ".pro" file

For the python-example read the README file in the python-example directory.
HOWEVER: before running the example program one has to create an
ADAPTED library.
The CMakeLists.txt file for creating such an adapted library is in the python
directory.


=============================================================================

Command-line Parameters for the C (C++) versions
-----------------------------------------------------------------------

The programs accept the following command line parameters:

	-B Band
selects the DAB band (default Band III),

	-M Mode
selects the DAB Mode (default Mode 1),

	-C the channel
the default is 11C, the channel I am listening to mostly,

	-P the program name
a prefix suffices. For e.g. "Classic FM" it suffices to give "Classic". However, when passing on a non-unique prefix (e.g. "Radio" for "Radio Maria" and "Radio Veronica") the software will select one arbitrarily. Note that letter case is IMPORTANT is the current version. The names of the programs in the ensemble being received in the selected channel will be printed during recognition.

Important: If no program names are found, or if no match can be made between the program name and the list of program names, the program has no other choice than to halt, what it does.

	-G the gain 
to be applied on the device, a value in the range from 1 .. 100.
The value will be translated to an acceptable value for the device. In case the gain is table driven, as in the case of a dabstick, a value of e.g. 75 is translated into the element on three quarters of the table (basic assumption is that the table elements are more or less linear). For e.g. the Airspy the values are mapped upon the range 0 .. 21 of the sensitivity slider. 
Note that when using the rtl_tcp interface, this does not hold. The sound
setting is passed on to the server.

	-W waiting time
the maximum time to wait for valid data.
If no valid data is found within this period, execution of the program will
stop.

Note that in example-2 the -W is replaced by TWO parameters, a -d xxx indicating
the maximum amount of time to reach time synchronization (which is implicit
in the other examples) and a -D xxx parameter indicating the maximum amount of
time to get the signal "ensemble found".

	-A the output channel (example 1 and 2 only)
again as with the program name, a prefix of the name suffices. As with the programs, the names of the sound channels identified will be printed. Note, however, that in Linux not all all names appearing on the namelist are useful,
some of them will just not work, a well known  issue with the combination portaudio/alsa under Linux. 
Important: If a name is selected for a channel that cannot be opened the program will try to open the default output device.

	-O filename or "-" (example 2 only) 
The PCM samples of the sound output are stored in the file <filename>. If "-"
is specified as filename the output is just written to stdout.
This output then can be made audible by some other program.


	-H hostname (example 2, 3 and 4 only)
If rtl_tcp is selected as input device, the -H option allows selection
of a hostname. Default is "127.0.0.1".

	-I port (example 2, 3, and 4 only)
If rtl_tcp is selected as input device, the -I option allows selection
of a port. Default is 1234.

For each of the parameters there is a default, i.e., if the command

	dab-cmdline-x
	
is given, the assumptions are 

1. the Mode is "1",
2. the band is "BAND III", that the channel selected is "11C",
3. the program we are looking for is "Classic FM", and 
4. the device to send the output to is "default". 
 
Note again, that the choice for the input device was fixed when creating the dab-library.

An example of a full specification of the command line is

	dab-cmdline -M 1 -B "BAND III" -C 12C -P "Radio 4" -G 80 -A default

The program - when started - will try to identify a DAB datastream in the selected channel (here channel 12C). If a valid DAB datastream is found, the names of the programs in the ensemble will be printed. If - after a given amount of time - no ensemble could be found, execution will halt, if the ensemble was identified, the data for the selected program (here Radio 4) will be decoded.

=========================================================================

simpleDAB
-------------------------------------------------------------------------

The simpleDAB directory contains the files for a simplified Qt GUI.
It binds to the dab-library, and can be created by qmake/make. No attempts
is made to create a CmakeLists.txt file, since the program is merely
an example to demonstrate the use of the library in a Qt context.

========================================================================

The API
-------------------------------------------------------------------------

The API specification, in dab-api.h, contains a specification of the
types for the callback functions and a specification for the real API functions.

===============================================================================

E X P E R I M E N T A L
--------------------------------------------------------------------------------

One of the issues still to be resolved is the handling of data. As an
experiment a callback function was added that is called from within the
tdc handler. In example-2 a simple TCP server was added, one that just
writes out packaged tdc frames.
The package structure is : an 8 byte header followed by the frame data.
The header starts with a -1 0 -1 0 pattern, followed by a two byte length,
followed by a zero, followed by a 0 for frametype 0 and 0xFF for frametype 1.
Install the server by adding "-DSERVER" to the cmake command line.

A simple "reader" (client), using qt is included in the sources.

===============================================================================

Copyrights

-----------------------------------------------------------------------------------
	
	Copyright (C)  2016, 2017
	Jan van Katwijk (J.vanKatwijk@gmail.com)
	Lazy Chair Programming

The dab-library software is made available under the GPL-2.0. The dab-library uses a number of GPL-ed libraries, all
rigfhts gratefully acknowledged.
All SDR-J software, among which dab-library is one - is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 	GNU General Public License for more details.

