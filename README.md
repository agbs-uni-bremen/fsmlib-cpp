# fsmlib-cpp
A C++ library containing algorithms for processing finite state machines and deriving test cases from FSMs

 0. Contributers and licence
 Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 Licensed under the EUPL V.1.1

 Please send questions and suggestions to 
        Jan Peleska (jp@informatik.uni-bremen.de) and/or 
        Gaël Dottel (dottel.gael@gmail.com).
 
 1. About this library
 This library provides algorithms for processing finite state machines (FSMs), such as, for example,
       - Minimising deterministic or nondeterministic state machines
       - Making nondeterministic state machines observable
       - Building the intersection of two state machines over the same input/output alphabets
 Moreover, the library contains "classical" test data generation algorithms, such as the W-Method and the Wp-Method.

 The underlying methods have been described in Part II of the lecture notes 
      Wen-ling Huang and Jan Peleska
      Test Automation - Foundations and Applications of Model-based Testing
      http://www.informatik.uni-bremen.de/agbs/jp/papers/test-automation-huang-peleska.pdf
 The material presented in this part of the lecture notes and the associated algorithms implemented in fsmlib-cpp is based on (sometimes quite old) results elaborated by other authors, these are all referenced in the lecture notes. 
 
 We do hope that the algorithms presented here are useful to researchers and practitioners working with FSMs and who are interested in model-based testing and complete testing theories and their application is practise.
 
 2. Structure of this repository
 This repository consists of source code contained in sub-folders of directory [src](src) this folder contains
    - several classes programmed in C++, with the main classes contained in folder [fsm](src/fsm);
      there the main classes are [Dfsm](src/fsm/Dfsm.h) and [Fsm](src/fsm/Fsm.h)
	- utility libraries in [trees](src/trees) and [sets](src/sets) as well as an example program in [main](src/main);
    - a Qt-based gui application in windows/;
 Some of the algorithms can be executed directly from the main program (if you build this, it is called fsm-main), in order to exemplify what the library functions can do. The library, however, contains also methods that are not exercised by fsm-main program. Therefore it is useful to browse all classes for other functionality that might also be useful for your applications.

Documentation can be found in folder doc/. File 'doc/doxyfile' can be used to create a class documentation using the doxygen tool. This will be quite useful to explore the contents of the library. If you do not wish to create the documentation yourself, you can unpack the zip-archive doc/fsmlib-cpp-doc.zip and open file html/index.html in your browser.
 
 3. How to build the code for different platforms
 
 We have prepared the code to be compiled and linked using the Cmake toolkit. For creating and using the main program fsm-main, a Qt-installation of Qt 5.7 is required. The compiler needs to support C++11. We have tested the code with these compilers: 
   gcc 5.4.0 under linux (ubuntu) - 64 bit
   gcc 4.8.3 under Linux (centos 7) - 64 bit
   MSVC 2015 under windows 10 - 64 and 32 bit
   clang-703.0.31 under MAC OS - 64 bit

 To compile under any of these platforms, follow these steps [some platform-specific parameters are specified below].
   - create a build directory, this will contain the results of the build process (library and main program, object files). For example, this directory could be called 'build' and located in the root directory of the repository.
   - change into this build directory
   - call cmake with the platform-sepcific commands described below 
   - this call creates the build files needed by your compiler (Makefiles, MSCV's files, etc.). Now compile the source code.
   - When you want to run the main program, change into sub-directory main which has been created by the
     the build process and run "fsm-main" (the files needed by the program should be in the same directory, if it's not, move them there) 


 3.1 Building the library for Linux
 
     For Linux platforms, the debug and release versions need to be built separately, using two different
     build directories. Command
	 
		cmake <relative path from debug build directory to the src-directory> -DCMAKE_BUILD_TYPE=Debug
     
                
     will create the makefiles for a debug version of the code. Command 
     
		cmake <relative path from release build directory to the src-directory> -DCMAKE_BUILD_TYPE=Release
		
     will create the makefiles for a release version of the code.
	 
	 
	 For building the graphical user interface as well, the above mentioned commands have to be changed to:
		
		cmake -Dgui=ON <relative path from debug build directory to the src-directory> -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/gcc_64/ -DCMAKE_BUILD_TYPE=Debug
		
	for debug mode, respectively
		
		cmake -Dgui=ON <relative path from release build directory to the src-directory> -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/gcc_64/ -DCMAKE_BUILD_TYPE=Release
		
	for release mode. 
		
 3.2 Building the library and main program for Mac OSX
 
     For Mac OSX platforms, the debug and release versions need to be built separately, using two different
     build directories (just as for Linux). Command
     
	    cmake <relative path from debug build directory to the src-directory> -DCMAKE_BUILD_TYPE=Debug
        
     will create the makefiles for a debug version of the code. Command
     
        cmake <relative path from release build directory to the src-directory> -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/clang_64/ -DCMAKE_BUILD_TYPE=Release
        
     will create the makefiles for a release version of the code.
	 
	 For building the graphical user interface use the following commands instead:
	 
	    cmake -Dgui=ON <relative path from debug build directory to the src-directory> -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/clang_64/ -DCMAKE_BUILD_TYPE=Debug
	 
	 for debug mode, respectively
	 
		cmake -Dgui=ON <relative path from release build directory to the src-directory> -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/clang_64/ -DCMAKE_BUILD_TYPE=Release
	 
	 for release mode.
     
 3.3 Building the library for Windows
 
     For Windows platforms, the debug and release versions will be built at the same time. Please note that you will have to add Qt's dll to your PATH in order to build the gui.
	 
	 Launch Cmake GUI.
	 
	 Specify the source code directory and the binaries directory.
	 
	 For building the graphical user interface the following steps are necessary:
	 
	 Enable the option `gui`
	 
	 Click Add Entry.
	 
	 Write this line for the Name. Command
	 
		CMAKE_PREFIX_PATH
		
	 Select PATH for the Type. For the value, put the absolute path to. Command
	 
		Qt/5.7/msvc2015_64
	 
	 Click Generate and verify if cmake is using the right compiler.
	 
	 Click Finish.
	 
	 will create all the files needed by MSVC. You can now open the project with the .sln file.
	 
	 
