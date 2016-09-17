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
 This repository consists of source code contained in sub-folders of directory src/; this folder contains
    - several classes programmed in C++, with the main classes contained in folder fsm/;
      there the main classes are Dfsm and Fsm
    - a Qt-based main program main.cpp in folder
    - a main window implementation MainWindow.h, .cpp in folder window
 Some of the algorithms can be executed directly from the main program (if you build this, it is called fsm-main), in order to exemplify what the library functions can do. The library, however, contains also methods that are not exercised by fsm-main program. Therefore it is useful to browse all classes for other functionality that might also be useful for your applications.

The repository also contains a file 'doxyfile' which can be used to create a class documentation using the doxygen tool. This will be quite useful to explore the contents of the library.
 
 3. How to build the code for different platforms
 We have prepared the code to be compiled and linked using the Cmake toolkit. For creating and using the main program fsm-main, a Qt-installation of Qt 5.7 is required. The compiler needs to support C++11. We have tested the code with these compilers: 
   gcc 5.4.0 under linux (ubuntu) - 64 bit
   gcc 4.8.3 under Linux (centos 7) - 64 bit
   MSVC 2015 under windows 10 - 64 and 32 bit
   clang-703.0.31 under MAC OS - 64 bit

 To compile under any of these platforms, follow these steps [some platform-specific parameters are specified below].
   - create a build directory, this will contain the results of the build process (library and main program, object files).
     For example, this directory could be called 'build' and located in the root directory of the repository.
   - change into this build directory
   - call cmake like this: 
        cmake <relative path from build directory to the src-directory> \
               -DCMAKE_PREFIX_PATH=<path to Qt 5.7, compiler-specific sub-directory>
   - this call creates the build files needed by your compiler (Makefiles, MSCV's files, etc.). Now compile the source code.
   - When you want to run the main program, change into sub-directory main which has been created by the
     the build process and run "fsm-main" (the files needed by the program should be in the same directory, if it's not, move them there) 


 3.1 Building the library and main program for Linux
     For Linux platforms, the debug and release versions need to be built separately, using two different
     build directories. Command
        cmake <relative path from debug build directory to the src-directory> \ 
                -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/gcc_64/ -DCMAKE_BUILD_TYPE=Debug
     will create the makefiles for a debug version of the code.
        cmake <relative path from debug build directory to the src-directory> \ 
                -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/gcc_64/ -DCMAKE_BUILD_TYPE=Release
     will create the makefiles for a release version of the code.
 3.2 Building the library and main program for Mac OSX
     For Mac OSX platforms, the debug and release versions need to be built separately, using two different
     build directories (just as for Linux). Command
        cmake <relative path from debug build directory to the src-directory> \ 
                -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/clang_64/ -DCMAKE_BUILD_TYPE=Debug
     will create the makefiles for a debug version of the code.
        cmake <relative path from debug build directory to the src-directory> \ 
                -DCMAKE_PREFIX_PATH=<absolute path to Qt>/Qt/5.7/clang_64/ -DCMAKE_BUILD_TYPE=Release
     will create the makefiles for a release version of the code.
 3.3 Building the library for Windows
     If you are compiling on windows with MSVC, the tool should copy the needed Qt DLLs for the 64 bit version into the right directory. Please note that you will need to do that by hand if you are compiling for 32 bit.
