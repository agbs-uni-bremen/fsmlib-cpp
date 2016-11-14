/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include <qapplication.h>

#include "window/MainWindow.h"

// TODO Bonus
//test the code (huge code and not so much testing done -> see test/test.cpp)
//doxygen comment in Fsm.h and FsmNode.h
//add button for action
//add shortchut
//system call without showing console
//zoom focused
//add zoom with wheel + ctrl
//optimisation

/**
\author Gaël Dottel, Jan Peleska, Cristoph Hilken
\version 1.0
\date 2016-08-05
*/

int main(int argc, char* argv [])
{
	QCoreApplication::addLibraryPath("./");
	QApplication app(argc, argv);

	//TODO This line is a work-around for MAC. It doesn't use the native menu bar
	app.setAttribute(Qt::AA_DontUseNativeMenuBar);

	MainWindow window;
	window.show();

	return app.exec();
}
