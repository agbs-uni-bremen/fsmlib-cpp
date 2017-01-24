/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_WINDOW_OPENFILEWINDOW_H_
#define FSM_WINDOW_OPENFILEWINDOW_H_

#include <memory>
#include <iostream>

#include <qfiledialog.h>
#include <qformlayout.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qmessagebox.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qwidget.h>

class OpenFileWindow : public QWidget
{
	Q_OBJECT
private:
	/**
	Name of the file to open, chosen by the user
	*/
	QString fileName;

	/**
	Layout used by the window
	*/
	QFormLayout layout;

	/**
	Label to explain that the following string is the name of the file to open
	*/
	QLabel fName;

	/**
	Name of the FSM to open, chosen by the user
	*/
	QLineEdit name;

	/**
	Spinbox to choose the number of nodes
	*/
	//QSpinBox nodes;

	/**
	Spinbox to choose the number of inputs
	*/
	//QSpinBox input;

	/**
	Spinbox to choose the number of outputs
	*/
	//QSpinBox output;

	/**
	Button to load the FSM and close this window
	*/
	QPushButton submit;
public:
	/**
	Create a window to choose which FSM to load into the main window
	*/
	OpenFileWindow();
	
	/**
	Gettter for the submit button (needed to connect it into the main window)
	@return The button
	*/
	QPushButton * getButton();
	
	/**
	Gettter for  the name of the file chosen by the user
	@return The file name
	*/
	std::string getFileName();
	
	/**
	Gettter for  the name of the fsm chosen by the user
	@return The name
	*/
	std::string getName();
	
	/**
	Gettter for  the number of nodes chosen by the user
	@return The number of nodes
	*/
	//int getMaxNodes();
	
	/**
	Gettter for  the number of input chosen by the user
	@return The number of input
	*/
	//int getMaxInput();
	
	/**
	Gettter for  the number of output chosen by the user
	@return The number of output
	*/
	//int getMaxOutput();
};
#endif //FSM_WINDOW_OPENFILEWINDOW_H_
