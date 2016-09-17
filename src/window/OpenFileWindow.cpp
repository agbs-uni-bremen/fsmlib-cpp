/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "window/OpenFileWindow.h"

OpenFileWindow::OpenFileWindow()
{
	fileName = QFileDialog::getOpenFileName(this, "Open a file", QString(), "Fsm (*.fsm)");
	fName.setText(fileName);
	submit.setText("Load the file");
	layout.addRow("File name", &fName);
	layout.addRow("Name of the FSM", &name);
	//layout.addRow("Number of nodes", &nodes);
	//layout.addRow("Maximal Input", &input);
	//layout.addRow("Maximal Output", &output);
	layout.addRow("", &submit);
	this->setLayout(&layout);
	this->show();
}

QPushButton * OpenFileWindow::getButton()
{
	return &submit;
}

std::string OpenFileWindow::getFileName()
{
	return fileName.toStdString();
}

std::string OpenFileWindow::getName()
{
	return name.text().toStdString();
}

/*int OpenFileWindow::getMaxNodes()
{
	return nodes.value();
}

int OpenFileWindow::getMaxInput()
{
	return input.value();
}

int OpenFileWindow::getMaxOutput()
{
	return output.value();
}*/
