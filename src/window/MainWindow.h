/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_WINDOW_MAINWINDOW_H_
#define FSM_WINDOW_MAINWINDOW_H_

#include <algorithm>
#include <vector>

#include <qapplication.h>
#include <qerrormessage.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qmainwindow.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpainter.h>

#include "window/ui_MainWindow.h"
#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/OutputTree.h"
#include "trees/TestSuite.h"
#include "trees/Tree.h"
#include "window/OpenFileWindow.h"

#define SCALE 0.8

namespace Ui
{
	class MainWindow;
}

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	/**
	Create the main window for the graphical interface. Every element of the
	UI (user Interface) is load from an UI file, editable with Qt designer
	@param parent The parent widget
	*/
	explicit MainWindow(QWidget * parent = 0);

	/**
	Delete the main window
	*/
	~MainWindow();
private:
	/**
	The current level of zoom for the image of the first tab (fsm graph)
	*/
	double scaleFactorFsm;

	/**
	The current level of zoom for the image of the second tab (fsm state cover graph)
	*/
	double scaleFactorStateCover;

	/**
	The current level of zoom for the image of the third tab (fsm transition cover graph)
	*/
	double scaleFactorTransitionCover;

	/**
	The UI file, loaded byt the constructor. It contain every element of the user interface
	*/
	Ui::MainWindow * ui;

	/**
	The window to open a new file (containing a new fsm)
	*/
	std::unique_ptr<OpenFileWindow> openFileWindow;

	/**
	The current FSM, the one which is displayed on the window
	*/
	std::shared_ptr<Fsm> currentFsm;

	/**
	The current index, corresponding to the current FSM
	*/
	int currentIndex;

	/**
	The list of all FSM
	*/
	std::vector<std::shared_ptr<Fsm>> fsms;

	/**
	The list of all test suite
	*/
	std::map<std::string, std::shared_ptr<TestSuite>> testSuites;

	/**
	Add a new FSM into the list
	@param newFsm The new FSM to be added into the list
	*/
	void addFsm(std::shared_ptr<Fsm> newFsm);

	/**
	Scale the image displayed into the qlabel, in the center of the screen. This method is used
	to zoom in and out WITHOUT losing a quality because it is alway loaded for the file
	@param factor The factor of zoom asked by the user (relative to the previous level of zoom)
	*/
	void scaleImage(const double & factor);

	/**
	Create the dot file for the current FSM, the state and transition cover of this FSM.
	After this step, it call dot to create a PNG file
	*/
	void createImage();

	/**
	This method is only here to make the dot system call and to warn the user
	if this system call had an issue. Please note that in order to use this
	program, you need to install dot and to add this executable to your
	PATH
	@param name The name of the file to be created
	*/
	void dotCall(const std::string & name);

	/**
	Update the view to show the current FSM into the UI. It actualise every characteristics in the
	panel at right and it load every image file in the center panel
	*/
	void updateView();

	/**
	Store several test cases into a .txt file
	@param fileName The name of the file in which you want to store the test cases
	@param testCases The test cases to store
	*/
	void storeTestCases(const std::string & fileName, const IOListContainer & testCases);

	/**
	Store ONE test case into a .txt file
	@param fileName The name of the file in which you want to store the test case
	@param testCase The test case to store
	*/
	void storeTestCase(const std::string & fileName, const std::vector<int> & testCase);

	/**
	Read several test cases from a .txt file
	@param fileName The name of the file from which you want to load the test cases
	*/
	IOListContainer readTestCases(const std::string & fileName);

	/**
	Read one test case from
	@param line The line containing the test case
	*/
	std::vector<int> readTestCase(const std::string & line);

	/**
	Store a test suite (several output tree) into a .txt file
	@param fileName The name of the file in which you want to store the test suite
	@param testSuite The test suite to store
	*/
	void storeTestSuite(const std::string & fileName, const TestSuite & testSuite);

	/**
	Store an output tree into a .txt file
	@param fileName The name of the file in which you want to store the output tree
	@param outputTree The output tree to store
	*/
	void storeOutputTree(const std::string & fileName, OutputTree & outputTree);
private slots:
	/**
	Change the current FSM displayed
	@param currentSelection The current FSM selectioned into the list
	*/
	void changeFsm(const QModelIndex & currentSelection);

	/**
	Create the windows in which you can load a new FSM from a .fsm file
	*/
	void openFile();

	/**
	Load a new FSM from a .fsm file, by using input given to the openFileWindow
	*/
	void loadFile();

	/**
	Save the current FSM into a .fsm file (and also her presentation layer into .in, .out and .state)
	*/
	void saveFile();

	/**
	Close the current FSM (unshow it from the left list)
	*/
	void closeFile();

	/**
	Transform the current FSM into an observable one and add this new FSM into the list
	*/
	void transformToObservable();

	/**
	Transform the current FSM into an minimised one and add this new FSM into the list
	*/
	void minimise();

	/**
	Calculate the caracterisation set of the current FSM and store it into a .txt file
	*/
	void calcCaracterisationSet();

	/**
	Intersect 2 FSM (the current one and one selected by the user) and add this new FSM into the list
	*/
	void intersection();

	/**
	Call the wpMethod for the current FSM (if the p chosen is big enough).
	Then, store the test case and the test suite into two .txt file
	*/
	void wpMethod();

	/**
	Call the wMethod for the current FSM (if the m chosen is big enough).
	Then, store the test case and the test suite into two .txt file
	*/
	void wMethod();

	/**
	Store a test suite into a .txt file
	*/
	void createTestSuite();

	/**
	Store an output tree into a .txt file
	*/
	void createOutputTree();

	/**
	Run a test suite equivalence between two test suite chosen by the user
	*/
	void runTestSuiteEquivalence();

	/**
	Run a test suite reduction between two test suite chosen by the user
	*/
	void runTestSuiteReduction();

	/**
	Zoom in into the center image
	*/
	void zoomIn();

	/**
	Zoom out into the center image
	*/
	void zoomOut();

	/**
	Rest the zoom of the center image to default value
	*/
	void zoomDefault();
};
#endif //FSM_WINDOW_MAINWINDOW_H_
