/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "window/MainWindow.h"

MainWindow::MainWindow(QWidget * parent) :
	QMainWindow(parent),
	scaleFactorFsm(1),
	scaleFactorStateCover(1),
	scaleFactorTransitionCover(1),
	ui(new Ui::MainWindow),
	//currentFsm(std::make_shared<Fsm> ("fsm.fsm", "fsm", 3, 1, 1, std::make_shared<FsmPresentationLayer>("fsm.in", "fsm.out", "fsm.state"))),
    currentFsm(std::make_shared<Fsm> ("fsm.fsm",std::make_shared<FsmPresentationLayer>("fsm.in", "fsm.out", "fsm.state"),"fsm")),
	currentIndex(0)
{
    
	ui->setupUi(this);
	fsms.push_back(currentFsm);
	ui->opened_list->addItem(QString::fromStdString(currentFsm->getName()));
	ui->opened_list->item(0)->setSelected(true);
	createImage();
	QObject::connect(ui->opened_list->selectionModel(), SIGNAL(currentChanged(QModelIndex, QModelIndex)), this, SLOT(changeFsm(QModelIndex)));
	QObject::connect(ui->action_Open, SIGNAL(triggered()), this, SLOT(openFile()));
	QObject::connect(ui->action_Save, SIGNAL(triggered()), this, SLOT(saveFile()));
	QObject::connect(ui->action_Close, SIGNAL(triggered()), this, SLOT(closeFile()));
	QObject::connect(ui->action_Quit, SIGNAL(triggered()), qApp, SLOT(quit()));

	QObject::connect(ui->actionTransform_to_Observable, SIGNAL(triggered()), this, SLOT(transformToObservable()));
	QObject::connect(ui->action_Minimise, SIGNAL(triggered()), this, SLOT(minimise()));
	QObject::connect(ui->actionCalculate_Characterisation_set, SIGNAL(triggered()), this, SLOT(calcCaracterisationSet()));
	QObject::connect(ui->action_Intersection_with_an_other_FSM, SIGNAL(triggered()), this, SLOT(intersection()));

	QObject::connect(ui->actionW_p_Method, SIGNAL(triggered()), this, SLOT(wpMethod()));
	QObject::connect (ui->action_W_Method, SIGNAL (triggered ()), this, SLOT (wMethod ()));
	QObject::connect(ui->action_Create_test_suite, SIGNAL(triggered()), this, SLOT(createTestSuite()));
	QObject::connect(ui->actionCreate_Output_tree, SIGNAL(triggered()), this, SLOT(createOutputTree()));
	QObject::connect(ui->actionRun_test_suite_Equivalence, SIGNAL(triggered()), this, SLOT(runTestSuiteEquivalence()));
	QObject::connect(ui->actionRun_test_suite_Reduction, SIGNAL(triggered()), this, SLOT(runTestSuiteReduction()));

	QObject::connect(ui->actionZoom_In, SIGNAL(triggered()), this, SLOT(zoomIn()));
	QObject::connect(ui->actionZoom_Out, SIGNAL(triggered()), this, SLOT(zoomOut()));
	QObject::connect(ui->action_Default_Zoom, SIGNAL(triggered()), this, SLOT(zoomDefault()));
    
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::addFsm(std::shared_ptr<Fsm> newFsm)
{
	currentIndex = static_cast<int> (fsms.size());
	currentFsm = newFsm;
	fsms.push_back(newFsm);
	ui->opened_list->addItem(QString::fromStdString(newFsm->getName()));
	ui->opened_list->item(currentIndex)->setSelected(true);
	createImage();
	updateView();
}

void MainWindow::scaleImage(const double & factor)
{
	double scale;
	QPixmap pixmap;
	switch (ui->graph->currentIndex())
	{
		case 0:
			scale = scaleFactorFsm *= factor;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + ".png");
			break;
		case 1:
			scale = scaleFactorStateCover *= factor;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + "_state_cover.png");
			break;
		case 2:
			scale = scaleFactorTransitionCover *= factor;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + "_transition_cover.png");
			break;
		default:
			return;
	}
	QLabel * label = ui->graph->currentWidget()->findChild<QLabel *>();
	label->setPixmap(pixmap.scaled(scale * pixmap.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
}

void MainWindow::createImage()
{
	currentFsm->toDot(currentFsm->getName());
	dotCall(currentFsm->getName());

	std::shared_ptr<Tree> scov = currentFsm->getStateCover();
	std::ofstream stateCover(currentFsm->getName() + "_state_cover.dot");
	scov->toDot(stateCover);
	stateCover.close();
	dotCall(currentFsm->getName() + "_state_cover");

	std::shared_ptr<Tree> tcov = currentFsm->getTransitionCover();
	std::ofstream transitionCover(currentFsm->getName() + "_transition_cover.dot");
	tcov->toDot(transitionCover);
	transitionCover.close();
	dotCall(currentFsm->getName() + "_transition_cover");
}

void MainWindow::dotCall(const std::string & name)
{
	std::string command = "dot -Tpng -o " + name + ".png" + " " + name + ".dot";
	if (system(command.c_str()))
	{
		QErrorMessage message;
		message.showMessage("There is an issue with dot system call. Please install dot (or verify your PATH)");
		message.exec();
	}
}

void MainWindow::updateView()
{
	ui->characteristics_name->setText(QString::fromStdString(currentFsm->getName()));
	ui->characteristics_maxNodes->setText(QString::number(currentFsm->getMaxNodes()));
	ui->characteristics_maxInput->setText(QString::number(currentFsm->getMaxInput()));
	ui->characteristics_maxOutput->setText(QString::number(currentFsm->getMaxOutput()));
	ui->characteristics_isCompletelyDefined->setText(currentFsm->isCompletelyDefined() ? "True" : "False");
	if (currentFsm->isObservable())
	{
		ui->characteristics_isObservable->setText("True");
		ui->actionTransform_to_Observable->setEnabled(false);
	}
	else
	{
		ui->characteristics_isObservable->setText("False");
		ui->actionTransform_to_Observable->setEnabled(true);
	}

	if (currentFsm->isMinimal() == True)
	{
		ui->characteristics_isMinimal->setText("True");
		ui->action_Minimise->setEnabled(false);
	}
	else
	{
		ui->characteristics_isMinimal->setText("???");
		ui->action_Minimise->setEnabled(true);
	}

	if (currentFsm->isDeterministic())
	{
		ui->characteristics_isDeterministic->setText("True");
		ui->action_W_Method->setEnabled(true);
	}
	else
	{
		ui->characteristics_isDeterministic->setText("False");
		ui->action_W_Method->setEnabled(false);

	}
	ui->graph_fsmLabel->setPixmap(QPixmap(QString::fromStdString(currentFsm->getName()) + ".png"));
	ui->graph_stateCoverLabel->setPixmap(QPixmap(QString::fromStdString(currentFsm->getName()) + "_state_cover.png"));
	ui->graph_transitionCoverLabel->setPixmap(QString::fromStdString(currentFsm->getName()) + "_transition_cover.png");

	QFileInfo fileCaracterisationSet(QString::fromStdString(currentFsm->getName() + "_caracterisation_set.txt"));
	if (fileCaracterisationSet.exists())
	{
		std::ifstream file(currentFsm->getName() + "_caracterisation_set.txt");
	}
}

void MainWindow::storeTestCases(const std::string & fileName, const IOListContainer & testCases)
{
	std::ofstream file(fileName);
	std::vector<std::vector<int>> tests = *testCases.getIOLists();
	for (unsigned int i = 0; i < tests.size (); ++ i)
	{
		std::vector<int> testCase = tests.at(i);
		for (unsigned int j = 0; j < testCase.size(); ++ j)
		{
			file << testCase.at(j);
			if (j != testCase.size() - 1)
			{
				file << ".";
			}
		}
		if (i != tests.size () - 1)
		{
			file << std::endl;
		}
	}
	file.close();
}

void MainWindow::storeTestCase(const std::string & fileName, const std::vector<int>& testCase)
{
	std::ofstream file(fileName);
	for (unsigned int j = 0; j < testCase.size(); ++ j)
	{
		file << testCase.at(j);
		if (j != testCase.size() - 1)
		{
			file << ".";
		}
	}
	file.close();
}

IOListContainer MainWindow::readTestCases(const std::string & fileName)
{
	std::ifstream file(fileName);
	std::shared_ptr<std::vector<std::vector<int>>> iolLst = std::make_shared<std::vector<std::vector<int>>> ();
	std::string line;
	while (getline(file, line))
	{
		iolLst->push_back(readTestCase(line));
	}
	file.close();
	return IOListContainer(iolLst, currentFsm->getPresentationLayer());
}

std::vector<int> MainWindow::readTestCase(const std::string & line)
{
	std::stringstream ss(line);
	std::string input;
	std::vector<int> inputTrace;
	while (getline(ss, input, '.'))
	{
		inputTrace.push_back(stoi(input));
	}
	return inputTrace;
}

void MainWindow::storeTestSuite(const std::string & fileName, const TestSuite & testSuite)
{
	std::shared_ptr<TestSuite> testSuitePtr = std::make_shared<TestSuite>();
	for (unsigned int i = 0; i < testSuite.size(); ++ i)
	{
		testSuitePtr->push_back(testSuite.at(i));
	}

	testSuites [fileName] = testSuitePtr;//insertion
	std::ofstream file(fileName);
	for (unsigned int i = 0; i < testSuite.size(); ++ i)
	{
		OutputTree tree = testSuite.at(i);
		tree.store(file);
		if (i != testSuite.size() - 1)
		{
			file << std::endl;
		}
	}
	file.close();
}

void MainWindow::storeOutputTree(const std::string & fileName, OutputTree & outputTree)
{
	std::shared_ptr<TestSuite> testSuitePtr = std::make_shared<TestSuite>();
	testSuitePtr->push_back(outputTree);

	testSuites [fileName] = testSuitePtr;//insertion
	std::ofstream file(fileName);
	outputTree.store(file);
	file.close();
}

void MainWindow::changeFsm(const QModelIndex & currentSelection)
{
	currentIndex = currentSelection.row();
	currentFsm = fsms.at(currentIndex);
	updateView();
}

void MainWindow::openFile()
{
	openFileWindow = std::unique_ptr<OpenFileWindow>(new OpenFileWindow());
	QObject::connect(openFileWindow->getButton(), SIGNAL(clicked()), this, SLOT(loadFile()));
}

void MainWindow::loadFile()
{
    /*
	std::string fileName = openFileWindow->getFileName();
	std::string name = openFileWindow->getName();
	int maxNodes = openFileWindow->getMaxNodes();
	int maxInput = openFileWindow->getMaxInput();
	int maxOutput = openFileWindow->getMaxOutput();
	QFileInfo info = QFileInfo(QString::fromStdString(fileName));
	std::string baseName = (info.path() + "/" + info.completeBaseName()).toStdString();
	std::shared_ptr<FsmPresentationLayer> presentationLayer =
        std::make_shared<FsmPresentationLayer>(baseName + ".in", baseName + ".out", baseName + ".state");
	addFsm(std::make_shared<Fsm>(fileName, name, maxNodes, maxInput, maxOutput, presentationLayer));
	openFileWindow.reset();
     */
    
    
    std::string fileName = openFileWindow->getFileName();
    std::string name = openFileWindow->getName();
    QFileInfo info = QFileInfo(QString::fromStdString(fileName));
    std::string baseName = (info.path() + "/" + info.completeBaseName()).toStdString();
    std::shared_ptr<FsmPresentationLayer> presentationLayer =
    std::make_shared<FsmPresentationLayer>(baseName + ".in", baseName + ".out", baseName + ".state");
    addFsm(std::make_shared<Fsm>(fileName, presentationLayer, name));
    openFileWindow.reset();
    
}

void MainWindow::saveFile()
{
	QString fileName = QFileDialog::getSaveFileName(this, "FSM Viewer", QString(), "Fsm (*.fsm)");
	std::ofstream fsmFile(fileName.toStdString());
	currentFsm->dumpFsm(fsmFile);
	fsmFile.close();

	QFileInfo info = QFileInfo(fileName);
	std::string baseName = (info.path() + "/" + info.completeBaseName()).toStdString();
	std::shared_ptr<FsmPresentationLayer> presentationLayer = currentFsm->getPresentationLayer();

	std::ofstream inputFile(baseName + ".in");
	presentationLayer->dumpIn(inputFile);
	inputFile.close();

	std::ofstream outputFile(baseName + ".out");
	presentationLayer->dumpOut(outputFile);
	outputFile.close();

	std::ofstream stateFile(baseName + ".state");
	presentationLayer->dumpState(stateFile);
	stateFile.close();
}

void MainWindow::closeFile()
{
	if (fsms.size() == 1)
	{
		QErrorMessage message;
		message.showMessage("You can't close the last FSM");
		message.exec();
		return;
	}
	int temp = currentIndex;
	delete ui->opened_list->currentItem();
	fsms.erase(fsms.begin() + temp);
	updateView();
}

void MainWindow::transformToObservable()
{
	addFsm(std::make_shared<Fsm>(currentFsm->transformToObservableFSM()));
}

void MainWindow::minimise()
{
	addFsm(std::make_shared<Fsm>(currentFsm->minimise()));
}

void MainWindow::calcCaracterisationSet()
{
	IOListContainer w = currentFsm->getCaracterisationSet();
	std::ofstream out(currentFsm->getName() + "_caracterisation_set.txt");
	out << w;
	out.close();
	QMessageBox::information(this, "FSM Viewer", "The file containing the result of the caracterisation set can be found in the current directory under the name of " + QString::fromStdString(currentFsm->getName()) + "_caracterisation_set.txt");
}

void MainWindow::intersection()
{
	QStringList items;
	for (std::shared_ptr<Fsm> fsm : fsms)
	{
		items << QString::fromStdString(fsm->getName());
	}

	bool ok;
	std::string item = QInputDialog::getItem(this, "FSM Viewer", "Select the other FSM", items, 0, false, &ok).toStdString();
	if (!ok)
	{
		return;
	}

	unsigned int i = 0;
	for (; i < fsms.size(); ++ i)
	{
		if (item == fsms.at(i)->getName())
		{
			break;
		}
	}

	std::shared_ptr<Fsm> otherFsm = fsms.at(i);
	std::shared_ptr<FsmPresentationLayer> presentationLayer1 = currentFsm->getPresentationLayer();
	std::shared_ptr<FsmPresentationLayer> presentationLayer2 = otherFsm->getPresentationLayer();
	if (presentationLayer1->compare(presentationLayer2) == false)
	{
		QErrorMessage message;
		message.showMessage("You cannot intersect these 2 FSM because they doesn't have the same input/output");
		message.exec();
		return;
	}

	addFsm(std::make_shared<Fsm>(currentFsm->intersect(*otherFsm)));
}

void MainWindow::wpMethod()
{
	int maxNodes = currentFsm->getMaxNodes();
	bool ok;
	int p = QInputDialog::getInt(this, "FSM Viewer", "Input the 'p' for Wp-method. Please remember that p must be higher or equal to the number of nodes (" + QString::number(maxNodes) + ")", maxNodes, maxNodes, 2147483647, 1, &ok);
	if (!ok)
	{
		return;
	}

	std::shared_ptr<FsmPresentationLayer> presentationLayer = currentFsm->getPresentationLayer();
	IOListContainer w = currentFsm->wpMethod(p);
	storeTestCases(currentFsm->getName() + "_test_cases.txt", w);

	TestSuite testSuite;
	for (std::vector<int> tc : *w.getIOLists())
	{
		testSuite.push_back(currentFsm->apply(InputTrace(tc, presentationLayer)));
	}
	storeTestSuite(currentFsm->getName() + "_test_suite.txt", testSuite);
	QMessageBox::information(this, "FSM Viewer", "The file containing the result of the Wp method (the test case and the test suite) can be found in the current directory under the name of " + QString::fromStdString(currentFsm->getName()) + "_test_cases.txt and " + QString::fromStdString(currentFsm->getName()) + "_test_suite.txt");
}

void MainWindow::wMethod()
{
	int maxNodes = currentFsm->getMaxNodes();
	bool ok;
	int m = QInputDialog::getInt(this, "FSM Viewer", "Input the 'm' for W-method. Please remember that m must be higher or equal to the number of nodes (" + QString::number(maxNodes) + ")", maxNodes, maxNodes, 2147483647, 1, &ok);
	if (!ok)
	{
		return;
	}

	Dfsm dfsm = *currentFsm;
	std::shared_ptr<FsmPresentationLayer> presentationLayer = currentFsm->getPresentationLayer();
	IOListContainer w = dfsm.wpMethod(m);
	storeTestCases(currentFsm->getName() + "_test_cases.txt", w);

	TestSuite testSuite;
	for (std::vector<int> tc : *w.getIOLists())
	{
		testSuite.push_back(currentFsm->apply(InputTrace(tc, presentationLayer)));
	}
	storeTestSuite(currentFsm->getName() + "_test_suite.txt", testSuite);
	QMessageBox::information(this, "FSM Viewer", "The file containing the result of the W method (the test case and the test suite) can be found in the current directory under the name of " + QString::fromStdString(currentFsm->getName()) + "_test_cases.txt and " + QString::fromStdString(currentFsm->getName()) + "_test_suite.txt");
}

void MainWindow::createTestSuite()
{
	QString fileName = QFileDialog::getOpenFileName(this, "FSM Viewer", QString(), "Test cases (*.txt)");
	std::shared_ptr<FsmPresentationLayer> presentationLayer = currentFsm->getPresentationLayer();
	IOListContainer w = readTestCases(fileName.toStdString());

	TestSuite testSuite;
	for (std::vector<int> tc : *w.getIOLists())
	{
		testSuite.push_back(currentFsm->apply(InputTrace(tc, presentationLayer)));
	}

	fileName = QFileDialog::getSaveFileName(this, "FSM Viewer", QString(), "Test suite (*.txt)");
	storeTestSuite(fileName.toStdString(), testSuite);
}

void MainWindow::createOutputTree()
{
	QString fileName = QFileDialog::getOpenFileName(this, "FSM Viewer", QString(), "Test case (*.txt)");
	std::shared_ptr<FsmPresentationLayer> presentationLayer = currentFsm->getPresentationLayer();
	InputTrace inputTrace = InputTrace(readTestCase (fileName.toStdString()), presentationLayer);
	OutputTree outputTree = currentFsm->apply(inputTrace);

	fileName = QFileDialog::getSaveFileName(this, "FSM Viewer", QString(), "Test suite (*.txt)");
	storeOutputTree(fileName.toStdString(), outputTree);
}

void MainWindow::runTestSuiteEquivalence()
{
	if (testSuites.empty())
	{
		QErrorMessage message;
		message.showMessage("There is no test suite created ! Please apply w(p) method before");
		message.exec();
		return;
	}

	QStringList items;
	for (auto it = testSuites.cbegin(); it != testSuites.cend(); ++ it)
	{
		items << QString::fromStdString(it->first);
	}

	bool ok;
	std::string item1 = QInputDialog::getItem(this, "FSM Viewer", "Select the first test suite", items, 0, false, &ok).toStdString();
	if (!ok)
	{
		return;
	}
	TestSuite testSuite = *testSuites.at(item1);

	std::string item2 = QInputDialog::getItem(this, "FSM Viewer", "Select the second test suite", items, 0, false, &ok).toStdString();
	if (!ok)
	{
		return;
	}
	TestSuite otherTestSuite = *testSuites.at(item2);
	
	if (testSuite.isEquivalentTo(otherTestSuite))
	{
		QMessageBox::information(this, "FSM Viewer", "The 2 test suite are equivalent");
	}
	else
	{
		QMessageBox::information(this, "FSM Viewer", "The 2 test suite are NOT equivalent");
	}
}

void MainWindow::runTestSuiteReduction()
{
	if (testSuites.empty())
	{
		QErrorMessage message;
		message.showMessage("There is no test suite created ! Please apply w(p) method before");
		message.exec();
		return;
	}

	QStringList items;
	for (auto it = testSuites.cbegin(); it != testSuites.cend(); ++ it)
	{
		items << QString::fromStdString(it->first);
	}

	bool ok;
	std::string item1 = QInputDialog::getItem(this, "FSM Viewer", "Select the first test suite", items, 0, false, &ok).toStdString();
	if (!ok)
	{
		return;
	}
	TestSuite testSuite = *testSuites.at(item1);

	std::string item2 = QInputDialog::getItem(this, "FSM Viewer", "Select the second test suite", items, 0, false, &ok).toStdString();
	if (!ok)
	{
		return;
	}
	TestSuite otherTestSuite = *testSuites.at(item2);
	
	if (testSuite.isReductionOf(otherTestSuite))
	{
		QMessageBox::information(this, "FSM Viewer", "The second test suite is a reduction of the first one");
	}
	else
	{
		QMessageBox::information(this, "FSM Viewer", "The second test suite is a NOT reduction of the first one");
	}
}

void MainWindow::zoomIn()
{
	scaleImage(1 / SCALE);
}

void MainWindow::zoomOut()
{
	scaleImage(SCALE);
}

void MainWindow::zoomDefault()
{
	QPixmap pixmap;
	switch (ui->graph->currentIndex())
	{
		case 0:
			scaleFactorFsm = 1;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + ".png");
			break;
		case 1:
			scaleFactorStateCover = 1;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + "_state_cover.png");
			break;
		case 2:
			scaleFactorTransitionCover = 1;
			pixmap = QPixmap(QString::fromStdString(currentFsm->getName()) + "_transition_cover.png");
			break;
		default:
			return;
	}
	QLabel * label = ui->graph->currentWidget()->findChild<QLabel *>();
	label->setPixmap(pixmap.scaled(pixmap.size()));
}
