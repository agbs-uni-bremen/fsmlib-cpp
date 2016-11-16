#include <iostream>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <interface/FsmPresentationLayer.h>
#include <fsm/Dfsm.h>
#include <fsm/Fsm.h>
#include <fsm/IOTrace.h>
#include <trees/IOListContainer.h>
#include <trees/OutputTree.h>


void test1() {
    
    std::cout << "TC-DFSM-0001 Show that Dfsm.applyDet() deals correctly with incomplete DFSMs "
    << std::endl;
    
    std::shared_ptr<FsmPresentationLayer> pl = std::make_shared<FsmPresentationLayer>();
    Dfsm d("m1.fsm",pl,"m1");
    d.toDot("m1");
 
    std::vector<int> inp;
    inp.push_back(1);
    inp.push_back(2);
    inp.push_back(1);
    inp.push_back(2);
    inp.push_back(2);


    InputTrace i(inp,pl);
    
    std::cout << "InputTrace = " << i << std::endl;

    
    IOTrace t = d.applyDet(i);
    
    std::cout << "IOTrace t = " << t << std::endl;
    
    
    inp.insert(inp.begin(),9);
    InputTrace j(inp,pl);
    IOTrace u = d.applyDet(j);
    std::cout << "IOTrace u = " << u << std::endl;

}

int main(int argc, char* argv[])
{
    
    test1();
    
    exit(0);
    
}


#if 0
/*
\file main.c This file contain the main function
\author Gaël Dottel
\version 0.1
\date 06 June 2016
*/

int main(int argc, char* argv[])
{
	std::string fsmName = argv[2];

	if (argc >= 6)
	{
		std::string inFile;
		std::string outFile;
		std::string stateFile;
		if (argc > 6)
		{
			inFile = argv[6];
		}
		if (argc > 7)
		{
			outFile = argv[7];
		}
		if (argc > 8)
		{
			stateFile = argv[8];
		}

		std::shared_ptr<FsmPresentationLayer> presentationLayer = std::make_shared<FsmPresentationLayer>(inFile, outFile, stateFile);
		Fsm fsm = Fsm(argv[1], fsmName, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), presentationLayer);
		fsm.toDot(fsmName);

		Fsm fsmObs = fsm.transformToObservableFSM();
		fsmObs.toDot(fsmObs.getName());

		Fsm fsmMin = fsm.minimise();
		fsmMin.toDot(fsmMin.getName());

		IOListContainer w = fsmMin.wpMethod(6);
		std::cout << "WP-Method Test Cases" << std::endl << "Number of Test Cases " << w.size() << std::endl << "Test Cases:" << std::endl << w << std::endl;

		int i = 0;
		for (std::vector<int> tc : *w.getIOLists())
		{
			OutputTree trOut = fsm.apply(InputTrace(tc, presentationLayer));
			std::cout << ++i << ". " << trOut << std::endl;
		}
		getchar();//in order to pause the program on windows*/
	}
	return EXIT_SUCCESS;
}
#endif
