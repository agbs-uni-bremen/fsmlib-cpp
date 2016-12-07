/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/Dfsm.h"
#include "fsm/DFSMTable.h"
#include "fsm/DFSMTableRow.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/FsmLabel.h"
#include "fsm/PkTable.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "trees/Tree.h"
#include "trees/IOListContainer.h"

void Dfsm::createAtRandom()
{
	srand((unsigned int) time(0));

	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		nodes [i] = std::make_shared<FsmNode>(i, presentationLayer);//insertion
	}

	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		std::shared_ptr<FsmNode> source = nodes.at(i);

		for (int input = 0; input <= maxInput; ++ input)
		{
			int nTarget = std::rand() % nodes.size();
			std::shared_ptr<FsmNode> target = nodes.at(nTarget);
			int output = std::rand() % (maxOutput + 1);
			FsmTransition transition = FsmTransition(source, target, FsmLabel(input, output, presentationLayer));
			source->addTransition(transition);
		}
	}
}

std::shared_ptr<DFSMTable> Dfsm::toDFSMTable() const
{
	std::shared_ptr<DFSMTable> tbl
          = std::make_shared<DFSMTable>(nodes.size(), maxInput, presentationLayer);

	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		if (nodes.at(i) == nullptr)
		{
			continue;
		}

		std::shared_ptr<DFSMTableRow> r = nodes.at(i)->getDFSMTableRow(maxInput);

		if (r == nullptr)
		{
			return nullptr;
		}
		tbl->setRow(i, r);
        
        
        
	}
    
	return tbl;
}

Dfsm::Dfsm(const std::string & fname, const std::string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Fsm(fname, fsmName, maxNodes, maxInput, maxOutput, presentationLayer)
{

}

Dfsm::Dfsm(const std::string& fname,
           const std::shared_ptr<FsmPresentationLayer> presentationLayer,
           const std::string & fsmName)
: Fsm(fname,presentationLayer,fsmName)
{
    
}

Dfsm::Dfsm(const std::string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Fsm(presentationLayer)
{
	name = fsmName;
	nodes.insert(nodes.end(), maxNodes, nullptr);
	initStateIdx = 0;
	this->maxInput = maxInput;
	this->maxOutput = maxOutput;
	currentParsedNode = nullptr;
	createAtRandom();
	std::ofstream out(getName() + ".txt");
	dumpFsm(out);
	out.close();
}

Dfsm::Dfsm(const std::string & fsmName, const int maxInput, const int maxOutput, const std::vector<std::shared_ptr<FsmNode>> lst, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Fsm(fsmName, maxInput, maxOutput, lst, presentationLayer)
{

}

Dfsm::Dfsm(const Fsm & fsm)
	: Fsm (fsm.getName(), fsm.getMaxInput(), fsm.getMaxOutput(), fsm.getNodes(), fsm.getPresentationLayer())
{
	initStateIdx = fsm.getInitStateIdx();;
	minimal = isMinimal();
	/*std::shared_ptr<FsmNode> currentParsedNode;
	std::vector<std::shared_ptr<OFSMTable>> ofsmTableLst;
	std::shared_ptr<Tree> characterisationSet;
	std::vector<std::shared_ptr<Tree>> stateIdentificationSets;*/

}

Dfsm Dfsm::minimise()
{
	dfsmTable = toDFSMTable();
	pktblLst.clear();
	std::shared_ptr<PkTable> p1 = dfsmTable->getP1Table();
	pktblLst.push_back(p1);
	std::shared_ptr<PkTable> pMin = p1;

	for (std::shared_ptr<PkTable> pk = p1->getPkPlusOneTable(); pk != nullptr; pk->getPkPlusOneTable())
	{
		pMin = pk;
		pktblLst.push_back(pk);
	}

	return pMin->toFsm(name);
}

void Dfsm::printTables() const
{
	std::ofstream file("tables.tex");
	if (dfsmTable != nullptr)
	{
		file << dfsmTable;
	}

	for (unsigned int i = 0; i < pktblLst.size(); ++ i)
	{
		file << pktblLst.at(i) << std::endl << std::endl;
	}
	file.close();
}

IOListContainer Dfsm::getCharacterisationSet()
{
	/*Create Pk-tables for the minimised FSM*/
	dfsmTable = toDFSMTable();
	pktblLst.clear();
	std::shared_ptr<PkTable> p1 = dfsmTable->getP1Table();
	pktblLst.push_back(p1);

	for (std::shared_ptr<PkTable> pk = p1->getPkPlusOneTable();
         pk != nullptr;
         pk = pk->getPkPlusOneTable())
	{
		pktblLst.push_back(pk);
	}

	/*Create an empty characterisation set as an empty InputTree instance*/
	std::shared_ptr<Tree> w = std::make_shared<Tree>(std::make_shared<TreeNode>(), presentationLayer);

	/*Loop over all non-equal pairs of states. If they are not already distinguished by 
	the input sequences contained in w, create a new input traces that distinguishes them
	and add it to w.*/
	for (unsigned int left = 0; left < nodes.size(); ++ left)
	{
		std::shared_ptr<FsmNode> leftNode = nodes.at(left);
		for (unsigned int right = left + 1; right < nodes.size(); ++ right)
		{
			std::shared_ptr<FsmNode> rightNode = nodes.at(right);

			if (leftNode->distinguished(rightNode, w) != nullptr)
			{
				continue;
			}

			/*We have to create a new input trace and add it to w, because
			leftNode and rightNode are not distinguished by the current
			input traces contained in w. This step is performed
			according to Gill's algorithm.*/
			InputTrace i = leftNode->calcDistinguishingTrace(rightNode, pktblLst, maxInput);
			std::shared_ptr<std::vector<std::vector<int>>> lli = std::make_shared<std::vector<std::vector<int>>>();
			lli->push_back(i.get());
			IOListContainer tcli = IOListContainer(lli, presentationLayer);
			w->addToRoot(tcli);
		}
	}

	/*Wrap list of lists by an IOListContainer instance*/
	IOListContainer tcl = w->getIOLists();
	return tcl;
}

IOTrace Dfsm::applyDet(const InputTrace & i)
{
	OutputTrace o = OutputTrace(presentationLayer);

	std::shared_ptr<FsmNode> currentNode = nodes.at(initStateIdx);
    
    // Apply input trace to FSM, as far as possible
    for ( int input : i.get() ) {
        if ( currentNode == nullptr ) break;
        currentNode = currentNode->apply(input, o);
    }

    // Handle the case where the very first input is not accpeted
    // by the incomplete DFSM, or even the initial node does not exist:
    // we return an empty IOTrace
	if (currentNode == nullptr && o.get().empty())
	{
		return IOTrace(InputTrace(presentationLayer),
                       OutputTrace(presentationLayer));
	}

    // Handle the case where only a prefix of the input trace
    // has been accepted by the incomplete DFSM: we return
    // an IOTrace whose input consist of this prefix, together
    // with the associated outputs already contained in o.
	if (currentNode == nullptr)
	{
        
        // Constant iterator to start of input trace.
        auto ifirst = i.cbegin();
        // Iterator pointing BEHIND the last input applied
        // @note The number of inputs processed so far equals o.size()
        auto ilast = ifirst + o.get().size();
        
        // Constant iterator to start of output trace.
        auto ofirst = o.cbegin();
        // Iterator pointing BEHIND last obtained output.
        auto olast = ofirst + o.get().size();
        
        return IOTrace(InputTrace(std::vector<int>(ifirst, ilast), presentationLayer),
                       OutputTrace(std::vector<int>(ofirst, olast), presentationLayer));
        
	}
    
    // The full input trace has been processed by the DFSM.
    // The associated outputs are contained in o.
	return IOTrace(InputTrace(i.get(), presentationLayer),
                   OutputTrace(o.get(), presentationLayer));
    
}

bool Dfsm::pass(const IOTrace & io)
{
	IOTrace myIO = applyDet(io.getInputTrace());
	return myIO.getOutputTrace() == io.getOutputTrace();
}

IOListContainer Dfsm::wMethod(const unsigned int m)
{
	if (m < nodes.size())
	{
		std::cout << "Illegal value " << m << " of m. Must be greater or equal " << nodes.size() << std::endl;
		exit(EXIT_FAILURE);
	}

	std::shared_ptr<Tree> iTree = getTransitionCover();

	std::cout << "Transition Cover = " << iTree->getIOLists() << std::endl;


	if (m > nodes.size())
	{
		IOListContainer inputEnum = IOListContainer(maxInput,
                                                    1,
                                                    m - static_cast<int> (nodes.size()),
                                                    presentationLayer);
		iTree->add(inputEnum);
	}

	IOListContainer iolc = getCharacterisationSet();
	std::cout << "Charset = " << iolc << std::endl;
	iTree->add(iolc);
	return iTree->getIOLists();
}

IOListContainer Dfsm::wpMethod(const int m)
{
	Fsm fMin = minimiseObservableFSM();
	return fMin.wpMethod(m);

}


IOListContainer Dfsm::tMethod()
{
    
    std::shared_ptr<Tree> iTree = getTransitionCover();
    
    return iTree->getIOLists();
    
}














