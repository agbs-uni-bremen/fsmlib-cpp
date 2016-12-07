/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/InputTrace.h"
#include "fsm/OFSMTable.h"
#include "sets/HittingSet.h"
#include "trees/TreeNode.h"
#include "trees/OutputTree.h"
#include "trees/Tree.h"
#include "trees/IOListContainer.h"
#include "trees/TestSuite.h"

std::shared_ptr<FsmNode> Fsm::newNode(const int id, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p)
{
	std::string nodeName = std::string("(" + p->first->getName() + std::to_string(p->first->getId()) + ","
									   + p->second->getName() + std::to_string(p->second->getId()) + ")");
	std::shared_ptr<FsmNode> n = std::make_shared<FsmNode>(id, nodeName, presentationLayer);
	n->setPair(p);
	return n;
}

bool Fsm::contains(const std::vector<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p)
{
	for (std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> pLst : lst)
	{
		if (*pLst == *p)
		{
			return true;
		}
	}
	return false;
}

bool Fsm::contains(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<FsmNode> n)
{
	for (std::shared_ptr<FsmNode> nLst : lst)
	{
		if (nLst->isDerivedFrom(n->getPair()))
		{
			return true;
		}
	}
	return false;
}

std::shared_ptr<FsmNode> Fsm::findp(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p)
{
	for (std::shared_ptr<FsmNode> nLst : lst)
	{
		if (nLst->isDerivedFrom(p))
		{
			return nLst;
		}
	}
	return nullptr;
}

void Fsm::parseLine(const std::string & line)
{
	std::stringstream ss(line);

	int source;
	int input;
	int output;
	int target;
	ss >> source;
	ss >> input;
	ss >> output;
	ss >> target;

	if (source < 0 || static_cast<int> (nodes.size()) <= source)
	{
		return;
	}
	if (target < 0 || static_cast<int> (nodes.size()) <= target)
	{
		return;
	}
	if (input < 0 || maxInput < input)
	{
		return;
	}
	if (output < 0 || maxOutput < output)
	{
		return;
	}

	/*First node number occurring in the file defines the initial state*/
	if (initStateIdx < 0)
	{
		initStateIdx = source;
	}

	if (currentParsedNode == nullptr)
	{
		currentParsedNode = std::make_shared<FsmNode>(source, name, presentationLayer);
		nodes [source] = currentParsedNode;
	}
	else if (currentParsedNode->getId() != source && nodes [source] == nullptr)
	{
		currentParsedNode = std::make_shared<FsmNode>(source, name, presentationLayer);
		nodes [source] = currentParsedNode;
	}
	else if (currentParsedNode->getId() != source)
	{
		currentParsedNode = nodes [source];
	}

	if (nodes [target] == nullptr)
	{
		nodes [target] = std::make_shared<FsmNode>(target, name, presentationLayer);
	}

	currentParsedNode->addTransition(FsmTransition(currentParsedNode, nodes [target], FsmLabel(input, output, presentationLayer)));
}

void Fsm::parseLineInitial (const std::string & line)
{
    std::stringstream ss(line);
    
    int source;
    int input;
    int output;
    int target;
    ss >> source;
    ss >> input;
    ss >> output;
    ss >> target;
    
    if ( source > maxState ) maxState = source;
    if ( target > maxState ) maxState = target;
    if ( input > maxInput ) maxInput = input;
    if ( output > maxOutput ) maxOutput = output;
    
}


void Fsm::readFsm(const std::string & fname)
{
    
    // Read the FSM file first to determine maxInput, maxOutput, maxState
    readFsmInitial(fname);
    
    // Create the node vector, but first with null-nodes only
    for ( int n = 0; n <= maxState; n++ ) {
        nodes.push_back(nullptr);
    }
    
    // Now read FSM file again to specify the FSM nodes and their transitions
    
	/* Mark that the initial state has not yet been determined
	(will be done in parseLine()) */
	initStateIdx = -1;
	std::ifstream inputFile(fname);
	if (inputFile.is_open())
	{
		std::string line;
		while (getline(inputFile, line))
		{
			parseLine(line);
		}
		inputFile.close();
        
	}
	else
	{
		std::cout << "Unable to open input file" << std::endl;
		exit(EXIT_FAILURE);
	}

}

void Fsm::readFsmInitial (const std::string & fname)
{
    

    initStateIdx = -1;
    std::ifstream inputFile (fname);
    if (inputFile.is_open())
    {
        std::string line;
        while (getline (inputFile, line))
        {
            parseLineInitial (line);
        }
        inputFile.close ();
    }
    else
    {
        std::cout << "Unable to open input file" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}


std::string Fsm::labelString(std::unordered_set<std::shared_ptr<FsmNode>>& lbl) const
{
	std::string s = "{ ";

	bool isFirst = true;
	for (std::shared_ptr<FsmNode> n : lbl)
	{
		if (!isFirst)
		{
			s += ",";
		}
		isFirst = false;
		s += n->getName() + "(" + std::to_string(n->getId()) + ")";
	}

	s += " }";
	return s;
}

Fsm::Fsm(const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: currentParsedNode(nullptr), initStateIdx(-1), characterisationSet(nullptr), presentationLayer(presentationLayer), minimal(Maybe)
{

}

Fsm::Fsm(const std::string& fname,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer,
         const std::string& fsmName)
: name(fsmName), currentParsedNode(nullptr),
  maxInput(-1), maxOutput(-1), maxState(-1),
  characterisationSet(nullptr), presentationLayer(presentationLayer), minimal(Maybe)
{
    readFsm(fname);
    
}

Fsm::Fsm(const std::string & fname,
         const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: name(fsmName), currentParsedNode(nullptr),
      maxInput(maxInput), maxOutput(maxOutput), maxState(maxNodes),
      characterisationSet(nullptr), presentationLayer(presentationLayer), minimal(Maybe)
{
    
    for (int i = 0; i < maxNodes; ++ i)
    {
        nodes.push_back (nullptr);
    }
    readFsm (fname);
    
}

Fsm::Fsm(const std::string & fsmName,
         const int maxInput,
         const int maxOutput,
         const std::vector<std::shared_ptr<FsmNode>> lst,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: name(fsmName), currentParsedNode(nullptr), maxInput(maxInput), maxOutput(maxOutput), initStateIdx(0), characterisationSet(nullptr), presentationLayer(presentationLayer), minimal(Maybe)
{
	nodes.insert(nodes.end(), lst.begin(), lst.end());
}

void Fsm::dumpFsm(std::ofstream & outputFile) const
{
	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		std::vector<FsmTransition> transitions = nodes.at(i)->getTransitions();
		for (unsigned int j = 0; j < transitions.size(); ++ j)
		{
			FsmTransition tr = transitions.at(j);
			outputFile << i << " " << tr.getLabel().getInput() << " " << tr.getLabel().getOutput() << " " << tr.getTarget()->getId();
			if (j < transitions.size() - 1 || i < nodes.size() - 1)
			{
				outputFile << std::endl;
			}
		}
	}
}

std::shared_ptr<FsmNode> Fsm::getInitialState() const
{
	return nodes.size() > 0 ? nodes.at(initStateIdx) : nullptr;
}

std::string Fsm::getName() const
{
	return name;
}

int Fsm::getMaxNodes() const
{
	return static_cast<int> (nodes.size());
}

int Fsm::getMaxInput() const
{
	return maxInput;
}

int Fsm::getMaxOutput() const
{
	return maxOutput;
}

std::vector<std::shared_ptr<FsmNode>> Fsm::getNodes() const
{
	return nodes;
}

std::shared_ptr<FsmPresentationLayer> Fsm::getPresentationLayer() const
{
	return presentationLayer;
}

int Fsm::getInitStateIdx() const
{
	return initStateIdx;
}

void Fsm::resetColor()
{
	for (std::shared_ptr<FsmNode> node : nodes)
	{
		node->setColor(FsmNode::white);
	}
}

void Fsm::toDot(const std::string & fname)
{
	std::ofstream out(fname + ".dot");
	out << *this;
	out.close();
}

Fsm Fsm::intersect(const Fsm & f)
{
	/*A list of node pairs which is used to control the breath-first search (BFS)*/
	std::vector<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>> nodeList;

	/*A list of new FSM states, each state created from a pair of this-nodes
	and f-nodes. At the end of this operation, the new FSM will be created from this list*/
	std::vector<std::shared_ptr<FsmNode>> fsmInterNodes;
	int id = 0;

	/*Initially, add the pair of initial this-node and f-node into the BFS list*/
	nodeList.push_back(std::make_shared<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>(getInitialState(), f.getInitialState()));

	/*This is the BFS loop, running over the (this,f)-node pairs*/
	while (!nodeList.empty())
	{
		/*Remove the head of the list and use p to refer to it
		p refers to the SOURCE node pair, from where all outgoing transitions
		are investigated in this loop cycle*/
		std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p = nodeList.front();
		nodeList.erase(nodeList.begin());

		/*current node of this FSM*/
		std::shared_ptr<FsmNode> myCurrentNode = p->first;

		/*current node of the f-FSM*/
		std::shared_ptr<FsmNode> theirCurrentNode = p->second;

		/*Do we already have an FSM state for the new FSM
		stored in fsmInterNodes, which is associated with the current pair p?*/
		std::shared_ptr<FsmNode> nSource = findp(fsmInterNodes, p);

		if (nSource == nullptr)
		{
			/*We create the new FSM state associated with p:
			nSource is created from the state pair (myCurrentNode,theirCurrentNode)
			which is identified by p.*/
			nSource = newNode(id ++, p);
			fsmInterNodes.push_back(nSource);
		}

		/*Mark this node: now all of its outgoing transitions are constructed*/
		nSource->setVisited();

		/*Loop over all transitions emanating from myCurrentNode*/
		for (FsmTransition tr : myCurrentNode->getTransitions())
		{
			/*Loop over all transitions emanating from theirCurrentNode*/
			for (FsmTransition trOther : theirCurrentNode->getTransitions())
			{
				/*If tr and trOther have identical labels, we can create a transition
				for the new FSM to be created. The transition has source node
				(myCurrentNode,theirCurrentNode), label tr.getLabel() which is the same as
				the label associated with the other transition, and target node
				(tr.getTarget(),trOther.getTarget()), which is the pair of the target nodes
				of each transition.*/
				if (tr.getLabel() == trOther.getLabel())
				{
					/*New target node represented as a pair (this-node,f-node)*/
					auto pTarget = std::make_shared<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>(tr.getTarget(), trOther.getTarget());

					/*If the target node does not yet exist in the list of state for the new FSM,
					then create it now*/
					std::shared_ptr<FsmNode> nTarget = findp(fsmInterNodes, pTarget);
					if (nTarget == nullptr)
					{
						nTarget = newNode(id ++, pTarget);
						fsmInterNodes.push_back(nTarget);
					}

					/*Add transition from nSource to nTarget*/
					FsmTransition newTr = FsmTransition(nSource, nTarget, tr.getLabel());
					nSource->addTransition(newTr);

					/*Conditions for insertion of the target pair into the nodeList:
					1. the target node corresponding to the pair has not yet been processed
					   (that is,  nTarget.hasBeenVisited() == false)
					2. The target pair is not already entered into the nodeList*/
					if (!(nTarget->hasBeenVisited() || contains(nodeList, pTarget)))
					{
						nodeList.push_back(pTarget);
					}
				}
			}
		}
	}

	return Fsm(f.getName(), maxInput, maxOutput, fsmInterNodes, presentationLayer);
}

std::shared_ptr<Tree> Fsm::getStateCover()
{
	resetColor();
	std::vector<std::shared_ptr<FsmNode>> bfsLst;
	std::unordered_map<std::shared_ptr<FsmNode>, std::shared_ptr<TreeNode>> f2t;

	std::shared_ptr<TreeNode> root = std::make_shared<TreeNode>();
	std::shared_ptr<Tree> scov = std::make_shared<Tree>(root, presentationLayer);

	std::shared_ptr<FsmNode> initState = getInitialState();
	initState->setColor(FsmNode::grey);
	bfsLst.push_back(initState);
	f2t [initState] = root;//insertion

	while (!bfsLst.empty())
	{
		std::shared_ptr<FsmNode> thisNode = bfsLst.front();
		bfsLst.erase(bfsLst.begin());
		std::shared_ptr<TreeNode> currentTreeNode = f2t.at(thisNode);

		for (int x = 0; x <= maxInput; ++ x)
		{
			for (std::shared_ptr<FsmNode> tgt : thisNode->after(x))
			{
				if (tgt->getColor() == FsmNode::white)
				{
					tgt->setColor(FsmNode::grey);
					bfsLst.push_back(tgt);
					std::shared_ptr<TreeNode> itn = currentTreeNode->add(x);
					f2t [tgt] = itn;//insertion
				}
			}
		}
		thisNode->setColor(FsmNode::black);
	}
	return scov;
}

std::shared_ptr<Tree> Fsm::getTransitionCover()
{
	std::shared_ptr<Tree> scov = getStateCover();
	resetColor();

	std::shared_ptr<std::vector<std::vector<int>>> tlst = std::make_shared<std::vector<std::vector<int>>>();

	for (int x = 0; x <= maxInput; ++ x)
	{
		std::vector<int> l;
		l.push_back(x);
		tlst->push_back(l);
	}

	IOListContainer tcl = IOListContainer(tlst, presentationLayer);

	scov->add(tcl);

	return scov;
}

OutputTree Fsm::apply(const InputTrace & itrc)
{
	return getInitialState()->apply(itrc);
}

Fsm Fsm::transformToObservableFSM() const
{
	std::vector<std::shared_ptr<FsmNode>> nodeLst;
	std::vector<std::shared_ptr<FsmNode>> bfsLst;
	std::unordered_map<std::shared_ptr<FsmNode>, std::unordered_set<std::shared_ptr<FsmNode>>> node2Label;
	std::unordered_set<std::shared_ptr<FsmNode>> theLabel;

	theLabel.insert(getInitialState());

	int id = 0;
	std::shared_ptr<FsmNode> q0 = std::make_shared<FsmNode>(id ++, labelString(theLabel), presentationLayer);
	nodeLst.push_back(q0);
	bfsLst.push_back(q0);
	node2Label [q0] = theLabel;//insertion

	while (!bfsLst.empty())
	{
		std::shared_ptr<FsmNode> q = bfsLst.front();
		bfsLst.erase(bfsLst.begin());

		q->setColor(FsmNode::black);

		for (int x = 0; x <= maxInput; ++ x)
		{
			for (int y = 0; y <= maxOutput; ++ y)
			{
				FsmLabel lbl = FsmLabel(x, y, presentationLayer);
				theLabel.clear();

				for (std::shared_ptr<FsmNode> n : node2Label.at(q))
				{
					for (FsmTransition tr : n->getTransitions())
					{
						if (tr.getLabel() == lbl)
						{
							theLabel.insert(tr.getTarget());
						}
					}
				}

				if (!theLabel.empty())
				{
					std::vector<std::pair<std::shared_ptr<FsmNode>, std::unordered_set<std::shared_ptr<FsmNode>>>> es;
					es.insert(es.end(), node2Label.begin(), node2Label.end());

					std::shared_ptr<FsmNode> tgtNode = nullptr;

					/*Use existing node if it has the same label*/
					for (std::pair<std::shared_ptr<FsmNode>, std::unordered_set<std::shared_ptr<FsmNode>>> entry : es)
					{
						if (entry.second == theLabel)
						{
							tgtNode = entry.first;
							break;
						}
					}

					/*We need to create a new node*/
					if (tgtNode == nullptr)
					{
						tgtNode = std::make_shared<FsmNode>(id ++, labelString(theLabel), presentationLayer);
						nodeLst.push_back(tgtNode);
						bfsLst.push_back(tgtNode);
						node2Label [tgtNode] = theLabel;//insertion
					}

					/*Create the transition from q to tgtNode*/
					FsmTransition trNew = FsmTransition(q, tgtNode, lbl);
					q->addTransition(trNew);
				}
			}
		}
	}
	return Fsm(name + "_O", maxInput, maxOutput, nodeLst, presentationLayer);
}

bool Fsm::isObservable() const
{
	for (std::shared_ptr<FsmNode> node : nodes)
	{
		if (!node->isObservable())
		{
			return false;
		}
	}
	return true;
}

Minimal Fsm::isMinimal() const
{
	return minimal;
}

Fsm Fsm::minimiseObservableFSM()
{
	/*Create new list to store all existing OFSMTables*/
	ofsmTableLst.clear();

	/*Create the initial OFSMTable representing the FSM,
	where all FSM states belong to the same class*/
	std::shared_ptr<OFSMTable> tbl = std::make_shared<OFSMTable>(nodes, maxInput, maxOutput, presentationLayer);

	/*Create all possible OFSMTables, each new one from its
	predecessor, and add them to the ofsmTableLst*/
	while (tbl != nullptr)
	{
		ofsmTableLst.push_back(tbl);
		tbl = tbl->next();
	}

	/*The last OFSMTable defined has classes corresponding to
	the minimised FSM to be constructed*/
	tbl = ofsmTableLst.back();

	/*Create the minimised FSM from tbl and return it*/
	Fsm fsm = tbl->toFsm(name + "_MIN");
	fsm.minimal = True;
	return fsm;
}

Fsm Fsm::minimise()
{
	if (!isObservable())
	{
		return transformToObservableFSM().minimiseObservableFSM();
	}

	return minimiseObservableFSM();
}

bool Fsm::isCharSet(const std::shared_ptr<Tree> w) const
{
	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		for (unsigned int j = i + 1; j < nodes.size(); ++ j)
		{
			if (nodes.at(i)->distinguished(nodes.at(j), w) == nullptr)
			{
				return false;
			}
		}
	}
	return true;
}

void Fsm::minimiseCharSet(const std::shared_ptr<Tree> w)
{
	IOListContainer wcnt = w->getIOLists();
	if (wcnt.size() <= 1)
	{
		return;
	}

	for (unsigned int i = 0; i < wcnt.getIOLists()->size(); ++ i)
	{
		IOListContainer wcntNew = IOListContainer(wcnt);
		wcnt.getIOLists()->erase(wcnt.getIOLists()->begin() + i);

		std::shared_ptr<Tree> itr = std::make_shared<Tree>(std::make_shared<TreeNode>(), presentationLayer);
		itr->addToRoot(wcntNew);
		if (isCharSet(itr))
		{
			if (itr->getIOLists().size() < characterisationSet->getIOLists().size())
			{
				characterisationSet = itr;
			}
		}
		minimiseCharSet(itr);
	}
}

IOListContainer Fsm::getCharacterisationSet()
{
	if (!isObservable())
	{
		std::cout << "This FSM is not observable - cannot calculate the charactersiation set." << std::endl;
		exit(EXIT_FAILURE);
	}

	/*Call minimisation algorithm again for creating the OFSM-Tables*/
	minimise();

	/*Create an empty characterisation set as an empty InputTree instance*/
	std::shared_ptr<Tree> w = std::make_shared<Tree>(std::make_shared<TreeNode>(), presentationLayer);

	/*Loop over all non-equal pairs of states.
	Calculate the state identification sets.*/
	for (unsigned int left = 0; left < nodes.size(); ++ left)
	{
		std::shared_ptr<FsmNode> leftNode = nodes.at(left);

		for (unsigned int right = left + 1; right < nodes.size(); ++ right)
		{
			std::shared_ptr<FsmNode> rightNode = nodes.at(right);

			/*Nothing to do if leftNode and rightNode are
			already distinguished by an element of w*/
			if (leftNode->distinguished(rightNode, w) != nullptr)
			{
				continue;
			}

			/*We have to create a new input trace and add it to w, because
			leftNode and rightNode are not distinguished by the current
			input traces contained in w. */
			InputTrace i = leftNode->calcDistinguishingTrace(rightNode, ofsmTableLst, maxInput, maxOutput);
			std::shared_ptr<std::vector<std::vector<int>>> lli = std::make_shared<std::vector<std::vector<int>>>();
			lli->push_back(i.get());
			IOListContainer tcli = IOListContainer(lli, presentationLayer);

			/*Insert this also into w*/
			w->addToRoot(tcli);
		}
	}

	/*Minimise and store characterisation set*/
	characterisationSet = w;
	minimiseCharSet(w);

	/*Wrap list of lists by an IOListContainer instance*/
	IOListContainer tcl = characterisationSet->getIOLists();

	std::cout << "W = " << tcl << std::endl;
	return tcl;
}

void Fsm::calcStateIdentificationSets()
{
	if (!isObservable())
	{
		std::cout << "This FSM is not observable - cannot calculate the charactersiation set." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (characterisationSet == nullptr)
	{
		std::cout << "Missing characterisation set - exit." << std::endl;
		exit(EXIT_FAILURE);
	}

	/*Create empty state identification sets for every FSM state*/
	stateIdentificationSets.clear();

	/*Identify W by integers 0..m*/
	IOListContainer wIC = characterisationSet->getIOLists();
	std::shared_ptr<std::vector<std::vector<int>>> wLst = wIC.getIOLists();

	/*wLst.get(0) is identified with Integer(0),
	wLst.get(1) is identified with Integer(1), ...*/

	std::vector<std::vector<std::unordered_set<int>>> z;
	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		z.push_back(std::vector<std::unordered_set<int>>());
		for (unsigned int j = 0; j < nodes.size(); ++ j)
		{
			z.at(i).push_back(std::unordered_set<int>());
		}
	}

	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		std::shared_ptr<FsmNode> iNode = nodes.at(i);

		for (unsigned int j = i + 1; j < nodes.size(); ++ j)
		{
			std::shared_ptr<FsmNode> jNode = nodes.at(j);

			for (unsigned int u = 0; u < wLst->size(); ++ u)
			{
				std::vector<int> thisTrace = wLst->at(u);

				if (iNode->distinguished(jNode, thisTrace))
				{
					z.at(i).at(j).insert(u);
					z.at(j).at(i).insert(u);
				}
			}
		}
	}

	for (unsigned int i = 0; i < nodes.size(); ++ i)
	{
		std::vector<std::unordered_set<int>> iLst;
		for (unsigned int j = 0; j < nodes.size(); ++ j)
		{
			if (i == j)
			{
				continue;
			}

			iLst.push_back(z.at(i).at(j));
		}

		/*Calculate minimal state identification set for
		FsmNode i*/
		HittingSet hs = HittingSet(iLst);
		std::unordered_set<int> h = hs.calcMinCardHittingSet();

		std::shared_ptr<Tree> iTree = std::make_shared<Tree>(std::make_shared<TreeNode>(), presentationLayer);
		for (int u : h)
		{
			std::vector<int> lli = wLst->at(u);
			std::shared_ptr<std::vector<std::vector<int>>> lllli = std::make_shared<std::vector<std::vector<int>>>();
			lllli->push_back(lli);
			iTree->addToRoot(IOListContainer(lllli, presentationLayer));
		}
		stateIdentificationSets.push_back(iTree);

		/*@debug*/
		for (unsigned int n = 0; n < stateIdentificationSets.size(); ++ n)
		{
			std::cout << "W(" << n << ") = " << stateIdentificationSets.at(n)->getTestCases() << std::endl;
		}
	}
}

void Fsm::appendStateIdentificationSets(const std::shared_ptr<Tree> Wp2) const
{
	IOListContainer cnt = Wp2->getIOLists();

	for (std::vector<int> lli : *cnt.getIOLists())
	{
		InputTrace itrc = InputTrace(lli, presentationLayer);

		/*Which are the target nodes reachable via input trace lli
		in this FSM?*/
		std::unordered_set<std::shared_ptr<FsmNode>> tgtNodes = getInitialState()->after(itrc);

		for (std::shared_ptr<FsmNode> n : tgtNodes)
		{
			int nodeId = n->getId();

			/*Get state identification set associated with n*/
			std::shared_ptr<Tree> wNodeId = stateIdentificationSets.at(nodeId);

			/*Append state identification set to Wp2 tree node
			reached after applying  itrc*/
			Wp2->addAfter(itrc, wNodeId->getIOLists());
		}
	}
}


IOListContainer Fsm::wMethod(const unsigned int m) {
    return IOListContainer(nullptr,nullptr);
}


IOListContainer Fsm::wpMethod(const int m)
{
	int mMinusN = static_cast<int> (m - nodes.size());
	if (mMinusN < 0)
	{
		mMinusN = 0;
	}

	std::shared_ptr<Tree> scov = getStateCover();
	std::cout << "State cover: " << scov->getTestCases() << std::endl;

	//added and not present in java
	std::ofstream stateCover(this->getName() + "_state_cover.dot");
	scov->toDot(stateCover);
	stateCover.close();

	std::shared_ptr<Tree> tcov = getTransitionCover();
	std::cout << "Transition cover: " << tcov->getTestCases() << std::endl;

	//added and not present in java
	std::ofstream transitionCover(this->getName() + "_transitionCover.dot");
	tcov->toDot(transitionCover);
	transitionCover.close();

	tcov->remove(scov);
	std::shared_ptr<Tree> r = tcov;
	std::cout << "R: " << r->getTestCases() << std::endl;

	IOListContainer w = getCharacterisationSet();
	std::cout << "Characterisation set: " << w << std::endl;

	calcStateIdentificationSets();

	std::shared_ptr<Tree> Wp1 = scov;
	if (mMinusN > 0)
	{
		IOListContainer inputEnum = IOListContainer(maxInput, 1, mMinusN, presentationLayer);
		Wp1->add(inputEnum);
	}
	Wp1->add(w);
	std::cout << "Wp1 = " << Wp1->getIOLists() << std::endl;

	std::shared_ptr<Tree> Wp2 = r;
	if (mMinusN > 0)
	{
		IOListContainer inputEnum = IOListContainer(maxInput, mMinusN, mMinusN, presentationLayer);
		Wp2->add(inputEnum);
	}
	appendStateIdentificationSets(Wp2);
	std::cout << "Wp2 = " << Wp2->getIOLists() << std::endl;

	Wp1->unionTree(Wp2);
	return Wp1->getIOLists();
}

TestSuite Fsm::createTestSuite(const IOListContainer & testCases)
{
	std::shared_ptr<std::vector<std::vector<int>>> tcLst = testCases.getIOLists();
	TestSuite theSuite;

	for (unsigned int i = 0; i < tcLst->size(); ++ i)
	{
		OutputTree ot = apply(InputTrace(tcLst->at(i), presentationLayer));
		theSuite.push_back(ot);
	}

	return theSuite;
}

bool Fsm::isCompletelyDefined() const
{
	bool cDefd = true;
	for (std::shared_ptr<FsmNode> nn : nodes)
	{
		for (int x = 0; x <= maxInput; ++ x)
		{
			bool found = false;
			for (FsmTransition tr : nn->getTransitions())
			{
				if (tr.getLabel().getInput() == x)
				{
					found = true;
					break;
				}
			}
			if (!found)
			{
				std::cout << "Incomplete FSM : for state " << nn->getName() << " " << nn->getId() << ", input " << x << " does not have a transition." << std::endl;
				cDefd = false;
			}
		}
	}
	return cDefd;
}

bool Fsm::isDeterministic() const
{
	for (std::shared_ptr<FsmNode> node : nodes)
	{
		if (!node->isDeterministic())
		{
			return false;
		}
	}
	return true;
}

void Fsm::setPresentationLayer(const std::shared_ptr<FsmPresentationLayer> ppresentationLayer)
{
	presentationLayer = ppresentationLayer;
}

std::ostream & operator<<(std::ostream & out, const Fsm & fsm)
{
	out << "digraph g {" << std::endl << std::endl << "node [shape = circle]" << std::endl << std::endl;
	for (int i = 0; i < static_cast<int> (fsm.nodes.size()); ++ i)
	{
		if (i == fsm.initStateIdx)
		{
			out << std::endl << "node [shape = doublecircle]" << std::endl;
		}

		if (fsm.nodes.at(i) == nullptr)
		{
			continue;
		}
		std::string nodeName = (fsm.nodes.at(i)->getName().empty()) ? "s" : fsm.nodes.at(i)->getName();
		out << i << "[label=\"" << nodeName << "(" << i << ")\"];" << std::endl;

		if (i == fsm.initStateIdx)
		{
			out << std::endl << "node [shape = circle]" << std::endl;
		}
	}

	for (std::shared_ptr<FsmNode> node : fsm.nodes)
	{
		if (node == nullptr)
		{
			continue;
		}
		out << *node;
	}
	out << std::endl << "}" << std::endl;
	return out;
}
