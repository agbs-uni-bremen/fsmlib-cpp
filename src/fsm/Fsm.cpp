/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <chrono>
#include <deque>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "fsm/IOTraceContainer.h"
#include "fsm/OFSMTable.h"
#include "fsm/RDistinguishability.h"
#include "sets/HittingSet.h"
#include "trees/TreeNode.h"
#include "trees/OutputTree.h"
#include "trees/Tree.h"
#include "trees/IOListContainer.h"
#include "trees/IOTreeContainer.h"
#include "trees/TestSuite.h"
#include "trees/InputOutputTree.h"
#include "trees/AdaptiveTreeNode.h"


using namespace std;
using namespace std::chrono;

shared_ptr<FsmNode> Fsm::newNode(const int id, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p)
{
    string nodeName = string("(" + p->first->getName() + to_string(p->first->getId()) + ","
                             + p->second->getName() + to_string(p->second->getId()) + ")");
    shared_ptr<FsmNode> n = make_shared<FsmNode>(id, nodeName, presentationLayer);
    n->setPair(p);
    return n;
}

bool Fsm::contains(const vector<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>>& lst, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p)
{
    for (shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> pLst : lst)
    {
        if (*pLst == *p)
        {
            return true;
        }
    }
    return false;
}

bool Fsm::contains(const vector<shared_ptr<FsmNode>>& lst, const shared_ptr<FsmNode> n)
{
    for (shared_ptr<FsmNode> nLst : lst)
    {
        if (nLst->isDerivedFrom(n->getPair()))
        {
            return true;
        }
    }
    return false;
}

shared_ptr<FsmNode> Fsm::findp(const vector<shared_ptr<FsmNode>>& lst, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p)
{
    for (shared_ptr<FsmNode> nLst : lst)
    {
        if (nLst->isDerivedFrom(p))
        {
            return nLst;
        }
    }
    return nullptr;
}

void Fsm::parseLine(const string & line)
{
    stringstream ss(line);
    
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
        currentParsedNode = make_shared<FsmNode>(source, name, presentationLayer);
        nodes[source] = currentParsedNode;
    }
    else if (currentParsedNode->getId() != source && nodes[source] == nullptr)
    {
        currentParsedNode = make_shared<FsmNode>(source, name, presentationLayer);
        nodes[source] = currentParsedNode;
    }
    else if (currentParsedNode->getId() != source)
    {
        currentParsedNode = nodes[source];
    }
    
    if (nodes[target] == nullptr)
    {
        nodes[target] = make_shared<FsmNode>(target, name, presentationLayer);
    }
    
    shared_ptr<FsmLabel> theLabel =
    make_shared<FsmLabel>(input, output, presentationLayer);
    currentParsedNode->addTransition(make_shared<FsmTransition>(currentParsedNode,
                                                                nodes[target],
                                                                theLabel));
}

void Fsm::parseLineInitial (const string & line)
{
    stringstream ss(line);
    
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


void Fsm::readFsm(const string & fname)
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
    ifstream inputFile(fname);
    if (inputFile.is_open())
    {
        string line;
        while (getline(inputFile, line))
        {
            parseLine(line);
        }
        inputFile.close();
        
    }
    else
    {
        cout << "Unable to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    
}

void Fsm::readFsmInitial (const string & fname)
{
    
    
    initStateIdx = -1;
    ifstream inputFile (fname);
    if (inputFile.is_open())
    {
        string line;
        while (getline (inputFile, line))
        {
            parseLineInitial (line);
        }
        inputFile.close ();
    }
    else
    {
        cout << "Unable to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    
}


string Fsm::labelString(unordered_set<shared_ptr<FsmNode>>& lbl) const
{
    string s = "{ ";
    
    bool isFirst = true;
    for (shared_ptr<FsmNode> n : lbl)
    {
        if (!isFirst)
        {
            s += ",";
        }
        isFirst = false;
        s += n->getName() + "(" + to_string(n->getId()) + ")";
    }
    
    s += " }";
    return s;
}

Fsm::Fsm() { }

Fsm::Fsm(const Fsm& other) {
    
    name = other.name;
    currentParsedNode = nullptr;
    maxInput = other.maxInput;
    maxOutput = other.maxOutput;
    maxState = other.maxState;
    initStateIdx = other.initStateIdx;
    characterisationSet = nullptr;
    minimal = other.minimal;
    presentationLayer = other.presentationLayer;
    
    for ( int n = 0; n <= maxState; n++ ) {
        nodes.push_back(make_shared<FsmNode>(n,name,presentationLayer));
    }
    
    // Now add transitions that correspond exactly to the transitions in
    // this FSM
    for ( int n = 0; n <= maxState; n++ ) {
        auto theNewFsmNodeSrc = nodes[n];
        auto theOldFsmNodeSrc = other.nodes[n];
        for ( auto tr : theOldFsmNodeSrc->getTransitions() ) {
            int tgtId = tr->getTarget()->getId();
            auto newLbl = make_shared<FsmLabel>(*tr->getLabel());
            shared_ptr<FsmTransition> newTr =
            make_shared<FsmTransition>(theNewFsmNodeSrc,nodes[tgtId],newLbl);
            theNewFsmNodeSrc->addTransition(newTr);
        }
    }
    
    // Mark the initial node
    nodes[initStateIdx]->markAsInitial();
    
}

Fsm::Fsm(const shared_ptr<FsmPresentationLayer> presentationLayer)
:
name(""),
currentParsedNode(nullptr),
maxInput(-1),
maxOutput(-1),
maxState(-1),
initStateIdx(-1),
characterisationSet(nullptr),
minimal(Maybe),
presentationLayer(presentationLayer)
{
    
}

Fsm::Fsm(const string& fname,
         const shared_ptr<FsmPresentationLayer> presentationLayer,
         const string& fsmName)
:
name(fsmName),
currentParsedNode(nullptr),
maxInput(-1),
maxOutput(-1),
maxState(-1),
characterisationSet(nullptr),
minimal(Maybe),
presentationLayer(presentationLayer)
{
    readFsm(fname);
    if ( initStateIdx >= 0 ) nodes[initStateIdx]->markAsInitial();
    
}

Fsm::Fsm(const string & fname,
         const string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const shared_ptr<FsmPresentationLayer> presentationLayer)
:
name(fsmName),
currentParsedNode(nullptr),
maxInput(maxInput),
maxOutput(maxOutput),
maxState(maxNodes),
characterisationSet(nullptr),
minimal(Maybe),
presentationLayer(presentationLayer)
{
    
    for (int i = 0; i < maxNodes; ++ i)
    {
        nodes.push_back (nullptr);
    }
    readFsm (fname);
    if ( initStateIdx >= 0 ) nodes[initStateIdx]->markAsInitial();

}

Fsm::Fsm(const string & fsmName,
         const int maxInput,
         const int maxOutput,
         const vector<shared_ptr<FsmNode>> lst,
         const shared_ptr<FsmPresentationLayer> presentationLayer)
:
name(fsmName),
currentParsedNode(nullptr),
maxInput(maxInput),
maxOutput(maxOutput),
maxState((int)(lst.size()-1)),
initStateIdx(0),
characterisationSet(nullptr),
minimal(Maybe),
presentationLayer(presentationLayer)
{
    nodes.insert(nodes.end(), lst.begin(), lst.end());
    // reset all nodes as 'white' and 'unvisited'
    
    for ( auto n : nodes ) {
        n->setColor(FsmNode::white);
        n->setUnvisited();
    }
    
    nodes[initStateIdx]->markAsInitial();
    
}

Fsm::Fsm(const string & fsmName,
         const int maxInput,
         const int maxOutput,
         const vector<shared_ptr<FsmNode>> lst,
         const int initStateIdx,
         const shared_ptr<FsmPresentationLayer> presentationLayer):
    name(fsmName),
    currentParsedNode(nullptr),
    maxInput(maxInput),
    maxOutput(maxOutput),
    maxState((int)(lst.size()-1)),
    initStateIdx(initStateIdx),
    characterisationSet(nullptr),
    minimal(Maybe),
    presentationLayer(presentationLayer)
    {
        nodes.insert(nodes.end(), lst.begin(), lst.end());
        // reset all nodes as 'white' and 'unvisited'

        for ( auto n : nodes ) {
            n->setColor(FsmNode::white);
            n->setUnvisited();
        }

        nodes[initStateIdx]->markAsInitial();

    }


void Fsm::dumpFsm(ofstream & outputFile) const
{
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        vector<shared_ptr<FsmTransition> > transitions = nodes.at(i)->getTransitions();
        for (unsigned int j = 0; j < transitions.size(); ++ j)
        {
            shared_ptr<FsmTransition> tr = transitions.at(j);
            outputFile << i << " "
            << tr->getLabel()->getInput()
            << " " << tr->getLabel()->getOutput()
            << " " << tr->getTarget()->getId();
            if (j < transitions.size() - 1 || i < nodes.size() - 1)
            {
                outputFile << endl;
            }
        }
    }
}

vector<shared_ptr<FsmNode>> Fsm::getDReachableStates()
{
    resetColor();
    deque<shared_ptr<FsmNode>> bfsLst;
    vector<shared_ptr<FsmNode>> nodes;
    map<shared_ptr<FsmNode>, shared_ptr<IOTrace>> paths;

    shared_ptr<FsmNode> initState = getInitialState();
    initState->setColor(FsmNode::grey);
    bfsLst.push_back(initState);
    nodes.push_back(initState);
    initState->setDReachable(IOTrace::getEmptyTrace(presentationLayer));
    paths.insert(make_pair(initState, IOTrace::getEmptyTrace(presentationLayer)));

    while (!bfsLst.empty())
    {
        shared_ptr<FsmNode> thisNode = bfsLst.front();
        bfsLst.pop_front();

        shared_ptr<IOTrace> thisNodePath;

        if (!thisNode->isInitial())
        {
            try
            {
                thisNodePath = paths.at(thisNode);
            }
            catch (out_of_range e)
            {
                // DO nothing.
            }
        }

        for (int x = 0; x <= maxInput; ++x)
        {
            vector<int> producedOutputs;
            vector<shared_ptr<FsmNode>> successorNodes = thisNode->after(x, producedOutputs);
            if (successorNodes.size() != 1)
            {
                continue;
            }
            shared_ptr<FsmNode> tgt = successorNodes.at(0);
            try
            {
                paths.at(tgt);
                // Path alrteady exists. Do nothing.
            }
            catch (out_of_range e)
            {
                // Create new path, since it doesn't exist.
                if (thisNodePath)
                {
                    shared_ptr<IOTrace> newPath(thisNodePath);
                    newPath->append(x, producedOutputs.at(0));
                    paths.insert(make_pair(tgt, newPath));
                }
                else
                {
                    InputTrace in = InputTrace({x}, presentationLayer);
                    OutputTrace out = OutputTrace({producedOutputs.at(0)}, presentationLayer);
                    shared_ptr<IOTrace> newPath = make_shared<IOTrace>(in, out);
                    paths.insert(make_pair(tgt, newPath));
                }
            }
            if (tgt->getColor() == FsmNode::white)
            {
                tgt->setColor(FsmNode::grey);
                bfsLst.push_back(tgt);
                nodes.push_back(tgt);
                tgt->setDReachable(paths.at(tgt));
            }
        }
    }
    resetColor();

    return nodes;
}


shared_ptr<FsmNode> Fsm::getInitialState() const
{
    return nodes.size() > 0 ? nodes.at(initStateIdx) : nullptr;
}

string Fsm::getName() const
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

vector<shared_ptr<FsmNode>> Fsm::getNodes() const
{
    return nodes;
}

shared_ptr<FsmPresentationLayer> Fsm::getPresentationLayer() const
{
    return presentationLayer;
}

int Fsm::getInitStateIdx() const
{
    return initStateIdx;
}

void Fsm::resetColor()
{
    for (auto node : nodes)
    {
        node->setColor(FsmNode::white);
    }
}

void Fsm::toDot(const string & fname)
{
    ofstream out(fname + ".dot");
    out << *this;
    out.close();
}

Fsm Fsm::intersect(const Fsm & f)
{
    /*A list of node pairs which is used to control the breath-first search (BFS)*/
    vector<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>> nodeList;
    
    /*A list of new FSM states, each state created from a pair of this-nodes
     and f-nodes. At the end of this operation, the new FSM will be created from this list*/
    vector<shared_ptr<FsmNode>> fsmInterNodes;
    int id = 0;
    
    /*Initially, add the pair of initial this-node and f-node into the BFS list*/
    nodeList.push_back(make_shared<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>(getInitialState(), f.getInitialState()));
    
    /*This is the BFS loop, running over the (this,f)-node pairs*/
    while (!nodeList.empty())
    {
        /*Remove the head of the list and use p to refer to it
         p refers to the SOURCE node pair, from where all outgoing transitions
         are investigated in this loop cycle*/
        shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p = nodeList.front();
        nodeList.erase(nodeList.begin());
        
        /*current node of this FSM*/
        shared_ptr<FsmNode> myCurrentNode = p->first;
        
        /*current node of the f-FSM*/
        shared_ptr<FsmNode> theirCurrentNode = p->second;
        
        /*Do we already have an FSM state for the new FSM
         stored in fsmInterNodes, which is associated with the current pair p?*/
        shared_ptr<FsmNode> nSource = findp(fsmInterNodes, p);
        
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
        for (auto tr : myCurrentNode->getTransitions())
        {
            /*Loop over all transitions emanating from theirCurrentNode*/
            for (auto trOther : theirCurrentNode->getTransitions())
            {
                /*If tr and trOther have identical labels, we can create a transition
                 for the new FSM to be created. The transition has source node
                 (myCurrentNode,theirCurrentNode), label tr.getLabel() which is the same as
                 the label associated with the other transition, and target node
                 (tr.getTarget(),trOther.getTarget()), which is the pair of the target nodes
                 of each transition.*/
                
                if (*tr->getLabel() == *trOther->getLabel())
                {
                    
                    /*New target node represented as a pair (this-node,f-node)*/
                    auto pTarget = make_shared<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>(tr->getTarget(), trOther->getTarget());
                    
                    /*If the target node does not yet exist in the list of state for the new FSM,
                     then create it now*/
                    shared_ptr<FsmNode> nTarget = findp(fsmInterNodes, pTarget);
                    if (nTarget == nullptr)
                    {
                        nTarget = newNode(id ++, pTarget);
                        fsmInterNodes.push_back(nTarget);
                    }
                    
                    /*Add transition from nSource to nTarget*/
                    auto newTr = make_shared<FsmTransition>(nSource,
                                                            nTarget,
                                                            tr->getLabel());
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

shared_ptr<Tree> Fsm::getDeterministicStateCover()
{
	resetColor();
	deque<shared_ptr<FsmNode>> bfsLst;
	unordered_map<shared_ptr<FsmNode>, shared_ptr<TreeNode>> f2t;

	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	shared_ptr<Tree> dscov = make_shared<Tree>(root, presentationLayer);

	shared_ptr<FsmNode> initState = getInitialState();
	initState->setColor(FsmNode::grey);
	bfsLst.push_back(initState);
	f2t[initState] = root;

	while (!bfsLst.empty())
	{
		shared_ptr<FsmNode> thisNode = bfsLst.front();
		bfsLst.pop_front();
		shared_ptr<TreeNode> currentTreeNode = f2t[thisNode];

		for (int x = 0; x <= maxInput; ++x)
		{
			vector<shared_ptr<FsmNode>> successorNodes = thisNode->after(x);
			if (successorNodes.size() != 1)
			{
				continue;
			}
			for (shared_ptr<FsmNode> tgt : successorNodes)
			{
				if (tgt->getColor() == FsmNode::white)
				{
					tgt->setColor(FsmNode::grey);
					shared_ptr<TreeNode> itn = currentTreeNode->add(x);
					bfsLst.push_back(tgt);
					f2t[tgt] = itn;
				}
			}
		}
		thisNode->setColor(FsmNode::black);
	}
	resetColor();
	return dscov;
}

shared_ptr<Tree> Fsm::getStateCover()
{
    resetColor();
    deque<shared_ptr<FsmNode>> bfsLst;
    unordered_map<shared_ptr<FsmNode>, shared_ptr<TreeNode>> f2t;
    
    shared_ptr<TreeNode> root = make_shared<TreeNode>();
    shared_ptr<Tree> scov = make_shared<Tree>(root, presentationLayer);
    
    shared_ptr<FsmNode> initState = getInitialState();
    initState->setColor(FsmNode::grey);
    bfsLst.push_back(initState);
    f2t[initState] = root;
    
    while (!bfsLst.empty())
    {
        shared_ptr<FsmNode> thisNode = bfsLst.front();
        bfsLst.pop_front();
        shared_ptr<TreeNode> currentTreeNode = f2t[thisNode];
        
        for (int x = 0; x <= maxInput; ++x)
        {
            for (shared_ptr<FsmNode> tgt : thisNode->after(x))
            {
                if (tgt->getColor() == FsmNode::white)
                {
                    tgt->setColor(FsmNode::grey);
                    shared_ptr<TreeNode> itn = currentTreeNode->add(x);
                    bfsLst.push_back(tgt);
                    f2t[tgt] = itn;
                }
            }
        }
        thisNode->setColor(FsmNode::black);
    }
    resetColor();
    return scov;
}

shared_ptr<Tree> Fsm::getTransitionCover()
{
    shared_ptr<Tree> scov = getStateCover();
    resetColor();
    
    shared_ptr<vector<vector<int>>> tlst = make_shared<vector<vector<int>>>();
    
    for (int x = 0; x <= maxInput; ++ x)
    {
        vector<int> l;
        l.push_back(x);
        tlst->push_back(l);
    }
    
    IOListContainer tcl = IOListContainer(tlst, presentationLayer);
    
    scov->add(tcl);
    
    return scov;
}

OutputTree Fsm::apply(const InputTrace & itrc, bool markAsVisited)
{
    return getInitialState()->apply(itrc,markAsVisited);
}

void Fsm::apply(const InputTrace& input, vector<shared_ptr<OutputTrace>>& producedOutputs, vector<shared_ptr<FsmNode>>& reachedNodes) const
{
    return getInitialState()->getPossibleOutputs(input, producedOutputs, reachedNodes);
}

Fsm Fsm::transformToObservableFSM() const
{
    vector<shared_ptr<FsmNode>> nodeLst;
    vector<shared_ptr<FsmNode>> bfsLst;
    unordered_map<shared_ptr<FsmNode>, unordered_set<shared_ptr<FsmNode>>> node2Label;
    unordered_set<shared_ptr<FsmNode>> theLabel;
    
    theLabel.insert(getInitialState());
    
    int id = 0;
    shared_ptr<FsmNode> q0 = make_shared<FsmNode>(id ++, labelString(theLabel), presentationLayer);
    nodeLst.push_back(q0);
    bfsLst.push_back(q0);
    node2Label[q0] = theLabel;
    
    while (!bfsLst.empty())
    {
        shared_ptr<FsmNode> q = bfsLst.front();
        bfsLst.erase(bfsLst.begin());
        
        q->setColor(FsmNode::black);
        
        for (int x = 0; x <= maxInput; ++ x)
        {
            for (int y = 0; y <= maxOutput; ++ y)
            {
                shared_ptr<FsmLabel> lbl =
                make_shared<FsmLabel>(x, y, presentationLayer);
                theLabel.clear();
                
                for (shared_ptr<FsmNode> n : node2Label.at(q))
                {
                    for (auto tr : n->getTransitions())
                    {
                        if (*tr->getLabel() == *lbl)
                        {
                            theLabel.insert(tr->getTarget());
                        }
                    }
                }
                
                if (!theLabel.empty())
                {
                    vector<pair<shared_ptr<FsmNode>, unordered_set<shared_ptr<FsmNode> > > > es;
                    es.insert(es.end(), node2Label.begin(), node2Label.end());
                    
                    shared_ptr<FsmNode> tgtNode = nullptr;
                    
                    /*Use existing node if it has the same label*/
                    for (pair<shared_ptr<FsmNode>, unordered_set<shared_ptr<FsmNode> > > entry : es)
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
                        tgtNode = make_shared<FsmNode>(id ++, labelString(theLabel), presentationLayer);
                        nodeLst.push_back(tgtNode);
                        bfsLst.push_back(tgtNode);
                        node2Label[tgtNode] = theLabel;
                    }
                    
                    /*Create the transition from q to tgtNode*/
                    auto trNew = make_shared<FsmTransition>(q, tgtNode, lbl);
                    q->addTransition(trNew);
                }
            }
        }
    }
    return Fsm(name + "_O", maxInput, maxOutput, nodeLst, presentationLayer);
}

bool Fsm::isObservable() const
{
    for (shared_ptr<FsmNode> node : nodes)
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

bool Fsm::isComplete() const
{
    return complete;
}

Fsm Fsm::minimiseObservableFSM()
{
    /*Create new list to store all existing OFSMTables*/
    ofsmTableLst.clear();
    
    /*Create the initial OFSMTable representing the FSM,
     where all FSM states belong to the same class*/
    shared_ptr<OFSMTable> tbl = make_shared<OFSMTable>(nodes, maxInput, maxOutput, presentationLayer);
    
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

    if (failState || errorState)
    {
        for (shared_ptr<FsmNode> node : fsm.nodes)
        {
            if (!fsm.failState && failState && tbl->getS2C().at(node->getId()) == tbl->getS2C().at(failState->getId()))
            {
                cout << "  Setting failState:" << node->getName() << "(" << node->getId() << ", " << node << ")" << endl;
                fsm.failState = node;
            }
            if (!fsm.errorState && errorState && tbl->getS2C().at(node->getId()) == tbl->getS2C().at(errorState->getId()))
            {
                cout << "  Setting errorState:" << node->getName() << endl;
                fsm.errorState = node;
            }
            if ((!failState || fsm.failState) && (!errorState || fsm.errorState))
            {
                break;
            }
        }
    }
    return fsm;
}

Fsm Fsm::minimise()
{
    
    vector<shared_ptr<FsmNode>> uNodes;
    removeUnreachableNodes(uNodes);
    
    if (!isObservable())
    {
        return transformToObservableFSM().minimiseObservableFSM();
    }
    
    return minimiseObservableFSM();
}

Fsm Fsm::makeComplete(CompleteMode mode)
{
    cout << "makeComplete():" << endl;
    vector<shared_ptr<FsmNode>> newNodes = nodes;
    bool addErrorState = false;
    bool newErrorState = false;
    shared_ptr<FsmNode> errorNode = errorState;
    if (mode == ErrorState)
    {
        if (!errorNode)
        {
            errorNode = make_shared<FsmNode>(getMaxNodes() + 1, "Error", presentationLayer);
            newErrorState = true;
        }
        presentationLayer->addState2String("Error");
        for (int x = 0; x <= maxInput; ++x)
        {
            shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::EPSILON_OUTPUT, presentationLayer);
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(errorNode, errorNode, label);
            errorNode->addTransition(transition);
        }
    }
    for (shared_ptr<FsmNode> node : newNodes)
    {
        cout << "  State " << node->getName() << ": " << endl;
        for (int x = 0; x <= maxInput; ++x)
        {
            cout << "    Input " << presentationLayer->getInId(x);
            if(!node->isPossibleInput(x))
            {
                cout << " is not defined." << endl;
                if (mode == ErrorState)
                {
                    addErrorState = true;
                    shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::ERROR_OUTPUT, presentationLayer);
                    shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(node, errorNode, label);
                    node->addTransition(transition);
                }
                else if (mode == SelfLoop)
                {
                    shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::EPSILON_OUTPUT, presentationLayer);
                    shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(node, node, label);
                    node->addTransition(transition);
                }
            }
            else
            {
                cout << " is defined." << endl;
            }
        }
    }
    if (mode == ErrorState && newErrorState && addErrorState)
    {
        newNodes.push_back(errorNode);
    }
    Fsm fsmComplete(name + "_COMPLETE", maxInput, maxOutput, newNodes, presentationLayer);
    fsmComplete.complete = true;
    if (mode == ErrorState)
    {
        fsmComplete.errorState = errorNode;
    }
    fsmComplete.failState = failState;
    return fsmComplete;
}

shared_ptr<FsmNode> Fsm::getErrorState()
{
    return errorState;
}

shared_ptr<FsmNode> Fsm::getFailState()
{
    return failState;
}

bool Fsm::isCharSet(const shared_ptr<Tree> w) const
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

void Fsm::minimiseCharSet(const shared_ptr<Tree> w)
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
        
        shared_ptr<Tree> itr = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
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
    
    // Do we already have a characterisation set ?
    if ( characterisationSet != nullptr ) {
        IOListContainer tcl = characterisationSet->getIOLists();
        return tcl;
    }
    
    
    // We have to calculate teh chracterisation set from scratch
    if (!isObservable())
    {
        cout << "This FSM is not observable - cannot calculate the charactersiation set." << endl;
        exit(EXIT_FAILURE);
    }
    
    /*Call minimisation algorithm again for creating the OFSM-Tables*/
    minimise();
    
    /*Create an empty characterisation set as an empty InputTree instance*/
    shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
    
    /*Loop over all non-equal pairs of states.
     Calculate the state identification sets.*/
    for (unsigned int left = 0; left < nodes.size(); ++ left)
    {
        shared_ptr<FsmNode> leftNode = nodes.at(left);
        
        for (unsigned int right = left + 1; right < nodes.size(); ++ right)
        {
            shared_ptr<FsmNode> rightNode = nodes.at(right);
            
            /*Nothing to do if leftNode and rightNode are
             already distinguished by an element of w*/
            if (leftNode->distinguished(rightNode, w) != nullptr)
            {
                continue;
            }
            
            /*We have to create a new input trace and add it to w, because
             leftNode and rightNode are not distinguished by the current
             input traces contained in w. */
            InputTrace i = leftNode->calcDistinguishingTrace(rightNode,
                                                             ofsmTableLst,
                                                             maxInput,
                                                             maxOutput);
            shared_ptr<vector<vector<int>>> lli = make_shared<vector<vector<int>>>();
            lli->push_back(i.get());
            IOListContainer tcli = IOListContainer(lli, presentationLayer);
            
            /*Insert this also into w*/
            w->addToRoot(tcli);
        }
    }
    
    /*Minimise and store characterisation set*/
    characterisationSet = w;
    //    minimiseCharSet(w);
    
    /*Wrap list of lists by an IOListContainer instance*/
    IOListContainer tcl = characterisationSet->getIOLists();
    
    return tcl;
}

vector<shared_ptr<OutputTrace>> Fsm::getOutputIntersection(shared_ptr<FsmNode> q1, shared_ptr<FsmNode> q2, int x) const
{
    vector<shared_ptr<OutputTrace>> possibleOutputs1 = q1->getPossibleOutputs(x);
    vector<shared_ptr<OutputTrace>> possibleOutputs2 = q2->getPossibleOutputs(x);
    vector<shared_ptr<OutputTrace>> intersection;
    for (shared_ptr<OutputTrace> a : possibleOutputs1)
    {
        for (shared_ptr<OutputTrace> b : possibleOutputs2)
        {
            if (*a == *b) {
                intersection.push_back(a);
            }
        }
    }
    return intersection;
}

void Fsm::calcROneDistinguishableStates()
{
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        for (size_t j = i + 1; j < nodes.size(); ++j)
        {
            nodes.at(i)->getRDistinguishability()->addNotDistinguishable(1, nodes.at(j));
        }
    }
    nodes.at(nodes.size() - 1)->getRDistinguishability()->addNotDistinguishable(1);

    for (size_t k = 0; k < nodes.size(); ++k)
    {
        shared_ptr<FsmNode> q1 = nodes.at(k);
        vector<shared_ptr<FsmNode>> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(1);

        for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
        {
            shared_ptr<FsmNode> q2 = *it;
            for (int x = 0; x <= maxInput; ++ x)
            {
                InputTrace input = InputTrace(vector<int>({x}), presentationLayer);
                vector<shared_ptr<OutputTrace>> intersection = getOutputIntersection(q1, q2, x);
                if (intersection.size() == 0)
                {
                    vector<shared_ptr<OutputTrace>> q1Output;
                    vector<shared_ptr<OutputTrace>> q2Output;
                    q1->getPossibleOutputs(x, q1Output);
                    q2->getPossibleOutputs(x, q2Output);
                    shared_ptr<AdaptiveTreeNode> q1Root = make_shared<AdaptiveTreeNode>(x);
                    shared_ptr<AdaptiveTreeNode> q2Root = make_shared<AdaptiveTreeNode>(x);

                    for (shared_ptr<OutputTrace> trace : q1Output)
                    {
                        shared_ptr<AdaptiveTreeNode> target = make_shared<AdaptiveTreeNode>();
                        shared_ptr<TreeEdge> edge = make_shared<TreeEdge>(trace->get()[0], target);
                        q1Root->add(edge);
                    }
                    for (shared_ptr<OutputTrace> trace : q2Output)
                    {
                        shared_ptr<AdaptiveTreeNode> target = make_shared<AdaptiveTreeNode>();
                        shared_ptr<TreeEdge> edge = make_shared<TreeEdge>(trace->get()[0], target);
                        q2Root->add(edge);
                    }

                    shared_ptr<InputOutputTree> q1Tree = make_shared<InputOutputTree>(q1Root, presentationLayer);
                    shared_ptr<InputOutputTree> q2Tree = make_shared<InputOutputTree>(q2Root, presentationLayer);
                    q1->getRDistinguishability()->addAdaptiveIOSequence(q2, q1Tree);
                    q2->getRDistinguishability()->addAdaptiveIOSequence(q1, q2Tree);

                    q1->getRDistinguishability()->addDistinguishable(1, q2);
                    q2->getRDistinguishability()->addDistinguishable(1, q1);

                    q1->getRDistinguishability()->removeNotDistinguishable(1, q2);
                    break;
                }
            }
        }
    }

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        for (size_t j = 0; j < nodes.size(); ++j)
        {
            try {
                shared_ptr<InputOutputTree> tree = nodes.at(i)->getRDistinguishability()->getAdaptiveIOSequence(nodes.at(j));
                if (tree->getRoot()->isLeaf())
                {
                    continue;
                }
                cout << "σ(" << nodes.at(i)->getName() << "," << nodes.at(j)->getName() << ") = ";
                cout << *tree << endl;
            } catch (std::out_of_range e) {
               // Do nothing.
            }

        }

    }
}

void Fsm::calcRDistinguishableStates()
{

    calcROneDistinguishableStates();

    size_t limit = nodes.size() * (nodes.size() - 1) / 2;
    for (size_t l = 2; l <= limit; ++l)
    {
        cout << "################ l = " << l << " ################" << endl;
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            nodes.at(k)->getRDistinguishability()->inheritDistinguishability(l);
        }
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            shared_ptr<FsmNode> q1 = nodes.at(k);
            cout << "q1 = " << q1->getName() << ":" << endl;
            vector<shared_ptr<FsmNode>> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(l);
            for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
            {
                shared_ptr<FsmNode> q2 = *it;
                cout << "  q2 = " << q2->getName() << ":" << endl;
                for (int x = 0; x <= maxInput; ++ x)
                {
                    vector<shared_ptr<OutputTrace>> intersection = getOutputIntersection(q1, q2, x);

                    shared_ptr<AdaptiveTreeNode> q1Root = make_shared<AdaptiveTreeNode>(x);
                    shared_ptr<AdaptiveTreeNode> q2Root = make_shared<AdaptiveTreeNode>(x);
                    vector<shared_ptr<TreeEdge>> q1Edges;
                    vector<shared_ptr<TreeEdge>> q2Edges;

                    bool isDistinguishable = true;
                    for (shared_ptr<OutputTrace> inter : intersection)
                    {
                        int y = inter->get()[0];
                        unordered_set<shared_ptr<FsmNode>> afterQ1 = q1->afterAsSet(x, y);
                        unordered_set<shared_ptr<FsmNode>> afterQ2 = q2->afterAsSet(x, y);
                        shared_ptr<FsmNode> afterNode1 = *afterQ1.begin();
                        shared_ptr<FsmNode> afterNode2 = *afterQ2.begin();
                        if (afterNode1 == afterNode2 || !afterNode1->getRDistinguishability()->isRDistinguishableWith(l - 1, afterNode2))
                        {

                            isDistinguishable = false;
                            break;
                        }
                        else
                        {
                            cout << "    x = " << presentationLayer->getInId(x) << ":\t"
                            << afterNode1->getName() << " != " << afterNode2->getName()
                            << " -> " << q1->getName() << " != " << q2->getName() << endl;

                            //shared_ptr<TreeNode> target1 = make_shared<TreeNode>();
                            shared_ptr<InputOutputTree> childTree1 = afterNode1->getRDistinguishability()->getAdaptiveIOSequence(afterNode2);
                            // TODO Fix
                            //      can't find linker symbol for virtual table for `TreeEdge' value
                            // messages when debugging.
                            // Put breakpoint at following line and debug.
                            cout << "      childIO1(" << afterNode1->getName() << "," << afterNode2->getName() << "): " << *childTree1 << endl;
                            shared_ptr<AdaptiveTreeNode> childNode1 = static_pointer_cast<AdaptiveTreeNode>(childTree1->getRoot());
                            shared_ptr<TreeEdge> edge1 = make_shared<TreeEdge>(y, childNode1);
                            q1Edges.push_back(edge1);

                            //shared_ptr<TreeNode> target2 = make_shared<TreeNode>();
                            shared_ptr<InputOutputTree> childTree2 = afterNode2->getRDistinguishability()->getAdaptiveIOSequence(afterNode1);
                            cout << "      childIO2(" << afterNode2->getName() << "," << afterNode1->getName() << "): " << *childTree2 << endl;
                            shared_ptr<AdaptiveTreeNode> childNode2 = static_pointer_cast<AdaptiveTreeNode>(childTree2->getRoot());
                            shared_ptr<TreeEdge> edge2 = make_shared<TreeEdge>(y, childNode2);
                            q2Edges.push_back(edge2);
                        }
                    }
                    if (isDistinguishable)
                    {
                        for (shared_ptr<TreeEdge> edge : q1Edges)
                        {
                            q1Root->add(edge);
                        }
                        for (shared_ptr<TreeEdge> edge : q2Edges)
                        {
                            q2Root->add(edge);
                        }
                        vector<shared_ptr<OutputTrace>> q1Outputs;
                        vector<shared_ptr<OutputTrace>> q2Outputs;
                        q1->getPossibleOutputs(x, q1Outputs);
                        q2->getPossibleOutputs(x, q2Outputs);

                        for (auto it = q1Outputs.begin(); it != q1Outputs.end(); ++it)
                        {
                            bool disjunct = true;
                            for (shared_ptr<OutputTrace> inter : intersection)
                            {
                                if (**it == *inter)
                                {
                                    disjunct = false;
                                    break;
                                }
                            }
                            if (disjunct)
                            {
                                shared_ptr<AdaptiveTreeNode> target = make_shared<AdaptiveTreeNode>();
                                shared_ptr<TreeEdge> edge = make_shared<TreeEdge>((*it)->get()[0], target);
                                q1Root->add(edge);
                            }
                        }

                        for (auto it = q2Outputs.begin(); it != q2Outputs.end(); ++it)
                        {
                            bool disjunct = true;
                            for (shared_ptr<OutputTrace> inter : intersection)
                            {
                                if (**it == *inter)
                                {
                                    disjunct = false;
                                    break;
                                }
                            }
                            if (disjunct)
                            {
                                shared_ptr<AdaptiveTreeNode> target = make_shared<AdaptiveTreeNode>();
                                shared_ptr<TreeEdge> edge = make_shared<TreeEdge>((*it)->get()[0], target);
                                q2Root->add(edge);
                            }
                        }

                        //InputTrace input = InputTrace(vector<int>({x}), presentationLayer);
                        shared_ptr<InputOutputTree> q1Tree = make_shared<InputOutputTree>(q1Root, presentationLayer);
                        shared_ptr<InputOutputTree> q2Tree = make_shared<InputOutputTree>(q2Root, presentationLayer);

                        cout << "    q1Tree: " << *q1Tree << endl;
                        cout << "    q2Tree: " << *q2Tree << endl;

                        q1->getRDistinguishability()->addAdaptiveIOSequence(q2, q1Tree);
                        q2->getRDistinguishability()->addAdaptiveIOSequence(q1, q2Tree);


                        q1->getRDistinguishability()->addDistinguishable(l, q2);
                        q2->getRDistinguishability()->addDistinguishable(l, q1);
                        q1->getRDistinguishability()->removeNotDistinguishable(l, q2);
                        q2->getRDistinguishability()->removeNotDistinguishable(l, q1);
                        break;
                    }
                }
            }
        }
    }
    for (auto node : nodes)
    {
        node->getRDistinguishability()->hasBeenCalculated(true);
    }
    for (auto node : nodes)
    {
        if (node->getRDistinguishability()->isNotRDistinguishable())
        {
            minimal = False;
            return;
        }
    }
    minimal = True;
    return;
}

IOListContainer Fsm::getRStateCharacterisationSet(shared_ptr<FsmNode> node) const
{
    if (!node->getRDistinguishability()->hasBeenCalculated())
    {
        cout << "r-characterisation sets haven't been calculated yet." << endl;
        exit(EXIT_FAILURE);
    }
    IOListContainer result = IOListContainer(presentationLayer);
    cout << "r-state characterisation set for " << node->getName() << ":\n";
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            cout << "r-characterisation sets haven't been calculated yet." << endl;
            exit(EXIT_FAILURE);
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            IOListContainer container = sequence->getInputLists();
            auto set = container.getIOLists();

            cout << "σ(" << node->getName() << "," << n->getName() << "): " << container << "\n";
            for (auto trace : *set)
            {
                result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
            }
        }
    }
    cout << "SCS(" << node->getName() << ") = " << result << endl;
    return result;
}

IOTreeContainer Fsm::getAdaptiveRStateCharacterisationSet(shared_ptr<FsmNode> node) const
{
    if (!node->getRDistinguishability()->hasBeenCalculated())
    {
        cout << "r-characterisation sets haven't been calculated yet." << endl;
        exit(EXIT_FAILURE);
    }
    IOTreeContainer result = IOTreeContainer(presentationLayer);
    cout << "Adaptive r-state characterisation set for " << node->getName() << endl;
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            cout << "r-characterisation sets haven't been calculated yet." << endl;
            exit(EXIT_FAILURE);
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            cout << "σ(" << node->getName() << "," << n->getName() << "): " << *sequence << endl;
            result.addUniqueRemovePrefixes(sequence);
        }
    }
    cout << "SCS(" << node->getName() << ") = " << result << endl;
    return result;
}

IOListContainer Fsm::getRCharacterisationSet() const
{
    IOListContainer result = IOListContainer(presentationLayer);
    for (shared_ptr<FsmNode> n : nodes)
    {
        IOListContainer container = getRStateCharacterisationSet(n);
        auto set = container.getIOLists();
        for (auto t : *set)
        {
            result.addUniqueRemovePrefixes(Trace(t, presentationLayer));
        }
    }
    return result;
}

IOTraceContainer Fsm::getPossibleIOTraces(shared_ptr<FsmNode> node,
                                          shared_ptr<InputOutputTree> tree,
                                          const bool cleanTrailingEmptyTraces) const
{
    cout << "(" << node->getName() << ") " << "getPossibleIOTraces()" << endl;
    cout << "(" << node->getName() << ") " << "  node: " << node->getName()  << endl;
    cout << "(" << node->getName() << ") " << "  tree: " << *tree  << endl;
    if (tree->isEmpty())
    {
        cout << "(" << node->getName() << ") " << "  tree is empty. returning." << endl;
        std::shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
        emptyTrace->setTargetNode(node);
        return IOTraceContainer(emptyTrace, presentationLayer);
    }
    IOTraceContainer result = IOTraceContainer(presentationLayer);

    for (int y = 0; y <= maxOutput; ++y)
    {
        shared_ptr<AdaptiveTreeNode> treeRoot = static_pointer_cast<AdaptiveTreeNode>(tree->getRoot());
        int x = treeRoot->getInput();
        bool isPossibleOutput = node->isPossibleOutput(x, y);
        cout << "(" << node->getName() << ") " << "  x: " << presentationLayer->getInId(x) << endl;
        cout << "(" << node->getName() << ") " << "  y: " << presentationLayer->getOutId(y) << endl;
        cout << "(" << node->getName() << ") " << "  isPossibleOutput: " << isPossibleOutput << endl;

        if (isPossibleOutput)
        {
            unordered_set<shared_ptr<FsmNode>> nextNodes = node->afterAsSet(x, y);
            if (nextNodes.size() != 1)
            {
                cerr << "The FSM does not seem to be observable." << endl;
                exit(EXIT_FAILURE);
            }
            shared_ptr<FsmNode> nextNode = *nextNodes.begin();

            if (!tree->isDefined(y))
            {
                IOTrace trace = IOTrace(x, y, nextNode, presentationLayer);
                cout << "(" << node->getName() << ") " << "  tree is NOT defined. Adding " << trace  << endl;
                result.addUnique(trace);
            }
            else if (tree->isDefined(y))
            {
                cout << "(" << node->getName() << ") " << "  tree is defined."  << endl;
                cout << "(" << node->getName() << ") " << "    nextNode: " << nextNode->getName()  << endl;
                shared_ptr<AdaptiveTreeNode> nextTreeNode = static_pointer_cast<AdaptiveTreeNode>(treeRoot->after(y));
                cout << "(" << node->getName() << ") " << "    nextTreeNode input: " << presentationLayer->getInId(nextTreeNode->getInput())  << endl;
                shared_ptr<InputOutputTree> nextTree = make_shared<InputOutputTree>(nextTreeNode, presentationLayer);
                cout << "(" << node->getName() << ") " << "    nextTree: " << *nextTree  << endl;
                IOTraceContainer iONext = getPossibleIOTraces(nextNode, nextTree);
                cout << "(" << node->getName() << ") " << "    iONext: " << iONext  << endl;
                IOTrace trace = IOTrace(x, y, nextNode, presentationLayer);
                if (iONext.isEmpty())
                {
                    iONext.add(trace);
                }
                else
                {
                    iONext.concatenateToFront(trace);
                }
                cout << "(" << node->getName() << ") " << "    iONext: " << iONext  << endl;
                cout << "(" << node->getName() << ") " << "    cleanTrailingEmptyTraces: " << cleanTrailingEmptyTraces  << endl;
                if (cleanTrailingEmptyTraces)
                {
                    shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
                    for (IOTrace& t : *iONext.getList())
                    {
                        cout << "(" << node->getName() << ") " << "    t.size(): " << t.size() << endl;
                        cout << "(" << node->getName() << ") " << "    isSuffix: " << t.isSuffix(*emptyTrace) << endl;
                        if (t.size() > 1 && t.isSuffix(*emptyTrace))
                        {
                            cout << "(" << node->getName() << ") " << "    REMOVING EMPTY SUFFIX from " << t << endl;
                            t = IOTrace(t, -1, t.getTargetNode());
                        }
                    }
                    cout << "(" << node->getName() << ") " << "    iONext: " << iONext  << endl;
                }
                result.addUnique(iONext);
            }
        }
        cout << "(" << node->getName() << ") " << "#####################################" << endl;
    }
    cout << "(" << node->getName() << ") " << "--- result: " << result << endl;
    return result;
}

IOTraceContainer Fsm::getPossibleIOTraces(shared_ptr<FsmNode> node,
                                          const IOTreeContainer& treeContainer,
                                          const bool cleanTrailingEmptyTraces) const
{
    IOTraceContainer result = IOTraceContainer(presentationLayer);;
    for (shared_ptr<InputOutputTree> tree : *treeContainer.getList())
    {
        IOTraceContainer container = getPossibleIOTraces(node, tree, cleanTrailingEmptyTraces);
        result.addUnique(container);
    }
    return result;
}

IOTraceContainer Fsm::bOmega(const IOTreeContainer& adaptiveTestCases, const IOTrace& trace) const
{
    IOTraceContainer result = IOTraceContainer(presentationLayer);

    shared_ptr<FsmNode> initialState = getInitialState();
    if (!initialState)
    {
        return result;
    }

    unordered_set<shared_ptr<FsmNode>> successorNodes = initialState->after(trace);
    if (successorNodes.size() == 0)
    {
        return result;
    }
    if (successorNodes.size() != 1)
    {
        cerr << "The FSM does not seem to be observable." << endl;
        exit(EXIT_FAILURE);
    }
    shared_ptr<FsmNode> successorNode = *successorNodes.begin();
    return getPossibleIOTraces(successorNode, adaptiveTestCases);
}

IOTraceContainer Fsm::bOmega(const IOTreeContainer& adaptiveTestCases, const vector<shared_ptr<InputTrace>>& inputTraces) const
{
    IOTraceContainer result = IOTraceContainer(presentationLayer);
    shared_ptr<FsmNode> initialState = getInitialState();
    if (!initialState)
    {
        return result;
    }
    for (shared_ptr<InputTrace> inputTrace : inputTraces)
    {
        vector<shared_ptr<OutputTrace>> producedOutputs;
        initialState->getPossibleOutputs(*inputTrace, producedOutputs);
        for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
        {
            IOTrace iOTrace = IOTrace(*inputTrace, *outputTrace);
            IOTraceContainer produced = bOmega(adaptiveTestCases, iOTrace);
            result.addUnique(produced);
        }
    }
    return result;
}

vector<IOTraceContainer> Fsm::getVPrime()
{
    vector<IOTraceContainer> result;

    shared_ptr<Tree> detStateCover = getDeterministicStateCover();
    IOListContainer testCases = detStateCover->getDeterministicTestCases();
    /* Get raw input sequences from the determinisitc state cover. */
    shared_ptr<vector<vector<int>>> testCasesRaw = testCases.getIOLists();

    vector<vector<shared_ptr<OutputTrace>>> allPossibleOutputTraces;
    size_t iterations = 1;

    /* Get all possible output traces generated by the determinisitc state cover. */
    for (size_t i = 0; i < testCasesRaw->size(); ++i)
    {
        vector<int> testCase = testCasesRaw->at(i);
        InputTrace input = InputTrace(testCase, presentationLayer);
        vector<shared_ptr<OutputTrace>> producedOutputs;
        vector<shared_ptr<FsmNode>> reached;
        getInitialState()->getPossibleOutputs(input, producedOutputs, reached);
        iterations *= producedOutputs.size();
        allPossibleOutputTraces.push_back(producedOutputs);
    }

    vector<vector<vector<int>>> ot;
    size_t repetitions = iterations / allPossibleOutputTraces.at(0).size();
    for (size_t i = 0; i < testCasesRaw->size(); ++i)
    {
        vector<vector<int>> t;
        size_t blocks = iterations / repetitions;
        size_t outputIndex = 0;

        for (size_t b = 0; b < blocks; ++b)
        {
            for (size_t idx = 0; idx < repetitions; ++idx)
            {
                t.push_back(allPossibleOutputTraces.at(i).at(outputIndex)->get());
            }
            outputIndex = (outputIndex + 1) % allPossibleOutputTraces.at(i).size();
        }

        if (i+1 < testCasesRaw->size())
        {
            repetitions /= allPossibleOutputTraces.at(i+1).size();
        }
        ot.push_back(t);
    }

    for (size_t j = 0; j < iterations; ++j)
    {
        IOTraceContainer container = IOTraceContainer(presentationLayer);
        for (size_t i = 0; i < testCasesRaw->size(); ++i)
        {
            InputTrace input = InputTrace(testCasesRaw->at(i), presentationLayer);
            OutputTrace output = OutputTrace(ot.at(i).at(j), presentationLayer);
            IOTrace iOTrace = IOTrace(input, output);
            container.add(iOTrace);
        }
        result.push_back(container);
    }
    return result;
}

IOTraceContainer Fsm::r(std::shared_ptr<FsmNode> node,
                   const IOTrace& base,
                   const IOTrace& suffix) const
{
//    cout << "node: " << node->getName() << endl;
//    cout << "base: " << base << endl;
//    cout << "suffix: " << suffix << endl;


    IOTraceContainer result = IOTraceContainer(presentationLayer);
    vector<IOTrace> prefs = suffix.getPrefixes();
    vector<IOTrace> prefixes;
    // Remove empty sequences from prefixes.
    for (IOTrace& prefix : prefs)
    {
        if (prefix.size() != 1 ||
                (prefix.getInputTrace().get().at(0) != -1 &&
                prefix.getOutputTrace().get().at(0) != -1))
        {
            prefixes.push_back(prefix);
        }
    }
//    cout << "prefixes:" << endl;
    for (auto p : prefixes)
    {
//        cout << "  " << p << endl;
    }

    for (IOTrace prefix : prefixes)
    {
//        cout << "prefix = " << prefix << endl;
        IOTrace baseCopy = base;
        baseCopy.append(prefix);
//        cout << "v = " << baseCopy << " reaches:\n  ";
        unordered_set<shared_ptr<FsmNode>> nodes = node->after(prefix.getInputTrace(), prefix.getOutputTrace());
        for (shared_ptr<FsmNode> n : nodes)
        {
//            cout << n->getName() << ", ";
            if (n == node)
            {
//                cout << "  adding " << baseCopy << " to result" << endl;
                result.add(baseCopy);
            }
        }
//        cout << endl;
    }

//    cout << "result: " << result << endl;

    return result;
}

IOTraceContainer Fsm::rPlus(std::shared_ptr<FsmNode> node,
                   const IOTrace& base,
                   const IOTrace& suffix,
                   const IOTraceContainer& vDoublePrime) const
{
    IOTraceContainer rResult = r(node, base, suffix);
//    cout << "rResult: " << rResult << endl;
    if (node->isDReachable() && vDoublePrime.contains(*node->getDReachTrace()))
    {
//        cout << "  Adding: " << *node->getDReachTrace() << endl;
        rResult.add(*node->getDReachTrace());
//        cout << "  rResult: " << rResult << endl;
    }
    return rResult;
}

size_t Fsm::lowerBound(const IOTrace& base,
          const IOTrace& suffix,
          const vector<shared_ptr<InputTrace>>& takenInputs,
          const vector<shared_ptr<FsmNode>>& states,
          const IOTreeContainer& adaptiveTestCases,
          const IOTraceContainer& vDoublePrime,
          const vector<shared_ptr<FsmNode>>& dReachableStates) const
{
    size_t result = 0;

    IOTraceContainer testTraces = bOmega(adaptiveTestCases, takenInputs);

    for (shared_ptr<FsmNode> state : states)
    {
        const IOTraceContainer& rResult = r(state, base, suffix);
        result += rResult.size();
        if(find(dReachableStates.begin(), dReachableStates.end(), state) != dReachableStates.end()) {
            ++result;
        }

        const IOTraceContainer rPlusResult = rPlus(state, base, suffix, vDoublePrime);
        for (IOTrace& trace : *rPlusResult.getList())
        {
            IOTraceContainer traces = bOmega(adaptiveTestCases, trace);
            testTraces.remove(traces);
        }
    }
    result += testTraces.size();
    return result;
}

bool Fsm::adaptiveStateCounting(const size_t m, IOTraceContainer& observedTraces)
{
    if (!isComplete())
    {
        cerr << "This FSM may not be completely specified.";
        exit(EXIT_FAILURE);
    }

    const shared_ptr<FsmNode> failState = getFailState();

    /**
     * Adaptive test cases for this FSM (Ω).
     */
    const IOTreeContainer& adaptiveTestCases = getAdaptiveRCharacterisationSet();
    cout << "adaptiveTestCases: " << adaptiveTestCases << endl;
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    cout << "adaptiveTestCases as input traces:\n" << adaptiveList << endl;
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = getMaximalSetsOfRDistinguishableStates();
    const vector<IOTraceContainer> vPrime = getVPrime();
    const vector<shared_ptr<FsmNode>>& dReachableStates = getDReachableStates();

    shared_ptr<Tree> detStateCoverTree = getDeterministicStateCover();
    IOListContainer detStateCoverRaw = detStateCoverTree->getDeterministicTestCases();
    cout << "detStateCoverRaw: " << detStateCoverRaw << endl;
    vector<shared_ptr<InputTrace>> detStateCover;
    for (vector<int> trace : *detStateCoverRaw.getIOLists())
    {
        detStateCover.push_back(make_shared<InputTrace>(trace, presentationLayer));
    }

    /**
     * T - set of input sequences that have been followed by Ω.
     */
    vector<shared_ptr<InputTrace>> t = detStateCover;
    /**
     * T_c - set of current elements of T: those that are being considered in the search
     * through state space. The elements in T_c are the maximal sequences considered that
     * do not meet the termination criterion.
     */
    vector<shared_ptr<InputTrace>> tC = detStateCover;
    while (tC.size() != 0)
    {
        cout << "tC:\n  ";
        for (auto w : tC)
        {
            cout << *w << ", ";
        }
        cout << endl;
        cout << "t:\n  ";
        for (auto w : t)
        {
            cout << *w << ", ";
        }
        cout << endl;
        map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;

        // Applying all input traces from T_c to this FSM.
        // All observed outputs are bein recorded.
        // If the FSM enters the error state, adaptive state counting terminates.
        for (shared_ptr<InputTrace> inputTrace : tC)
        {

            cout << "  inputTrace: " << *inputTrace << endl;
            /**
             * Holds the produced output traces for the current input trace.
             */
            vector<shared_ptr<OutputTrace>> producedOutputs;
            /**
             * Holds the reached nodes for the current input trace.
             */
            vector<shared_ptr<FsmNode>> reachedNodes;

            apply(*inputTrace, producedOutputs, reachedNodes);
            cout << "    producedOutputs:\n      ";
            for (size_t i = 0; i < producedOutputs.size(); ++i)
            {
                cout << *producedOutputs.at(i);
                if (i != producedOutputs.size() - 1)
                {
                    cout << ", ";
                }
            }
            cout << endl;
            cout << "    reachedNodes:\n      ";
            for (size_t i = 0; i < reachedNodes.size(); ++i)
            {
                cout << reachedNodes.at(i)->getName();
                if (i != reachedNodes.size() - 1)
                {
                    cout << ", ";
                }
            }
            cout << endl;
            observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputs));

            if(failState && std::find(reachedNodes.begin(), reachedNodes.end(), failState) != reachedNodes.end()) {
                // FSM entered the error state.
                cout << "  Failure observed:" << endl;
                cout << "    Input Trace: " << *inputTrace << endl;
                cout << "    Produced Outputs:\n      ";
                for (size_t i = 0; i < producedOutputs.size(); ++i)
                {
                    cout << *producedOutputs.at(i);
                    if (i != producedOutputs.size() - 1)
                    {
                        cout << ", ";
                    }
                    // Adding observed traces for simple input traces only in case of error.
                    // When no error is being observed, this traces are bein used later when
                    // concatenating them with the adaptive test cases.
                    IOTrace iOTrace(*inputTrace, *producedOutputs.at(i), reachedNodes.at(i));
                    observedTraces.addUnique(iOTrace);
                }
                cout << endl;
                cout << "    Reached nodes:\n      ";
                for (size_t i = 0; i < reachedNodes.size(); ++i)
                {
                    cout << reachedNodes.at(i)->getName();
                    if (i != reachedNodes.size() - 1)
                    {
                        cout << ", ";
                    }
                }
                cout << endl;
                return false;
            }

            if (producedOutputs.size() != reachedNodes.size())
            {
                cerr << "Number of produced outputs and number of reached nodes do not match.";
                exit(EXIT_FAILURE);
            }

            // Applying adaptive test cases to every node reached by the current input trace.
            for (size_t i = 0; i< producedOutputs.size(); ++i)
            {
                shared_ptr<FsmNode> node = reachedNodes.at(i);
                shared_ptr<OutputTrace> outputTrace = producedOutputs.at(i);
                cout << "----------------- Getting adaptive traces -----------------" << endl;
                IOTraceContainer observedAdaptiveTraces = getPossibleIOTraces(node, adaptiveTestCases);
                cout << "  observedAdaptiveTraces: " << observedAdaptiveTraces << endl;
                cout << "  concatenating: " << *inputTrace << "/" << *outputTrace << endl;
                observedAdaptiveTraces.concatenateToFront(*inputTrace, *outputTrace);
                cout << "  observedAdaptiveTraces after concatenation to front: " << observedAdaptiveTraces << endl;
                observedTraces.addUnique(observedAdaptiveTraces);
                for (IOTrace& trace : *observedAdaptiveTraces.getList())
                {
                    if(failState && trace.getTargetNode() == failState) {
                        // FSM entered the error state.
                        cout << "  Failure observed:" << endl;
                        cout << "    Input Trace: " << *inputTrace << endl;
                        cout << "    Observed adaptive traces:" << endl;
                        cout << observedAdaptiveTraces << endl;
                        return false;
                    }
                }
            }
        }

        bool cancel = false;
        vector<shared_ptr<InputTrace>> newT = t;
        vector<shared_ptr<InputTrace>> newTC;;
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                if (cancel)
                {
                    break;
                }
                IOTrace currentTrace(*inputTrace, *outputTrace);
                for (const IOTraceContainer& vDoublePrime : vPrime)
                {
                    if (cancel)
                    {
                        break;
                    }
                    IOTrace* maxPrefix = nullptr;
                    for (IOTrace& iOTrace : *vDoublePrime.getList())
                    {
                        if (currentTrace.isPrefix(iOTrace) && (!maxPrefix || iOTrace.size() > maxPrefix->size()))
                        {
                            maxPrefix = &iOTrace;
                        }
                    }
                    if (!maxPrefix)
                    {
                        continue;
                    }
                    for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
                    {
                        if (cancel)
                        {
                            break;
                        }
                        for (size_t i = 0; i < rDistStates.size() - 1; ++i)
                        {
                            if (cancel)
                            {
                                break;
                            }
                            shared_ptr<FsmNode> s1 = rDistStates.at(i);
                            for (size_t j = i + 1; j < rDistStates.size(); ++j)
                            {
                                shared_ptr<FsmNode> s2 = rDistStates.at(j);
                                if (s1 == s2)
                                {
                                    continue;
                                }
                                const IOTraceContainer& s1RPlus = rPlus(s1, *maxPrefix, currentTrace, vDoublePrime);
                                const IOTraceContainer& s2RPlus = rPlus(s2, *maxPrefix, currentTrace, vDoublePrime);

                                unordered_set<shared_ptr<FsmNode>> reached1;
                                unordered_set<shared_ptr<FsmNode>> reached2;

                                for (IOTrace& trace : *s1RPlus.getList())
                                {
                                    unordered_set<shared_ptr<FsmNode>> reached = getInitialState()->after(trace);
                                    reached1.insert(reached.begin(), reached.end());

                                }
                                for (IOTrace& trace : *s2RPlus.getList())
                                {
                                    unordered_set<shared_ptr<FsmNode>> reached = getInitialState()->after(trace);
                                    reached2.insert(reached.begin(), reached.end());
                                }

                                vector<shared_ptr<FsmNode>> reached1V(reached1.begin(), reached1.end());
                                vector<shared_ptr<FsmNode>> reached2V(reached2.begin(), reached2.end());

                                if (rDistinguishesAllStates(reached1V, reached2V, adaptiveTestCases))
                                {
                                    if (lowerBound(*maxPrefix, currentTrace, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates) > m)
                                    {
                                        cancel = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (!cancel)
            {
                // Expanding sequences.
                for (int x = 0; x <= maxInput; ++x)
                {
                    shared_ptr<InputTrace> concat  = make_shared<InputTrace>(*inputTrace);
                    concat->add(x);
                    if (!InputTrace::contains(t, *concat) && !InputTrace::contains(newTC, *concat))
                    {
                        newTC.push_back(concat);
                    }
                    if (!InputTrace::contains(newT, *concat))
                    {
                        newT.push_back(concat);
                    }
                }
            }
        }
        tC = newTC;
        t = newT;
    }
    cout << "  RESULT: " << observedTraces << endl;
    return true;
}

bool Fsm::rDistinguishesAllStates(std::vector<std::shared_ptr<FsmNode>>& nodesA,
                            std::vector<std::shared_ptr<FsmNode>>& nodesB,
                            const IOTreeContainer& adaptiveTestCases) const
{
    for (size_t i = 0; i < nodesA.size(); ++i)
    {
        shared_ptr<FsmNode> nodeA = nodesA.at(i);
        for (size_t j = i + 1; j < nodesB.size(); ++j)
        {
            shared_ptr<FsmNode> nodeB = nodesB.at(j);
            if (nodeA == nodeB)
            {
                cout << nodeA->getName() << " == " << nodeB->getName() << ". This should not happen." << endl;
                return false;
            }
            if (!rDistinguishes(nodeA, nodeB, adaptiveTestCases))
            {
                return false;
            }
        }
    }
    return true;
}

bool Fsm::rDistinguishes(shared_ptr<FsmNode> nodeA,
                         shared_ptr<FsmNode> nodeB,
                         const IOTreeContainer& adaptiveTestCases) const
{
    for (shared_ptr<InputOutputTree> tree : *adaptiveTestCases.getList())
    {
        if (rDistinguishes(nodeA, nodeB, tree))
        {
            return true;
        }
    }
    return false;
}

bool Fsm::rDistinguishes(shared_ptr<FsmNode> nodeA,
                         shared_ptr<FsmNode> nodeB,
                         shared_ptr<InputOutputTree> adaptiveTestCase) const
{
    IOTraceContainer containerA = getPossibleIOTraces(nodeA, adaptiveTestCase);
    IOTraceContainer containerB = getPossibleIOTraces(nodeB, adaptiveTestCase);
    vector<OutputTrace> outputsA = containerA.getOutputTraces();
    vector<OutputTrace> outputsB = containerB.getOutputTraces();

    for (OutputTrace& traceA: outputsA)
    {
        for (OutputTrace& traceB : outputsB)
        {
            if (traceA == traceB)
            {
                return false;
            }
        }
    }
    return true;
}


IOTreeContainer Fsm::getAdaptiveRCharacterisationSet() const
{
    IOTreeContainer result = IOTreeContainer(presentationLayer);
    for (shared_ptr<FsmNode> n : nodes)
    {
        IOTreeContainer container = getAdaptiveRStateCharacterisationSet(n);
        for (auto tree : *container.getList())
        {
            result.addUniqueRemovePrefixes(tree);
        }
    }
    return result;
}

vector<vector<shared_ptr<FsmNode>>> Fsm::getMaximalSetsOfRDistinguishableStates() const
{
    vector<vector<shared_ptr<FsmNode>>> result;
    for (shared_ptr<FsmNode> node : nodes)
    {
        bool skip = false;
        for (auto set : result)
        {
            for (auto n : set)
            {
                if (n == node)
                {
                    skip = true;
                    break;
                }
            }
            if (skip)
            {
                break;
            }
        }
        if (skip)
        {
            continue;
        }
        vector<shared_ptr<FsmNode>> set = {node};
        for (shared_ptr<FsmNode> n : nodes)
        {
            if (node == n)
            {
                continue;
            }
            if (n->getRDistinguishability()->isRDistinguishableWith(set))
            {
                set.push_back(n);
            }
        }
        result.push_back(set);
    }

    return result;
}

void Fsm::calcStateIdentificationSets()
{
    if (!isObservable())
    {
        cout << "This FSM is not observable - cannot calculate the charactersiation set." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (characterisationSet == nullptr)
    {
        cout << "Missing characterisation set - exit." << endl;
        exit(EXIT_FAILURE);
    }
    
    /*Create empty state identification sets for every FSM state*/
    stateIdentificationSets.clear();
    
    /*Identify W by integers 0..m*/
    IOListContainer wIC = characterisationSet->getIOLists();
    shared_ptr<vector<vector<int>>> wLst = wIC.getIOLists();
    
    /*wLst.get(0) is identified with Integer(0),
     wLst.get(1) is identified with Integer(1), ...*/
    
    vector<vector<unordered_set<int>>> z;
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        z.push_back(vector<unordered_set<int>>());
        for (unsigned int j = 0; j < nodes.size(); ++ j)
        {
            z.at(i).push_back(unordered_set<int>());
        }
    }
    
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        shared_ptr<FsmNode> iNode = nodes.at(i);
        
        for (unsigned int j = i + 1; j < nodes.size(); ++ j)
        {
            shared_ptr<FsmNode> jNode = nodes.at(j);
            
            for (unsigned int u = 0; u < wLst->size(); ++ u)
            {
                vector<int> thisTrace = wLst->at(u);
                
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
        vector<unordered_set<int>> iLst;
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
        unordered_set<int> h = hs.calcMinCardHittingSet();
        
        shared_ptr<Tree> iTree = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
        for (int u : h)
        {
            vector<int> lli = wLst->at(u);
            shared_ptr<vector<vector<int>>> lllli = make_shared<vector<vector<int>>>();
            lllli->push_back(lli);
            iTree->addToRoot(IOListContainer(lllli, presentationLayer));
        }
        stateIdentificationSets.push_back(iTree);
        
    }
    
#if 0
    for (unsigned int n = 0; n < stateIdentificationSets.size(); ++ n)
    {
        cout << "W(" << n << ") = " << stateIdentificationSets.at(n)->getTestCases() << endl;
    }
#endif
    
}



void Fsm::calcStateIdentificationSetsFast()
{
    if (!isObservable())
    {
        cout << "This FSM is not observable - cannot calculate the charactersiation set." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (characterisationSet == nullptr)
    {
        cout << "Missing characterisation set - exit." << endl;
        exit(EXIT_FAILURE);
    }
    
    /*Create empty state identification sets for every FSM state*/
    stateIdentificationSets.clear();
    
    /*Identify W by integers 0..m*/
    IOListContainer wIC = characterisationSet->getIOLists();
    shared_ptr<vector<vector<int>>> wLst = wIC.getIOLists();
    
    // Matrix indexed over nodes
    vector< vector<int> > distinguish;
    
    // Every node is associated with an IOListContainer
    // containing its distinguishing traces
    vector< shared_ptr<IOListContainer> > node2iolc;
    
    for (size_t i = 0; i < size(); ++ i) {
        node2iolc.push_back(make_shared<IOListContainer>(presentationLayer));
        vector<int> v;
        distinguish.push_back(v);
        for (size_t j = 0; j < size(); j++ ) {
            distinguish.at(i).push_back(-1);
        }
    }
    
    for (size_t i = 0; i < size(); ++ i) {
        
        int traceIdx = 0;
        for (auto trc : *wLst) {
            
            bool complete = true;
            for ( size_t j = i+1; j < size(); j++ ) {
                if ( distinguish.at(i).at(j) == -1 ) {
                    
                    if (nodes[i]->distinguished(nodes[j], trc))
                    {
                        distinguish.at(i).at(j) = traceIdx;
                        distinguish.at(j).at(i) = traceIdx;
                        Trace tr(trc,presentationLayer);
                        node2iolc.at(i)->add(tr);
                        node2iolc.at(j)->add(tr);
                    }
                    else {
                        complete = false;
                    }
                    
                }
            }
            
            if ( complete ) break;
            traceIdx++;
        
        }
    }
    
    
    for (size_t i = 0; i < size(); ++ i) {
        shared_ptr<Tree> iTree = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
        iTree->addToRoot(*node2iolc.at(i));
        stateIdentificationSets.push_back(iTree);
    }
    
#if 0
    for (unsigned int n = 0; n < stateIdentificationSets.size(); ++ n)
    {
        cout << "W(" << n << ") = " << stateIdentificationSets.at(n)->getTestCases() << endl;
    }
#endif
    
}


void Fsm::appendStateIdentificationSets(const shared_ptr<Tree> Wp2) const
{
    IOListContainer cnt = Wp2->getIOLists();
    
    for (vector<int> lli : *cnt.getIOLists())
    {
        InputTrace itrc = InputTrace(lli, presentationLayer);
        
        /*Which are the target nodes reachable via input trace lli
         in this FSM?*/
        unordered_set<shared_ptr<FsmNode>> tgtNodes = getInitialState()->after(itrc);
        
        for (shared_ptr<FsmNode> n : tgtNodes)
        {
            int nodeId = n->getId();
            
            /*Get state identification set associated with n*/
            shared_ptr<Tree> wNodeId = stateIdentificationSets.at(nodeId);
            
            /*Append state identification set to Wp2 tree node
             reached after applying  itrc*/
            Wp2->addAfter(itrc, wNodeId->getIOLists());
        }
    }
}


IOListContainer Fsm::wMethod(const unsigned int numAddStates) {
    
    Fsm fo = transformToObservableFSM();
    Fsm fom = fo.minimise();
    
    return fom.wMethodOnMinimisedFsm(numAddStates);
}


IOListContainer Fsm::wMethodOnMinimisedFsm(const unsigned int numAddStates) {
    
    shared_ptr<Tree> iTree = getTransitionCover();
    
    if ( numAddStates > 0 ) {
        IOListContainer inputEnum = IOListContainer(maxInput,
                                                    1,
                                                    (int)numAddStates,
                                                    presentationLayer);
        iTree->add(inputEnum);
    }
    
    
    IOListContainer w = getCharacterisationSet();
    iTree->add(w);
    
    return iTree->getIOLists();
    
}

IOListContainer Fsm::wpMethod(const unsigned int numAddStates)
{
    
    shared_ptr<Tree> scov = getStateCover();
    
    shared_ptr<Tree> tcov = getTransitionCover();
    
    tcov->remove(scov);
    shared_ptr<Tree> r = tcov;
    
    IOListContainer w = getCharacterisationSet();
        
    calcStateIdentificationSetsFast();
    
    shared_ptr<Tree> Wp1 = scov;
    if (numAddStates > 0)
    {
        IOListContainer inputEnum = IOListContainer(maxInput, 1,
                                                    (int)numAddStates,
                                                    presentationLayer);
        
        Wp1->add(inputEnum);
    }
    Wp1->add(w);
    
    shared_ptr<Tree> Wp2 = r;
    if (numAddStates > 0)
    {
        IOListContainer inputEnum = IOListContainer(maxInput,
                                                    (int)numAddStates,
                                                    (int)numAddStates,
                                                    presentationLayer);
        
        Wp2->add(inputEnum);
    }
    appendStateIdentificationSets(Wp2);

    Wp1->unionTree(Wp2);
    return Wp1->getIOLists();
}


IOListContainer Fsm::hsiMethod(const unsigned int numAddStates)
{

    if (!isObservable())
    {
        cout << "This FSM is not observable - cannot calculate the harmonized state identification set." << endl;
        exit(EXIT_FAILURE);
    }
    
    IOListContainer wSet = getCharacterisationSet();

    shared_ptr<Tree> scov = getStateCover();

    /* V.(Inputs from length 1 to m-n+1) */
    shared_ptr<Tree> hsi = scov;
    IOListContainer inputEnum = IOListContainer(maxInput,
                                                    1,
                                                    (int)numAddStates + 1,
                                                    presentationLayer);
    hsi->add(inputEnum);

    /* initialize HWi trees */
    std::vector<shared_ptr<Tree>> hwiTrees;
    for (unsigned i = 0; i < nodes.size(); i++)
    {
        shared_ptr<TreeNode> root = make_shared<TreeNode>();
        shared_ptr<Tree> emptyTree = make_shared<Tree>(root, presentationLayer);
        hwiTrees.push_back(emptyTree);
    }

    /* Create harmonised state identification set for every FSM state.
     * For each pair of nodes i and j get one element
     * of the characterisation set that distinguishes the two nodes.
     * Add the distinguishing sequence to both HWi and HWj.
     */
    for (unsigned i = 0; i < nodes.size()-1; i++)
    {
        shared_ptr<FsmNode> node1 = nodes[i];
        for (unsigned j = i+1; j < nodes.size(); j++)
        {
            shared_ptr<FsmNode>node2 = nodes[j];
            bool distinguished = false;
            for (auto iolst : *wSet.getIOLists())
            {
                if (node1->distinguished(node2, iolst)){
                    distinguished = true;
                    hwiTrees[i]->addToRoot(iolst);
                    hwiTrees[j]->addToRoot(iolst);
                    break;
                }
            }
            if (!distinguished) {
                cout << "[ERR] Found inconsistency when applying HSI-Method: FSM not minimal." << endl;
            }
        }
    }

    /* Append harmonised state identification sets */
    IOListContainer cnt = hsi->getIOLists();
    for (auto lli : *cnt.getIOLists())
    {
        InputTrace itrc = InputTrace(lli, presentationLayer);

        /*Which are the target nodes reachable via input trace lli
         in this FSM?*/
        unordered_set<shared_ptr<FsmNode>> tgtNodes = getInitialState()->after(itrc);

        for (auto n : tgtNodes)
        {
            int nodeId = n->getId();

            /* harmonised state identification set associated with n*/
            shared_ptr<Tree> hwNodeId = hwiTrees[nodeId];

            /* Append harmonised state identification set to hsi tree node
               reached after applying itrc */
            hsi->addAfter(itrc,hwNodeId->getIOLists());
        }
    }

    return hsi->getIOLists();
}

TestSuite Fsm::createTestSuite(const IOListContainer & testCases)
{
    shared_ptr<vector<vector<int>>> tcLst = testCases.getIOLists();
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
    for (shared_ptr<FsmNode> nn : nodes)
    {
        for (int x = 0; x <= maxInput; ++ x)
        {
            bool found = false;
            for (auto tr : nn->getTransitions())
            {
                if (tr->getLabel()->getInput() == x)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                cout << "Incomplete FSM : for state " << nn->getName() << " " << nn->getId() << ", input " << x << " does not have a transition." << endl;
                cDefd = false;
            }
        }
    }
    return cDefd;
}

bool Fsm::isDeterministic() const
{
    for (shared_ptr<FsmNode> node : nodes)
    {
        if (!node->isDeterministic())
        {
            return false;
        }
    }
    return true;
}

void Fsm::setPresentationLayer(const shared_ptr<FsmPresentationLayer> ppresentationLayer)
{
    presentationLayer = ppresentationLayer;
}

ostream & operator<<(ostream & out, const Fsm & fsm)
{
    out << "digraph g {" << endl << endl << "node [shape = circle]" << endl << endl;
    for (int i = 0; i < static_cast<int> (fsm.nodes.size()); ++ i)
    {
        if (i == fsm.initStateIdx)
        {
            out << endl << "node [shape = doublecircle]" << endl;
        }
        
        if (fsm.nodes.at(i) == nullptr)
        {
            continue;
        }
        string nodeName = (fsm.nodes.at(i)->getName().empty()) ? "s" : fsm.nodes.at(i)->getName();
        out << i << "[label=\"" << nodeName << "(" << fsm.nodes.at(i)->getId() << ")\"];" << endl;
        
        if (i == fsm.initStateIdx)
        {
            out << endl << "node [shape = ellipse]" << endl;
        }
    }
    
    for (shared_ptr<FsmNode> node : fsm.nodes)
    {
        if (node != nullptr)
        {
            out << *node;
        }
    }
    out << endl << "}" << endl;
    return out;
}


unsigned int Fsm::getRandomSeed() {
    
    return static_cast<unsigned int>
    (high_resolution_clock::now().time_since_epoch().count());
    
}


shared_ptr<Fsm>
Fsm::createRandomFsm(const string & fsmName,
                     const int maxInput,
                     const int maxOutput,
                     const int maxState,
                     const shared_ptr<FsmPresentationLayer> pl,
                     const unsigned seed) {
    
    // Initialisation of random number generation
    if ( seed == 0 ) {
        srand(getRandomSeed());
    }
    else {
        srand(seed);
    }
    
    // Produce the nodes and put them into a vector.
    // All nodes are marked 'white' by the costructor - this is now
    // used to mark unreachable states which have to be made reachable
    vector<shared_ptr<FsmNode> > lst;
    for ( int n = 0; n <= maxState; n++ ) {
        lst.push_back(make_shared<FsmNode>(n,fsmName,pl));
    }
    
    // At index 0 of the vector, the initial state is store, and
    // this is reachable
    lst[0]->setColor(FsmNode::black);
    
    // We create transitions by starting from black (reachable) nodes
    // and trying to reach at least one white node from there.
    deque< shared_ptr<FsmNode> > bfsq;
    bfsq.push_back(lst[0]);
    
    while ( not bfsq.empty() ) {
        
        shared_ptr<FsmNode> srcNode = bfsq.front();
        bfsq.pop_front();
        
        // Generation part 1.
        // Select an uncovered node at random
        int whiteNodeIndex = rand() % (maxState+1);
        shared_ptr<FsmNode> whiteNode = nullptr;
        shared_ptr<FsmNode> startNode = lst[whiteNodeIndex];
        shared_ptr<FsmNode> thisNode = startNode;
        
        do {
            
            if ( thisNode->getColor() == FsmNode::white ) {
                whiteNode = thisNode;
            }
            else {
                whiteNodeIndex = (whiteNodeIndex + 1) % (maxState+1);
                thisNode = lst[whiteNodeIndex];
            }
            
        } while ( whiteNode == nullptr and thisNode != startNode );
        
        // Link srcNode by random transition to thisNode
        // and mark thisNode as black. Also insert into BFS queue
        int x0 = -1;
        int y0 = -1;
        
        if ( whiteNode != nullptr ) {
            x0 = rand() % (maxInput+1);
            y0 = rand() % (maxOutput+1);
            auto theTrans =
            make_shared<FsmTransition>(srcNode,whiteNode,
                                       make_shared<FsmLabel>(x0,y0,pl));
            // Add transition to adjacency list of the source node
            srcNode->addTransition(theTrans);
            thisNode->setColor(FsmNode::black);
            bfsq.push_back(thisNode);
        }
        
        // Generation part 2.
        // Random transition generation.
        // At least one transition for every input, with
        // arbitrary target nodes.
        for ( int x = 0; x <= maxInput; x++ ) {
            // If x equals x0 produced already above,
            // we may skip it at random
            if ( x == x0 and (rand() % 2) ) continue;
            
            // How many transitions do we want for input x?
            // We construct at most 2 of these transitions
            int numTrans = rand() % 2;
            for ( int t = 0; t <= numTrans; t++ ) {
                // Which output do we want?
                int y = rand() % (maxOutput+1);
                // Which target node?
                int tgtNodeId = rand() % (maxState+1);
                auto tgtNode = lst[tgtNodeId];
                if ( tgtNode->getColor() == FsmNode::white ) {
                    tgtNode->setColor(FsmNode::black);
                    bfsq.push_back(tgtNode);
                }
                auto theTrans =
                make_shared<FsmTransition>(srcNode,tgtNode,
                                           make_shared<FsmLabel>(x,y,pl));
                // Add transition to adjacency list of the source node
                srcNode->addTransition(theTrans);
            }
            
        }
        
    }
    
    return make_shared<Fsm>(fsmName,maxInput,maxOutput,lst,pl);
    
}

shared_ptr<Fsm> Fsm::createProductMachine(shared_ptr<Fsm> reference, shared_ptr<Fsm> iut, const string& fsmName)
{
    shared_ptr<FsmPresentationLayer> pl = FsmPresentationLayer::mergeAlphabets(reference->getPresentationLayer(), iut->getPresentationLayer());
    const size_t failOutput = static_cast<size_t>(pl->addOut2String("fail"));
    size_t maxInput = pl->getIn2String().size() - 1;
    size_t maxOutput = pl->getOut2String().size() - 1;
    vector<string> state2String;
    vector<shared_ptr<FsmNode>> nodes;
    map<shared_ptr<FsmNode>, pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> productNodesToOriginalNodes;
    map<pair<int, int>, shared_ptr<FsmNode>> originalNodeIdsToProductNode;
    int productNodeId = 0;
    for(shared_ptr<FsmNode> nodeA : reference->getNodes())
    {
        for(shared_ptr<FsmNode> nodeB : iut->getNodes())
        {
            shared_ptr<FsmNode> productNode = make_shared<FsmNode>(productNodeId, fsmName, pl);
            productNodesToOriginalNodes.insert(make_pair(productNode, make_pair(nodeA, nodeB)));
            originalNodeIdsToProductNode.insert(make_pair(make_pair(nodeA->getId(), nodeB->getId()), productNode));
            string nodeName = "(" + nodeA->getName() + "," + nodeB->getName() + ")";
            state2String.push_back(nodeName);
            nodes.push_back(productNode);
            productNodeId++;
        }
    }

    shared_ptr<FsmNode> failState = make_shared<FsmNode>(productNodeId, fsmName, pl);

    for (size_t x = 0; x <= maxInput; ++x)
    {
        shared_ptr<FsmLabel> failLabel = make_shared<FsmLabel>(x, failOutput, pl);
        shared_ptr<FsmTransition> failTransition = make_shared<FsmTransition>(
                    failState,
                    failState,
                    failLabel);
        failState->addTransition(failTransition);
    }

    for(auto it = productNodesToOriginalNodes.begin(); it != productNodesToOriginalNodes.end(); ++it)
    {
        shared_ptr<FsmNode> productNode = it->first;
        shared_ptr<FsmNode> nodeA = it->second.first;
        shared_ptr<FsmNode> nodeB = it->second.second;
        for (size_t x = 0; x <= maxInput; ++x)
        {
            for (size_t y = 0; y <= maxOutput; ++y)
            {
                if (y == failOutput)
                {
                    continue;
                }
                shared_ptr<FsmLabel> productLabel = make_shared<FsmLabel>(x, y, pl);
                unordered_set<shared_ptr<FsmNode>> targetsA = nodeA->afterAsSet(static_cast<int>(x), static_cast<int>(y));
                unordered_set<shared_ptr<FsmNode>> targetsB = nodeB->afterAsSet(static_cast<int>(x), static_cast<int>(y));
                if (targetsB.size() > 0 && targetsA.size() == 0)
                {
                    // The transition with input x and output y is defined for the IUT, but it isn't for the reference.
                    // Therefore, the product machine will lead to the fail state.
                    shared_ptr<FsmTransition> productTransition = make_shared<FsmTransition>(
                                productNode,
                                failState,
                                productLabel);
                    productNode->addTransition(productTransition);
                    continue;
                }
                for (shared_ptr<FsmNode> targetA : targetsA)
                {
                    for (shared_ptr<FsmNode> targetB : targetsB)
                    {
                        shared_ptr<FsmNode> productTarget = originalNodeIdsToProductNode.at(make_pair(targetA->getId(), targetB->getId()));
                        shared_ptr<FsmTransition> productTransition = make_shared<FsmTransition>(
                                    productNode,
                                    productTarget,
                                    productLabel);
                        productNode->addTransition(productTransition);
                        //shared_ptr<FsmNode> originInProduct =
                    }
                }
            }
        }
    }

    int initStateIdx = originalNodeIdsToProductNode.at(
                make_pair(reference->getInitialState()->getId(),iut->getInitialState()->getId()))->getId();

    state2String.push_back("Fail");
    nodes.push_back(failState);
    pl->setState2String(state2String);
    shared_ptr<Fsm> result = make_shared<Fsm>(fsmName, maxInput, maxOutput, nodes, initStateIdx, pl);
    result->failState = failState;
    return result;
}



shared_ptr<Fsm> Fsm::createMutant(const std::string & fsmName,
                                  const size_t numOutputFaults,
                                  const size_t numTransitionFaults){
    
    srand(getRandomSeed());
    
    // Create new nodes for the mutant.
    vector<shared_ptr<FsmNode> > lst;
    for ( int n = 0; n <= maxState; n++ ) {
        lst.push_back(make_shared<FsmNode>(n,fsmName,presentationLayer));
    }
    
    // Now add transitions that correspond exactly to the transitions in
    // this FSM
    for ( int n = 0; n <= maxState; n++ ) {
        auto theNewFsmNodeSrc = lst[n];
        auto theOldFsmNodeSrc = nodes[n];
        for ( auto tr : theOldFsmNodeSrc->getTransitions() ) {
            int tgtId = tr->getTarget()->getId();
            auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
            shared_ptr<FsmTransition> newTr =
            make_shared<FsmTransition>(theNewFsmNodeSrc,lst[tgtId],newLbl);
            theNewFsmNodeSrc->addTransition(newTr);
        }
    }
    
    // Now add transition faults to the new machine
    for ( size_t tf = 0; tf < numTransitionFaults; tf++ ) {
        int srcNodeId = rand() % (maxState+1);
        int newTgtNodeId = rand() % (maxState+1);
        int trNo = rand() % lst[srcNodeId]->getTransitions().size();
        auto tr = lst[srcNodeId]->getTransitions()[trNo];
        if ( tr->getTarget()->getId() == newTgtNodeId ) {
            newTgtNodeId = (newTgtNodeId+1) % (maxState+1);
        }
        lst[srcNodeId]->getTransitions()[trNo]->setTarget(lst[newTgtNodeId]);
    }
    
    // Now add output faults to the new machine
    for (size_t of = 0; of < numOutputFaults; of++ ) {
        
        int srcNodeId = rand() % (maxState+1);
        int trNo = rand() % lst[srcNodeId]->getTransitions().size();
        auto tr = lst[srcNodeId]->getTransitions()[trNo];
        int theInput = tr->getLabel()->getInput();
        int newOutVal = rand() % (maxOutput+1);
        int originalNewOutVal = rand() % (maxOutput+1);
        bool newOutValOk;
        
        // We don't want to modify this transition in such a way
        // that another one with the same label and the same
        // source/target nodes already exists.
        do {
            
            newOutValOk = true;
            
            for ( auto trOther : lst[srcNodeId]->getTransitions() ) {
                if ( tr == trOther ) continue;
                if ( trOther->getTarget()->getId() != tr->getTarget()->getId() )
                    continue;
                if ( trOther->getLabel()->getInput() != theInput ) continue;
                if ( trOther->getLabel()->getOutput() == newOutVal ) {
                    newOutValOk = false;
                }
            }
            
            if ( not newOutValOk ) {
                newOutVal = (newOutVal+1) % (maxOutput+1);
            }
            
        } while ( (not newOutValOk) and (originalNewOutVal != newOutVal) );
        
        if ( newOutValOk ) {
            
            auto newLbl = make_shared<FsmLabel>(tr->getLabel()->getInput(),
                                                
                                                newOutVal,
                                                presentationLayer);
            
            tr->setLabel(newLbl);
        }
    }
    
    return make_shared<Fsm>(fsmName,maxInput,maxOutput,lst,presentationLayer);
    
}



std::vector< std::unordered_set<int> >
Fsm::getEquivalentInputsFromPrimeMachine() {
    
    vector< std::unordered_set<int> > v;
    
    shared_ptr<OFSMTable> ot =
        make_shared<OFSMTable>(nodes,
                               maxInput,
                               maxOutput,
                               presentationLayer);
    
    // mark all inputs as non-equivalent
    vector<bool> equivalentToSmallerInput;
    for ( int x = 0; x <= maxInput; x++ )
        equivalentToSmallerInput.push_back(false);
    
    // Check inputs for equivalence
    for ( int x1 = 0; x1 <= maxInput; x1++ ) {
        
        if ( equivalentToSmallerInput[x1] ) continue;
        
        unordered_set<int> classOfX1;
        classOfX1.insert(x1);
        
        for ( int x2 = x1 + 1; x2 <= maxInput; x2++ ) {
            
            bool x2EquivX1 = true;
            
            // To check whether x1 is equivalent to x2,
            // loop over all outputs y and compare OFSM
            // table columns x1/y and x2/y
            for ( int y = 0; y <= maxOutput; y++ ) {
                if ( not ot->compareColumns(x1,y,x2,y) ) {
                    x2EquivX1 = false;
                    break;
                }
            }
            
            if ( x2EquivX1 ) {
                equivalentToSmallerInput[x2] = true;
                classOfX1.insert(x2);
            }
            
        }
        
        v.push_back(classOfX1);

    }
    
    return v;
    
}

std::vector< std::unordered_set<int> > Fsm::getEquivalentInputs() {
    
    if ( minimal != True ) {
        return minimise().getEquivalentInputsFromPrimeMachine();
    }
    else {
        return getEquivalentInputsFromPrimeMachine();
    }
    
}


void Fsm::accept(FsmVisitor& v) {
    
    deque< shared_ptr<FsmNode> > bfsq;
    
    resetColor();
    
    v.visit(*this);
    
    bfsq.push_back(nodes[initStateIdx]);
    
    while ( not bfsq.empty() ) {
        shared_ptr<FsmNode> theNode = bfsq.front();
        bfsq.pop_front();
        v.setNew(true);
        theNode->accept(v,bfsq);
    }
    
    
}



bool Fsm::removeUnreachableNodes(std::vector<shared_ptr<FsmNode>>& unreachableNodes) {
    
    vector<shared_ptr<FsmNode>> newNodes;
    FsmVisitor v;
    
    // Mark all reachable nodes as 'visited'
    accept(v);
    
    // When removing nodes from the FSM, the node ids of all remaining nodes
    // have to be adapted, in order to match the index in the list of nodes.
    // This is necessary, because during minimisation with OFSM tables or
    // Pk-tables, the algorithms rely on the range of row numbers being
    // identical to the range of node ids of the reachable nodes.
    int subtractFromId = 0;
    for ( auto n : nodes ) {
        if ( not n->hasBeenVisited() ) {
            unreachableNodes.push_back(n);
            presentationLayer->removeState2String(n->getId() - subtractFromId);
            ++subtractFromId;
        }
        else {
            n->setId(n->getId() - subtractFromId);
            newNodes.push_back(n);
        }
    }
    
    nodes = newNodes;
    
    return (unreachableNodes.size() > 0);
}


