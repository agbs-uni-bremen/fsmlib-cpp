/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <chrono>
#include <deque>
#include <algorithm>
#include <regex>

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
#include "logging/easylogging++.h"
#include "logging/Logging.h"

using namespace std;
using namespace std::chrono;

shared_ptr<FsmNode> Fsm::newNode(const int id, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p, const shared_ptr<FsmPresentationLayer>& pl)
{
    string nodeName = string("(" + p->first->getName() + to_string(p->first->getId()) + ","
                             + p->second->getName() + to_string(p->second->getId()) + ")");
    shared_ptr<FsmNode> n = make_shared<FsmNode>(id, nodeName, pl);
    n->setPair(p);
    return n;
}

bool Fsm::contains(const vector<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>>& lst, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p)
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

bool Fsm::contains(const vector<shared_ptr<FsmNode>>& lst, const shared_ptr<FsmNode>& n)
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

shared_ptr<FsmNode> Fsm::findp(const vector<shared_ptr<FsmNode>>& lst, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p)
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
        LOG(FATAL) << "Unable to open input file";
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
        LOG(FATAL) << "Unable to open input file";
    }
    
}

void Fsm::readFsmFromDot (const string & fname)
{

    maxInput = 0;
    maxOutput = 0;
    maxState = 0;
    initStateIdx = -1;
    int nodeIdCount = 0;
    map<int,shared_ptr<FsmNode>> existingNodes;
    initStateIdx = -1;
    presentationLayer = make_shared<FsmPresentationLayer>();

    ifstream inputFile (fname);
    if (inputFile.is_open())
    {
        string line;
        regex regInitial("\\s*node\\s*\\[\\s*shape\\s*=\\s*doublecircle\\s*\\]\\s*");
        regex regNode("\\s*(\\d)\\s*\\[\\s*label=\"(.*)\"\\s*\\]\\s*;");
        bool nextIsInitial = false;
        bool initialSet = false;
        while (getline (inputFile, line))
        {
            if (!initialSet && !nextIsInitial && regex_match(line, regInitial))
            {
                nextIsInitial = true;
                continue;
            }
            cmatch matches;
            regex_match(line.c_str(), matches, regNode);
            LOG(INFO) << "line:";
            LOG(INFO) << "  " << line;
            LOG(INFO) << "matches:";
            for (unsigned i=0; i<matches.size(); ++i) {
                LOG(INFO) << "  " << matches[i];
            }
            if (matches.size() == 3)
            {
                int nodeId = stoi(matches[1]);
                string nodeName = matches[2];
                if(existingNodes.find(nodeId) != existingNodes.end()) {
                    LOG(FATAL) << "Error while parsing dot file. The node id " << nodeId << "has been assigned more than once.";
                }
                presentationLayer->addState2String(nodeName);
                shared_ptr<FsmNode> node = make_shared<FsmNode>(nodeIdCount++, presentationLayer);
                nodes.push_back(node);
                existingNodes.insert(make_pair(nodeId, node));
                maxState++;
                if (nextIsInitial)
                {
                    initStateIdx = node->getId();
                    node->markAsInitial();
                    nextIsInitial = false;
                    initialSet = true;
                }
            }
        }
        inputFile.close ();
    }

    inputFile.open(fname);
    if (inputFile.is_open())
    {
        string line;
        regex regTransition("\\s*(\\d)\\s*->\\s*(\\d)\\s*\\[\\s*label=\"(.+)/(.+)\"\\s*\\]\\s*;");
        int inputCount = 0;
        int outputCount = 0;
        map<string, int> inputs;
        map<string, int> outputs;
        while (getline (inputFile, line))
        {
            cmatch matches;
            regex_match(line.c_str(), matches, regTransition);
            if (matches.size() == 5)
            {
                LOG(INFO) << matches[1] << " -" << matches[3] << "/" << matches[4] << "-> " << matches[2];
                int sourceId = stoi(matches[1]);
                int targetId = stoi(matches[2]);
                string input = matches[3];
                string output = matches[4];
                int in;
                int out;

                if(inputs.find(input) == inputs.end()) {
                    maxInput++;
                    in = inputCount++;
                    presentationLayer->addIn2String(input);
                    inputs.insert(make_pair(input, in));
                }
                else
                {
                    in = inputs.at(input);
                }
                if(outputs.find(output) == outputs.end()) {
                    maxOutput++;
                    out = outputCount++;
                    presentationLayer->addOut2String(output);
                    outputs.insert(make_pair(output, out));
                }
                else
                {
                    out = outputs.at(output);
                }

                shared_ptr<FsmLabel> label = make_shared<FsmLabel>(in, out, presentationLayer);
                shared_ptr<FsmTransition> trans = make_shared<FsmTransition>(existingNodes.at(sourceId), existingNodes.at(targetId), label);
                existingNodes.at(sourceId)->addTransition(trans);
            }
        }
    }
    else
    {
        LOG(FATAL) << "Unable to open input file";
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

Fsm::Fsm(const Fsm& other):
name(other.name), currentParsedNode(nullptr), characterisationSet(nullptr),
presentationLayer(other.presentationLayer) {
    
    maxInput = other.maxInput;
    maxOutput = other.maxOutput;
    maxState = other.maxState;
    failOutput = other.failOutput;
    initStateIdx = other.initStateIdx;
    minimal = other.minimal;
    
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

Fsm::Fsm(const shared_ptr<FsmPresentationLayer>& presentationLayer)
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
         const shared_ptr<FsmPresentationLayer>& presentationLayer,
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

Fsm::Fsm(const std::string& dotFileName,
    const std::string& fsmName)
{
    readFsmFromDot(dotFileName);
}

Fsm::Fsm(const string & fname,
         const string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const shared_ptr<FsmPresentationLayer>& presentationLayer)
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
         const vector<shared_ptr<FsmNode>>& lst,
         const shared_ptr<FsmPresentationLayer>& presentationLayer)
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
         const vector<shared_ptr<FsmNode>>& lst,
         const int initStateIdx,
         const shared_ptr<FsmPresentationLayer>& presentationLayer):
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

vector<shared_ptr<FsmNode>> Fsm::getDReachableStates(vector<std::shared_ptr<InputTrace>>& detStateCover)
{
    TIMED_FUNC(timerObj);
    VLOG(2) << "getDReachableStates()";
    resetColor();
    deque<shared_ptr<FsmNode>> bfsLst;
    vector<shared_ptr<FsmNode>> nodes;
    map<shared_ptr<FsmNode>, shared_ptr<IOTrace>> paths;

    shared_ptr<FsmNode> initState = getInitialState();
    initState->setColor(FsmNode::grey);
    bfsLst.push_back(initState);
    nodes.push_back(initState);
    shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
    detStateCover.push_back(make_shared<InputTrace>(FsmLabel::EPSILON, presentationLayer));
    emptyTrace->setTargetNode(initState);
    initState->setDReachable(emptyTrace);
    paths.insert(make_pair(initState, IOTrace::getEmptyTrace(presentationLayer)));

    while (!bfsLst.empty())
    {
        shared_ptr<FsmNode> thisNode = bfsLst.front();
        bfsLst.pop_front();
        VLOG(2) << "thisNode: " << thisNode->getName();

        shared_ptr<IOTrace> thisNodePath;

        if (!thisNode->isInitial())
        {
            try
            {
                thisNodePath = paths.at(thisNode);
                VLOG(2) << "thisNodePath: " << *thisNodePath;
            }
            catch (out_of_range e)
            {
                // DO nothing.
            }
        }

        for (int x = 0; x <= maxInput; ++x)
        {
            VLOG(2) << "x: " << presentationLayer->getInId(x);
            vector<int> producedOutputs;
            vector<shared_ptr<FsmNode>> successorNodes = thisNode->after(x, producedOutputs);
            VLOG(2) << "successorNodes:";
            for (auto n : successorNodes)
            {
                VLOG(2) << "  " << n->getName();
            }
            VLOG(2) << "producedOutputs:";
            for (auto n : producedOutputs)
            {
                VLOG(2) << "  " << presentationLayer->getOutId(n);
            }
            if (successorNodes.size() > 1)
            {
                bool skip = false;
                const shared_ptr<FsmNode>& n = successorNodes.at(0);
                for (size_t i = 1; i < successorNodes.size(); ++i)
                {
                    const shared_ptr<FsmNode>& other  = successorNodes.at(i);
                    if (n != other)
                    {
                        VLOG(2) << "Skipping.";
                        skip = true;
                        break;
                    }
                }
                if (skip)
                {
                    continue;
                }
            }
            shared_ptr<FsmNode> tgt = successorNodes.at(0);
            VLOG(2) << "tgt:" << tgt->getName();
            try
            {
                paths.at(tgt);
                // Path already exists. Do nothing.
                VLOG(2) << "Path already exists. Do nothing.";
            }
            catch (out_of_range e)
            {
                // Create new path, since it doesn't exist.
                shared_ptr<IOTrace> newPath;
                if (thisNodePath)
                {
                    newPath = make_shared<IOTrace>(*thisNodePath);
                    newPath->append(x, producedOutputs.at(0));
                    VLOG(2) << "newPath (appended): " << *newPath;
                }
                else
                {
                    InputTrace in = InputTrace({x}, presentationLayer);
                    OutputTrace out = OutputTrace({producedOutputs.at(0)}, presentationLayer);
                    newPath = make_shared<IOTrace>(in, out);
                    VLOG(2) << "newPath (new): " << *newPath;
                }
                newPath->setTargetNode(tgt);
                paths.insert(make_pair(tgt, newPath));
            }
            if (tgt->getColor() == FsmNode::white)
            {
                VLOG(2) << "Target color is white. Setting grey, adding node, setting d-reach path.";
                tgt->setColor(FsmNode::grey);
                bfsLst.push_back(tgt);
                nodes.push_back(tgt);
                detStateCover.push_back(make_shared<InputTrace>(paths.at(tgt)->getInputTrace()));
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

    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
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
            nSource = newNode(id ++, p, pl);
            string nodeName = string("(" + p->first->getName() + "," + p->second->getName() + ")");
            pl->addState2String(nodeName);
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
                        nTarget = newNode(id ++, pTarget, pl);
                        string nodeName = string("(" + pTarget->first->getName() + "," + pTarget->second->getName() + ")");
                        pl->addState2String(nodeName);
                        fsmInterNodes.push_back(nTarget);
                    }
                    
                    /*Add transition from nSource to nTarget*/
                    auto newTr = make_shared<FsmTransition>(nSource,
                                                            nTarget,
                                                            tr->getLabel());

                    pl->addIn2String(newTr->getLabel()->getInput(), presentationLayer->getInId(newTr->getLabel()->getInput()));
                    pl->addOut2String(newTr->getLabel()->getOutput(), presentationLayer->getOutId(newTr->getLabel()->getOutput()));

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

    return Fsm(f.getName(), maxInput, maxOutput, fsmInterNodes, pl);
}

shared_ptr<Tree> Fsm::getDeterministicStateCover()
{
    VLOG(2) << "getDeterministicStateCover()";
    TIMED_FUNC(timerObj);
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
        VLOG(2) << "thisNode: " << thisNode->getName();
		bfsLst.pop_front();
		shared_ptr<TreeNode> currentTreeNode = f2t[thisNode];

		for (int x = 0; x <= maxInput; ++x)
		{
            VLOG(2) << "x: " << x;
			vector<shared_ptr<FsmNode>> successorNodes = thisNode->after(x);
            VLOG(2) << "successorNodes: ";
            for (auto n : successorNodes)
            {
                VLOG(2) << "  " + n->getName();
            }
            //BUG It is possible to have the same successor node with different input
            // and different outputs.
			if (successorNodes.size() != 1)
			{
				continue;
			}
            shared_ptr<FsmNode>& tgt = successorNodes.at(0);
            if (tgt->getColor() == FsmNode::white)
            {
                tgt->setColor(FsmNode::grey);
                shared_ptr<TreeNode> itn = currentTreeNode->add(x);
                bfsLst.push_back(tgt);
                f2t[tgt] = itn;
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
    TIMED_FUNC(timerObj);
    return getInitialState()->getPossibleOutputs(input, producedOutputs, reachedNodes);
}

Fsm Fsm::transformToObservableFSM() const
{
    VLOG(1) << "transformToObservableFSM()";
    vector<shared_ptr<FsmNode>> nodeLst;
    vector<shared_ptr<FsmNode>> bfsLst;
    unordered_map<shared_ptr<FsmNode>, unordered_set<shared_ptr<FsmNode>>> node2Label;
    unordered_set<shared_ptr<FsmNode>> theLabel;
    
    theLabel.insert(getInitialState());

    vector<string> obsState2String;
    shared_ptr<FsmPresentationLayer> obsPl =
    make_shared<FsmPresentationLayer>(presentationLayer->getIn2String(),
                                      presentationLayer->getOut2String(),
                                      obsState2String);
    
    int id = 0;
    string nodeName = labelString(theLabel);
    shared_ptr<FsmNode> q0 = make_shared<FsmNode>(id ++, nodeName, obsPl);
    CVLOG(2, logging::fsmConversion) << "Initial state: " << q0->getName();
    nodeLst.push_back(q0);
    bfsLst.push_back(q0);
    node2Label[q0] = theLabel;
    obsPl->addState2String(nodeName);
    
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
                make_shared<FsmLabel>(x, y, obsPl);
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
                            CVLOG(4, logging::fsmConversion) << "Use existing node: " << tgtNode->getName() << ", id: " << tgtNode->getId();
                            break;
                        }
                    }
                    
                    /*We need to create a new node*/
                    if (tgtNode == nullptr)
                    {
                        nodeName = labelString(theLabel);
                        tgtNode = make_shared<FsmNode>(id ++, nodeName, obsPl);
                        CVLOG(2, logging::fsmConversion) << "New node: " << tgtNode->getName() << ", id: " << tgtNode->getId() << ", label: " << labelString(theLabel);
                        nodeLst.push_back(tgtNode);
                        bfsLst.push_back(tgtNode);
                        node2Label[tgtNode] = theLabel;
                        obsPl->addState2String(nodeName);
                    }
                    
                    /*Create the transition from q to tgtNode*/
                    auto trNew = make_shared<FsmTransition>(q, tgtNode, lbl);
                    q->addTransition(trNew);
                }
            }
        }
    }
    Fsm obsFsm(name + "_O", maxInput, maxOutput, nodeLst, obsPl);
    obsFsm.failOutput = failOutput;
    return obsFsm;
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
    fsm.failOutput = failOutput;
    return fsm;
}

Fsm Fsm::minimise()
{
    TIMED_FUNC(timerObj);
    vector<shared_ptr<FsmNode>> uNodes;
    removeUnreachableNodes(uNodes);
    
    if (!isObservable())
    {
        LOG(INFO) << "Fsm is not observable. Converting.";
        return transformToObservableFSM().minimiseObservableFSM();
    }
    
    return minimiseObservableFSM();
}

Fsm Fsm::makeComplete(CompleteMode mode)
{
    TIMED_FUNC(timerObj);
    CLOG(DEBUG, logging::fsmConversion) << "makeComplete():";
    vector<shared_ptr<FsmNode>> newNodes = nodes;
    bool addErrorState = false;
    shared_ptr<FsmNode> errorNode;
    if (mode == ErrorState)
    {
        errorNode = make_shared<FsmNode>(getMaxNodes() + 1, "Error", presentationLayer);
        presentationLayer->addState2String("Error");
        for (int x = 0; x <= maxInput; ++x)
        {
            shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::EPSILON, presentationLayer);
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(errorNode, errorNode, label);
            errorNode->addTransition(transition);
        }
    }
    for (shared_ptr<FsmNode> node : newNodes)
    {
        CVLOG(2, logging::fsmConversion) << "  State " << node->getName() << ": ";
        for (int x = 0; x <= maxInput; ++x)
        {
            std::stringstream ss;
            ss << "    Input " << presentationLayer->getInId(x);
            if(!node->isPossibleInput(x))
            {
                ss << " is not defined.";
                if (mode == ErrorState)
                {
                    addErrorState = true;
                    shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::ERROR_OUTPUT, presentationLayer);
                    shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(node, errorNode, label);
                    node->addTransition(transition);
                }
                else if (mode == SelfLoop)
                {
                    shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, FsmLabel::EPSILON, presentationLayer);
                    shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(node, node, label);
                    node->addTransition(transition);
                }
            }
            else
            {
                ss << " is defined.";
            }
            CVLOG(2, logging::fsmConversion) << ss.str();
        }
    }
    if (mode == ErrorState && addErrorState)
    {
        newNodes.push_back(errorNode);
    }
    Fsm fsmComplete(name + "_COMPLETE", maxInput, maxOutput, newNodes, presentationLayer);
    fsmComplete.complete = true;
    fsmComplete.failOutput = failOutput;
    return fsmComplete;
}

int Fsm::getFailOutput() const
{
    return failOutput;
}

bool Fsm::isCharSet(const shared_ptr<Tree>& w) const
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

void Fsm::minimiseCharSet(const shared_ptr<Tree>& w)
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
        LOG(FATAL) << "This FSM is not observable - cannot calculate the charactersiation set.";
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
        nodes.at(i)->getRDistinguishability()->initRDistinguishable(1);
        for (size_t j = i + 1; j < nodes.size(); ++j)
        {
            nodes.at(i)->getRDistinguishability()->addNotRDistinguishable(1, nodes.at(j));
        }
    }
    nodes.at(nodes.size() - 1)->getRDistinguishability()->addNotRDistinguishable(1);

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

                    q1->getRDistinguishability()->addRDistinguishable(1, q2);
                    q2->getRDistinguishability()->addRDistinguishable(1, q1);

                    q1->getRDistinguishability()->removeNotRDistinguishable(1, q2);
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
                std::stringstream ss;
                ss << "σ(" << nodes.at(i)->getName() << "," << nodes.at(j)->getName() << ") = " << *tree;
                VLOG(2) << ss.str();
            } catch (std::out_of_range e) {
               // Do nothing.
            }

        }

    }
}

void Fsm::calcRDistinguishableStates()
{
    TIMED_FUNC(timerObj);
    VLOG(2) << "calcRDistinguishableStates():";
    calcROneDistinguishableStates();

    size_t limit = nodes.size() * (nodes.size() - 1) / 2;
    bool allRDistinguishable = false;
    bool newDistinguishabilityCalculated = true;
    size_t maxL = 0;
    for (size_t l = 2; !allRDistinguishable && newDistinguishabilityCalculated && l <= limit; ++l)
    {
        maxL = l;
        VLOG(2) << "################ l = " << l << " (max " << limit << ") ################";
        allRDistinguishable = true;
        newDistinguishabilityCalculated = false;
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            nodes.at(k)->getRDistinguishability()->inheritDistinguishability(l);
        }
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            shared_ptr<FsmNode> q1 = nodes.at(k);
            VLOG(3) << "q1 = " << q1->getName() << ":";
            vector<shared_ptr<FsmNode>> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(l);
            for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
            {
                // There are still nodes that can not be r-distuinguisehd from each other. Do one more iteration.
                allRDistinguishable = false;
                shared_ptr<FsmNode> q2 = *it;
                VLOG(3) << "  q2 = " << q2->getName() << ":";
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
                            VLOG(3) << "    x = " << presentationLayer->getInId(x) << ":    "
                            << afterNode1->getName() << " != " << afterNode2->getName()
                            << "  ->  " << q1->getName() << " != " << q2->getName();

                            //shared_ptr<TreeNode> target1 = make_shared<TreeNode>();
                            shared_ptr<InputOutputTree> childTree1 = afterNode1->getRDistinguishability()->getAdaptiveIOSequence(afterNode2);
                            // TODO Fix
                            //      can't find linker symbol for virtual table for `TreeEdge' value
                            // messages when debugging.
                            // Put breakpoint at following line and debug.
                            stringstream ss;
                            ss << "      childIO1(" << afterNode1->getName() << "," << afterNode2->getName() << "): " << *childTree1;
                            VLOG(3) << ss.str();
                            ss.str(std::string());

                            shared_ptr<AdaptiveTreeNode> childNode1 = static_pointer_cast<AdaptiveTreeNode>(childTree1->getRoot());
                            shared_ptr<TreeEdge> edge1 = make_shared<TreeEdge>(y, childNode1);
                            q1Edges.push_back(edge1);

                            //shared_ptr<TreeNode> target2 = make_shared<TreeNode>();
                            shared_ptr<InputOutputTree> childTree2 = afterNode2->getRDistinguishability()->getAdaptiveIOSequence(afterNode1);
                            ss << "      childIO2(" << afterNode2->getName() << "," << afterNode1->getName() << "): " << *childTree2 << endl;
                            VLOG(3) << ss.str();
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

                        stringstream ss;
                        ss << "    q1Tree: " << *q1Tree;
                        VLOG(2) << ss.str();
                        ss.str(std::string());
                        ss << "    q2Tree: " << *q2Tree;
                        VLOG(2) << ss.str();

                        q1->getRDistinguishability()->addAdaptiveIOSequence(q2, q1Tree);
                        q2->getRDistinguishability()->addAdaptiveIOSequence(q1, q2Tree);


                        q1->getRDistinguishability()->addRDistinguishable(l, q2);
                        q2->getRDistinguishability()->addRDistinguishable(l, q1);
                        q1->getRDistinguishability()->removeNotRDistinguishable(l, q2);
                        q2->getRDistinguishability()->removeNotRDistinguishable(l, q1);
                        newDistinguishabilityCalculated = true;
                        break;
                    }
                }
            }
        }
    }
    // Deducing non-r-distinguishability from r-distinguishability.
    for (auto node : nodes)
    {
        for (size_t l = 1; l <= maxL; ++l)
        {
            node->getRDistinguishability()->addNotRDistinguishable(l);
            vector<shared_ptr<FsmNode>> dist = node->getRDistinguishability()->getRDistinguishableWith(l);
            for (auto n : nodes)
            {
                if(node != n && find(dist.begin(), dist.end(), n) == dist.end()) {
                    node->getRDistinguishability()->addNotRDistinguishable(l, n);
                }
            }
        }
        // Setting flag for every node.
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
        LOG(FATAL) << "r-characterisation sets haven't been calculated yet.";
    }
    IOListContainer result = IOListContainer(presentationLayer);
    VLOG(1) << "r-state characterisation set for " << node->getName();
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            LOG(FATAL) << "r-characterisation sets haven't been calculated yet.";
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            IOListContainer container = sequence->getInputLists();
            auto set = container.getIOLists();

            VLOG(2) << "σ(" << node->getName() << "," << n->getName() << "): " << container;
            for (auto trace : *set)
            {
                result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
            }
        }
        else
        {
            VLOG(2) << "Nodes " << node->getName() << " and " << n->getName() << " are not r-distinguishable.";
        }
    }
    VLOG(1) << "SCS(" << node->getName() << ") = " << result;
    return result;
}

IOTreeContainer Fsm::getAdaptiveRStateCharacterisationSet(shared_ptr<FsmNode> node) const
{
    if (!node->getRDistinguishability()->hasBeenCalculated())
    {
        LOG(FATAL) << "r-characterisation sets haven't been calculated yet.";
    }
    IOTreeContainer result = IOTreeContainer(presentationLayer);
    VLOG(1) << "Adaptive r-state characterisation set for " << node->getName();
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            LOG(FATAL) << "r-characterisation sets haven't been calculated yet.";
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            VLOG(2) << "σ(" << node->getName() << "," << n->getName() << "): " << sequence->str();
            result.addUniqueRemovePrefixes(sequence);
        }
    }
    VLOG(1) << "SCS(" << node->getName() << ") = " << result;
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


void Fsm::addPossibleIOTraces(shared_ptr<FsmNode> node,
                              shared_ptr<InputOutputTree> tree,
                              IOTraceContainer& iOTraceContainer,
                              const bool cleanTrailingEmptyTraces) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    VLOG(2) << "(" << node->getName() << ") " << "getPossibleIOTraces()";
    VLOG(2) << "(" << node->getName() << ") " << "  node: " << node->getName();
    VLOG(2) << "(" << node->getName() << ") " << "  tree: " << tree->str() ;
    if (tree->isEmpty())
    {
        VLOG(2)  << "(" << node->getName() << ") " << "  tree is empty. returning.";
        std::shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
        emptyTrace->setTargetNode(node);
        return;
    }

    for (int y = 0; y <= maxOutput; ++y)
    {
        shared_ptr<AdaptiveTreeNode> treeRoot = static_pointer_cast<AdaptiveTreeNode>(tree->getRoot());
        int x = treeRoot->getInput();
        bool isPossibleOutput = node->isPossibleOutput(x, y);
        VLOG(2)  << "(" << node->getName() << ") " << "  x: " << presentationLayer->getInId(x);
        VLOG(2)  << "(" << node->getName() << ") " << "  y: " << presentationLayer->getOutId(y);
        VLOG(2)  << "(" << node->getName() << ") " << "  isPossibleOutput: " << isPossibleOutput;

        if (isPossibleOutput)
        {
            unordered_set<shared_ptr<FsmNode>> nextNodes = node->afterAsSet(x, y);
            if (nextNodes.size() != 1)
            {
                LOG(FATAL)  << "The FSM does not seem to be observable.";
            }
            shared_ptr<FsmNode> nextNode = *nextNodes.begin();

            if (!tree->isDefined(y))
            {
                shared_ptr<IOTrace> trace = make_shared<IOTrace>(x, y, nextNode, presentationLayer);
                VLOG(2)  << "(" << node->getName() << ") " << "  tree is NOT defined. Adding " << *trace;
                iOTraceContainer.add(trace);
            }
            else if (tree->isDefined(y))
            {
                VLOG(2)  << "(" << node->getName() << ") " << "  tree is defined.";
                VLOG(2)  << "(" << node->getName() << ") " << "    nextNode: " << nextNode->getName();
                shared_ptr<AdaptiveTreeNode> nextTreeNode = static_pointer_cast<AdaptiveTreeNode>(treeRoot->after(y));
                VLOG(2)  << "(" << node->getName() << ") " << "    nextTreeNode input: " << presentationLayer->getInId(nextTreeNode->getInput());
                shared_ptr<InputOutputTree> nextTree = make_shared<InputOutputTree>(nextTreeNode, presentationLayer);
                VLOG(2) << "(" << node->getName() << ") " << "    nextTree: " << nextTree->str();
                VLOG(2) << "++ ENTERING RECURSION.";
                IOTraceContainer iONext;
                addPossibleIOTraces(nextNode, nextTree, iONext);
                VLOG(2) << "-- LEAVING RECURSION.";
                VLOG(2)  << "(" << node->getName() << ") " << "    iONext: " << iONext;
                shared_ptr<IOTrace> trace = make_shared<IOTrace>(x, y, nextNode, presentationLayer);
                VLOG(2) << "trace: " << trace;
                if (iONext.isEmpty())
                {
                    iONext.add(trace);
                }
                else
                {
                    iONext.concatenateToFront(trace);
                }
                VLOG(2)  << "(" << node->getName() << ") " << "    iONext: " << iONext;
                VLOG(2)  << "(" << node->getName() << ") " << "    cleanTrailingEmptyTraces: " << cleanTrailingEmptyTraces;
                if (cleanTrailingEmptyTraces)
                {
                    shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
                    //for (IOTrace& t : *iONext.getList())
                    for (auto traceIt = iONext.begin(); traceIt != iONext.end(); ++traceIt)
                    {
                        shared_ptr<const IOTrace> t = *traceIt;
                        VLOG(2)  << "(" << node->getName() << ") " << "    t.size(): " << t->size();
                        VLOG(2)  << "(" << node->getName() << ") " << "    isSuffix: " << t->isSuffix(*emptyTrace);
                        if (t->size() > 1 && t->isSuffix(*emptyTrace))
                        {
                            VLOG(2)  << "(" << node->getName() << ") " << "    REMOVING EMPTY SUFFIX from " << *t;
                            t = make_shared<IOTrace>(*t, -1, t->getTargetNode());
                        }
                    }
                    VLOG(2)  << "(" << node->getName() << ") " << "    iONext: " << iONext;
                }
                VLOG(2)  << "Adding " << iONext << " to result.";
                iOTraceContainer.add(iONext);
            }
        }
        VLOG(2)  << "(" << node->getName() << ") " << "#####################################";
    }
    VLOG(2)  << "(" << node->getName() << ") " << "--- result: " << iOTraceContainer;
}

void Fsm::addPossibleIOTraces(std::shared_ptr<FsmNode> node,
                         const IOTreeContainer& treeContainer,
                         IOTraceContainer& iOTraceContainer,
                         const bool cleanTrailingEmptyTraces) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(2));
    for (shared_ptr<InputOutputTree> tree : *treeContainer.getList())
    {
        addPossibleIOTraces(node, tree, iOTraceContainer, cleanTrailingEmptyTraces);
    }
}

bool Fsm::hasFailure() const
{
    for (const shared_ptr<FsmNode>& node : nodes)
    {

        shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> pair = node->getPair();
        if (pair == nullptr)
        {
            LOG(FATAL) << "This FSM does not seem to be a valid intersection.";
        }

        const shared_ptr<FsmNode>& specNode = pair->first;
        const shared_ptr<FsmNode>& otherNode = pair->second;

        for (const shared_ptr<FsmTransition>& otherTrans : otherNode->getTransitions())
        {
            bool foundTransition = false;
            for (const shared_ptr<FsmTransition>& specTrans : specNode->getTransitions())
            {
                if (*otherTrans->getLabel() == *specTrans->getLabel())
                {
                    foundTransition = true;
                    break;
                }
            }
            if (!foundTransition)
            {
                VLOG(1) << "The IUT has transition " << otherTrans->str() << " in state " << otherNode->getName() << " but it is missing in the specification's state " << specNode->getName() << ".";
                return true;
            }
        }
    }
    return false;
}

IOTraceContainer Fsm::bOmega(const IOTreeContainer& adaptiveTestCases, const IOTrace& trace) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    VLOG(7) << "bOmega() - adaptiveTestCases.size: " << adaptiveTestCases.size() << ", trace.size(): " << trace.size();
    IOTraceContainer result;

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
        LOG(FATAL) << "The FSM does not seem to be observable.";
    }
    shared_ptr<FsmNode> successorNode = *successorNodes.begin();
    VLOG(2) << "bOmega successorNode with " << trace << ": " << successorNode->getName();
    addPossibleIOTraces(successorNode, adaptiveTestCases, result);
    return result;
}

vector<IOTraceContainer> Fsm::bOmega(const IOTreeContainer& adaptiveTestCases, const vector<shared_ptr<InputTrace>>& inputTraces) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(6));
    VLOG(6) << "bOmega() - adaptiveTestCases.size: " << adaptiveTestCases.size() << ", inputTraces.size(): " << inputTraces.size();
    vector<IOTraceContainer> result;
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
            VLOG(2) << "produced bOmega with " << iOTrace << ": " << produced;
            IOTraceContainer::addUnique(result, produced);
        }
    }
    return result;
}

vector<IOTraceContainer> Fsm::getVPrime(const vector<shared_ptr<InputTrace>>& detStateCover)
{
    VLOG(1) << "getVPrime()";
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(4));
    vector<IOTraceContainer> result;

    VLOG(1) << "detStateCover:";
    for (auto i : detStateCover)
    {
        VLOG(1) << "  " << *i;
    }

    vector<vector<shared_ptr<OutputTrace>>> allPossibleOutputTraces;
    size_t iterations = 1;

    /* Get all possible output traces generated by the determinisitc state cover. */
    VLOG(1) << "allPossibleOutputTraces:";
    for (size_t i = 0; i < detStateCover.size(); ++i)
    {
        const shared_ptr<InputTrace>& input = detStateCover.at(i);
        //InputTrace input = InputTrace(testCase, presentationLayer);
        vector<shared_ptr<OutputTrace>> producedOutputs;
        vector<shared_ptr<FsmNode>> reached;
        getInitialState()->getPossibleOutputs(*input, producedOutputs, reached);
        iterations *= producedOutputs.size();
        for (auto v : producedOutputs)
        {
            VLOG(1) << *v;
        }
        allPossibleOutputTraces.push_back(producedOutputs);
        VLOG(1) << "--------------";
    }

    VLOG(1) << "iterations: " << iterations;

    vector<vector<vector<int>>> ot;
    size_t repetitions = iterations / allPossibleOutputTraces.at(0).size();
    for (size_t i = 0; i < detStateCover.size(); ++i)
    {
        vector<vector<int>> t;
        size_t blocks = iterations / repetitions;
        size_t outputIndex = 0;
        VLOG(1) << "repetitions: " << repetitions;
        VLOG(1) << "blocks: " << blocks;
        VLOG(1) << "outputIndex: " << outputIndex;

        for (size_t b = 0; b < blocks; ++b)
        {
            for (size_t idx = 0; idx < repetitions; ++idx)
            {
                t.push_back(allPossibleOutputTraces.at(i).at(outputIndex)->get());
            }
            outputIndex = (outputIndex + 1) % allPossibleOutputTraces.at(i).size();
        }

        if (i+1 < detStateCover.size())
        {
            repetitions /= allPossibleOutputTraces.at(i+1).size();
        }
        ot.push_back(t);
    }

    for (size_t j = 0; j < iterations; ++j)
    {
        IOTraceContainer container = IOTraceContainer();
        for (size_t i = 0; i < detStateCover.size(); ++i)
        {
            const shared_ptr<InputTrace>& input = detStateCover.at(i);
            OutputTrace output = OutputTrace(ot.at(i).at(j), presentationLayer);
            shared_ptr<IOTrace> iOTrace = make_shared<IOTrace>(*input, output);
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
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    VLOG(3) << "r():";
    VLOG(3) << "node: " << node->getName();
    VLOG(3) << "base: " << base;
    VLOG(3) << "suffix: " << suffix;


    IOTraceContainer result = IOTraceContainer();
    vector<IOTrace> prefs = suffix.getPrefixes();
    vector<IOTrace> prefixes;
    // Remove empty sequences from prefixes.
    for (IOTrace& prefix : prefs)
    {
        if (prefix.size() != 1 ||
                (prefix.getInputTrace().get().at(0) != FsmLabel::EPSILON &&
                prefix.getOutputTrace().get().at(0) != FsmLabel::EPSILON))
        {
            prefixes.push_back(prefix);
        }
    }
    VLOG(3) << "prefixes:";
    for (auto p : prefixes)
    {
        VLOG(3) << "  " << p;
    }

    for (IOTrace prefix : prefixes)
    {
        VLOG(3) << "prefix = " << prefix;
        shared_ptr<IOTrace> baseCopy = make_shared<IOTrace>(base);
        baseCopy->append(prefix);
        VLOG(3) << "v = " << baseCopy << " reaches:";
        unordered_set<shared_ptr<FsmNode>> nodes = getInitialState()->after(baseCopy->getInputTrace(), baseCopy->getOutputTrace());
        for (shared_ptr<FsmNode> n : nodes)
        {
            if (n == node)
            {
                VLOG(3) << "  " << n->getName() << " (adding " << *baseCopy << " to result), ";
                result.add(baseCopy);
            }
            else
            {
                VLOG(3) << "  " <<  n->getName() << ", ";
            }
        }
    }

    VLOG(3) << "result: " << result;

    return result;
}

IOTraceContainer Fsm::rPlus(std::shared_ptr<FsmNode> node,
                   const IOTrace& base,
                   const IOTrace& suffix,
                   const IOTraceContainer& vDoublePrime) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    VLOG(2) << "rPlus()";
    VLOG(2) << "node: " << node->getName();
    VLOG(2) << "base: " << base;
    VLOG(2) << "suffix: " << suffix;
    IOTraceContainer rResult = r(node, base, suffix);
    VLOG(2) << "rResult: " << rResult;
    if (node->isDReachable())
    {
        IOTraceCont::const_iterator vDoublePrimeElement = vDoublePrime.get(node->getDReachTrace()->getInputTrace());
        if (vDoublePrime.cend() != vDoublePrimeElement)
        {
            VLOG(2) << "  Adding: " << **vDoublePrimeElement;
            rResult.add(*vDoublePrimeElement);
            VLOG(2) << "  rPlusResult: " << rResult;
        }
    }
    return rResult;
}

bool Fsm::exceedsBound(const size_t m,
                       const IOTrace& base,
                       const IOTrace& suffix,
                       const vector<shared_ptr<InputTrace>>& takenInputs,
                       const vector<shared_ptr<FsmNode>>& states,
                       const IOTreeContainer& adaptiveTestCases,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec,
                       const Fsm& iut,
                       const bool useErroneousImplementation)
{
    if (useErroneousImplementation)
    {
        size_t lB = Fsm::lowerBound(base, suffix, takenInputs, states, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut);
        VLOG(1) << "lB: " << lB;
        return lB > m;
    }
    else
    {
        VLOG(1) << "lowerBound()";
        VLOG(1) << "base: " << base;
        VLOG(1) << "suffix: " << suffix;
        VLOG(1) << "takenInputs:";
        for (auto i : takenInputs)
        {
            VLOG(1) << "  " << *i;
        }
        VLOG(1) << "states:";
        for (auto s : states)
        {
            VLOG(1) << "  " << s->getName();
        }
        VLOG(1) << "adaptiveTestCases: " << adaptiveTestCases;
        VLOG(1) << "vDoublePrime: " << vDoublePrime;
        VLOG(1) << "dReachableStates: ";
        for (auto s : dReachableStates)
        {
            VLOG(1) << "  " << s->getName();
        }
        unordered_set<shared_ptr<FsmNode>> baseSuccessors = spec.getInitialState()->after(base);
        if (baseSuccessors.size() != 1)
        {
            LOG(FATAL) << "The Specification does not seem to be observable.";
        }
        const shared_ptr<FsmNode>& node = *baseSuccessors.begin();

        vector<shared_ptr<FsmNode>> traversedStates = node->getTraversedStates(suffix);
        size_t traversedCount = 0;
        size_t numDReach = 0;
        for (const shared_ptr<FsmNode>& n : states)
        {
            if (n->isDReachable())
            {
                ++numDReach;
            }
        }
        for (const shared_ptr<FsmNode>& n : traversedStates)
        {
            if(find(states.begin(), states.end(), n) != states.end())
            {
                ++traversedCount;
            }
        }
        size_t limit = m - numDReach + 1;
        VLOG(1) << "lBlimit " << limit;
        VLOG(1) << "traversedCount: " << traversedCount;
        return traversedCount > limit;
    }
}

size_t Fsm::lowerBound(const IOTrace& base,
                       const IOTrace& suffix,
                       const vector<shared_ptr<InputTrace>>& takenInputs,
                       const vector<shared_ptr<FsmNode>>& states,
                       const IOTreeContainer& adaptiveTestCases,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec,
                       const Fsm& iut)
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(5));
    VLOG(1) << "lowerBound()";
    VLOG(1) << "base: " << base;
    VLOG(1) << "suffix: " << suffix;
    VLOG(1) << "takenInputs:";
    for (auto i : takenInputs)
    {
        VLOG(1) << "  " << *i;
    }
    VLOG(1) << "states:";
    for (auto s : states)
    {
        VLOG(1) << "  " << s->getName();
    }
    VLOG(1) << "adaptiveTestCases: " << adaptiveTestCases;
    VLOG(1) << "vDoublePrime: " << vDoublePrime;
    VLOG(1) << "dReachableStates: ";
    for (auto s : dReachableStates)
    {
        VLOG(1) << "  " << s->getName();
    }
    size_t result = 0;
    VLOG(1) << "lb result: " << result;

    vector<IOTraceContainer> testTraces = iut.bOmega(adaptiveTestCases, takenInputs);
    VLOG(1) << "bOmega testTraces:";
    for (const auto& cont : testTraces)
    {
        VLOG(1) << "  " << cont;
    }

    for (shared_ptr<FsmNode> state : states)
    {
        const IOTraceContainer& rResult = spec.r(state, base, suffix);
        VLOG(1) << "--- state: " << state->getName();
        VLOG(1) << "rResult: " << rResult;
        result += rResult.size();
        VLOG(1) << "lb result: " << result;
        if(find(dReachableStates.begin(), dReachableStates.end(), state) != dReachableStates.end()) {
            ++result;
            VLOG(1) << "State " << state->getName() << " is d-reachable. Incrementing.";
            VLOG(1) << "lb result: " << result;
        }

        //TODO no need to calculate rResult again, as it has already been calculated above.
        const IOTraceContainer rPlusResult = spec.rPlus(state, base, suffix, vDoublePrime);
        VLOG(1) << "rPlusResult: " << rPlusResult;
        for (auto traceIt = rPlusResult.cbegin(); traceIt != rPlusResult.cend(); ++traceIt)
        {
            const shared_ptr<const IOTrace>& trace = *traceIt;
            IOTraceContainer traces = iut.bOmega(adaptiveTestCases, *trace);
            VLOG(1) << "Removing " << traces << " from testTraces.";

            IOTraceContainer::remove(testTraces, traces);

            VLOG(1) << "testTraces:";
            for (const auto& cont : testTraces)
            {
                VLOG(1) << "  " << cont;
            }
        }
    }
    result += testTraces.size();
    VLOG(1) << "lowerBound() result: " << result;
    return result;
}

bool Fsm::adaptiveStateCounting(Fsm& spec, Fsm& iut, const size_t m, IOTraceContainer& observedTraces, bool useErroneousImplementation)
{
    VLOG(1)<< "adaptiveStateCounting()";
    if (spec.isMinimal() != True)
    {
        LOG(FATAL) << "Please ensure to minimize the specification before starting adaptive state counting.";
    }
    if (iut.isMinimal() != True)
    {
        LOG(FATAL) << "Please ensure to minimize the specification before starting adaptive state counting.";
    }
#ifdef ENABLE_DEBUG_MACRO

    const string dotPrefix = "../../../resources/adaptive-test/" + spec.getName() + "-";

#endif
    spec.calcRDistinguishableStates();


    TIMED_FUNC(timerObj);
    observedTraces.clear();
    LOG(INFO) << "m: " << m;
    /**
     * Adaptive test cases for the product FSM (Ω).
     */
    const IOTreeContainer& adaptiveTestCases = spec.getAdaptiveRCharacterisationSet();
    VLOG(1) << "adaptiveTestCases: " << adaptiveTestCases;
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    VLOG(1) << "adaptiveTestCases as input traces:";
    VLOG(1) << adaptiveList;
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = spec.getMaximalSetsOfRDistinguishableStates();
    VLOG(1) << "maximalSetsOfRDistinguishableStates:";
    for (auto v : maximalSetsOfRDistinguishableStates)
    {
        stringstream ss;
        ss << "{";
        for (auto e : v)
        {
            ss << e->getName() << ", ";
        }
        ss << "}";
        VLOG(1) << ss.str();
    }

    vector<shared_ptr<InputTrace>> detStateCover;
    const vector<shared_ptr<FsmNode>>& dReachableStates = spec.getDReachableStates(detStateCover);

    VLOG(1) << "dReachableStates:";
    for (auto s : dReachableStates)
    {
        VLOG(1) << s->getName();
    }
    VLOG(1) << "detStateCover:";
    for (auto t : detStateCover)
    {
        VLOG(1) << *t;
    }

    const vector<IOTraceContainer> vPrime = iut.getVPrime(detStateCover);
    VLOG(1) << "vPrime:";
    for (auto v : vPrime)
    {
        VLOG(1) << v;
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
        stringstream ss;
        ss << "tC: ";
        for (auto w : tC)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str();
        ss.str(std::string());
        ss << "t: ";
        for (auto w : t)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str();
        ss.str(std::string());
        VLOG(1) << "adaptiveTestCases as input traces:";
        VLOG(1) << adaptiveList;
        map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;

        // Applying all input traces from T_c to this FSM.
        // All observed outputs are bein recorded.
        // If the FSM observes a failure, adaptive state counting terminates.
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            TIMED_SCOPE(timerBlkObj, "apply inputTrace");
            VLOG(1) << "############################################################";
            VLOG(1) << "  inputTrace: " << *inputTrace;
            /**
             * Holds the produced output traces for the current input trace.
             */
            vector<shared_ptr<OutputTrace>> producedOutputsSpec;
            vector<shared_ptr<OutputTrace>> producedOutputsIut;
            /**
             * Holds the reached nodes for the current input trace.
             */
            vector<shared_ptr<FsmNode>> reachedNodesSpec;
            vector<shared_ptr<FsmNode>> reachedNodesIut;

            spec.apply(*inputTrace, producedOutputsSpec, reachedNodesSpec);
            iut.apply(*inputTrace, producedOutputsIut, reachedNodesIut);
            ss << "    producedOutputs spec: ";
            for (size_t i = 0; i < producedOutputsSpec.size(); ++i)
            {
                ss << *producedOutputsSpec.at(i);
                if (i != producedOutputsSpec.size() - 1)
                {
                    ss << ", ";
                }
            }
            VLOG(1) << ss.str();
            ss.str(std::string());
            ss << "    producedOutputs IUT: ";
            for (size_t i = 0; i < producedOutputsIut.size(); ++i)
            {
                ss << *producedOutputsIut.at(i);
                if (i != producedOutputsIut.size() - 1)
                {
                    ss << ", ";
                }
            }
            VLOG(1) << ss.str();
            ss.str(std::string());
            ss << "    reachedNodes spec: ";
            for (size_t i = 0; i < reachedNodesSpec.size(); ++i)
            {
                ss << reachedNodesSpec.at(i)->getName();
                if (i != reachedNodesSpec.size() - 1)
                {
                    ss << ", ";
                }
            }
            VLOG(1) << ss.str();
            ss.str(std::string());
            ss << "    reachedNodes IUT: ";
            for (size_t i = 0; i < reachedNodesIut.size(); ++i)
            {
                ss << reachedNodesIut.at(i)->getName();
                if (i != reachedNodesIut.size() - 1)
                {
                    ss << ", ";
                }
            }
            VLOG(1) << ss.str();
            ss.str(std::string());
            PERFORMANCE_CHECKPOINT_WITH_ID(timerBlkObj, "apply inputTrace before insertion");
            observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputsIut));
            PERFORMANCE_CHECKPOINT_WITH_ID(timerBlkObj, "apply inputTrace after insertion");

            //Chek if the IUT has produced any output that can not be produced by the specification.
            for (size_t i = 0; i < producedOutputsIut.size(); ++i)
            {
                TIMED_SCOPE_IF(timerBlkObj, "Check produced output for failure", VLOG_IS_ON(1));
                const shared_ptr<OutputTrace>& outIut = producedOutputsIut.at(i);
                bool allowed = false;
                for (size_t j = 0; j < producedOutputsSpec.size(); ++j)
                {
                    const shared_ptr<OutputTrace>& outSpec = producedOutputsSpec.at(j);
                    if (*outIut == *outSpec)
                    {
                        allowed = true;
                        // Applying adaptive test cases to every node reached by the current input/output trace.
                        VLOG(1) << "----------------- Getting adaptive traces -----------------";
                        IOTraceContainer observedAdaptiveTracesIut;
                        IOTraceContainer observedAdaptiveTracesSpec;
                        const shared_ptr<FsmNode>& nodeIut = reachedNodesIut.at(i);
                        const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

                        iut.addPossibleIOTraces(nodeIut, adaptiveTestCases, observedAdaptiveTracesIut);
                        spec.addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

                        VLOG(1) << "  observedAdaptiveTracesIut: " << observedAdaptiveTracesIut;
                        VLOG(1) << "  observedAdaptiveTracesSpec: " << observedAdaptiveTracesSpec;

                        bool failure = false;
                        for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
                        {
                            const shared_ptr<const IOTrace>& trace = *traceIt;
                            if (!observedAdaptiveTracesSpec.contains(trace))
                            {
                                LOG(INFO) << "  Specification does not contain " << *trace;
                                failure = true;
                                break;
                            }
                        }
        //                PERFORMANCE_CHECKPOINT_WITH_ID(timerBlkObj, "after observedAdaptiveTracesIut loop");
                        VLOG(1) << "  concatenating: " << *inputTrace << "/" << *outIut;
                        observedAdaptiveTracesIut.concatenateToFront(inputTrace, outIut);
                        VLOG(1) << "  observedAdaptiveTraces after concatenation to front: " << observedAdaptiveTracesIut;
                        observedTraces.add(observedAdaptiveTracesIut);
                        if (failure)
                        {
                            // IUT produced an output that can not be produced by the specification.
                            LOG(INFO) << "  Failure observed:";
                            LOG(INFO) << "    Input Trace: " << *inputTrace;
                            LOG(INFO) << "    Observed adaptive traces:";
                            LOG(INFO) << observedAdaptiveTracesIut;
                            VLOG(1) << "IUT is not a reduction of the specification.";
                            return false;
                        }
                        break;
                    }
                }
                if (!allowed)
                {
                    // IUT produced an output that can not be produced by the specification.
                    LOG(INFO) << "  Failure observed:";
                    LOG(INFO) << "    Input Trace: " << *inputTrace;
                    ss << "    Produced Outputs Iut: ";
                    for (size_t i = 0; i < producedOutputsIut.size(); ++i)
                    {
                        ss << *producedOutputsIut.at(i);
                        if (i != producedOutputsIut.size() - 1)
                        {
                            ss << ", ";
                        }
                        // Adding observed traces for simple input traces only in case of failure.
                        // When no error is being observed, this traces are bein used later when
                        // concatenating them with the adaptive test cases.
                        shared_ptr<IOTrace> iOTrace = make_shared<IOTrace>(*inputTrace, *producedOutputsIut.at(i), reachedNodesIut.at(i));
                        observedTraces.add(iOTrace);
                    }
                    LOG(INFO) << ss.str();
                    ss.str(std::string());
                    ss << "    Produced Outputs Spec: ";
                    for (size_t i = 0; i < producedOutputsSpec.size(); ++i)
                    {
                        ss << *producedOutputsSpec.at(i);
                        if (i != producedOutputsSpec.size() - 1)
                        {
                            ss << ", ";
                        }
                    }
                    ss << "    Reached nodes: ";
                    for (size_t i = 0; i < reachedNodesIut.size(); ++i)
                    {
                        ss << reachedNodesIut.at(i)->getName();
                        if (i != reachedNodesIut.size() - 1)
                        {
                            ss << ", ";
                        }
                    }
                    LOG(INFO) << ss.str();
                    ss.str(std::string());
                    VLOG(1) << "IUT is not a reduction of the specification.";
                    return false;
                }
            }

            if (producedOutputsIut.size() != reachedNodesIut.size())
            {
                cerr << "Number of produced outputs and number of reached nodes do not match.";
                exit(EXIT_FAILURE);
            }
        }

        long numberToCheck = 0;
        VLOG(1) << "observedOutputsTCElements:";
        for (auto e : observedOutputsTCElements)
        {
            VLOG(1) << "  " << *e.first << ":";
            for (auto o : e.second)
            {
                VLOG(1) << "    " << *o;
                ++numberToCheck;
            }
        }
        VLOG(1) << "Number of input/output combinations: " << numberToCheck;
        bool discardInputTrace = false;
        vector<shared_ptr<InputTrace>> newT = t;
        vector<shared_ptr<InputTrace>> newTC;
        long inputTraceCount = 0;
        size_t numberInputTraces = tC.size();
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            discardInputTrace = false;
            TIMED_SCOPE(timerBlkObj, "Check input trace");
            VLOG(1) << "check inputTrace: " << *inputTrace << " (" << ++inputTraceCount << " of " << numberInputTraces << ")";
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
            VLOG(1) << "producedOutputs:";
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                VLOG(1) << "  " << *outputTrace;
            }
            long outputTraceCount = 0;
            size_t numberOutputTraces = producedOutputs.size();
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                TIMED_SCOPE_IF(timerBlkObj, "Check output trace", VLOG_IS_ON(1));
                if (discardInputTrace)
                {
                    break;
                }
                VLOG(1) << "outputTrace: " << *outputTrace << " (" << ++outputTraceCount << " of " << numberOutputTraces << ")";
                IOTrace currentTrace(*inputTrace, *outputTrace);
                VLOG(1) << "currentTrace (x_1/y_1): " << currentTrace;
                bool isLastOutputTrace = outputTrace == producedOutputs.back();
                bool maxPrefixFound = false;
                bool discardVDoublePrime = false;
                bool inputTraceMeetsCriteria = false;
                bool outputTraceMeetsCriteria = false;
                for (const IOTraceContainer& vDoublePrime : vPrime)
                {
                    TIMED_SCOPE_IF(timerBlkObj, "Check vDoublePrime", VLOG_IS_ON(1));
                    if (discardInputTrace || outputTraceMeetsCriteria || inputTraceMeetsCriteria)
                    {
                        break;
                    }
                    VLOG(1) << "vDoublePrime: " << vDoublePrime;
                    VLOG(1) << "Looking for the maximum prefix of " << currentTrace;

                    maxPrefixFound = false;
                    discardVDoublePrime = false;

                    bool isFirstVDoublePrime = vDoublePrime == vPrime.at(0);
                    if (!isFirstVDoublePrime && !useErroneousImplementation)
                    {
                        VLOG(1) << "Discarding vDoublePrime. Since using correct impllementation, discarding input trace.";
                        discardInputTrace = true;
                        break;
                    }
                    bool isLastVDoublePrime = vDoublePrime == vPrime.back();
                    shared_ptr<const IOTrace> maxPrefix = nullptr;
                    IOTrace suffix(InputTrace(spec.presentationLayer), OutputTrace(spec.presentationLayer));

                    if (useErroneousImplementation)
                    {
                        for (auto traceIt = vDoublePrime.cbegin(); traceIt != vDoublePrime.cend(); ++traceIt)
                        {
                            const shared_ptr<const IOTrace>& iOTrace = *traceIt;
                            TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-2", VLOG_IS_ON(4));
                            if (currentTrace.isPrefix(*iOTrace, false, true) &&
                                    ( !maxPrefix || maxPrefix->isEmptyTrace() || iOTrace->size() > maxPrefix->size()))
                            {
                                maxPrefix = iOTrace;
                                maxPrefixFound = true;
                            }
                        }
                        if (!maxPrefix)
                        {
                            LOG(ERROR) << "No maxPrefix (v/v'). This should not happen.";
                            continue;
                        }
                        suffix = currentTrace.getSuffix(*maxPrefix);
                    }
                    else
                    {

                        size_t prefixLength = 0;
                        for (const shared_ptr<FsmNode>& node : dReachableStates)
                        {

                            const InputTrace& dReachInput = node->getDReachTrace()->getInputTrace();
                            if (inputTrace->isPrefix(dReachInput) &&
                                    (dReachInput.size() > prefixLength ||
                                     (!dReachInput.isEmptyTrace() && dReachInput.size() >= prefixLength)))
                            {
                                prefixLength = dReachInput.size();
                                maxPrefix = node->getDReachTrace();
                                maxPrefixFound = true;
                                VLOG(1) << "new prefixLength: : " << prefixLength << " (" << dReachInput << ")";
                            }
                        }
                        suffix = IOTrace(currentTrace.removeLeadingEpsilons(), prefixLength, true);
                    }

                    VLOG(1) << "maxPrefix (v/v'): " << *maxPrefix;
                    VLOG(1) << "suffix (x/y): " << suffix;
                    bool discardSet;
                    for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
                    {
                        if (discardVDoublePrime || outputTraceMeetsCriteria || inputTraceMeetsCriteria)
                        {
                            break;
                        }
                        discardSet = false;
                        bool isLastSet = rDistStates == maximalSetsOfRDistinguishableStates.back();
                        VLOG(1) << "rDistStates:";
                        for (auto r : rDistStates)
                        {
                            VLOG(1) << "  " << r->getName();
                        }
                        TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3", VLOG_IS_ON(4));
                        for (size_t i = 0; i < rDistStates.size() - 1; ++i)
                        {
                            TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3-1", VLOG_IS_ON(5));
                            shared_ptr<FsmNode> s1 = rDistStates.at(i);
                            VLOG(1) << "############################################################";
                            VLOG(1) << "s1:" << s1->getName();
                            for (size_t j = i + 1; j < rDistStates.size(); ++j)
                            {
                                TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3-1-1", VLOG_IS_ON(6));
                                shared_ptr<FsmNode> s2 = rDistStates.at(j);
                                VLOG(1) << "----------------------------------------------------------";
                                VLOG(1) << "s2:" << s2->getName();
                                if (s1 == s2)
                                {
                                    continue;
                                }
                                const IOTraceContainer& s1RPlus = spec.rPlus(s1, *maxPrefix, suffix, vDoublePrime);
                                const IOTraceContainer& s2RPlus = spec.rPlus(s2, *maxPrefix, suffix, vDoublePrime);

                                VLOG(1) << "s1RPlus:" << s1RPlus;
                                VLOG(1) << "s2RPlus:" << s2RPlus;

                                //TODO reached nodes should already be in the traces. No need to calculate them again.
                                unordered_set<shared_ptr<FsmNode>> reached1;
                                unordered_set<shared_ptr<FsmNode>> reached2;

                                for (auto traceIt = s1RPlus.cbegin(); traceIt != s1RPlus.cend(); ++traceIt)
                                {
                                    const shared_ptr<const IOTrace>& trace = *traceIt;
                                    TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3-1-1-1", VLOG_IS_ON(7));
                                    unordered_set<shared_ptr<FsmNode>> reached = iut.getInitialState()->after(*trace);
                                    reached1.insert(reached.begin(), reached.end());
                                }
                                for (auto traceIt = s2RPlus.cbegin(); traceIt != s2RPlus.cend(); ++traceIt)
                                {
                                    const shared_ptr<const IOTrace>& trace = *traceIt;
                                    TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3-1-1-2", VLOG_IS_ON(7));
                                    unordered_set<shared_ptr<FsmNode>> reached = iut.getInitialState()->after(*trace);
                                    reached2.insert(reached.begin(), reached.end());
                                }

                                vector<shared_ptr<FsmNode>> reached1V(reached1.begin(), reached1.end());
                                vector<shared_ptr<FsmNode>> reached2V(reached2.begin(), reached2.end());

                                VLOG(1) << "reached1V:";
                                for (auto r : reached1V)
                                {
                                    VLOG(1) << "  " << r->getName();
                                }
                                VLOG(1) << "reached2V:";
                                for (auto r : reached2V)
                                {
                                    VLOG(1) << "  " << r->getName();
                                }

                                if (iut.distinguishesAllStates(reached1V, reached2V, adaptiveTestCases))
                                {
                                    VLOG(1) << "Omega ( " << adaptiveList << ") distinguishes every state of the IUT reached by an IO sequence"
                                            << " from rPlus(s1, v/v', x/y, V'') (" << s1RPlus << ") from every state of the IUT reached by an IO sequence"
                                            << " from rPlus(s2, v/v', x/y, V'') (" << s2RPlus << ").";
                                }
                                else
                                {
                                    discardSet = true;
                                    VLOG(1) << "rDistStates:";
                                    for (auto r : rDistStates)
                                    {
                                        VLOG(1) << "  " << r->getName();
                                    }
                                    VLOG(1) << "Does not r-distinguish all states.";
                                    if (isLastSet)
                                    {
                                        VLOG(1) << "isLastSet, discarding vDoublePrime: " << vDoublePrime;
                                        discardVDoublePrime = true;
                                    }
                                    break;
                                }
                            }
                            if (discardSet)
                            {
                                break;
                            }
                        }
                        if (discardVDoublePrime)
                        {
                            if (isLastVDoublePrime)
                            {
                                VLOG(1) << "isLastVDoublePrime, discarding input trace: " << *inputTrace;
                                discardInputTrace = true;
                            }
                            break;
                        }
                        if (!discardSet)
                        {
                            //size_t lB = Fsm::lowerBound(*maxPrefix, suffix, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut);
                            //VLOG(1) << "lB: " << lB;
                            bool exceedsBound = Fsm::exceedsBound(m, *maxPrefix, suffix, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut, useErroneousImplementation);
                            VLOG(1) << "exceedsBound: " << exceedsBound;
                            if (exceedsBound)
                            {
                                VLOG(1) << "Exceeded lower bound. Output trace " << *outputTrace << " meets criteria.";
                                outputTraceMeetsCriteria = true;
                                if (isLastOutputTrace)
                                {
                                    VLOG(1) << "Input trace " << *inputTrace << " meets all criteria.";
                                    inputTraceMeetsCriteria = true;
                                }
                            }
                            else
                            {
                                VLOG(1) << "Lower bound not exceeded.";
                                if (isLastSet)
                                {
                                    discardVDoublePrime = true;
                                    VLOG(1) << "isLastSet. Discarding vDoublePrime: " << vDoublePrime ;
                                    if (isLastVDoublePrime)
                                    {
                                        VLOG(1) << "isLastVDoublePrime, discarding input trace: " << *inputTrace;
                                        discardInputTrace = true;
                                    }
                                    break;
                                }
                                VLOG(1) << "Going to skip rDistStates:";
                                for (auto r : rDistStates)
                                {
                                    VLOG(1) << "  " << r->getName();
                                }
                                continue;
                            }
                        }
                    }
                }
                if (!maxPrefixFound)
                {
                    discardInputTrace = true;
                    VLOG(1) << "No max prefix found, discarding input trace: " << *inputTrace;
                }
            }

            if (discardInputTrace)
            {
                // Keeping current input trace in T_C
                VLOG(1) << "Keeping " << *inputTrace << " in T_C.";
                if (!InputTrace::contains(newTC, *inputTrace))
                {
                    newTC.push_back(inputTrace);
                }
                // Next input trace.
                continue;
            }
            else
            {
                VLOG(1) << "Removing " << *inputTrace << " from T_C.";
            }
        }
        ss << "newTC: ";
        for (auto w : newTC)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str();
        ss.str(std::string());
        // Expanding sequences.
        vector<shared_ptr<InputTrace>> expandedTC;
        for (int x = 0; x <= spec.maxInput; ++x)
        {
            for (shared_ptr<InputTrace>& inputTrace : newTC)
            {
                TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-expansion", VLOG_IS_ON(2));

                shared_ptr<InputTrace> concat = make_shared<InputTrace>(*inputTrace);
                concat->add(x);
                if (!InputTrace::contains(t, *concat) && !InputTrace::contains(expandedTC, *concat))
                {
                    expandedTC.push_back(concat);
                }
                if (!InputTrace::contains(newT, *concat))
                {
                    newT.push_back(concat);
                }
            }
        }
        ss << "expandedTC: ";
        for (auto w : expandedTC)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str();
        ss.str(std::string());
        ss << "newT: ";
        for (auto w : newT)
        {
            ss << *w << ", ";
        }
        tC = expandedTC;
        t = newT;
    }
    VLOG(1) << "  RESULT: " << observedTraces;
    VLOG(1) << "IUT is a reduction of the specification.";
    return true;
}

bool Fsm::rDistinguishesAllStates(std::vector<std::shared_ptr<FsmNode>>& nodesA,
                            std::vector<std::shared_ptr<FsmNode>>& nodesB,
                            const IOTreeContainer& adaptiveTestCases) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    for (size_t i = 0; i < nodesA.size(); ++i)
    {
        shared_ptr<FsmNode> nodeA = nodesA.at(i);
        for (size_t j = i + 1; j < nodesB.size(); ++j)
        {
            shared_ptr<FsmNode> nodeB = nodesB.at(j);
            if (nodeA == nodeB)
            {
                LOG(DEBUG) << nodeA->getName() << " == " << nodeB->getName();
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

bool Fsm::distinguishesAllStates(std::vector<std::shared_ptr<FsmNode>>& nodesA,
                            std::vector<std::shared_ptr<FsmNode>>& nodesB,
                            const IOTreeContainer& adaptiveTestCases) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    for (size_t i = 0; i < nodesA.size(); ++i)
    {
        shared_ptr<FsmNode> nodeA = nodesA.at(i);
        for (size_t j = i + 1; j < nodesB.size(); ++j)
        {
            shared_ptr<FsmNode> nodeB = nodesB.at(j);
            if (nodeA == nodeB)
            {
                LOG(DEBUG) << nodeA->getName() << " == " << nodeB->getName();
                return false;
            }
            if (!distinguishes(nodeA, nodeB, adaptiveTestCases))
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

bool Fsm::distinguishes(shared_ptr<FsmNode> nodeA,
                         shared_ptr<FsmNode> nodeB,
                         const IOTreeContainer& adaptiveTestCases) const
{
    for (shared_ptr<InputOutputTree> tree : *adaptiveTestCases.getList())
    {
        if (distinguishes(nodeA, nodeB, tree))
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
    IOTraceContainer containerA;
    addPossibleIOTraces(nodeA, adaptiveTestCase, containerA);
    IOTraceContainer containerB;
    addPossibleIOTraces(nodeB, adaptiveTestCase, containerB);
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

bool Fsm::distinguishes(shared_ptr<FsmNode> nodeA,
                         shared_ptr<FsmNode> nodeB,
                         shared_ptr<InputOutputTree> adaptiveTestCase) const
{
    IOTraceContainer containerA;
    addPossibleIOTraces(nodeA, adaptiveTestCase, containerA);
    IOTraceContainer containerB;
    addPossibleIOTraces(nodeB, adaptiveTestCase, containerB);

    return containerA != containerB;
}


IOTreeContainer Fsm::getAdaptiveRCharacterisationSet() const
{
    TIMED_FUNC(timerObj);
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
    TIMED_FUNC(timerObj);
    VLOG(1) << "getMaximalSetsOfRDistinguishableStates()";
    vector<vector<shared_ptr<FsmNode>>> result;
    result.reserve(static_cast<size_t>(getMaxNodes()));
    for (shared_ptr<FsmNode> node : nodes)
    {
        PERFORMANCE_CHECKPOINT(timerObj);
        VLOG(3) << "Looking for node " << node->getName();
        bool skip = false;
        for (vector<shared_ptr<FsmNode>>& set : result)
        {
            for (shared_ptr<FsmNode> n : set)
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
            VLOG(3) << "Skipping node " << node->getName();
            continue;
        }
        vector<shared_ptr<FsmNode>> set = {node};
        set.reserve(static_cast<size_t>(getMaxNodes()));
        VLOG(2) << "Creating set for node " << node->getName();
        for (shared_ptr<FsmNode> n : nodes)
        {
            if (node == n)
            {
                continue;
            }
            TIMED_SCOPE_IF(timerBlkObj, "Inserting node", VLOG_IS_ON(3));
            if (n->getRDistinguishability()->isRDistinguishableWith(set))
            {
                set.push_back(n);
            }
        }
        VLOG(2) << "Set size: " << set.size();
        set.resize(set.size());
        result.push_back(set);
    }
    VLOG(2) << "result size: " << result.size();
    result.resize(result.size());
    return result;
}

void Fsm::calcStateIdentificationSets()
{
    if (!isObservable())
    {
        LOG(FATAL) << "This FSM is not observable - cannot calculate the charactersiation set.";
    }
    
    if (characterisationSet == nullptr)
    {
        LOG(FATAL) << "Missing characterisation set - exit.";
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
        LOG(FATAL) << "This FSM is not observable - cannot calculate the charactersiation set.";
    }
    
    if (characterisationSet == nullptr)
    {
        LOG(FATAL) << "Missing characterisation set - exit.";
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


void Fsm::appendStateIdentificationSets(const shared_ptr<Tree>& Wp2) const
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
        LOG(FATAL) << "This FSM is not observable - cannot calculate the harmonized state identification set.";
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
                LOG(ERROR)  << "[ERR] Found inconsistency when applying HSI-Method: FSM not minimal." << endl;
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
                LOG(INFO) << "Incomplete FSM : for state " << nn->getName() << " " << nn->getId() << ", input " << x << " does not have a transition." << endl;
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

void Fsm::setPresentationLayer(const shared_ptr<FsmPresentationLayer>& ppresentationLayer)
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
                     const shared_ptr<FsmPresentationLayer>& pl,
                     const bool observable,
                     const unsigned seed) {
    TIMED_FUNC(timerObj);
    // Initialisation of random number generation
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG(DEBUG) << "createRandomFsm seed: " << s;
    }
    else {
        srand(seed);
        LOG(DEBUG) << "createRandomFsm seed: " << seed;
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
            if (observable && srcNode->hasTransition(x0, y0))
            {
                continue;
            }
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
                if (observable && srcNode->hasTransition(x, y))
                {
                    continue;
                }
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

shared_ptr<Fsm> Fsm::createMutant(const std::string & fsmName,
                                  const size_t numOutputFaults,
                                  const size_t numTransitionFaults,
                                  const unsigned seed){
    TIMED_FUNC(timerObj);
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG(DEBUG) << "createMutant seed: " << s;
    }
    else {
        srand(seed);
        LOG(DEBUG) << "createMutant seed: " << seed;
    }

    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>(*presentationLayer);

    vector<shared_ptr<FsmTransition>> cantTouchThis;
    vector<shared_ptr<FsmTransition>> transitions;
    vector<int> srcNodeIds;
    vector<int> srcNodeIdsCpy;
    vector<int> tgtNodeIdsCpy;
    
    // Create new nodes for the mutant.
    vector<shared_ptr<FsmNode> > lst;
    for ( int n = 0; n <= maxState; n++ ) {
        lst.push_back(make_shared<FsmNode>(n,fsmName,pl));
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
            srcNodeIds.push_back(n);
        }
    }
    
    // Now add transition faults to the new machine
    size_t createdTransitionFaults = 0;
    bool addedFault = true;
    while (createdTransitionFaults < numTransitionFaults && addedFault)
    {
        addedFault = false;
        srcNodeIdsCpy = srcNodeIds;
        while (srcNodeIdsCpy.size() > 0)
        {
            if (addedFault)
            {
                break;
            }
            std::vector<int>::iterator srcNodeIt = srcNodeIdsCpy.begin() + (rand() % srcNodeIdsCpy.size());
            size_t srcNodeId = static_cast<size_t>(*srcNodeIt);
            LOG(INFO) << "srcNodeId: " << srcNodeId;
            srcNodeIdsCpy.erase(srcNodeIt);

            tgtNodeIdsCpy = srcNodeIds;

            while (tgtNodeIdsCpy.size() > 0)
            {
                if (addedFault)
                {
                    break;
                }
                std::vector<int>::iterator tgtNodeIt = tgtNodeIdsCpy.begin() + (rand() % tgtNodeIdsCpy.size());
                size_t newTgtNodeId = static_cast<size_t>(*tgtNodeIt);
                LOG(INFO) << "  newTgtNodeId: " << newTgtNodeId;
                tgtNodeIdsCpy.erase(tgtNodeIt);


                transitions = lst[srcNodeId]->getTransitions();
                while (transitions.size() > 0)
                {
                    std::vector<shared_ptr<FsmTransition>>::iterator transitionIt = transitions.begin() + (rand() % transitions.size());
                    shared_ptr<FsmTransition> tr = *transitionIt;
                    LOG(INFO) << "    tr: " << tr->str();
                    transitions.erase(transitionIt);

                    if (find(cantTouchThis.begin(), cantTouchThis.end(), tr) != cantTouchThis.end())
                    {
                        LOG(INFO) << "(Transition fault) Won't touch transition " << tr->str();
                        continue;
                    }

                    if (tr->getTarget()->getId() == static_cast<int>(newTgtNodeId)) {
                        continue;
                    }
                    LOG(INFO) << "Adding transition fault:";
                    LOG(INFO) << "  Old transition: " << tr->str();
                    tr->setTarget(lst[newTgtNodeId]);
                    LOG(INFO) << "  New transition: " << tr->str();
                    ++createdTransitionFaults;
                    cantTouchThis.push_back(tr);
                    addedFault = true;
                    break;
                }
            }

        }
    }

    if (createdTransitionFaults < numTransitionFaults)
    {
        LOG(ERROR) << "Could not create all requested transition faults.";
        throw std::logic_error("Could not create all requested transition faults.");
    }
    
    // Now add output faults to the new machine
    size_t createdOutputFaults = 0;

    addedFault = true;
    while (createdOutputFaults < numOutputFaults && addedFault)
    {
        addedFault = false;
        srcNodeIdsCpy = srcNodeIds;
        while (srcNodeIdsCpy.size() > 0)
        {
            if (addedFault)
            {
                break;
            }
            std::vector<int>::iterator srcNodeIt = srcNodeIdsCpy.begin() + (rand() % srcNodeIdsCpy.size());
            size_t srcNodeId = static_cast<size_t>(*srcNodeIt);
            srcNodeIdsCpy.erase(srcNodeIt);

            transitions = lst[srcNodeId]->getTransitions();
            while (transitions.size() > 0)
            {
                std::vector<shared_ptr<FsmTransition>>::iterator transitionIt = transitions.begin() + (rand() % transitions.size());
                shared_ptr<FsmTransition> tr = *transitionIt;
                transitions.erase(transitionIt);

                if (find(cantTouchThis.begin(), cantTouchThis.end(), tr) != cantTouchThis.end())
                {
                    LOG(INFO) << "(Output fault) Won't touch transition " << tr->str();
                    continue;
                }

                int theInput = tr->getLabel()->getInput();
                int newOutVal = rand() % (maxOutput+1);
                int originalNewOutVal = newOutVal;
                bool newOutValOk;

                // We don't want to modify this transition in such a way
                // that another one with the same label and the same
                // source/target nodes already exists.
                do {

                    newOutValOk = true;

                    //BUG This is always true in the first iteration.
                    if (newOutVal == originalNewOutVal)
                    {
                        newOutValOk = false;
                    }
                    else if (lst[srcNodeId]->getTransitions().size() == 1)
                    {
                        newOutValOk = false;
                    }
                    else
                    {
                        for ( auto trOther : lst[srcNodeId]->getTransitions() ) {
                            if ( tr == trOther ) continue;
                            if ( trOther->getTarget()->getId() != tr->getTarget()->getId() )
                                continue;
                            if ( trOther->getLabel()->getInput() != theInput ) continue;
                            if ( trOther->getLabel()->getOutput() == newOutVal ) {
                                newOutValOk = false;
                                break;
                            }
                        }
                    }

                    if ( not newOutValOk ) {
                        newOutVal = (newOutVal+1) % (maxOutput+1);
                    }

                } while ( (not newOutValOk) and (originalNewOutVal != newOutVal) );

                if ( newOutValOk ) {

                    auto newLbl = make_shared<FsmLabel>(tr->getLabel()->getInput(),
                                                        newOutVal,
                                                        pl);
                    LOG(INFO) << "Adding output fault:";
                    LOG(INFO) << "  Old transition: " << tr->str();
                    tr->setLabel(newLbl);
                    LOG(INFO) << "  New transition: " << tr->str();
                    ++createdOutputFaults;
                    cantTouchThis.push_back(tr);
                    addedFault = true;
                    break;
                }
            }
        }
    }
    if (createdOutputFaults < numOutputFaults)
    {
        LOG(ERROR) << "Could not create all requested output faults.";
        throw std::logic_error("Could not create all requested output faults.");
    }
    
    shared_ptr<Fsm> result = make_shared<Fsm>(fsmName,maxInput,maxOutput,lst,pl);
    result->failOutput = failOutput;
    return result;
    
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
    VLOG(2) << "removeUnreachableNodes()";
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

    map<int,string> oldNames;
    for ( auto n : nodes ) {
        oldNames.insert(make_pair(n->getId(), n->getName()));
    }

    for ( auto n : nodes ) {
        if ( not n->hasBeenVisited() ) {
            VLOG(1) << "Removing node " << oldNames.at(n->getId()) << " (" << n->getId() << ", " << n << ").";
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


