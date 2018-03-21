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
#include "fsm/VPrimeLazy.h"
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

too_many_transition_faults::too_many_transition_faults(const std::string& msg): runtime_error(msg)
{

}

too_many_output_faults::too_many_output_faults(const std::string& msg): runtime_error(msg)
{

}

unexpected_reduction::unexpected_reduction(const std::string& msg): runtime_error(msg)
{

}

reduction_not_possible::reduction_not_possible(const std::string& msg): runtime_error(msg)
{

}

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

void Fsm::readFsmFromDot (const string & fname, const string name)
{

    maxInput = -1;
    maxOutput = -1;
    maxState = 0;
    initStateIdx = -1;
    int nodeIdCount = 0;
    map<int,shared_ptr<FsmNode>> existingNodes;
    initStateIdx = -1;
    presentationLayer = make_shared<FsmPresentationLayer>();

    if (name.empty())
    {
        // Get FSM name from file name.
        regex regFileName(".*\\/(.*?)(?:\\.\\w*)?");

        cmatch matches;
        regex_match(fname.c_str(), matches, regFileName);

        if (matches.size() > 1)
        {
            this->name = matches[1];
        }
    }

    regex regTransition("\\s*(\\d)\\s*->\\s*(\\d)\\s*\\[\\s*label=\"(.+)/(.+)\"\\s*\\]\\s*;");
    // Get maxInput and maxOutput value
    set<string> parsedInputs;
    set<string> parsedOutputs;
    map<string, int> inputStringsToIndex;
    map<string, int> outputStringsToIndex;
    ifstream inputFile (fname);
    if (inputFile.is_open())
    {
        string line;
        while (getline (inputFile, line))
        {
            cmatch matches;
            regex_match(line.c_str(), matches, regTransition);
            if (matches.size() == 5)
            {
                string in = matches[3];
                string out = matches[4];

                if (in == to_string(FsmLabel::EPSILON) || out == to_string(FsmLabel::EPSILON))
                {
                    LOG(FATAL) << "The emty input is not being supported as input or output.";
                }

                parsedInputs.insert(in);
                parsedOutputs.insert(out);
            }
        }
        maxInput = static_cast<int>(parsedInputs.size() - 1);
        maxOutput = static_cast<int>(parsedOutputs.size() - 1);
        inputFile.close();
    }
    else
    {
        LOG(FATAL) << "Unable to open input file '" << fname << "'";
    }

    VLOG(1) << "maxInput: " << maxInput;
    VLOG(1) << "maxOutput: " << maxOutput;

    // Fill presentation layer with inputs
    int i = 0;
    for (string in : parsedInputs)
    {
        presentationLayer->addIn2String(in);
        inputStringsToIndex.insert(make_pair(in, i++));
    }
    // Fill presentation layer with outputs
    i = 0;
    for (string out : parsedOutputs)
    {
        presentationLayer->addOut2String(out);
        outputStringsToIndex.insert(make_pair(out, i++));
    }

    // Getting all nodes
    inputFile.open(fname);
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
            VLOG(2) << "line:";
            VLOG(2) << "  " << line;
            VLOG(2) << "matches:";
            for (unsigned i=0; i<matches.size(); ++i) {
                VLOG(2) << "  " << matches[i];
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
                VLOG(1) << "Found state " << node->getName() << ".";
                if (nextIsInitial)
                {
                    VLOG(1) << "State " << node->getName() << " is initial state.";
                    initStateIdx = node->getId();
                    node->markAsInitial();
                    nextIsInitial = false;
                    initialSet = true;
                }
            }
        }
        inputFile.close();
    }
    else
    {
        LOG(FATAL) << "Unable to open input file '" << fname << "'";
    }

    inputFile.open(fname);
    if (inputFile.is_open())
    {
        string line;
        while (getline (inputFile, line))
        {
            cmatch matches;
            regex_match(line.c_str(), matches, regTransition);
            if (matches.size() == 5)
            {
                VLOG(1) << "Transition: " << matches[1] << " -- (" << matches[3] << "/" << matches[4] << ") --> " << matches[2];
                int sourceId = stoi(matches[1]);
                int targetId = stoi(matches[2]);
                string input = matches[3];
                string output = matches[4];

                int in = inputStringsToIndex.at(input);
                int out = outputStringsToIndex.at(output);

                shared_ptr<FsmLabel> label = make_shared<FsmLabel>(in, out, presentationLayer);
                shared_ptr<FsmTransition> trans = make_shared<FsmTransition>(existingNodes.at(sourceId), existingNodes.at(targetId), label);
                existingNodes.at(sourceId)->addTransition(trans);
            }
        }
        inputFile.close();
    }
    else
    {
        LOG(FATAL) << "Unable to open input file '" << fname << "'";
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

Fsm::Fsm(const Fsm& other): Fsm(other, other.name, other.presentationLayer)
{

}

Fsm::Fsm(const Fsm& other,
    const std::string& fsmName,
    const std::shared_ptr<FsmPresentationLayer>& presentationLayer):
    name(fsmName), currentParsedNode(nullptr), characterisationSet(nullptr),
    presentationLayer(presentationLayer)
{
    maxInput = other.maxInput;
    maxOutput = other.maxOutput;
    maxState = other.maxState;
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
    name = fsmName;
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

vector<shared_ptr<FsmNode>> Fsm::calcDReachableStates(InputTraceSet& detStateCover)
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
    detStateCover.insert(make_shared<InputTrace>(FsmLabel::EPSILON, presentationLayer));
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

            if (successorNodes.size() == 0 || producedOutputs.size() == 0)
            {
                // The Fsm isn't completely specified;
                continue;
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
                detStateCover.insert(make_shared<InputTrace>(paths.at(tgt)->getInputTrace()));
                tgt->setDReachable(paths.at(tgt));
            }
        }
    }
    resetColor();
    dReachableStates = nodes;
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


int Fsm::getNumberOfPossibleTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    return (maxInput + 1) * (maxOutput + 1) * static_cast<int>(nodePool.size());
}

float Fsm::getDegreeOfCompleteness(const int& minus, vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "getDegreeOfCompleteness()";
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    int numberOfTransitionsFound = getNumberOfDifferentInputTransitions(nodePool) - minus;
    float numberOfTransitionsPossible = (maxInput + 1.0f) * nodePool.size();
    VLOG(2) << "  numberOfTransitionsFound: " << numberOfTransitionsFound;
    VLOG(2) << "  numberOfTransitionsPossible: " << numberOfTransitionsPossible;
    return numberOfTransitionsFound / numberOfTransitionsPossible;
}

int Fsm::getNumberOfNotDefinedDeterministicTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "getNumberOfNotDefinedDeterministicTransitions()";
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    int numIn = maxInput + 1;

    int result = 0;

    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        for (int i = 0; i < numIn; ++i)
        {
            if (!n->hasTransition(i))
            {
                ++result;
            }
        }
    }
    VLOG(2) << "  result: " << result;
    return result;
}

int Fsm::getNumberOfNonDeterministicTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "getNumberOfDeterministicTransitions()";
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    int result = 0;
    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        const vector<shared_ptr<FsmTransition>>& transitions = n->getTransitions();

        unordered_map<int, int> inputOccurences;
        for (const shared_ptr<FsmTransition>& t : transitions)
        {
            inputOccurences[t->getLabel()->getInput()]++;
        }

        for (const shared_ptr<FsmTransition>& t : transitions)
        {
            if (inputOccurences.at(t->getLabel()->getInput()) > 1)
            {
                ++result;
            }
        }
    }
    VLOG(2) << "  result: " << result;
    return result;
}

int Fsm::getNumberOfTotalTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "getNumberOfTotalTransitions()";

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    int result = 0;
    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        result += n->getTransitions().size();
    }
    VLOG(2) << "  result: " << result;
    return result;
}

vector<shared_ptr<FsmTransition>> Fsm::getNonDeterministicTransitions() const
{
    vector<shared_ptr<FsmTransition>> result;
    for (const shared_ptr<FsmNode>& n : nodes)
    {
        vector<shared_ptr<FsmTransition>> r = n->getNonDeterminisitcTransitions();
        result.insert(result.end(), r.begin(), r.end());
    }
    return result;
}

float Fsm::getDegreeOfNonDeterminism(const int& diff, vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "calcDegreeOfNondeterminism()";
    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    float totalTransitions = 0;
    float numberNonDeterministicTransitions = diff;

    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        const vector<shared_ptr<FsmTransition>>& transitions = n->getTransitions();
        totalTransitions += transitions.size();

        unordered_map<int, int> inputOccurences;
        for (const shared_ptr<FsmTransition>& t : transitions)
        {
            inputOccurences[t->getLabel()->getInput()]++;
        }

        for (const shared_ptr<FsmTransition>& t : transitions)
        {
            if (inputOccurences.at(t->getLabel()->getInput()) > 1)
            {
                ++numberNonDeterministicTransitions;
            }
        }
    }
    float result = 0;
    if (totalTransitions > 0)
    {
        result = numberNonDeterministicTransitions / totalTransitions;
    }
    VLOG(2) << "  totalTransitions: " << totalTransitions;
    VLOG(2) << "  numberNonDeterministicTransitions: " << numberNonDeterministicTransitions;
    VLOG(2) << "  result: " << result;
    return result;
}

int Fsm::getNumberOfDifferentInputTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    int result = 0;
    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        unordered_set<int> differentInputs;
        for (const shared_ptr<FsmTransition>& t : n->getTransitions())
        {
            differentInputs.insert(t->getLabel()->getInput());
        }
        result += differentInputs.size();
    }
    return result;
}

Fsm Fsm::intersect(const Fsm & f, string name)
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
    
    bool isInitialState = true;

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
            // Adding the trace that reaches the new state.
            if (isInitialState)
            {
                isInitialState = false;
                nSource->setReachTrace(IOTrace::getEmptyTrace(pl));
            }
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
                        // Adding the trace that reaches the new state.
                        shared_ptr<IOTrace> nSourceReachTrace = nSource->getReachTrace();
                        shared_ptr<IOTrace> nTargetReachTrace = make_shared<IOTrace>(*tr->getLabel()->toIOTrace());
                        nTargetReachTrace->prepend(*nSourceReachTrace);
                        nTarget->setReachTrace(nTargetReachTrace);

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

    if (name.empty())
    {
        name = this->name + "x" + f.getName();
    }

    return Fsm(name, maxInput, maxOutput, fsmInterNodes, pl);
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

Fsm Fsm::transformToObservableFSM(const string& nameSuffix) const
{
    TIMED_FUNC(timerObj);
    VLOG(1) << "transformToObservableFSM()";
    CLOG(INFO, logging::globalLogger) << "Transforming to observable: " << getName();

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
    Fsm obsFsm(name + nameSuffix, maxInput, maxOutput, nodeLst, obsPl);
    return obsFsm;
}

bool Fsm::isObservable() const
{
    TIMED_FUNC(timerObj);
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

Fsm Fsm::minimiseObservableFSM(bool storeOFSMTables, const string& nameSuffix, bool prependFsmName)
{
    TIMED_FUNC(timerObj);
    if (storeOFSMTables)
    {
        /*Create new list to store all existing OFSMTables*/
        ofsmTableLst.clear();
    }
    
    /*Create the initial OFSMTable representing the FSM,
     where all FSM states belong to the same class*/
    shared_ptr<OFSMTable> tbl = make_shared<OFSMTable>(nodes, maxInput, maxOutput, presentationLayer);
    
    /*Create all possible OFSMTables, each new one from its
     predecessor, and add them to the ofsmTableLst*/
    while (tbl != nullptr)
    {
        if (storeOFSMTables)
        {
            ofsmTableLst.push_back(tbl);
        }
        const shared_ptr<OFSMTable>& nextTbl = tbl->next();
        if (nextTbl == nullptr)
        {
            break;
        }
        else
        {
            tbl = nextTbl;
        }
    }

    /*Create the minimised FSM from the last OFSMTable defined and return it*/
    Fsm fsm = tbl->toFsm(name + nameSuffix, prependFsmName);
    fsm.minimal = True;
    return fsm;
}

Fsm Fsm::minimise(bool storeOFSMTables, const string& nameSuffixMin, const string& nameSuffixObs, bool prependFsmName)
{
    VLOG(1) << "minimise()";
    TIMED_FUNC(timerObj);
    vector<shared_ptr<FsmNode>> uNodes;
    removeUnreachableNodes(uNodes);
    
    if (!isObservable())
    {
        LOG(INFO) << "Fsm is not observable. Converting.";
        return transformToObservableFSM(nameSuffixObs)
                .minimiseObservableFSM(storeOFSMTables, nameSuffixMin, prependFsmName);
    }
    
    return minimiseObservableFSM(storeOFSMTables, nameSuffixMin, prependFsmName);
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

shared_ptr<FsmNode> Fsm::getNode(int id) const
{
    for(auto n : nodes)
    {
        if (id == n->getId())
        {
            return n;
        }
    }
    return nullptr;
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
        vector<int> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(1);

        for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
        {
            int q2Id = *it;
            for (int x = 0; x <= maxInput; ++ x)
            {
                InputTrace input = InputTrace(vector<int>({x}), presentationLayer);
                shared_ptr<FsmNode> q2 = getNode(q2Id);
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
            vector<int> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(l);
            for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
            {
                // There are still nodes that can not be r-distuinguisehd from each other. Do one more iteration.
                allRDistinguishable = false;
                int q2Id = *it;
                shared_ptr<FsmNode> q2 = getNode(q2Id);
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
            vector<int> dist = node->getRDistinguishability()->getRDistinguishableWith(l);
            for (auto n : nodes)
            {
                if(node != n && find(dist.begin(), dist.end(), n->getId()) == dist.end()) {
                    node->getRDistinguishability()->addNotRDistinguishable(l, n);
                }
            }
        }
        // Setting flag for every node.
        node->getRDistinguishability()->hasBeenCalculated(true);
    }
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

    VLOG(2) << "Node " << node->getName() << " is being distuingished from: ";
    for (const shared_ptr<FsmNode>& n1 : nodes) {
        VLOG(2) << "  " << n1->getName() << ": " << rDistinguishes(node, n1, result);
    }

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

    VLOG(2) << "Node " << node->getName() << " is being r-distuingished from: ";
    for (const shared_ptr<FsmNode>& n1 : nodes) {
        VLOG(2) << "  " << n1->getName() << ": " << rDistinguishes(node, n1, result);
    }

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
                const shared_ptr<const IOTrace>& trace = make_shared<const IOTrace>(x, y, nextNode, presentationLayer);
                VLOG(2)  << "(" << node->getName() << ") " << "  tree is NOT defined. Adding " << *trace;
                iOTraceContainer.add(trace);
            }
            else if (tree->isDefined(y))
            {
                VLOG(2)  << "(" << node->getName() << ") " << "  tree is defined.";
                VLOG(2)  << "(" << node->getName() << ") " << "    nextNode: " << nextNode->getName();
                shared_ptr<AdaptiveTreeNode> nextTreeNode = static_pointer_cast<AdaptiveTreeNode>(treeRoot->after(y));
                shared_ptr<InputOutputTree> nextTree = make_shared<InputOutputTree>(nextTreeNode, presentationLayer);
                VLOG(2) << "(" << node->getName() << ") " << "    nextTree: " << nextTree->str();
                VLOG(2) << "++ ENTERING RECURSION.";
                IOTraceContainer iONext;
                addPossibleIOTraces(nextNode, nextTree, iONext);
                VLOG(2) << "-- LEAVING RECURSION.";
                VLOG(2)  << "(" << node->getName() << ") " << "    iONext: " << iONext;
                const shared_ptr<const IOTrace>& trace = make_shared<const IOTrace>(x, y, nextNode, presentationLayer);
                VLOG(2) << "trace: " << *trace;
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
                LOG(INFO) << "The IUT has transition " << otherTrans->str() << " in state " << otherNode->getName()
                          << " but it is missing in the specification's state " << specNode->getName() << ".";
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
    if (adaptiveTestCases.size() == 0)
    {
        return result;
    }

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

void Fsm::bOmega(const IOTreeContainer& adaptiveTestCases,
                 const InputTraceSet& inputTraces,
                 unordered_set<IOTraceContainer>& result) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(6));
    VLOG(6) << "bOmega() - adaptiveTestCases.size: " << adaptiveTestCases.size() << ", inputTraces.size(): " << inputTraces.size();
    if (adaptiveTestCases.size() == 0)
    {
        return;
    }
    shared_ptr<FsmNode> initialState = getInitialState();
    if (!initialState)
    {
        return;
    }
    for (const shared_ptr<InputTrace>& inputTrace : inputTraces)
    {
        vector<shared_ptr<OutputTrace>> producedOutputs;
        initialState->getPossibleOutputs(*inputTrace, producedOutputs);
        for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
        {
            IOTrace iOTrace = IOTrace(*inputTrace, *outputTrace);
            IOTraceContainer produced = bOmega(adaptiveTestCases, iOTrace);
            VLOG(1) << "produced bOmega with " << iOTrace << ": " << produced;
            result.insert(produced);
        }
    }
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

    for (const IOTrace& prefix : prefixes)
    {
        VLOG(3) << "prefix = " << prefix;
        const shared_ptr<const IOTrace>& baseExtension = make_shared<const IOTrace>(base, prefix);
        VLOG(3) << "v = " << baseExtension << " reaches:";
        unordered_set<shared_ptr<FsmNode>> nodes = getInitialState()->after(baseExtension->getInputTrace(), baseExtension->getOutputTrace());
        for (shared_ptr<FsmNode> n : nodes)
        {
            if (n == node)
            {
                VLOG(3) << "  " << n->getName() << " (adding " << *baseExtension << " to result), ";
                result.add(baseExtension);
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
                            const IOTraceContainer& vDoublePrime,
                            const bool onlyPlusPortion) const
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(7));
    VLOG(2) << "rPlus()";
    VLOG(2) << "node: " << node->getName();
    VLOG(2) << "base: " << base;
    VLOG(2) << "suffix: " << suffix;
    IOTraceContainer rResult;
    if (!onlyPlusPortion)
    {
        rResult = r(node, base, suffix);
    }
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
                       const vector<shared_ptr<FsmNode>>& states,
                       const IOTreeContainer& adaptiveTestCases,
                       unordered_set<IOTraceContainer> bOmegaT,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec,
                       const Fsm& iut)
{
    size_t lB = Fsm::lowerBound(base, suffix, states, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
    VLOG(1) << "lB: " << lB;
    return lB > m;
}

size_t Fsm::lowerBound(const IOTrace& base,
                       const IOTrace& suffix,
                       const vector<shared_ptr<FsmNode>>& states,
                       const IOTreeContainer& adaptiveTestCases,
                       unordered_set<IOTraceContainer> bOmegaT,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec,
                       const Fsm& iut)
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(5));
    VLOG(1) << "lowerBound()";
    VLOG(1) << "base: " << base;
    VLOG(1) << "suffix: " << suffix;
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

    VLOG(1) << "bOmegaT:";
    for (const auto& cont : bOmegaT)
    {
        VLOG(1) << "  " << cont;
    }

    for (shared_ptr<FsmNode> state : states)
    {
        const IOTraceContainer& rResult = spec.r(state, base, suffix);
        VLOG(1) << "--- state: " << state->getName();
        VLOG(1) << "rResult(" << state->getName() << ", " << base << ", " << suffix << "): " << rResult;
        result += rResult.size();
        VLOG(1) << "lb result: " << result;
        if(find(dReachableStates.begin(), dReachableStates.end(), state) != dReachableStates.end()) {
            ++result;
            VLOG(1) << "State " << state->getName() << " is d-reachable. Incrementing.";
            VLOG(1) << "lb result: " << result;
        }

        IOTraceContainer rPlusResult = spec.rPlus(state, base, suffix, vDoublePrime, true);
        rPlusResult.add(rResult);
        VLOG(1) << "rPlusResult: " << rPlusResult;
        for (auto traceIt = rPlusResult.cbegin(); traceIt != rPlusResult.cend(); ++traceIt)
        {
            const shared_ptr<const IOTrace>& trace = *traceIt;
            IOTraceContainer traces = iut.bOmega(adaptiveTestCases, *trace);
            VLOG(1) << "Removing " << traces << " from testTraces.";

            IOTraceContainer::remove(bOmegaT, traces);

            VLOG(1) << "testTraces:";
            for (const auto& cont : bOmegaT)
            {
                VLOG(1) << "  " << cont;
            }
        }
    }
    VLOG(1) << "bOmegaT size: " << bOmegaT.size();
    VLOG(1) << "bOmegaT:";
    for (const auto& cont : bOmegaT)
    {
        VLOG(1) << "  " << cont;
    }
    result += bOmegaT.size();
    VLOG(1) << "lowerBound() result: " << result;
    return result;
}

bool Fsm::adaptiveStateCounting(Fsm& spec, Fsm& iut, const size_t m,
                                IOTraceContainer& observedTraces,
                                shared_ptr<IOTrace>& failTrace,
                                int& iterations)
{
    VLOG(1)<< "adaptiveStateCounting()";
    if (spec.isMinimal() != True)
    {
        LOG(FATAL) << "Please ensure to minimize the specification before starting adaptive state counting.";
    }
    if (iut.isMinimal() != True)
    {
        LOG(FATAL) << "Please ensure to minimize the IUT before starting adaptive state counting.";
    }
#ifdef ENABLE_DEBUG_MACRO

    const string dotPrefix = "../../../resources/adaptive-test/" + spec.getName() + "-";

#endif
    spec.calcRDistinguishableStates();
    IOListContainer rCharacterisationSet = spec.getRCharacterisationSet();
    VLOG(1) << "Spec rCharacterisationSet:" << rCharacterisationSet;

    TIMED_FUNC(timerObj);
    observedTraces.clear();
    LOG(INFO) << "m: " << m;
    /**
     * Adaptive test cases (Ω) for the specification FSM.
     */
    const IOTreeContainer& adaptiveTestCases = spec.getAdaptiveRCharacterisationSet();
    LOG(INFO) << "adaptiveTestCases: " << adaptiveTestCases;
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    LOG(INFO) << "adaptiveTestCases as input traces:";
    LOG(INFO) << adaptiveList;
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = spec.getMaximalSetsOfRDistinguishableStates();
    LOG(INFO) << "maximalSetsOfRDistinguishableStates:";
    for (auto v : maximalSetsOfRDistinguishableStates)
    {
        stringstream ss;
        ss << "{";
        for (auto e : v)
        {
            ss << e->getName() << ", ";
        }
        ss << "}";
        LOG(INFO) << ss.str();
    }

    InputTraceSet detStateCover;
    const vector<shared_ptr<FsmNode>>& dReachableStates = spec.calcDReachableStates(detStateCover);

    LOG(INFO) << "dReachableStates:";
    for (auto s : dReachableStates)
    {
        LOG(INFO) << s->getName();
    }
    LOG(INFO) << "detStateCover:";
    for (auto t : detStateCover)
    {
        LOG(INFO) << *t;
    }

    VPrimeLazy vPrimeLazy(detStateCover, iut);

    /**
     * T - set of input sequences that have been followed by Ω.
     */
    InputTraceSet t = detStateCover;
    /**
     * Holds all B_Ω(T) for the current t.
     */
    unordered_set<IOTraceContainer> bOmegaT;
    iut.bOmega(adaptiveTestCases, t, bOmegaT);
    /**
     * T_c - set of current elements of T: those that are being considered in the search
     * through state space. The elements in T_c are the maximal sequences considered that
     * do not meet the termination criterion.
     */
    InputTraceSet tC = detStateCover;
    iterations = 0;
    while (tC.size() != 0)
    {
        ++iterations;
        stringstream ss;
#ifdef ENABLE_DEBUG_MACRO
        ss << "tC: ";
        for (auto w : tC)
        {
            ss << *w << ", ";
        }
        LOG(INFO) << ss.str();
        ss.str(std::string());
        ss << "t: ";
        for (auto w : t)
        {
            ss << *w << ", ";
        }
        LOG(INFO) << ss.str();
        ss.str(std::string());
#endif
        VLOG(1) << "adaptiveTestCases as input traces:";
        VLOG(1) << adaptiveList;
        map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;
        size_t numberInputTraces = tC.size();
        size_t inputTraceCount = 0;
        // Applying all input traces from T_c to this FSM.
        // All observed outputs are bein recorded.
        // If the FSM observes a failure, adaptive state counting terminates.
        for (const shared_ptr<InputTrace>& inputTrace : tC)
        {
            TIMED_SCOPE(timerBlkObj, "apply inputTrace");
            VLOG(1) << "############################################################";
            VLOG(1) << "  Applying inputTrace " << ++inputTraceCount << " of " << numberInputTraces << ": " << *inputTrace;
            /**
             * Hold the produced output traces for the current input trace.
             */
            vector<shared_ptr<OutputTrace>> producedOutputsSpec;
            vector<shared_ptr<OutputTrace>> producedOutputsIut;
            /**
             * Hold the reached nodes for the current input trace.
             */
            vector<shared_ptr<FsmNode>> reachedNodesSpec;
            vector<shared_ptr<FsmNode>> reachedNodesIut;

            spec.apply(*inputTrace, producedOutputsSpec, reachedNodesSpec);
            iut.apply(*inputTrace, producedOutputsIut, reachedNodesIut);
#ifdef ENABLE_DEBUG_MACRO
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
#endif
            observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputsIut));

            for (const shared_ptr<OutputTrace>& oTrace : producedOutputsIut)
            {
                observedTraces.add(make_shared<const IOTrace>(*inputTrace, *oTrace));
            }

            VLOG(1) << "Checking produced outputs for failures";
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
                        // No need to apply adaptive test cases, if there are no adaptive test cases.
                        if (adaptiveTestCases.size() > 0)
                        {
                            // Applying adaptive test cases to every node reached by the current input/output trace.
                            VLOG(1) << "----------------- Getting adaptive traces -----------------";
                            IOTraceContainer observedAdaptiveTracesIut;
                            IOTraceContainer observedAdaptiveTracesSpec;
                            const shared_ptr<FsmNode>& nodeIut = reachedNodesIut.at(i);
                            const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

                            iut.addPossibleIOTraces(nodeIut, adaptiveTestCases, observedAdaptiveTracesIut);
                            spec.addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

                            VLOG(1) << "  observedAdaptiveTracesIut (" << nodeIut->getName() << "): " << observedAdaptiveTracesIut;
                            VLOG(1) << "  observedAdaptiveTracesSpec (" << nodeSpec->getName() << "): " << observedAdaptiveTracesSpec;

                            bool failure = false;
                            for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
                            {
                                const shared_ptr<const IOTrace>& trace = *traceIt;
                                if (!observedAdaptiveTracesSpec.contains(trace))
                                {
                                    LOG(INFO) << "  Specification does not contain " << *trace;
                                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                                    IOTrace traceCopy = IOTrace(*trace);
                                    failTrace->append(traceCopy);
                                    LOG(INFO) << "failTrace: " << *failTrace;
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
                        }
                        // No failure observed, IUT output is allowed by specification.
                        // No need to search through remaining specification outputs.
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
                    }
                    LOG(INFO) << ss.str();
                    ss.str(std::string());
#ifdef ENABLE_DEBUG_MACRO
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
#endif
                    VLOG(1) << "Specification does not produce output " << *outIut << ".";
                    VLOG(1) << "IUT is not a reduction of the specification.";
                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                    LOG(INFO) << "failTrace: " << *failTrace;
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
        InputTraceSet newT = t;
        InputTraceSet newTC;
        inputTraceCount = 0;
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            bool inputTraceMeetsCriteria = true;
            TIMED_SCOPE(timerBlkObj, "Check input trace");
            LOG(INFO) << "check inputTrace: " << *inputTrace << " (" << ++inputTraceCount << " of " << numberInputTraces << ")";
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
            VLOG(1) << "producedOutputs:";
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                VLOG(1) << "  " << *outputTrace;
            }
            long outputTraceCount = 0;
            size_t numberOutputTraces = producedOutputs.size();

            shared_ptr<const InputTrace> maxInputPrefixInV = nullptr;
            for (const shared_ptr<InputTrace>& detStateTransition : detStateCover)
            {
                TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-find-prefix-in-v", VLOG_IS_ON(4));
                if (inputTrace->isPrefix(*detStateTransition, false, true) &&
                        ( !maxInputPrefixInV || maxInputPrefixInV->isEmptyTrace() || detStateTransition->size() > maxInputPrefixInV->size()))
                {
                    maxInputPrefixInV = detStateTransition;
                }
            }
            if (!maxInputPrefixInV)
            {
                LOG(FATAL) << "No prefix for input trace " << *inputTrace << " found in V. This should not happen.";
            }
            VLOG(1) << "maxInputPrefixInV: " << *maxInputPrefixInV;

            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                TIMED_SCOPE_IF(timerBlkObj, "Check output trace", VLOG_IS_ON(1));
                if (!inputTraceMeetsCriteria)
                {
                    break;
                }
                LOG(INFO) << "outputTrace: " << *outputTrace << " (" << ++outputTraceCount << " of " << numberOutputTraces << ")";
                IOTrace currentTrace(*inputTrace, *outputTrace);
                VLOG(1) << "currentTrace (x_1/y_1): " << currentTrace;
                bool outputTraceMeetsCriteria = false;
                vPrimeLazy.reset();

                VLOG(1) << "maxInputPrefixInV.size(): " << maxInputPrefixInV->size();
                shared_ptr<const IOTrace> maxIOPrefixInV = make_shared<const IOTrace>(*static_pointer_cast<const Trace>(maxInputPrefixInV),
                                                                                      *outputTrace->getPrefix(maxInputPrefixInV->size(), true));
                VLOG(1) << "maxIOPrefixInV (v/v'): " << *maxIOPrefixInV;
                IOTrace suffix(InputTrace(spec.presentationLayer), OutputTrace(spec.presentationLayer));
                suffix = currentTrace.getSuffix(*maxIOPrefixInV);
                VLOG(1) << "suffix (x/y): " << suffix;

                VLOG(1) << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext();
                while (vPrimeLazy.hasNext())
                {
                    const IOTraceContainer& vDoublePrime = vPrimeLazy.getNext();
                    TIMED_SCOPE_IF(timerBlkObj, "Check vDoublePrime", VLOG_IS_ON(1));
                    if (outputTraceMeetsCriteria)
                    {
                        break;
                    }
                    VLOG(1) << "vDoublePrime: " << vDoublePrime;

                    if (!vDoublePrime.contains(maxIOPrefixInV))
                    {
                        VLOG(1) << "vDoublePrime does not contain prefix " << *maxIOPrefixInV << ". Skipping.";
                        VLOG(1) << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext();
                        continue;
                    }
                    for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
                    {
                        VLOG(1) << "rDistStates:";
                        for (auto r : rDistStates)
                        {
                            VLOG(1) << "  " << r->getName();
                        }
                        TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-loop-2-1-3", VLOG_IS_ON(4));
                         //size_t lB = Fsm::lowerBound(*maxPrefix, suffix, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut);
                        //VLOG(1) << "lB: " << lB;
                        bool exceedsBound = Fsm::exceedsBound(m, *maxIOPrefixInV, suffix, rDistStates, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
                        VLOG(1) << "exceedsBound: " << exceedsBound;
                        if (exceedsBound)
                        {
                            VLOG(1) << "Exceeded lower bound. Output trace " << *outputTrace << " meets criteria.";
                            outputTraceMeetsCriteria = true;
                            break;
                        }
                    }
                }
                if (outputTraceMeetsCriteria == false)
                {
                    inputTraceMeetsCriteria = false;
                }
            }

            if (!inputTraceMeetsCriteria)
            {
                // Keeping current input trace in T_C
                VLOG(1) << "Keeping " << *inputTrace << " in T_C.";
                newTC.insert(inputTrace);
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
        VLOG(1) << ss.str() << endl;
        ss.str(std::string());
        // Expanding sequences.
        InputTraceSet expandedTC;
        InputTraceSet tracesAddedToT;
        LOG(INFO) << "Expanding input sequences.";
        for (int x = 0; x <= spec.maxInput; ++x)
        {
            for (const shared_ptr<InputTrace>& inputTrace : newTC)
            {
                TIMED_SCOPE_IF(timerBlkObj, "adaptiveStateCounting-expansion", VLOG_IS_ON(2));

                shared_ptr<InputTrace> concat;
                if (inputTrace->isEmptyTrace())
                {
                    concat = make_shared<InputTrace>(inputTrace->getPresentationLayer());
                }
                else
                {
                    concat = make_shared<InputTrace>(*inputTrace);
                }

                concat->add(x);
                if (!InputTrace::contains(t, concat))
                {
                    expandedTC.insert(concat);
                }

                if (newT.insert(concat).second)
                {
                    tracesAddedToT.insert(concat);
                }
            }
        }
        LOG(INFO) << "Finished expansion.";
        iut.bOmega(adaptiveTestCases, tracesAddedToT, bOmegaT);
        LOG(INFO) << "Finished calculating bOmega.";

        ss << "expandedTC: ";
        for (auto w : expandedTC)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str() << endl;
        ss.str(std::string());
        ss << "newT: ";
        for (auto w : newT)
        {
            ss << *w << ", ";
        }
        VLOG(1) << ss.str() << endl;
        ss.str(std::string());
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


bool Fsm::rDistinguishes(std::shared_ptr<FsmNode> nodeA,
                         std::shared_ptr<FsmNode> nodeB,
                         const IOListContainer& testCases) const
{
    for (vector<int> list : *testCases.getIOLists())
    {
        if (rDistinguishes(nodeA, nodeB, list))
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

bool Fsm::rDistinguishes(shared_ptr<FsmNode> nodeA,
                         shared_ptr<FsmNode> nodeB,
                         const vector<int>& list) const
{
    const InputTrace input = InputTrace(list, presentationLayer);

    vector<shared_ptr<OutputTrace>> outA;
    vector<shared_ptr<OutputTrace>> outB;

    nodeA->getPossibleOutputs(input, outA);
    nodeB->getPossibleOutputs(input, outB);

    for (const shared_ptr<OutputTrace>& traceA: outA)
    {
        for (const shared_ptr<OutputTrace>& traceB: outB)
        {
            if (*traceA == *traceB)
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
    VLOG(2) << "getMaximalSetsOfRDistinguishableStates()";
    vector<vector<shared_ptr<FsmNode>>> result;
    result.reserve(static_cast<size_t>(getMaxNodes()));
    for (shared_ptr<FsmNode> node : nodes)
    {
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
                LOG(INFO) << "Incomplete FSM : for state " << nn->getName() << " (" << nn->getId() << "), input " << x << " does not have a transition." << endl;
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

shared_ptr<Fsm> Fsm::createRandomFsm(const std::string& fsmName,
                                     const int& maxInput,
                                     const int& maxOutput,
                                     const int& maxState,
                                     const std::shared_ptr<FsmPresentationLayer>& pl,
                                     const float& degreeOfCompleteness,
                                     const float& maxDegreeOfNonDeterminism,
                                     const bool& forceNonDeterminism,
                                     const bool& minimal,
                                     const bool& observable,
                                     const unsigned& seed)
{
    VLOG(1) << "**createRandomFsm()";
    VLOG(1) << "maxInput: " << maxInput;
    VLOG(1) << "maxOutput: " << maxOutput;
    VLOG(1) << "maxState: " << maxState;
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

    int numIn = maxInput + 1;
    int numOut = maxOutput + 1;
    int numStates = maxState + 1;
    const bool degreeOfCompletenessRequired = degreeOfCompleteness > 0;

    VLOG(1) << "numIn: " << numIn;
    VLOG(1) << "numOut: " << numOut;
    VLOG(1) << "numStates: " << numStates;
    VLOG(1) << "degreeOfCompleteness: " << degreeOfCompleteness;
    VLOG(1) << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism;
    VLOG(1) << "forceNonDeterminism: " << boolalpha << forceNonDeterminism;

    if (forceNonDeterminism && numOut < 2)
    {
        LOG(FATAL) << "Can not create non-determinism with less than two output symbols.";
    }

    // Produce the nodes and put them into a vector.
    vector<shared_ptr<FsmNode>> createdNodes;
    vector<shared_ptr<FsmNode>> reachedNodes;
    vector<shared_ptr<FsmNode>> unReachedNodes;
    for (int n = 0; n < numStates; ++n) {
        shared_ptr<FsmNode> node = make_shared<FsmNode>(n, fsmName, pl);
        createdNodes.push_back(node);
        unReachedNodes.push_back(node);
    }

    shared_ptr<Fsm> fsm = make_shared<Fsm>(fsmName, maxInput, maxOutput, createdNodes, pl);

    // At index 0 of the vector, the initial state is stored, and
    // this is reachable
    reachedNodes.push_back(fsm->nodes.at(0));
    unReachedNodes.erase(unReachedNodes.begin());

    // Connecting all nodes.
    while (reachedNodes.size() != fsm->nodes.size())
    {

        const shared_ptr<FsmNode>& targetNode = unReachedNodes.back();
        unReachedNodes.pop_back();

        shared_ptr<FsmNode> srcNode;
        shared_ptr<FsmLabel> label;

        fsm->selectRandomNodeAndCreateLabel(reachedNodes, maxDegreeOfNonDeterminism, false, observable, srcNode, label);

        // We could not find a source node or a valid label.
        if (!srcNode || !label)
        {
            LOG(FATAL) << "createRandomFsm(): Could not create requested number of transitions. This shouldn't happen.";
        }

        shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
        srcNode->addTransition(transition);

        reachedNodes.push_back(targetNode);
        VLOG(1) << "Created transition " << transition->str();
    }
    VLOG(2) << "Connected all nodes.";

    fsm->addRandomTransitions(maxDegreeOfNonDeterminism, false, observable, 1.0f);

    if (degreeOfCompletenessRequired)
    {
        VLOG(2) << "Creating or removing transitions to comply with the given degree of completeness.";
        fsm->meetDegreeOfCompleteness(degreeOfCompleteness, maxDegreeOfNonDeterminism, observable);
    }

    if(forceNonDeterminism && fsm->getNumberOfNonDeterministicTransitions() < 1)
    {
        fsm->addRandomTransitions(maxDegreeOfNonDeterminism, true, observable, 1.0f);
    }

    if (minimal)
    {
        VLOG(2) << "Fsm has to be minimal. Minimizing. Num states: " << fsm->size();
        Fsm fsmMin = fsm->minimise(false, "", "", false);
        fsmMin.presentationLayer = pl;
        VLOG(2) << "Num states after minimizing: " << fsmMin.size();

        float degreeOfCompletenessMin = fsmMin.getDegreeOfCompleteness();
        size_t numStatesMin = fsmMin.size();

        bool metDegreeOfCompleteness = fsmMin.doesMeetDegreeOfCompleteness(degreeOfCompleteness);
        bool metNumberOfStates = (numStatesMin == static_cast<size_t>(numStates));
        bool metNonDeterminism = (forceNonDeterminism && fsm->getNumberOfNonDeterministicTransitions() > 0);

        int retryCount = 0;
        VLOG(2) << "metDegreeOfCompleteness: " << std::boolalpha << metDegreeOfCompleteness;
        VLOG(2) << "metnumberOfStates: " << std::boolalpha << metNumberOfStates;
        VLOG(2) << "metNonDeterminism: " << std::boolalpha << metNonDeterminism;
        while (!metDegreeOfCompleteness || !metNumberOfStates || !metNonDeterminism)
        {
            ++retryCount;
            if (retryCount < 100)
            {
                LOG(WARNING) << "Could not create the requested FSM. Trying new seed.";
                const unsigned int newSeed = static_cast<unsigned int>(rand());
                return createRandomFsm(fsmName, maxInput, maxOutput, maxState, pl, observable, newSeed);

            }
            VLOG(2) << "FSM does not meet all criteria yet.";
            if (!metNumberOfStates)
            {
                VLOG(1) << "Minimal FSM does not contain requested number of states: "
                        << numStatesMin << " < " << numStates;
                vector<shared_ptr<FsmNode>> newNodes;
                fsmMin.meetNumberOfStates(maxState, maxDegreeOfNonDeterminism, observable, newNodes);
                fsmMin.addRandomTransitions(maxDegreeOfNonDeterminism, false, observable, 1.0f, newNodes);
            }
            else if (!metDegreeOfCompleteness)
            {
                VLOG(1) << "Minimal FSM does not meet degree of completeness: "
                        << degreeOfCompletenessMin << " != " << degreeOfCompleteness;
                fsmMin.meetDegreeOfCompleteness(degreeOfCompleteness, maxDegreeOfNonDeterminism, observable);
            }
            else if (!metNonDeterminism)
            {
                VLOG(1) << "Minimal FSM is not non-deterministic.";
                fsm->addRandomTransitions(maxDegreeOfNonDeterminism, true, observable, 1.0f);
            }

            fsmMin = fsmMin.minimise(false, "", "", false);
            fsmMin.presentationLayer = pl;
            degreeOfCompletenessMin = fsmMin.getDegreeOfCompleteness();
            numStatesMin = fsmMin.size();
            metDegreeOfCompleteness = fsmMin.doesMeetDegreeOfCompleteness(degreeOfCompleteness);
            metNumberOfStates = (numStatesMin == static_cast<size_t>(numStates));
            metNonDeterminism = (forceNonDeterminism && fsm->getNumberOfNonDeterministicTransitions() > 0);
            VLOG(2) << "metDegreeOfCompleteness: " << std::boolalpha << metDegreeOfCompleteness;
            VLOG(2) << "metnumberOfStates: " << std::boolalpha << metNumberOfStates;
            VLOG(2) << "metNonDeterminism: " << std::boolalpha << metNonDeterminism;
        }
        fsm = make_shared<Fsm>(fsmMin);
    }

    VLOG(1) << "Created FSM with " << fsm->size() << " states and degreeOfCompleteness: "
            << fsm->getDegreeOfCompleteness();

    return fsm;
}


shared_ptr<Fsm> Fsm::createMutant(const std::string & fsmName,
                                  const int numOutputFaults,
                                  const int numTransitionFaults,
                                  const bool keepObservability,
                                  const unsigned seed,
                                  const shared_ptr<FsmPresentationLayer>& pLayer){
    TIMED_FUNC(timerObj);

    if (keepObservability && !isObservable())
    {
        LOG(FATAL) << "Can not keep an FSM observable that is not already observable.";
    }

    if (numOutputFaults > 0 && maxOutput < 1)
    {
        throw too_many_output_faults("Can not create output faults on FSMs with output alphabet size < 2");
    }

    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG(DEBUG) << "createMutant seed: " << s;
    }
    else {
        srand(seed);
        LOG(DEBUG) << "createMutant seed: " << seed;
    }

    LOG(DEBUG) << "numOutputFaults: " << numOutputFaults;
    LOG(DEBUG) << "numTransitionFaults: " << numTransitionFaults;

    shared_ptr<FsmPresentationLayer> pl;
    if (pLayer == nullptr)
    {
        pl = make_shared<FsmPresentationLayer>(*presentationLayer);
    }
    else
    {
        pl = make_shared<FsmPresentationLayer>(*pLayer);
    }

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
    int createdTransitionFaults = 0;
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
            VLOG(2) << "srcNodeId: " << srcNodeId;
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
                VLOG(2) << "  newTgtNodeId: " << newTgtNodeId;
                tgtNodeIdsCpy.erase(tgtNodeIt);


                transitions = lst[srcNodeId]->getTransitions();
                while (transitions.size() > 0)
                {
                    std::vector<shared_ptr<FsmTransition>>::iterator transitionIt = transitions.begin() + (rand() % transitions.size());
                    shared_ptr<FsmTransition> tr = *transitionIt;
                    VLOG(2) << "    tr: " << tr->str();
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
        LOG(INFO) << "Could not create all requested transition faults.";
        throw too_many_transition_faults("Could not create all requested transition faults.");
    }
    
    // Now add output faults to the new machine
    int createdOutputFaults = 0;

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
                int originalOutVal = tr->getLabel()->getOutput();

                int newOutVal = rand() % (maxOutput+1);
                int originalNewOutVal = newOutVal;

                if (newOutVal == originalOutVal)
                {
                    newOutVal = (newOutVal+1) % (maxOutput+1);
                }


                bool newOutValOk;

                // We don't want to modify this transition in such a way
                // that another one with the same label and the same
                // source/target nodes already exists.
                do {

                    newOutValOk = true;

                    const vector<shared_ptr<FsmTransition>>& transitions = lst[srcNodeId]->getTransitions();
                    if (transitions.size() > 1)
                    {
                        for (auto it = transitions.begin(); it != transitions.end(); ++it) {
                            const shared_ptr<FsmTransition>& trOther = *it;
                            if (tr == trOther)
                            {
                                continue;
                            }
                            if (!keepObservability && trOther->getTarget()->getId() != tr->getTarget()->getId())
                            {
                                continue;
                            }
                            if (trOther->getLabel()->getInput() != theInput)
                            {
                                continue;
                            }
                            if (trOther->getLabel()->getOutput() == newOutVal) {
                                newOutValOk = false;
                                break;
                            }
                        }
                    }

                    if ( not newOutValOk ) {
                        newOutVal = (newOutVal+1) % (maxOutput+1);
                        // Checking if we already have tried all output values.
                        if (newOutVal != originalNewOutVal &&
                                // We want to create a faulty transition, not a clone.
                                newOutVal == originalOutVal)
                        {
                            newOutVal = (newOutVal+1) % (maxOutput+1);
                        }
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
        LOG(INFO) << "Could not create all requested output faults.";
        throw too_many_output_faults("Could not create all requested output faults.");
    }
    
    shared_ptr<Fsm> result = make_shared<Fsm>(fsmName,maxInput,maxOutput,lst,pl);
    return result;
    
}

shared_ptr<Fsm> Fsm::createReduction(const string& fsmName,
                                     const bool& force,
                                     int& removedTransitions,
                                     const unsigned seed,
                                     const std::shared_ptr<FsmPresentationLayer>& pLayer) const
{
    VLOG(1) << "**createReduction()";
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG(DEBUG) << "createReduction seed: " << s;
    }
    else {
        srand(seed);
        LOG(DEBUG) << "createReduction seed: " << seed;
    }

    VLOG(2) << "Fsm:";
    VLOG(2) << *this;

    shared_ptr<FsmPresentationLayer> pl;
    if (pLayer == nullptr)
    {
        pl = make_shared<FsmPresentationLayer>(*presentationLayer);
    }
    else
    {
        pl = make_shared<FsmPresentationLayer>(*pLayer);
    }

    removedTransitions = 0;

    shared_ptr<Fsm> red = make_shared<Fsm>(*this, fsmName, pl);

    vector<shared_ptr<FsmTransition>> nonDetTransitions = red->getNonDeterministicTransitions();

    if (nonDetTransitions.empty() && force)
    {
        VLOG(1) << "Could not create reduction.";
        throw reduction_not_possible("There are no deterministic transitions.");
    }

    bool keepGoing = true;
    while (!nonDetTransitions.empty() && keepGoing)
    {
        size_t idx = static_cast<size_t>(rand()) % nonDetTransitions.size();
        const shared_ptr<FsmTransition>& transition = nonDetTransitions.at(idx);
        VLOG(2) << "Removing transition " << transition->str();
        transition->getSource()->removeTransition(transition);
        ++removedTransitions;

        nonDetTransitions = red->getNonDeterministicTransitions();
        size_t size = nonDetTransitions.size();
        int mod = 10 * static_cast<int>(ceil(size * size / 2.0f));
        VLOG(2) << "size: " << size;
        VLOG(2) << "mod: " << mod;
        keepGoing = (mod == 0) ? false : (rand() % mod) > 7;
        VLOG(2) << "keepGoing: " << boolalpha << keepGoing;
    }

    return red;

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

bool Fsm::moreTransitionsPossible(const float& maxDegreeOfNonDeterminism,
                             const bool& onlyNonDeterministic,
                             vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "moreTransitionsPossible()";
    const float newDegreeofNonDet = getDegreeOfNonDeterminism(1, nodePool);
    const int notDefDet = getNumberOfNotDefinedDeterministicTransitions();
    const int transPossible = getNumberOfPossibleTransitions(nodePool);
    const int totalDefined = getNumberOfTotalTransitions(nodePool);

    VLOG(2) << "newDegreeofNonDet: " << newDegreeofNonDet;
    VLOG(2) << "notDefDet: " << notDefDet;
    VLOG(2) << "transPossible: " << transPossible;
    VLOG(2) << "totalDefined: " << totalDefined;

    if (totalDefined >= transPossible)
    {
        return false;
    }

    if (newDegreeofNonDet >= maxDegreeOfNonDeterminism && notDefDet <= 0)
    {
        return false;
    }
    else if (newDegreeofNonDet >= maxDegreeOfNonDeterminism && notDefDet > 0 && !onlyNonDeterministic)
    {
        return true;
    }
    else if (newDegreeofNonDet >= maxDegreeOfNonDeterminism && notDefDet > 0 && onlyNonDeterministic)
    {
        return false;
    }
    else if (newDegreeofNonDet < maxDegreeOfNonDeterminism)
    {
        return true;
    }
    // This should never be reached.
    return false;
}

void Fsm::addRandomTransitions(const float& maxDegreeOfNonDeterminism,
                               const bool& onlyNonDeterministic,
                               const bool& observable,
                               const float& factor,
                               vector<shared_ptr<FsmNode>> nodePool)
{
    VLOG(1) << "**addRandomTransitions()";
    VLOG(2) << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism;
    VLOG(2) << "onlyNonDeterministic: " << onlyNonDeterministic;
    VLOG(2) << "observable: " << observable;
    VLOG(2) << "factor: " << factor;

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    if (VLOG_IS_ON(2))
    {
        VLOG(2) << "Add random transitions for nodes";
        for (const shared_ptr<FsmNode>& n : nodePool)
        {
            VLOG(2) << "  " << n->getName() << " (" << n << ")";
        }
    }

    const int numStates = static_cast<int>(nodePool.size());
    int numberOfTransitionsCreated = 0;
    bool keepGoing = moreTransitionsPossible(maxDegreeOfNonDeterminism, onlyNonDeterministic, nodePool);
    bool impossible = false;

    while (keepGoing && !impossible)
    {
        VLOG(2) << "Allowed target nodes:";
        vector<shared_ptr<FsmNode>> allowedTargetNodes;
        for (const shared_ptr<FsmNode>& n : nodes)
        {
            VLOG(2) << "  " << n->getName() << " (" << n << ")";
            allowedTargetNodes.push_back(n);
        }

        shared_ptr<FsmNode> targetNode;
        shared_ptr<FsmNode> srcNode;
        shared_ptr<FsmLabel> label;

        while (!allowedTargetNodes.empty() && (!srcNode || !label))
        {
            size_t targetNodeIndex = static_cast<size_t>(rand()) % (allowedTargetNodes.size());
            targetNode = allowedTargetNodes.at(targetNodeIndex);

            VLOG(2) << "Trying to create transition to target node " << targetNode->getName();

            selectRandomNodeAndCreateLabel(nodePool, maxDegreeOfNonDeterminism, onlyNonDeterministic, observable, srcNode, label);

            // We could not find a source node or a valid label.
            if (!srcNode || !label)
            {
                VLOG(2) << "Could not create transition to target node " << targetNode->getName() << ". Trying next node.";
                allowedTargetNodes.erase(allowedTargetNodes.begin()
                                         + static_cast<vector<shared_ptr<FsmNode>>::difference_type>(targetNodeIndex));
                continue;

            }
        }

        if (!srcNode || !label)
        {
            LOG(ERROR) << "Could not create requested number of transitions.";
            impossible = true;
        }
        else
        {
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);
            ++numberOfTransitionsCreated;
            VLOG(1) << "Created transition " << transition->str();
            VLOG(2) << "numberOfTransitionsCreated: " << numberOfTransitionsCreated;
        }

        keepGoing = moreTransitionsPossible(maxDegreeOfNonDeterminism, onlyNonDeterministic, nodePool);
        if (keepGoing)
        {
            float observableFactor = (observable) ? 1.0f : 1.75f;
            keepGoing = (rand() % static_cast<int>(round((10.0f * numStates * observableFactor * factor)))) >= numStates * 2;
        }
        VLOG(2) << "keepGoing: " << boolalpha << keepGoing;
    }

}

bool Fsm::meetDegreeOfCompleteness(const float& degreeOfCompleteness,
                                   const float& maxDegreeOfNonDeterminism,
                                   const bool& observable,
                                   vector<shared_ptr<FsmNode>> nodePool)
{
    VLOG(1) << "**meetDegreeOfCompleteness()";

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    float actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
    VLOG(2) << "actualDegreeOfCompleteness: " << actualDegreeOfCompleteness;

    bool metRequirement = false;
    if (actualDegreeOfCompleteness < degreeOfCompleteness)
    {
        VLOG(2) << "Degree of completeness: " << actualDegreeOfCompleteness << " < " << degreeOfCompleteness;
        VLOG(2) << "Going to add transitions.";
        while (actualDegreeOfCompleteness < degreeOfCompleteness)
        {
            size_t targetNodeIndex = static_cast<size_t>(rand()) % (nodes.size());
            const shared_ptr<FsmNode>& targetNode = nodes.at(targetNodeIndex);

            shared_ptr<FsmNode> srcNode;
            shared_ptr<FsmLabel> label;

            selectRandomNodeAndCreateLabel(nodePool, maxDegreeOfNonDeterminism, false, observable, srcNode, label);

            // We could not find a source node or a valid label.
            if (!srcNode || !label)
            {
                LOG(FATAL) << "meetDegreeOfCompleteness(): Could not create requested number of transitions. This shouldn't happen.";
            }

            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);

            VLOG(1) << "Created transition " << transition->str();
            actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
        }
        metRequirement = true;
    }
    else if (actualDegreeOfCompleteness > degreeOfCompleteness)
    {
        VLOG(2) << "Degree of completeness: " << actualDegreeOfCompleteness << " > " << degreeOfCompleteness;
        VLOG(2) << "Going to remove transitions.";
        while (actualDegreeOfCompleteness >= degreeOfCompleteness && !metRequirement)
        {

            if (doesMeetDegreeOfCompleteness(degreeOfCompleteness, nodePool))
            {
                VLOG(1) << "Won't remove any transition, as degree of completeness would be too small afterwards.";
                metRequirement = true;
                break;
            }

            // We have to find a transition that would affect the degree of completeness
            // when removing it. Therefore we have to find transitions per node, that are
            // deterministic.
            // If there isn't any, we have to remove several non-deterministic transitions
            // from the same node.

            // Select random node
            vector<shared_ptr<FsmNode>> selectFrom = nodePool;
            size_t nodeIdx ;
            shared_ptr<FsmNode> node;
            bool removedTransition = false;
            while (!selectFrom.empty() && !removedTransition)
            {
                nodeIdx = static_cast<size_t>(rand()) % selectFrom.size();
                node = selectFrom.at(nodeIdx);
                VLOG(2) << "Selected node " << node->getName();
                vector<shared_ptr<FsmTransition>> detTrans = node->getDeterminisitcTransitions();
                VLOG(2) << "Found " << detTrans.size() << " deterministic transitions.";
                if (!detTrans.empty())
                {
                    size_t transIndex = static_cast<size_t>(rand()) % detTrans.size();
                    const shared_ptr<FsmTransition>& tr = detTrans.at(transIndex);
                    VLOG(2) << "Removing transition " << tr->str();
                    if (!node->removeTransition(tr))
                    {
                        LOG(FATAL) << "Could not remove transition " << tr->str() << " from node " << node->getName();
                    }
                    detTrans.erase(detTrans.begin()
                                   + static_cast<vector<shared_ptr<FsmTransition>>::difference_type>(transIndex));
                    removedTransition = true;
                }
                else
                {
                    // No deterministic transition found. Trying next node.
                    VLOG(2) << "No deterministic transitions found. Trying next node.";
                    selectFrom.erase(selectFrom.begin()
                                     + static_cast<vector<shared_ptr<FsmTransition>>::difference_type>(nodeIdx));
                }

            }

            if (removedTransition)
            {
                actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
                continue;
            }

            VLOG(2) << "Could not find any node with deterministic transitions. "
                    << "Going to remove several transitions instead.";

            selectFrom = nodePool;
            removedTransition = false;
            while (!selectFrom.empty() && !removedTransition)
            {
                nodeIdx = static_cast<size_t>(rand()) % selectFrom.size();
                node = selectFrom.at(nodeIdx);
                VLOG(2) << "Selected node " << node->getName();
                vector<shared_ptr<FsmTransition>> transitions = node->getTransitions();
                VLOG(2) << "Found " << transitions.size() << " transitions.";
                if (!transitions.empty())
                {
                    const shared_ptr<FsmTransition>& trans = transitions.at(static_cast<size_t>(rand()) % transitions.size());
                    int input = trans->getLabel()->getInput();
                    VLOG(2) << "Removing all transitions with input " << input;

                    vector<shared_ptr<FsmTransition>> keepTrans;

                    for (const shared_ptr<FsmTransition>& t : transitions)
                    {
                        if (t->getLabel()->getInput() != input)
                        {
                            keepTrans.push_back(t);
                        }
                        else
                        {
                            VLOG(2) << "Removing transition " << t->str();
                        }
                    }
                    node->setTransitions(keepTrans);
                    removedTransition = true;
                }
                else
                {
                    // No deterministic transition found. Trying next node.
                    VLOG(2) << "No transitions found. Trying next node.";
                    selectFrom.erase(selectFrom.begin()
                                     + static_cast<vector<shared_ptr<FsmTransition>>::difference_type>(nodeIdx));
                }

            }

            if (removedTransition)
            {
                actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
                continue;
            }
            else
            {
                LOG(ERROR) << "Could not comply with the required degree of completeness.";
                break;
            }

        }
        if (!metRequirement)
        {
            LOG(ERROR) << "Could not comply with the required degree of completeness.";
        }
    }
    VLOG(2) << "Finished with new degree of completeness: " << actualDegreeOfCompleteness;
    return metRequirement;
}

bool Fsm::doesMeetDegreeOfCompleteness(const float& degreeOfCompleteness, vector<shared_ptr<FsmNode>> nodePool) const
{
    VLOG(2) << "doesMeetDegreeOfCompleteness()";

    if (degreeOfCompleteness <= 0)
    {
        return true;
    }

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    const float current = getDegreeOfCompleteness(0, nodePool);
    if(VLOG_IS_ON(2))
    {
        VLOG(2) << "Current degree of completeness: " << current;
    }

    bool yes;
    if (current >= degreeOfCompleteness)
    {
        float newDegreeOfCompleteness = getDegreeOfCompleteness(1, nodePool);
        yes = current >= degreeOfCompleteness && newDegreeOfCompleteness < degreeOfCompleteness;
        if (newDegreeOfCompleteness < degreeOfCompleteness)
        {
            VLOG(2) << "Fsm does meet degree of completeness, as removing one transition would "
                    << "reduce degree to " << newDegreeOfCompleteness << " < " << degreeOfCompleteness;
        }
    }
    else
    {
        VLOG(2) << "Fsm does not meet degree of completeness.";
        yes = false;
    }
    return yes;
}

void Fsm::meetNumberOfStates(const int& maxState,
                             const float& maxDegreeOfNonDeterminism,
                             const bool& observable,
                             vector<shared_ptr<FsmNode>>& createdNodes)
{
    VLOG(1) << "**meetNumberOfStates()";
    VLOG(2) << "maxState: " << maxState;

    int numIn = maxInput + 1;
    int numOut = maxOutput + 1;
    int numStates = maxState + 1;

    VLOG(2) << "numIn: " << numIn;
    VLOG(2) << "numOut: " << numOut;
    VLOG(2) << "numStates: " << numStates;

    int currentNumberNodes = static_cast<int>(size());
    int missingStates = numStates - currentNumberNodes;
    VLOG(2) << "missingStates: " << missingStates;
    VLOG(2) << "currentNumberNodes: " << currentNumberNodes;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        VLOG(2) << "  Node at index " << i << " has ID " << nodes.at(i)->getId();
    }

    // Produce the nodes and put them into a vector.
    vector<shared_ptr<FsmNode>> unReachedNodes;
    unReachedNodes.reserve(static_cast<size_t>(missingStates));
    int lowestId = currentNumberNodes;
    int highestId = lowestId + missingStates - 1;
    for (int n = highestId; n >=lowestId; --n) {
        shared_ptr<FsmNode> node = make_shared<FsmNode>(n, name, presentationLayer);
        VLOG(2) << "Created node " << node->getName() << " with id " << n << " (" << node << ")";
        unReachedNodes.push_back(node);
        createdNodes.push_back(node);
    }

    // Connecting all nodes.
    while (unReachedNodes.size() > 0)
    {
        const shared_ptr<FsmNode>& targetNode = unReachedNodes.back();
        VLOG(2) << "targetNode: " << targetNode->getName();

        shared_ptr<FsmNode> srcNode;
        shared_ptr<FsmLabel> label;

        selectRandomNodeAndCreateLabel(nodes, maxDegreeOfNonDeterminism, false, observable, srcNode, label);

        // We could not find a source node or a valid label.
        if (!srcNode || !label)
        {
            if (srcNode)
            {
                VLOG(1) << "Could not create requested number of transitions.";
                VLOG(1) << "Going to change the target of an existing one instead";
                const vector<shared_ptr<FsmTransition>>& transitions = srcNode->getTransitions();
                shared_ptr<FsmTransition> transition = transitions.at(static_cast<size_t>(rand()) % transitions.size());
                VLOG(2) << "Selected transition: " << transition->str();
                VLOG(2) << "Replacing target node " << transition->getTarget()->getName() << " with node " << targetNode->getName();
                transition->setTarget(targetNode);
                nodes.push_back(targetNode);
                VLOG(2) << "Modified transition: " << transition->str();
                unReachedNodes.pop_back();
            }
            else
            {
                LOG(FATAL) << "meetNumberOfStates(): Could not create requested number of transitions. This shouldn't happen.";
            }
        }
        else
        {
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);
            nodes.push_back(targetNode);
            unReachedNodes.pop_back();
            VLOG(1) << "Created transition " << transition->str();
        }
    }

    VLOG(2) << "Checking node IDs:";

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        VLOG(2) << "Node at index " << i << " has ID " << nodes.at(i)->getId();
        if (i != static_cast<size_t>(nodes.at(i)->getId()))
        {
            LOG(FATAL) << "Node at index " << i << " has ID " << nodes.at(i)->getId() << ". "
                       << "This is an invalid internal state and should not happen!.";
        }
    }

    VLOG(2) << "Connected all nodes.";
}

shared_ptr<FsmLabel> Fsm::createRandomLabel(const shared_ptr<FsmNode>& srcNode,
                                            const float& maxDegreeOfNonDeterminism,
                                            const bool& onlyNonDeterministic,
                                            const bool& observable) const
{
    VLOG(2) << "createRandomLabel()";

    const int numIn = maxInput + 1;
    const int numOut = maxOutput + 1;

    const bool couldAddMoreNonDet = getDegreeOfNonDeterminism(1, nodes) <= maxDegreeOfNonDeterminism;
    VLOG(2) << "couldAddMoreNonDet: " << boolalpha << couldAddMoreNonDet;

    shared_ptr<FsmLabel> label;

    if (numIn == 0 || numOut == 0)
    {
        // We can't create a label without an input or an output.
        VLOG(2) << "No input or output allowed.";
        return label;
    }

    if (onlyNonDeterministic && maxDegreeOfNonDeterminism <= 0)
    {
        LOG(FATAL) << "Invalid choice of parameters.";
    }

    // Find valid input and output values.
    vector<int> allowedInputs;
    vector<int> allowedOutputs;
    bool impossible = false;
    while (!label && !impossible)
    {
        if (!couldAddMoreNonDet)
        {
            if (onlyNonDeterministic)
            {
                LOG(FATAL) << "Requested to create only non-determinisitc lacels, "
                           << "but the degree of non-determinism is already too high.";
            }
            // Allow only inputs that are not defined in the source node,
            // but allow every output.
            // Observability isn't an issue in this case, since we choose only
            // inputs that are not yet defined.
            VLOG(2) << "Use only inputs that are not yet defined in node " << srcNode->getName();
            allowedInputs = srcNode->getNotDefinedInputs(maxInput);

            if (allowedInputs.empty())
            {
                VLOG(2) << "No input allowed. Impossible to create label.";
                // The source node has no input left under the given circumstances.
                impossible = true;
                break;
            }

            // Allow every output.
            VLOG(2) << "All outputs allowed:";
            for (int o = 0; o < numOut; ++o)
            {
                VLOG(2) << "  " << presentationLayer->getOutId(static_cast<unsigned int>(o));
                allowedOutputs.push_back(o);
            }
            int input = allowedInputs.at(static_cast<size_t>(rand()) % allowedInputs.size());
            int output = allowedOutputs.at(static_cast<size_t>(rand()) % allowedOutputs.size());
            VLOG(2) << "Selected input: " << presentationLayer->getInId(static_cast<unsigned int>(input));
            VLOG(2) << "Selected output: " << presentationLayer->getOutId(static_cast<unsigned int>(output));
            label = make_shared<FsmLabel>(input, output, presentationLayer);
        }
        else
        {
            // We can still create non-deterministic transitions.
            VLOG(2) << "Non-determinsism allowed. Allowed inputs:";
            for (int i = 0; i < numIn; ++i)
            {
                // Check if we have to create non-deterministic transitions only.
                if (!onlyNonDeterministic || srcNode->hasTransition(i))
                {
                    VLOG(2) << "  " << presentationLayer->getInId(static_cast<unsigned int>(i));
                    allowedInputs.push_back(i);
                }
            }
            if (observable)
            {
                // But we have to stay observable. Therefore we have to pick an
                // input and see, if there is any non-defined output left for
                // that input.
                VLOG(2) << "Fsm has to be observable.";
                while (!allowedInputs.empty())
                {
                    VLOG(2) << "Still inputs left.";
                    size_t inputIndex = static_cast<size_t>(rand()) % allowedInputs.size();
                    int input = allowedInputs.at(inputIndex);
                    VLOG(2) << "Getting allowed outputs for input "
                            << presentationLayer->getInId(static_cast<unsigned int>(input));
                    allowedOutputs = srcNode->getNotDefinedOutputs(input, maxOutput);

                    if (allowedOutputs.empty())
                    {
                        // There are no more outputs left. Remove selected input from
                        // allowed inputs and try again.
                        VLOG(2) << "No outputs allowed for the given input. Trying next input.";
                        allowedInputs.erase(allowedInputs.begin()
                                            + static_cast<vector<shared_ptr<int>>::difference_type>(inputIndex));
                        continue;
                    }
                    else
                    {
                        int output = allowedOutputs.at(static_cast<size_t>(rand()) % allowedOutputs.size());
                        VLOG(2) << "Selected output: " << presentationLayer->getOutId(static_cast<unsigned int>(output));
                        label = make_shared<FsmLabel>(input, output, presentationLayer);
                        break;
                    }
                }
                if (allowedInputs.empty())
                {
                    // The source node has no input left under the given circumstances.
                    VLOG(2) << "No input allowed. Impossible to create label.";
                    impossible = true;
                    break;
                }
            }
            else
            {
                // We can have non-determinism and the FSM does not have to
                // be observable. Any output is allowed.
                VLOG(2) << "No need for observability. All outputs allowed:";
                for (int o = 0; o < numOut; ++o)
                {
                    VLOG(2) << "  " << presentationLayer->getOutId(static_cast<unsigned int>(o));
                    allowedOutputs.push_back(o);
                }
                int input = allowedInputs.at(static_cast<size_t>(rand()) % allowedInputs.size());
                int output = allowedOutputs.at(static_cast<size_t>(rand()) % allowedOutputs.size());
                label = make_shared<FsmLabel>(input, output, presentationLayer);
            }
        }
    }
    if (label)
    {
        VLOG(2) << "Created label: " << *label;
    }
    else
    {
        VLOG(2) << "Could not create a label.";
    }
    return label;
}

void Fsm::selectRandomNodeAndCreateLabel(
        const vector<shared_ptr<FsmNode>> srcNodePool,
        const float& maxDegreeOfNonDeterminism,
        const bool& onlyNonDeterministic,
        const bool& observable,
        std::shared_ptr<FsmNode>& node,
        std::shared_ptr<FsmLabel>& label) const
{
    VLOG(2) << "selectRandomNodeAndCreateLabel()";

    node = nullptr;
    label = nullptr;

    // Select a reached node at random.
    VLOG(2) << "Trying to find a source node. Allowed:";
    vector<shared_ptr<FsmNode>> allowedSourceNodes;
    for (const shared_ptr<FsmNode>& n : srcNodePool)
    {
        if (!onlyNonDeterministic || !n->getTransitions().empty())
        {
            allowedSourceNodes.push_back(n);
            VLOG(2) << "  " << n->getName();
        }
    }
    while ((!node || !label) && !allowedSourceNodes.empty())
    {
        size_t srcNodeIndex = static_cast<size_t>(rand()) % (allowedSourceNodes.size());
        node = allowedSourceNodes.at(srcNodeIndex);
        VLOG(2) << "Trying node " << node->getName();
        label = createRandomLabel(node, maxDegreeOfNonDeterminism, onlyNonDeterministic, observable);
        if (!label)
        {
            // No label found. We have to try another source node
            VLOG(2) << "Could not find a label. Discarding node " << node->getName();
            allowedSourceNodes.erase(allowedSourceNodes.begin()
                                     + static_cast<vector<shared_ptr<FsmNode>>::difference_type>(srcNodeIndex));
        }
    }
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
    VLOG(1) << "removeUnreachableNodes()";
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


