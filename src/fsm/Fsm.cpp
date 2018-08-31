/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <chrono>
#include <deque>
#include <algorithm>
#include <regex>
#include <math.h>

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
#include "utils/Logger.hpp"

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

shared_ptr<FsmNode> Fsm::newNode(const int id, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p,
                                 const shared_ptr<FsmPresentationLayer>& pl)
{
    string nodeName = string("(" + p->first->getName() + to_string(p->first->getId()) + ","
                             + p->second->getName() + to_string(p->second->getId()) + ")");
    shared_ptr<FsmNode> n = make_shared<FsmNode>(id, nodeName, pl);
    n->setPair(p);
    return n;
}

bool Fsm::contains(const deque<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>>& lst,
                   const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p)
{
    for (auto pLst : lst)
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

shared_ptr<FsmNode> Fsm::findp(const vector<shared_ptr<FsmNode>>& lst,
                               const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>& p)
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
        std::cerr << "Unable to open input file";
        throw "Unable to open input file";
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
        std::cerr << "Unable to open input file";
        throw "Unable to open input file";
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
                    stringstream ss;
ss << "The emty input is not being supported as input or output.";
std::cerr << ss.str();
throw ss.str();
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
        stringstream ss;
ss << "Unable to open input file '" << fname << "'";
std::cerr << ss.str();
throw ss.str();
    }

    LOG("VERBOSE_1") << "maxInput: " << maxInput << std::endl;
    LOG("VERBOSE_1") << "maxOutput: " << maxOutput << std::endl;

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
            LOG("VERBOSE_2") << "line:" << std::endl;
            LOG("VERBOSE_2") << "  " << line << std::endl;
            LOG("VERBOSE_2") << "matches:" << std::endl;
            for (unsigned i=0; i<matches.size(); ++i) {
                LOG("VERBOSE_2") << "  " << matches[i] << std::endl;
            }
            if (matches.size() == 3)
            {
                int nodeId = stoi(matches[1]);
                string nodeName = matches[2];
                if(existingNodes.find(nodeId) != existingNodes.end()) {
                    stringstream ss;
ss << "Error while parsing dot file. The node id " << nodeId << "has been assigned more than once.";
std::cerr << ss.str();
throw ss.str();
                }
                presentationLayer->addState2String(nodeName);
                shared_ptr<FsmNode> node = make_shared<FsmNode>(nodeIdCount++, presentationLayer);
                nodes.push_back(node);
                existingNodes.insert(make_pair(nodeId, node));
                maxState++;
                LOG("VERBOSE_1") << "Found state " << node->getName() << "." << std::endl;
                if (nextIsInitial)
                {
                    LOG("VERBOSE_1") << "State " << node->getName() << " is initial state." << std::endl;
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
        stringstream ss;
ss << "Unable to open input file '" << fname << "'";
std::cerr << ss.str();
throw ss.str();
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
                LOG("VERBOSE_1") << "Transition: " << matches[1] << " -- (" << matches[3] << "/" << matches[4] << ") --> " << matches[2] << std::endl;
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
        stringstream ss;
ss << "Unable to open input file '" << fname << "'";
std::cerr << ss.str();
throw ss.str();
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
    LOG("VERBOSE_2") << "getDReachableStates()" << std::endl;
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
        LOG("VERBOSE_2") << "thisNode: " << thisNode->getName() << std::endl;

        shared_ptr<IOTrace> thisNodePath;

        if (!thisNode->isInitial())
        {
            try
            {
                thisNodePath = paths.at(thisNode);
                LOG("VERBOSE_2") << "thisNodePath: " << *thisNodePath << std::endl;
            }
            catch (out_of_range &e)
            {
                // DO nothing.
            }
        }

        for (int x = 0; x <= maxInput; ++x)
        {
            LOG("VERBOSE_2") << "x: " << presentationLayer->getInId(x) << std::endl;
            vector<int> producedOutputs;
            vector<shared_ptr<FsmNode>> successorNodes = thisNode->after(x, producedOutputs);
            LOG("VERBOSE_2") << "successorNodes:" << std::endl;
            for (auto n : successorNodes)
            {
                LOG("VERBOSE_2") << "  " << n->getName() << std::endl;
            }
            LOG("VERBOSE_2") << "producedOutputs:" << std::endl;
            for (auto n : producedOutputs)
            {
                LOG("VERBOSE_2") << "  " << presentationLayer->getOutId(n) << std::endl;
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
                        LOG("VERBOSE_2") << "Skipping." << std::endl;
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
            LOG("VERBOSE_2") << "tgt:" << tgt->getName() << std::endl;
            try
            {
                paths.at(tgt);
                // Path already exists. Do nothing.
                LOG("VERBOSE_2") << "Path already exists. Do nothing." << std::endl;
            }
            catch (out_of_range &e)
            {
                // Create new path, since it doesn't exist.
                shared_ptr<IOTrace> newPath;
                if (thisNodePath)
                {
                    newPath = make_shared<IOTrace>(*thisNodePath);
                    newPath->append(x, producedOutputs.at(0));
                    LOG("VERBOSE_2") << "newPath (appended): " << *newPath << std::endl;
                }
                else
                {
                    InputTrace in = InputTrace(x, presentationLayer);
                    OutputTrace out = OutputTrace({producedOutputs.at(0)}, presentationLayer);
                    newPath = make_shared<IOTrace>(in, out);
                    LOG("VERBOSE_2") << "newPath (new): " << *newPath << std::endl;
                }
                newPath->setTargetNode(tgt);
                paths.insert(make_pair(tgt, newPath));
            }
            if (tgt->getColor() == FsmNode::white)
            {
                LOG("VERBOSE_2") << "Target color is white. Setting grey, adding node, setting d-reach path." << std::endl;
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
    LOG("VERBOSE_2") << "getDegreeOfCompleteness()" << std::endl;
    if (nodePool.empty())
    {
        nodePool = nodes;
    }
    int numberOfTransitionsFound = getNumberOfDifferentInputTransitions(nodePool) - minus;
    float numberOfTransitionsPossible = (maxInput + 1.0f) * nodePool.size();
    LOG("VERBOSE_2") << "  numberOfTransitionsFound: " << numberOfTransitionsFound << std::endl;
    LOG("VERBOSE_2") << "  numberOfTransitionsPossible: " << numberOfTransitionsPossible << std::endl;
    return numberOfTransitionsFound / numberOfTransitionsPossible;
}

int Fsm::getNumberOfNotDefinedDeterministicTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    LOG("VERBOSE_2") << "getNumberOfNotDefinedDeterministicTransitions()" << std::endl;
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
    LOG("VERBOSE_2") << "  result: " << result << std::endl;
    return result;
}

int Fsm::getNumberOfNonDeterministicTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    LOG("VERBOSE_2") << "getNumberOfDeterministicTransitions()" << std::endl;
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
    LOG("VERBOSE_2") << "  result: " << result << std::endl;
    return result;
}

int Fsm::getNumberOfTotalTransitions(vector<shared_ptr<FsmNode>> nodePool) const
{
    LOG("VERBOSE_2") << "getNumberOfTotalTransitions()" << std::endl;

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    int result = 0;
    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        result += n->getTransitions().size();
    }
    LOG("VERBOSE_2") << "  result: " << result << std::endl;
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
    LOG("VERBOSE_2") << "calcDegreeOfNondeterminism()" << std::endl;
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
    LOG("VERBOSE_2") << "  totalTransitions: " << totalTransitions << std::endl;
    LOG("VERBOSE_2") << "  numberNonDeterministicTransitions: " << numberNonDeterministicTransitions << std::endl;
    LOG("VERBOSE_2") << "  result: " << result << std::endl;
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
    // A list of node pairs which is used to
    // control the breath-first search (BFS)
    deque<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>> nodeList;
    
    // A list of new FSM states, each state created from a pair of
    // this-nodes and f-nodes. At the end of this operation,
    // the new FSM will be created from this list.
    vector<shared_ptr<FsmNode>> fsmInterNodes;
    int id = 0;
    
    // Initially, add the pair of initial this-node and f-node
    // into the BFS list.
    nodeList.push_back(make_shared<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>(getInitialState(), f.getInitialState()));

    bool isInitialState = true;

    // We need a new presentation layer. It has the same inputs and
    // outputs as this Fsm, but the state names will be pairs of
    // state names from this FSM and f
    vector<std::string> stateNames;
    shared_ptr<FsmPresentationLayer> newPl =
    make_shared<FsmPresentationLayer>(presentationLayer->getIn2String(),
                                      presentationLayer->getOut2String(),
                                      stateNames);
    
    // This is the BFS loop, running over the (this,f)-node pairs
    while (!nodeList.empty())
    {
        // Remove the head of the list and use p to refer to it
        // p refers to the SOURCE node pair, from where all
        // outgoing transitions are investigated in this loop cycle
        shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>
            p = nodeList.front();
        nodeList.pop_front();
        
        // current node of this FSM
        shared_ptr<FsmNode> myCurrentNode = p->first;
        
        // current node of the f-FSM
        shared_ptr<FsmNode> theirCurrentNode = p->second;
        
        // Do we already have an FSM state for the new FSM
        // stored in fsmInterNodes, which is associated
        // with the current pair p?
        shared_ptr<FsmNode> nSource = findp(fsmInterNodes, p);
        
        if (nSource == nullptr)
        {
            // Set the node name as pair of the individual node names
            string newNodeName("(" + myCurrentNode->getName() + "," +
                               theirCurrentNode->getName() + ")");
            
            // Register node name in new presentation layer
            newPl->addState2String(newNodeName);
            
            // We create the new FSM state associated with p:
            // nSource is created from the state
            // pair(myCurrentNode,theirCurrentNode)
            // which is identified by p.
            nSource = newNode(id++, p, newPl);
            // Adding the trace that reaches the new state.
            if (isInitialState)
            {
                isInitialState = false;
                nSource->setReachTrace(IOTrace::getEmptyTrace(newPl));
            }
            fsmInterNodes.push_back(nSource);
            
        }
        
        // Mark this node: now all of its outgoing transitions are constructed
        nSource->setVisited();
        
        // Loop over all transitions emanating from myCurrentNode
        for (auto tr : myCurrentNode->getTransitions())
        {
            // Loop over all transitions emanating from theirCurrentNode
            for (auto trOther : theirCurrentNode->getTransitions())
            {
                /* If tr and trOther have identical labels, we can create a
                   transition for the new FSM to be created.
                   The transition has source node (myCurrentNode,theirCurrentNode)
                   and label tr.getLabel(), which is the same as
                   the label associated with the other transition,
                   and target node (tr.getTarget(),trOther.getTarget()),
                   which is the pair of the target nodes
                   of each transition.*/
                
                if (*tr->getLabel() == *trOther->getLabel())
                {
                    // New target node represented as a pair (this-node,f-node)
                    auto pTarget = make_shared<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>(tr->getTarget(), trOther->getTarget());
                    
                    // If the target node does not yet exist in the list
                    // of state for the new FSM, then create it now
                    shared_ptr<FsmNode> nTarget = findp(fsmInterNodes, pTarget);
                    if (nTarget == nullptr)
                    {
                        // Set the node name as pair of the individual node names
                        string newNodeName("(" + tr->getTarget()->getName() +
                                           "," +
                                           trOther->getTarget()->getName() + ")");
                        
                        // Register node name in new presentation layer
                        newPl->addState2String(newNodeName);
                        
                        nTarget = newNode(id++, pTarget, newPl);

                        // Adding the trace that reaches the new state.
                        shared_ptr<IOTrace> nSourceReachTrace = nSource->getReachTrace();
                        shared_ptr<IOTrace> nTargetReachTrace = make_shared<IOTrace>(*tr->getLabel()->toIOTrace());
                        nTargetReachTrace->prepend(*nSourceReachTrace);
                        nTarget->setReachTrace(nTargetReachTrace);

                        fsmInterNodes.push_back(nTarget);
                    }
                    
                    // Add transition from nSource to nTarget
                    auto newTr = make_shared<FsmTransition>(nSource,
                                                            nTarget,
                                                            tr->getLabel());

                    nSource->addTransition(newTr);
                    
                    /* Conditions for insertion of the target pair
                       into the nodeList:
                     1. the target node corresponding to the pair
                        has not yet been processed
                        (that is,  nTarget.hasBeenVisited() == false)
                     2. The target pair is not already entered into the nodeList
                     */
                    if (not (nTarget->hasBeenVisited() or
                             contains(nodeList, pTarget)))
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

    return Fsm(name, maxInput, maxOutput, fsmInterNodes, newPl);

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

Fsm Fsm::transformToObservableFSM(const string& nameSuffix) const
{
    // List to be filled with the new states to be created
    // for the observable FSM
    vector<shared_ptr<FsmNode>> nodeLst;
    
    // Breadth first search list, containing the
    // new FSM nodes, still to be processed by the algorithm
    vector<shared_ptr<FsmNode>> bfsLst;
    
    // Map a newly created node to the set of nodes from the
    // original FSM, comprised in the new FSM state
    unordered_map<shared_ptr<FsmNode>, unordered_set<shared_ptr<FsmNode>>> node2NodeLabel;
    
    // Set of nodes from the original FSM, comprised in the
    // current new node
    unordered_set<shared_ptr<FsmNode>> theNodeLabel;
    
    // For the first step of the algorithm, the initial
    // state of the original FSM is the only state comprised
    // in the initial state of the new FSM
    theNodeLabel.insert(getInitialState());

    
    // Create a new presentation layer which has
    // the same names for the inputs and outputs as
    // the old presentation layer, but still an EMPTY vector
    // of node names.
    vector<string> obsState2String;
    shared_ptr<FsmPresentationLayer> obsPl =
    make_shared<FsmPresentationLayer>(presentationLayer->getIn2String(),
                                      presentationLayer->getOut2String(),
                                      obsState2String);
    
    
    // id to be taken for the next state of the new
    // observable FSM to be created.
    int id = 0;
    
    // The initial state of the new FSM is labelled with
    // the set containing just the initial state of the old FSM
    string nodeName = labelString(theNodeLabel);
    shared_ptr<FsmNode> q0 = make_shared<FsmNode>(id++, nodeName, obsPl);
    nodeLst.push_back(q0);
    bfsLst.push_back(q0);
    node2NodeLabel[q0] = theNodeLabel;
    
    // The node label is added to the presentation layer,
    // as name of the initial state
    obsPl->addState2String(nodeName);
    
    // Loop while there is at least one node in the BFS list.
    // Initially, the list contains just the new initial state
    // of the new FSM.
    while (!bfsLst.empty())
    {
        // Pop the first node from the list
        shared_ptr<FsmNode> q = bfsLst.front();
        bfsLst.erase(bfsLst.begin());
        
        // Nested loop over all input/output labels that
        // might be associated with one or more outgoing transitions
        for (int x = 0; x <= maxInput; ++ x)
        {
            for (int y = 0; y <= maxOutput; ++ y)
            {
                // This is the transition label currently processed
                shared_ptr<FsmLabel> lbl =
                make_shared<FsmLabel>(x, y, obsPl);
                
                // Clear the set of node labels that may
                // serve as node name for the target node to be
                // created (or already exists). This target node
                // comprises all nodes of the original FSM that
                // can be reached from q under a transition labelled
                // with lbl
                theNodeLabel.clear();
                
                // Loop over all nodes of the original FSM which
                // are comprised by q
                for (shared_ptr<FsmNode> n : node2NodeLabel.at(q))
                {
                    // For each node comprised by q, check
                    // its outgoing transitions, whether they are
                    // labelled by lbl
                    for (auto tr : n->getTransitions())
                    {
                        if (*tr->getLabel() == *lbl)
                        {
                            // If so, insert the target node
                            // into the node label set
                            theNodeLabel.insert(tr->getTarget());
                        }
                    }
                }
                
                // Process only non-empty label sets.
                // An empty label set means that no transition labelled
                // by lbl exists.
                if (!theNodeLabel.empty())
                {
                    
                    shared_ptr<FsmNode> tgtNode = nullptr;
                    
                    // Loop over the node-2-label map and check
                    // if there exists already a node labelled
                    // by theNodeLabel. If this is the case,
                    // it will become the target node reached from q
                    // under lbl
                    for ( auto u : node2NodeLabel ) {
                        
                        if ( u.second == theNodeLabel ) {
                            tgtNode = u.first;
                            break;
                        }
                        
                    }
                    
                    if (tgtNode == nullptr)
                    {
                        // We need to create a new target node, to be reached
                        // from q under lbl
                        nodeName = labelString(theNodeLabel);
                        tgtNode = make_shared<FsmNode>(id++, nodeName, obsPl);
                        nodeLst.push_back(tgtNode);
                        bfsLst.push_back(tgtNode);
                        node2NodeLabel[tgtNode] = theNodeLabel;
                        obsPl->addState2String(nodeName);
                    }
                    
                    // Create the transition from q to tgtNode
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

void Fsm::calcOFSMTables() {
    
    ofsmTableLst.clear();
    
    // Create the initial OFSMTable representing the FSM,
    //  where all FSM states belong to the same class
    shared_ptr<OFSMTable> tbl = make_shared<OFSMTable>(nodes, maxInput, maxOutput, presentationLayer);
    
    // Create all possible OFSMTables, each new one from its
    // predecessor, and add them to the ofsmTableLst
    while (tbl != nullptr)
    {
        ofsmTableLst.push_back(tbl);
        tbl = tbl->next();
    }

}

Fsm Fsm::minimiseObservableFSM(const std::string& nameSuffix, bool prependFsmName)
{
    calcOFSMTables();

    // The last OFSMTable defined has classes corresponding to
    // the minimised FSM to be constructed*/
    shared_ptr<OFSMTable> tbl = ofsmTableLst.back();
    
    // Create the minimised FSM from tbl and return it
    Fsm fsm = tbl->toFsm(name + nameSuffix, prependFsmName);

    fsm.minimal = True;
    return fsm;
}

Fsm Fsm::minimise(const string& nameSuffixMin, const string& nameSuffixObs, bool prependFsmName)
{
    LOG("VERBOSE_1") << "minimise()" << std::endl;
    vector<shared_ptr<FsmNode>> uNodes;
    removeUnreachableNodes(uNodes);
    
    if (!isObservable())
    {
        LOG("INFO") << "Fsm is not observable. Converting." << std::endl;
        return transformToObservableFSM(nameSuffixObs)
                .minimiseObservableFSM(nameSuffixMin, prependFsmName);
    }
    
    return minimiseObservableFSM(nameSuffixMin, prependFsmName);
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
        stringstream ss;
ss << "This FSM is not observable - cannot calculate the charactersiation set.";
std::cerr << ss.str();
throw ss.str();
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
                LOG("VERBOSE_2") << ss.str() << std::endl;
            } catch (std::out_of_range &e) {
               // Do nothing.
            }

        }

    }
}

void Fsm::calcRDistinguishableStates()
{
    LOG("VERBOSE_2") << "calcRDistinguishableStates():" << std::endl;
    calcROneDistinguishableStates();

    size_t limit = nodes.size() * (nodes.size() - 1) / 2;
    bool allRDistinguishable = false;
    bool newDistinguishabilityCalculated = true;
    size_t maxL = 0;
    for (size_t l = 2; !allRDistinguishable && newDistinguishabilityCalculated && l <= limit; ++l)
    {
        maxL = l;
        LOG("VERBOSE_2") << "################ l = " << l << " (max " << limit << ") ################" << std::endl;
        allRDistinguishable = true;
        newDistinguishabilityCalculated = false;
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            nodes.at(k)->getRDistinguishability()->inheritDistinguishability(l);
        }
        for (size_t k = 0; k < nodes.size(); ++k)
        {
            shared_ptr<FsmNode> q1 = nodes.at(k);
            LOG("VERBOSE_3") << "q1 = " << q1->getName() << ":" << std::endl;
            vector<int> notROneDist = q1->getRDistinguishability()->getNotRDistinguishableWith(l);
            for (auto it = notROneDist.begin(); it != notROneDist.end(); ++it)
            {
                // There are still nodes that can not be r-distuinguisehd from each other. Do one more iteration.
                allRDistinguishable = false;
                int q2Id = *it;
                shared_ptr<FsmNode> q2 = getNode(q2Id);
                LOG("VERBOSE_3") << "  q2 = " << q2->getName() << ":" << std::endl;
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
                            LOG("VERBOSE_3") << "    x = " << presentationLayer->getInId(x) << ":    "
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
                            LOG("VERBOSE_3") << ss.str() << std::endl;
                            ss.str(std::string());

                            shared_ptr<AdaptiveTreeNode> childNode1 = static_pointer_cast<AdaptiveTreeNode>(childTree1->getRoot());
                            shared_ptr<TreeEdge> edge1 = make_shared<TreeEdge>(y, childNode1);
                            q1Edges.push_back(edge1);

                            //shared_ptr<TreeNode> target2 = make_shared<TreeNode>();
                            shared_ptr<InputOutputTree> childTree2 = afterNode2->getRDistinguishability()->getAdaptiveIOSequence(afterNode1);
                            ss << "      childIO2(" << afterNode2->getName() << "," << afterNode1->getName() << "): " << *childTree2 << endl;
                            LOG("VERBOSE_3") << ss.str() << std::endl;
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
                        LOG("VERBOSE_2") << ss.str() << std::endl;
                        ss.str(std::string());
                        ss << "    q2Tree: " << *q2Tree;
                        LOG("VERBOSE_2") << ss.str() << std::endl;

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
        stringstream ss;
ss << "r-characterisation sets haven't been calculated yet.";
std::cerr << ss.str();
throw ss.str();
    }
    IOListContainer result = IOListContainer(presentationLayer);
    LOG("VERBOSE_1") << "r-state characterisation set for " << node->getName() << std::endl;
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            stringstream ss;
ss << "r-characterisation sets haven't been calculated yet.";
std::cerr << ss.str();
throw ss.str();
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            IOListContainer container = sequence->getInputLists();
            auto set = container.getIOLists();

            LOG("VERBOSE_2") << "σ(" << node->getName() << "," << n->getName() << "): " << container << std::endl;
            for (auto trace : *set)
            {
                result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
            }
        }
        else
        {
            LOG("VERBOSE_2") << "Nodes " << node->getName() << " and " << n->getName() << " are not r-distinguishable." << std::endl;
        }
    }
    LOG("VERBOSE_1") << "SCS(" << node->getName() << ") = " << result << std::endl;

    LOG("VERBOSE_2") << "Node " << node->getName() << " is being distuingished from: " << std::endl;
    for (const shared_ptr<FsmNode>& n1 : nodes) {
        LOG("VERBOSE_2") << "  " << n1->getName() << ": " << rDistinguishes(node, n1, result) << std::endl;
    }

    return result;
}

IOTreeContainer Fsm::getAdaptiveRStateCharacterisationSet(shared_ptr<FsmNode> node) const
{
    if (!node->getRDistinguishability()->hasBeenCalculated())
    {
        stringstream ss;
ss << "r-characterisation sets haven't been calculated yet.";
std::cerr << ss.str();
throw ss.str();
    }
    IOTreeContainer result = IOTreeContainer(presentationLayer);
    LOG("VERBOSE_1") << "Adaptive r-state characterisation set for " << node->getName() << std::endl;
    for (shared_ptr<FsmNode> n : nodes)
    {
        if (n == node)
        {
            continue;
        }
        if (!n->getRDistinguishability()->hasBeenCalculated())
        {
            stringstream ss;
ss << "r-characterisation sets haven't been calculated yet.";
std::cerr << ss.str();
throw ss.str();
        }
        shared_ptr<InputOutputTree> sequence = node->getRDistinguishability()->getAdaptiveIOSequence(n);
        if (!sequence->isEmpty())
        {
            LOG("VERBOSE_2") << "σ(" << node->getName() << "," << n->getName() << "): " << sequence->str() << std::endl;
            result.addUniqueRemovePrefixes(sequence);
        }
    }
    LOG("VERBOSE_1") << "SCS(" << node->getName() << ") = " << result << std::endl;

    LOG("VERBOSE_2") << "Node " << node->getName() << " is being r-distuingished from: " << std::endl;
    for (const shared_ptr<FsmNode>& n1 : nodes) {
        LOG("VERBOSE_2") << "  " << n1->getName() << ": " << rDistinguishes(node, n1, result) << std::endl;
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
    LOG("VERBOSE_2") << "(" << node->getName() << ") " << "getPossibleIOTraces()" << std::endl;
    LOG("VERBOSE_2") << "(" << node->getName() << ") " << "  node: " << node->getName() << std::endl;
    LOG("VERBOSE_2") << "(" << node->getName() << ") " << "  tree: " << tree->str()  << std::endl;
    if (tree->isEmpty())
    {
        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  tree is empty. returning." << std::endl;
        std::shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
        emptyTrace->setTargetNode(node);
        return;
    }

    for (int y = 0; y <= maxOutput; ++y)
    {
        shared_ptr<AdaptiveTreeNode> treeRoot = static_pointer_cast<AdaptiveTreeNode>(tree->getRoot());
        int x = treeRoot->getInput();
        bool isPossibleOutput = node->isPossibleOutput(x, y);
        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  x: " << presentationLayer->getInId(x) << std::endl;
        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  y: " << presentationLayer->getOutId(y) << std::endl;
        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  isPossibleOutput: " << isPossibleOutput << std::endl;

        if (isPossibleOutput)
        {
            unordered_set<shared_ptr<FsmNode>> nextNodes = node->afterAsSet(x, y);
            if (nextNodes.size() != 1)
            {
                stringstream ss;
ss << "The FSM does not seem to be observable.";
std::cerr << ss.str();
throw ss.str();
            }
            shared_ptr<FsmNode> nextNode = *nextNodes.begin();

            if (!tree->isDefined(y))
            {
                const shared_ptr<const IOTrace>& trace = make_shared<const IOTrace>(x, y, nextNode, presentationLayer);
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  tree is NOT defined. Adding " << *trace << std::endl;
                iOTraceContainer.add(trace);
            }
            else if (tree->isDefined(y))
            {
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "  tree is defined." << std::endl;
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    nextNode: " << nextNode->getName() << std::endl;
                shared_ptr<AdaptiveTreeNode> nextTreeNode = static_pointer_cast<AdaptiveTreeNode>(treeRoot->after(y));
                shared_ptr<InputOutputTree> nextTree = make_shared<InputOutputTree>(nextTreeNode, presentationLayer);
                LOG("VERBOSE_2") << "(" << node->getName() << ") " << "    nextTree: " << nextTree->str() << std::endl;
                LOG("VERBOSE_2") << "++ ENTERING RECURSION." << std::endl;
                IOTraceContainer iONext;
                addPossibleIOTraces(nextNode, nextTree, iONext);
                LOG("VERBOSE_2") << "-- LEAVING RECURSION." << std::endl;
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    iONext: " << iONext << std::endl;
                const shared_ptr<const IOTrace>& trace = make_shared<const IOTrace>(x, y, nextNode, presentationLayer);
                LOG("VERBOSE_2") << "trace: " << *trace << std::endl;
                if (iONext.isEmpty())
                {
                    iONext.add(trace);
                }
                else
                {
                    iONext.concatenateToFront(trace);
                }
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    iONext: " << iONext << std::endl;
                LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    cleanTrailingEmptyTraces: " << cleanTrailingEmptyTraces << std::endl;
                if (cleanTrailingEmptyTraces)
                {
                    shared_ptr<IOTrace> emptyTrace = IOTrace::getEmptyTrace(presentationLayer);
                    //for (IOTrace& t : *iONext.getList())
                    for (auto traceIt = iONext.begin(); traceIt != iONext.end(); ++traceIt)
                    {
                        shared_ptr<const IOTrace> t = *traceIt;
                        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    t.size(): " << t->size() << std::endl;
                        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    isSuffix: " << t->isSuffix(*emptyTrace) << std::endl;
                        if (t->size() > 1 && t->isSuffix(*emptyTrace))
                        {
                            LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    REMOVING EMPTY SUFFIX from " << *t << std::endl;
                            t = make_shared<IOTrace>(*t, -1, t->getTargetNode());
                        }
                    }
                    LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "    iONext: " << iONext << std::endl;
                }
                LOG("VERBOSE_2")  << "Adding " << iONext << " to result." << std::endl;
                iOTraceContainer.add(iONext);
            }
        }
        LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "#####################################" << std::endl;
    }
    LOG("VERBOSE_2")  << "(" << node->getName() << ") " << "--- result: " << iOTraceContainer << std::endl;
}

void Fsm::addPossibleIOTraces(std::shared_ptr<FsmNode> node,
                         const IOTreeContainer& treeContainer,
                         IOTraceContainer& iOTraceContainer,
                         const bool cleanTrailingEmptyTraces) const
{
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
            stringstream ss;
ss << "This FSM does not seem to be a valid intersection.";
std::cerr << ss.str();
throw ss.str();
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
                LOG("INFO") << "The IUT has transition " << otherTrans->str() << " in state " << otherNode->getName()
                          << " but it is missing in the specification's state " << specNode->getName() << ".";
                return true;
            }
        }
    }
    return false;
}

IOTraceContainer Fsm::bOmega(const IOTreeContainer& adaptiveTestCases, const IOTrace& trace) const
{
    LOG("VERBOSE_7") << "bOmega() - adaptiveTestCases.size: " << adaptiveTestCases.size() << ", trace.size(): " << trace.size() << std::endl;
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
        stringstream ss;
ss << "The FSM does not seem to be observable.";
std::cerr << ss.str();
throw ss.str();
    }
    shared_ptr<FsmNode> successorNode = *successorNodes.begin();
    LOG("VERBOSE_2") << "bOmega successorNode with " << trace << ": " << successorNode->getName() << std::endl;
    addPossibleIOTraces(successorNode, adaptiveTestCases, result);
    return result;
}

void Fsm::bOmega(const IOTreeContainer& adaptiveTestCases,
                 const InputTraceSet& inputTraces,
                 unordered_set<IOTraceContainer>& result) const
{
    LOG("VERBOSE_6") << "bOmega() - adaptiveTestCases.size: " << adaptiveTestCases.size() << ", inputTraces.size(): " << inputTraces.size() << std::endl;
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
            LOG("VERBOSE_1") << "produced bOmega with " << iOTrace << ": " << produced << std::endl;
            result.insert(produced);
        }
    }
}

IOTraceContainer Fsm::r(std::shared_ptr<FsmNode> node,
                   const IOTrace& base,
                   const IOTrace& suffix) const
{
    LOG("VERBOSE_3") << "r():" << std::endl;
    LOG("VERBOSE_3") << "node: " << node->getName() << std::endl;
    LOG("VERBOSE_3") << "base: " << base << std::endl;
    LOG("VERBOSE_3") << "suffix: " << suffix << std::endl;


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
    LOG("VERBOSE_3") << "prefixes:" << std::endl;
    for (auto p : prefixes)
    {
        LOG("VERBOSE_3") << "  " << p << std::endl;
    }

    for (const IOTrace& prefix : prefixes)
    {
        LOG("VERBOSE_3") << "prefix = " << prefix << std::endl;
        const shared_ptr<const IOTrace>& baseExtension = make_shared<const IOTrace>(base, prefix);
        LOG("VERBOSE_3") << "v = " << baseExtension << " reaches:" << std::endl;
        unordered_set<shared_ptr<FsmNode>> nodes = getInitialState()->after(baseExtension->getInputTrace(), baseExtension->getOutputTrace());
        for (shared_ptr<FsmNode> n : nodes)
        {
            if (n == node)
            {
                LOG("VERBOSE_3") << "  " << n->getName() << " (adding " << *baseExtension << " to result), " << std::endl;
                result.add(baseExtension);
            }
            else
            {
                LOG("VERBOSE_3") << "  " <<  n->getName() << ", " << std::endl;
            }
        }
    }

    LOG("VERBOSE_3") << "result: " << result << std::endl;

    return result;
}

IOTraceContainer Fsm::rPlus(std::shared_ptr<FsmNode> node,
                            const IOTrace& base,
                            const IOTrace& suffix,
                            const IOTraceContainer& vDoublePrime,
                            const bool onlyPlusPortion) const
{
    LOG("VERBOSE_2") << "rPlus()" << std::endl;
    LOG("VERBOSE_2") << "node: " << node->getName() << std::endl;
    LOG("VERBOSE_2") << "base: " << base << std::endl;
    LOG("VERBOSE_2") << "suffix: " << suffix << std::endl;
    IOTraceContainer rResult;
    if (!onlyPlusPortion)
    {
        rResult = r(node, base, suffix);
    }
    LOG("VERBOSE_2") << "rResult: " << rResult << std::endl;
    if (node->isDReachable())
    {
        IOTraceCont::const_iterator vDoublePrimeElement = vDoublePrime.get(node->getDReachTrace()->getInputTrace());
        if (vDoublePrime.cend() != vDoublePrimeElement)
        {
            LOG("VERBOSE_2") << "  Adding: " << **vDoublePrimeElement << std::endl;
            rResult.add(*vDoublePrimeElement);
            LOG("VERBOSE_2") << "  rPlusResult: " << rResult << std::endl;
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
    LOG("VERBOSE_1") << "lB: " << lB << std::endl;
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
    LOG("VERBOSE_1") << "lowerBound()" << std::endl;
    LOG("VERBOSE_1") << "base: " << base << std::endl;
    LOG("VERBOSE_1") << "suffix: " << suffix << std::endl;
    LOG("VERBOSE_1") << "states:" << std::endl;
    for (auto s : states)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    LOG("VERBOSE_1") << "adaptiveTestCases: " << adaptiveTestCases << std::endl;
    LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;
    LOG("VERBOSE_1") << "dReachableStates: " << std::endl;
    for (auto s : dReachableStates)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    size_t result = 0;
    LOG("VERBOSE_1") << "lb result: " << result << std::endl;

    LOG("VERBOSE_1") << "bOmegaT:" << std::endl;
    for (const auto& cont : bOmegaT)
    {
        LOG("VERBOSE_1") << "  " << cont << std::endl;
    }

    for (shared_ptr<FsmNode> state : states)
    {
        const IOTraceContainer& rResult = spec.r(state, base, suffix);
        LOG("VERBOSE_1") << "--- state: " << state->getName() << std::endl;
        LOG("VERBOSE_1") << "rResult(" << state->getName() << ", " << base << ", " << suffix << "): " << rResult << std::endl;
        result += rResult.size();
        LOG("VERBOSE_1") << "lb result: " << result << std::endl;
        if(find(dReachableStates.begin(), dReachableStates.end(), state) != dReachableStates.end()) {
            ++result;
            LOG("VERBOSE_1") << "State " << state->getName() << " is d-reachable. Incrementing." << std::endl;
            LOG("VERBOSE_1") << "lb result: " << result << std::endl;
        }

        IOTraceContainer rPlusResult = spec.rPlus(state, base, suffix, vDoublePrime, true);
        rPlusResult.add(rResult);
        LOG("VERBOSE_1") << "rPlusResult: " << rPlusResult << std::endl;
        for (auto traceIt = rPlusResult.cbegin(); traceIt != rPlusResult.cend(); ++traceIt)
        {
            const shared_ptr<const IOTrace>& trace = *traceIt;
            IOTraceContainer traces = iut.bOmega(adaptiveTestCases, *trace);
            LOG("VERBOSE_1") << "Removing " << traces << " from testTraces." << std::endl;

            IOTraceContainer::remove(bOmegaT, traces);

            LOG("VERBOSE_1") << "testTraces:" << std::endl;
            for (const auto& cont : bOmegaT)
            {
                LOG("VERBOSE_1") << "  " << cont << std::endl;
            }
        }
    }
    LOG("VERBOSE_1") << "bOmegaT size: " << bOmegaT.size() << std::endl;
    LOG("VERBOSE_1") << "bOmegaT:" << std::endl;
    for (const auto& cont : bOmegaT)
    {
        LOG("VERBOSE_1") << "  " << cont << std::endl;
    }
    result += bOmegaT.size();
    LOG("VERBOSE_1") << "lowerBound() result: " << result << std::endl;
    return result;
}

bool Fsm::adaptiveStateCounting(Fsm& spec, Fsm& iut, const size_t m,
                                IOTraceContainer& observedTraces,
                                shared_ptr<IOTrace>& failTrace,
                                int& iterations)
{
    LOG("VERBOSE_1")<< "adaptiveStateCounting()" << std::endl;
    if (spec.isMinimal() != True)
    {
        stringstream ss;
ss << "Please ensure to minimize the specification before starting adaptive state counting.";
std::cerr << ss.str();
throw ss.str();
    }
    if (iut.isMinimal() != True)
    {
        stringstream ss;
ss << "Please ensure to minimize the IUT before starting adaptive state counting.";
std::cerr << ss.str();
throw ss.str();
    }
#ifdef ENABLE_DEBUG_MACRO

    const string dotPrefix = "../../../resources/adaptive-test/" + spec.getName() + "-";

#endif
    spec.calcRDistinguishableStates();
    IOListContainer rCharacterisationSet = spec.getRCharacterisationSet();
    LOG("VERBOSE_1") << "Spec rCharacterisationSet:" << rCharacterisationSet << std::endl;

    observedTraces.clear();
    LOG("INFO") << "m: " << m << std::endl;
    /**
     * Adaptive test cases (Ω) for the specification FSM.
     */
    const IOTreeContainer& adaptiveTestCases = spec.getAdaptiveRCharacterisationSet();
    LOG("INFO") << "adaptiveTestCases: " << adaptiveTestCases << std::endl;
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    LOG("INFO") << "adaptiveTestCases as input traces:" << std::endl;
    LOG("INFO") << adaptiveList << std::endl;
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = spec.getMaximalSetsOfRDistinguishableStates();
    LOG("VERBOSE_1") << "maximalSetsOfRDistinguishableStates:" << std::endl;
    for (auto v : maximalSetsOfRDistinguishableStates)
    {
        stringstream ss;
        ss << "{";
        for (auto e : v)
        {
            ss << e->getName() << ", ";
        }
        ss << "}";
        LOG("VERBOSE_1") << ss.str() << std::endl;
    }

    InputTraceSet detStateCover;
    const vector<shared_ptr<FsmNode>>& dReachableStates = spec.calcDReachableStates(detStateCover);

    LOG("INFO") << "dReachableStates:" << std::endl;
    for (auto s : dReachableStates)
    {
        LOG("INFO") << s->getName() << std::endl;
    }
    LOG("INFO") << "detStateCover:" << std::endl;
    for (auto t : detStateCover)
    {
        LOG("INFO") << *t << std::endl;
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
        LOG("VERBOSE_1") << ss.str() << std::endl;
        ss.str(std::string());
        ss << "t: ";
        for (auto w : t)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << std::endl;
        ss.str(std::string());
#endif
        LOG("VERBOSE_1") << "adaptiveTestCases as input traces:" << std::endl;
        LOG("VERBOSE_1") << adaptiveList << std::endl;
        map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;
        size_t numberInputTraces = tC.size();
        size_t inputTraceCount = 0;
        // Applying all input traces from T_c to this FSM.
        // All observed outputs are bein recorded.
        // If the FSM observes a failure, adaptive state counting terminates.
        for (const shared_ptr<InputTrace>& inputTrace : tC)
        {
            LOG("VERBOSE_1") << "############################################################" << std::endl;
            LOG("VERBOSE_1") << "  Applying inputTrace " << ++inputTraceCount << " of " << numberInputTraces << ": " << *inputTrace << std::endl;
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
            LOG("VERBOSE_1") << ss.str() << std::endl;
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
            LOG("VERBOSE_1") << ss.str() << std::endl;
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
            LOG("VERBOSE_1") << ss.str() << std::endl;
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
            LOG("VERBOSE_1") << ss.str() << std::endl;
            ss.str(std::string());
#endif
            observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputsIut));

            for (const shared_ptr<OutputTrace>& oTrace : producedOutputsIut)
            {
                observedTraces.add(make_shared<const IOTrace>(*inputTrace, *oTrace));
            }

            LOG("VERBOSE_1") << "Checking produced outputs for failures" << std::endl;
            //Chek if the IUT has produced any output that can not be produced by the specification.
            for (size_t i = 0; i < producedOutputsIut.size(); ++i)
            {
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
                            LOG("VERBOSE_1") << "----------------- Getting adaptive traces -----------------" << std::endl;
                            IOTraceContainer observedAdaptiveTracesIut;
                            IOTraceContainer observedAdaptiveTracesSpec;
                            const shared_ptr<FsmNode>& nodeIut = reachedNodesIut.at(i);
                            const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

                            iut.addPossibleIOTraces(nodeIut, adaptiveTestCases, observedAdaptiveTracesIut);
                            spec.addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

                            LOG("VERBOSE_1") << "  observedAdaptiveTracesIut (" << nodeIut->getName() << "): " << observedAdaptiveTracesIut << std::endl;
                            LOG("VERBOSE_1") << "  observedAdaptiveTracesSpec (" << nodeSpec->getName() << "): " << observedAdaptiveTracesSpec << std::endl;

                            bool failure = false;
                            for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
                            {
                                const shared_ptr<const IOTrace>& trace = *traceIt;
                                if (!observedAdaptiveTracesSpec.contains(trace))
                                {
                                    LOG("INFO") << "  Specification does not contain " << *trace << std::endl;
                                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                                    IOTrace traceCopy = IOTrace(*trace);
                                    failTrace->append(traceCopy);
                                    LOG("INFO") << "failTrace: " << *failTrace << std::endl;
                                    failure = true;
                                    break;
                                }
                            }
            //                PERFORMANCE_CHECKPOINT_WITH_ID(timerBlkObj, "after observedAdaptiveTracesIut loop");
                            LOG("VERBOSE_1") << "  concatenating: " << *inputTrace << "/" << *outIut << std::endl;
                            observedAdaptiveTracesIut.concatenateToFront(inputTrace, outIut);
                            LOG("VERBOSE_1") << "  observedAdaptiveTraces after concatenation to front: " << observedAdaptiveTracesIut << std::endl;
                            observedTraces.add(observedAdaptiveTracesIut);
                            if (failure)
                            {
                                // IUT produced an output that can not be produced by the specification.
                                LOG("INFO") << "  Failure observed:" << std::endl;
                                LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
                                LOG("INFO") << "    Observed adaptive traces:" << std::endl;
                                LOG("INFO") << observedAdaptiveTracesIut << std::endl;
                                LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
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
                    LOG("INFO") << "  Failure observed:" << std::endl;
                    LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
                    ss << "    Produced Outputs Iut: ";
                    for (size_t i = 0; i < producedOutputsIut.size(); ++i)
                    {
                        ss << *producedOutputsIut.at(i);
                        if (i != producedOutputsIut.size() - 1)
                        {
                            ss << ", ";
                        }
                    }
                    LOG("INFO") << ss.str() << std::endl;
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
                    LOG("VERBOSE_1") << ss.str() << std::endl;
                    ss.str(std::string());
#endif
                    LOG("VERBOSE_1") << "Specification does not produce output " << *outIut << "." << std::endl;
                    LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                    LOG("INFO") << "failTrace: " << *failTrace << std::endl;
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
        LOG("VERBOSE_1") << "observedOutputsTCElements:" << std::endl;
        for (auto e : observedOutputsTCElements)
        {
            LOG("VERBOSE_1") << "  " << *e.first << ":" << std::endl;
            for (auto o : e.second)
            {
                LOG("VERBOSE_1") << "    " << *o << std::endl;
                ++numberToCheck;
            }
        }
        LOG("VERBOSE_1") << "Number of input/output combinations: " << numberToCheck << std::endl;
        InputTraceSet newT = t;
        InputTraceSet newTC;
        inputTraceCount = 0;
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            bool inputTraceMeetsCriteria = true;
            LOG("INFO") << "check inputTrace: " << *inputTrace << " (" << ++inputTraceCount << " of " << numberInputTraces << ")" << std::endl;
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
            LOG("VERBOSE_1") << "producedOutputs:" << std::endl;
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                LOG("VERBOSE_1") << "  " << *outputTrace << std::endl;
            }
            long outputTraceCount = 0;
            size_t numberOutputTraces = producedOutputs.size();

            shared_ptr<const InputTrace> maxInputPrefixInV = nullptr;
            for (const shared_ptr<InputTrace>& detStateTransition : detStateCover)
            {
                if (inputTrace->isPrefix(*detStateTransition, false, true) &&
                        ( !maxInputPrefixInV || maxInputPrefixInV->isEmptyTrace() || detStateTransition->size() > maxInputPrefixInV->size()))
                {
                    maxInputPrefixInV = detStateTransition;
                }
            }
            if (!maxInputPrefixInV)
            {
                stringstream ss;
ss << "No prefix for input trace " << *inputTrace << " found in V. This should not happen.";
std::cerr << ss.str();
throw ss.str();
            }
            LOG("VERBOSE_1") << "maxInputPrefixInV: " << *maxInputPrefixInV << std::endl;

            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                if (!inputTraceMeetsCriteria)
                {
                    break;
                }
                LOG("INFO") << "outputTrace: " << *outputTrace << " (" << ++outputTraceCount << " of " << numberOutputTraces << ")" << std::endl;
                IOTrace currentTrace(*inputTrace, *outputTrace);
                LOG("VERBOSE_1") << "currentTrace (x_1/y_1): " << currentTrace << std::endl;
                bool outputTraceMeetsCriteria = false;
                vPrimeLazy.reset();

                LOG("VERBOSE_1") << "maxInputPrefixInV.size(): " << maxInputPrefixInV->size() << std::endl;
                shared_ptr<const IOTrace> maxIOPrefixInV = make_shared<const IOTrace>(*static_pointer_cast<const Trace>(maxInputPrefixInV),
                                                                                      *outputTrace->getPrefix(maxInputPrefixInV->size(), true));
                LOG("VERBOSE_1") << "maxIOPrefixInV (v/v'): " << *maxIOPrefixInV << std::endl;
                IOTrace suffix(InputTrace(spec.presentationLayer), OutputTrace(spec.presentationLayer));
                suffix = currentTrace.getSuffix(*maxIOPrefixInV);
                LOG("VERBOSE_1") << "suffix (x/y): " << suffix << std::endl;

                LOG("VERBOSE_1") << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext() << std::endl;
                while (vPrimeLazy.hasNext())
                {
                    const IOTraceContainer& vDoublePrime = vPrimeLazy.getNext();
                    if (outputTraceMeetsCriteria)
                    {
                        break;
                    }
                    LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;

                    if (!vDoublePrime.contains(maxIOPrefixInV))
                    {
                        LOG("VERBOSE_1") << "vDoublePrime does not contain prefix " << *maxIOPrefixInV << ". Skipping." << std::endl;
                        LOG("VERBOSE_1") << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext() << std::endl;
                        continue;
                    }
                    for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
                    {
                        LOG("VERBOSE_1") << "rDistStates:" << std::endl;
                        for (auto r : rDistStates)
                        {
                            LOG("VERBOSE_1") << "  " << r->getName() << std::endl;
                        }
                         //size_t lB = Fsm::lowerBound(*maxPrefix, suffix, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut);
                        //LOG("VERBOSE_1") << "lB: " << lB << std::endl;
                        bool exceedsBound = Fsm::exceedsBound(m, *maxIOPrefixInV, suffix, rDistStates, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
                        LOG("VERBOSE_1") << "exceedsBound: " << exceedsBound << std::endl;
                        if (exceedsBound)
                        {
                            LOG("VERBOSE_1") << "Exceeded lower bound. Output trace " << *outputTrace << " meets criteria." << std::endl;
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
                LOG("VERBOSE_1") << "Keeping " << *inputTrace << " in T_C." << std::endl;
                newTC.insert(inputTrace);
                // Next input trace.
                continue;
            }
            else
            {
                LOG("VERBOSE_1") << "Removing " << *inputTrace << " from T_C." << std::endl;
            }
        }
        ss << "newTC: ";
        for (auto w : newTC)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        // Expanding sequences.
        InputTraceSet expandedTC;
        InputTraceSet tracesAddedToT;
        LOG("INFO") << "Expanding input sequences." << std::endl;
        for (int x = 0; x <= spec.maxInput; ++x)
        {
            for (const shared_ptr<InputTrace>& inputTrace : newTC)
            {

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
        LOG("INFO") << "Finished expansion." << std::endl;
        iut.bOmega(adaptiveTestCases, tracesAddedToT, bOmegaT);
        LOG("INFO") << "Finished calculating bOmega." << std::endl;

        ss << "expandedTC: ";
        for (auto w : expandedTC)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        ss << "newT: ";
        for (auto w : newT)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        tC = expandedTC;
        t = newT;
    }
    LOG("VERBOSE_1") << "  RESULT: " << observedTraces << std::endl;
    LOG("VERBOSE_1") << "IUT is a reduction of the specification." << std::endl;
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
                LOG("DEBUG") << nodeA->getName() << " == " << nodeB->getName() << std::endl;
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
    for (size_t i = 0; i < nodesA.size(); ++i)
    {
        shared_ptr<FsmNode> nodeA = nodesA.at(i);
        for (size_t j = i + 1; j < nodesB.size(); ++j)
        {
            shared_ptr<FsmNode> nodeB = nodesB.at(j);
            if (nodeA == nodeB)
            {
                LOG("DEBUG") << nodeA->getName() << " == " << nodeB->getName() << std::endl;
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
    LOG("VERBOSE_2") << "getMaximalSetsOfRDistinguishableStates()" << std::endl;
    vector<vector<shared_ptr<FsmNode>>> result;
    result.reserve(static_cast<size_t>(getMaxNodes()));
    for (shared_ptr<FsmNode> node : nodes)
    {
        LOG("VERBOSE_3") << "Looking for node " << node->getName() << std::endl;
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
            LOG("VERBOSE_3") << "Skipping node " << node->getName() << std::endl;
            continue;
        }
        vector<shared_ptr<FsmNode>> set = {node};
        set.reserve(static_cast<size_t>(getMaxNodes()));
        LOG("VERBOSE_2") << "Creating set for node " << node->getName() << std::endl;
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
        LOG("VERBOSE_2") << "Set size: " << set.size() << std::endl;
        set.resize(set.size());
        result.push_back(set);
    }
    LOG("VERBOSE_2") << "result size: " << result.size() << std::endl;
    result.resize(result.size());
    return result;
}

void Fsm::calcStateIdentificationSets()
{
    if (!isObservable())
    {
        stringstream ss;
ss << "This FSM is not observable - cannot calculate the charactersiation set.";
std::cerr << ss.str();
throw ss.str();
    }
    
    if (characterisationSet == nullptr)
    {
        stringstream ss;
ss << "Missing characterisation set - exit.";
std::cerr << ss.str();
throw ss.str();
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
        stringstream ss;
ss << "This FSM is not observable - cannot calculate the charactersiation set.";
std::cerr << ss.str();
throw ss.str();
    }
    
    if (characterisationSet == nullptr)
    {
        stringstream ss;
ss << "Missing characterisation set - exit.";
std::cerr << ss.str();
throw ss.str();
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
        stringstream ss;
ss << "This FSM is not observable - cannot calculate the harmonized state identification set.";
std::cerr << ss.str();
throw ss.str();
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
                LOG("ERROR")  << "[ERR] Found inconsistency when applying HSI-Method: FSM not minimal." << endl << std::endl;
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
                LOG("INFO") << "Incomplete FSM : for state " << nn->getName() << " (" << nn->getId() << "), input " << x << " does not have a transition." << endl << std::endl;
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
    // Initialisation of random number generation
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG("DEBUG") << "createRandomFsm seed: " << s << std::endl;
    }
    else {
        srand(seed);
        LOG("DEBUG") << "createRandomFsm seed: " << seed << std::endl;
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
    LOG("VERBOSE_1") << "**createRandomFsm()" << std::endl;
    LOG("VERBOSE_1") << "maxInput: " << maxInput << std::endl;
    LOG("VERBOSE_1") << "maxOutput: " << maxOutput << std::endl;
    LOG("VERBOSE_1") << "maxState: " << maxState << std::endl;
    // Initialisation of random number generation
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG("DEBUG") << "createRandomFsm seed: " << s << std::endl;
    }
    else {
        srand(seed);
        LOG("DEBUG") << "createRandomFsm seed: " << seed << std::endl;
    }

    int numIn = maxInput + 1;
    int numOut = maxOutput + 1;
    int numStates = maxState + 1;
    const bool degreeOfCompletenessRequired = degreeOfCompleteness > 0;

    LOG("VERBOSE_1") << "numIn: " << numIn << std::endl;
    LOG("VERBOSE_1") << "numOut: " << numOut << std::endl;
    LOG("VERBOSE_1") << "numStates: " << numStates << std::endl;
    LOG("VERBOSE_1") << "degreeOfCompleteness: " << degreeOfCompleteness << std::endl;
    LOG("VERBOSE_1") << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism << std::endl;
    LOG("VERBOSE_1") << "forceNonDeterminism: " << boolalpha << forceNonDeterminism << std::endl;

    if (forceNonDeterminism && numOut < 2)
    {
        stringstream ss;
ss << "Can not create non-determinism with less than two output symbols.";
std::cerr << ss.str();
throw ss.str();
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
            stringstream ss;
ss << "createRandomFsm(): Could not create requested number of transitions. This shouldn't happen.";
std::cerr << ss.str();
throw ss.str();
        }

        shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
        srcNode->addTransition(transition);

        reachedNodes.push_back(targetNode);
        LOG("VERBOSE_1") << "Created transition " << transition->str() << std::endl;
    }
    LOG("VERBOSE_2") << "Connected all nodes." << std::endl;

    fsm->addRandomTransitions(maxDegreeOfNonDeterminism, false, observable, 1.0f);

    if (degreeOfCompletenessRequired)
    {
        LOG("VERBOSE_2") << "Creating or removing transitions to comply with the given degree of completeness." << std::endl;
        fsm->meetDegreeOfCompleteness(degreeOfCompleteness, maxDegreeOfNonDeterminism, observable);
    }

    if(forceNonDeterminism && fsm->getNumberOfNonDeterministicTransitions() < 1)
    {
        fsm->addRandomTransitions(maxDegreeOfNonDeterminism, true, observable, 1.0f);
    }

    if (minimal)
    {
        LOG("VERBOSE_2") << "Fsm has to be minimal. Minimizing. Num states: " << fsm->size() << std::endl;
        LOG("VERBOSE_2") << *fsm << std::endl;
        Fsm fsmMin = fsm->minimise("", "", false);
        fsmMin.presentationLayer = pl;
        LOG("VERBOSE_2") << "Num states after minimizing: " << fsmMin.size() << std::endl;

        float degreeOfCompletenessMin = fsmMin.getDegreeOfCompleteness();
        size_t numStatesMin = fsmMin.size();

        bool metDegreeOfCompleteness = fsmMin.doesMeetDegreeOfCompleteness(degreeOfCompleteness);
        bool metNumberOfStates = (numStatesMin == static_cast<size_t>(numStates));
        bool metNonDeterminism = (!forceNonDeterminism || fsm->getNumberOfNonDeterministicTransitions() > 0);

        int retryCount = 0;
        LOG("VERBOSE_2") << "metDegreeOfCompleteness: " << std::boolalpha << metDegreeOfCompleteness << std::endl;
        LOG("VERBOSE_2") << "metnumberOfStates: " << std::boolalpha << metNumberOfStates << std::endl;
        LOG("VERBOSE_2") << "metNonDeterminism: " << std::boolalpha << metNonDeterminism << std::endl;
        while (!metDegreeOfCompleteness || !metNumberOfStates || !metNonDeterminism)
        {
            ++retryCount;
            if (retryCount > 15)
            {
                LOG("WARNING") << "Could not create the requested FSM. Trying new seed." << std::endl;
                const unsigned int newSeed = static_cast<unsigned int>(rand());
                return createRandomFsm(fsmName, maxInput, maxOutput, maxState, pl, degreeOfCompleteness, maxDegreeOfNonDeterminism,
                                       forceNonDeterminism, minimal, observable, newSeed);

            }
            LOG("VERBOSE_2") << "FSM does not meet all criteria yet:" << std::endl;
            LOG("VERBOSE_2") << fsmMin << std::endl;
            if (!metNumberOfStates)
            {
                LOG("VERBOSE_1") << "Minimal FSM does not contain requested number of states: "
                        << numStatesMin << " < " << numStates;
                vector<shared_ptr<FsmNode>> newNodes;
                fsmMin.meetNumberOfStates(maxState, maxDegreeOfNonDeterminism, observable, newNodes);
                fsmMin.addRandomTransitions(maxDegreeOfNonDeterminism, false, observable, 1.0f, newNodes);
            }
            else if (!metDegreeOfCompleteness)
            {
                LOG("VERBOSE_1") << "Minimal FSM does not meet degree of completeness: "
                        << degreeOfCompletenessMin << " != " << degreeOfCompleteness;
                fsmMin.meetDegreeOfCompleteness(degreeOfCompleteness, maxDegreeOfNonDeterminism, observable);
            }
            else if (!metNonDeterminism)
            {
                LOG("VERBOSE_1") << "Minimal FSM is not non-deterministic." << std::endl;
                fsm->addRandomTransitions(maxDegreeOfNonDeterminism, true, observable, 1.0f);
            }

            fsmMin = fsmMin.minimise("", "", false);
            fsmMin.presentationLayer = pl;
            degreeOfCompletenessMin = fsmMin.getDegreeOfCompleteness();
            numStatesMin = fsmMin.size();
            metDegreeOfCompleteness = fsmMin.doesMeetDegreeOfCompleteness(degreeOfCompleteness);
            metNumberOfStates = (numStatesMin == static_cast<size_t>(numStates));
            metNonDeterminism = (!forceNonDeterminism || fsm->getNumberOfNonDeterministicTransitions() > 0);
            LOG("VERBOSE_2") << "metDegreeOfCompleteness: " << std::boolalpha << metDegreeOfCompleteness << std::endl;
            LOG("VERBOSE_2") << "metnumberOfStates: " << std::boolalpha << metNumberOfStates << std::endl;
            LOG("VERBOSE_2") << "metNonDeterminism: " << std::boolalpha << metNonDeterminism << std::endl;
        }
        fsm = make_shared<Fsm>(fsmMin);
    }

    LOG("VERBOSE_1") << "Created FSM with " << fsm->size() << " states and degreeOfCompleteness: "
            << fsm->getDegreeOfCompleteness();

    return fsm;
}


shared_ptr<Fsm> Fsm::createMutant(const std::string & fsmName,
                                  const int numOutputFaults,
                                  const int numTransitionFaults,
                                  const bool keepObservability,
                                  const unsigned seed,
                                  const shared_ptr<FsmPresentationLayer>& pLayer){
    if (keepObservability && !isObservable())
    {
        stringstream ss;
ss << "Can not keep an FSM observable that is not already observable.";
std::cerr << ss.str();
throw ss.str();
    }

    if (numOutputFaults > 0 && maxOutput < 1)
    {
        throw too_many_output_faults("Can not create output faults on FSMs with output alphabet size < 2");
    }

    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG("DEBUG") << "createMutant seed: " << s << std::endl;
    }
    else {
        srand(seed);
        LOG("DEBUG") << "createMutant seed: " << seed << std::endl;
    }

    LOG("DEBUG") << "numOutputFaults: " << numOutputFaults << std::endl;
    LOG("DEBUG") << "numTransitionFaults: " << numTransitionFaults << std::endl;

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
            LOG("VERBOSE_2") << "srcNodeId: " << srcNodeId << std::endl;
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
                LOG("VERBOSE_2") << "  newTgtNodeId: " << newTgtNodeId << std::endl;
                tgtNodeIdsCpy.erase(tgtNodeIt);


                transitions = lst[srcNodeId]->getTransitions();
                while (transitions.size() > 0)
                {
                    std::vector<shared_ptr<FsmTransition>>::iterator transitionIt = transitions.begin() + (rand() % transitions.size());
                    shared_ptr<FsmTransition> tr = *transitionIt;
                    LOG("VERBOSE_2") << "    tr: " << tr->str() << std::endl;
                    transitions.erase(transitionIt);

                    if (find(cantTouchThis.begin(), cantTouchThis.end(), tr) != cantTouchThis.end())
                    {
                        LOG("INFO") << "(Transition fault) Won't touch transition " << tr->str() << std::endl;
                        continue;
                    }

                    if (tr->getTarget()->getId() == static_cast<int>(newTgtNodeId)) {
                        continue;
                    }
                    LOG("INFO") << "Adding transition fault:" << std::endl;
                    LOG("INFO") << "  Old transition: " << tr->str() << std::endl;
                    tr->setTarget(lst[newTgtNodeId]);
                    LOG("INFO") << "  New transition: " << tr->str() << std::endl;
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
        LOG("INFO") << "Could not create all requested transition faults." << std::endl;
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
                    LOG("INFO") << "(Output fault) Won't touch transition " << tr->str() << std::endl;
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
                    LOG("INFO") << "Adding output fault:" << std::endl;
                    LOG("INFO") << "  Old transition: " << tr->str() << std::endl;
                    tr->setLabel(newLbl);
                    LOG("INFO") << "  New transition: " << tr->str() << std::endl;
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
        LOG("INFO") << "Could not create all requested output faults." << std::endl;
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
    LOG("VERBOSE_1") << "**createReduction()" << std::endl;
    if ( seed == 0 ) {
        unsigned int s = getRandomSeed();
        srand(s);
        LOG("DEBUG") << "createReduction seed: " << s << std::endl;
    }
    else {
        srand(seed);
        LOG("DEBUG") << "createReduction seed: " << seed << std::endl;
    }

    LOG("VERBOSE_2") << "Fsm:" << std::endl;
    LOG("VERBOSE_2") << *this << std::endl;

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
        LOG("VERBOSE_1") << "Could not create reduction." << std::endl;
        throw reduction_not_possible("There are no deterministic transitions.");
    }

    bool keepGoing = true;
    while (!nonDetTransitions.empty() && keepGoing)
    {
        size_t idx = static_cast<size_t>(rand()) % nonDetTransitions.size();
        const shared_ptr<FsmTransition>& transition = nonDetTransitions.at(idx);
        LOG("VERBOSE_2") << "Removing transition " << transition->str() << std::endl;
        transition->getSource()->removeTransition(transition);
        ++removedTransitions;

        nonDetTransitions = red->getNonDeterministicTransitions();
        size_t size = nonDetTransitions.size();
        int mod = 10 * static_cast<int>(ceil(size * size / 2.0f));
        LOG("VERBOSE_2") << "size: " << size << std::endl;
        LOG("VERBOSE_2") << "mod: " << mod << std::endl;
        keepGoing = (mod == 0) ? false : (rand() % mod) > 7;
        LOG("VERBOSE_2") << "keepGoing: " << boolalpha << keepGoing << std::endl;
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
    LOG("VERBOSE_2") << "moreTransitionsPossible()" << std::endl;
    const float newDegreeofNonDet = getDegreeOfNonDeterminism(1, nodePool);
    const int notDefDet = getNumberOfNotDefinedDeterministicTransitions();
    const int transPossible = getNumberOfPossibleTransitions(nodePool);
    const int totalDefined = getNumberOfTotalTransitions(nodePool);

    LOG("VERBOSE_2") << "newDegreeofNonDet: " << newDegreeofNonDet << std::endl;
    LOG("VERBOSE_2") << "notDefDet: " << notDefDet << std::endl;
    LOG("VERBOSE_2") << "transPossible: " << transPossible << std::endl;
    LOG("VERBOSE_2") << "totalDefined: " << totalDefined << std::endl;

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
    LOG("VERBOSE_1") << "**addRandomTransitions()" << std::endl;
    LOG("VERBOSE_2") << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism << std::endl;
    LOG("VERBOSE_2") << "onlyNonDeterministic: " << onlyNonDeterministic << std::endl;
    LOG("VERBOSE_2") << "observable: " << observable << std::endl;
    LOG("VERBOSE_2") << "factor: " << factor << std::endl;

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    LOG("VERBOSE_2") << "Add random transitions for nodes" << std::endl;
    for (const shared_ptr<FsmNode>& n : nodePool)
    {
        LOG("VERBOSE_2") << "  " << n->getName() << " (" << n << ")" << std::endl;
    }

    const int numStates = static_cast<int>(nodePool.size());
    int numberOfTransitionsCreated = 0;
    bool keepGoing = moreTransitionsPossible(maxDegreeOfNonDeterminism, onlyNonDeterministic, nodePool);
    bool impossible = false;

    while (keepGoing && !impossible)
    {
        LOG("VERBOSE_2") << "Allowed target nodes:" << std::endl;
        vector<shared_ptr<FsmNode>> allowedTargetNodes;
        for (const shared_ptr<FsmNode>& n : nodes)
        {
            LOG("VERBOSE_2") << "  " << n->getName() << " (" << n << ")" << std::endl;
            allowedTargetNodes.push_back(n);
        }

        shared_ptr<FsmNode> targetNode;
        shared_ptr<FsmNode> srcNode;
        shared_ptr<FsmLabel> label;

        while (!allowedTargetNodes.empty() && (!srcNode || !label))
        {
            size_t targetNodeIndex = static_cast<size_t>(rand()) % (allowedTargetNodes.size());
            targetNode = allowedTargetNodes.at(targetNodeIndex);

            LOG("VERBOSE_2") << "Trying to create transition to target node " << targetNode->getName() << std::endl;

            selectRandomNodeAndCreateLabel(nodePool, maxDegreeOfNonDeterminism, onlyNonDeterministic, observable, srcNode, label);

            // We could not find a source node or a valid label.
            if (!srcNode || !label)
            {
                LOG("VERBOSE_2") << "Could not create transition to target node " << targetNode->getName() << ". Trying next node." << std::endl;
                allowedTargetNodes.erase(allowedTargetNodes.begin()
                                         + static_cast<vector<shared_ptr<FsmNode>>::difference_type>(targetNodeIndex));
                continue;

            }
        }

        if (!srcNode || !label)
        {
            LOG("ERROR") << "Could not create transition." << std::endl;
            impossible = true;
        }
        else
        {
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);
            ++numberOfTransitionsCreated;
            LOG("VERBOSE_1") << "Created transition " << transition->str() << std::endl;
            LOG("VERBOSE_2") << "numberOfTransitionsCreated: " << numberOfTransitionsCreated << std::endl;
        }

        keepGoing = !impossible && moreTransitionsPossible(maxDegreeOfNonDeterminism, onlyNonDeterministic, nodePool);
        if (keepGoing)
        {
            float observableFactor = (observable) ? 1.0f : 1.75f;
            keepGoing = (rand() % static_cast<int>(round((10.0f * numStates * observableFactor * factor)))) >= numStates * 2;
        }
        LOG("VERBOSE_2") << "keepGoing: " << boolalpha << keepGoing << std::endl;
    }

}

bool Fsm::meetDegreeOfCompleteness(const float& degreeOfCompleteness,
                                   const float& maxDegreeOfNonDeterminism,
                                   const bool& observable,
                                   vector<shared_ptr<FsmNode>> nodePool)
{
    LOG("VERBOSE_1") << "**meetDegreeOfCompleteness()" << std::endl;

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    float actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
    LOG("VERBOSE_2") << "actualDegreeOfCompleteness: " << actualDegreeOfCompleteness << std::endl;

    bool metRequirement = false;
    if (actualDegreeOfCompleteness < degreeOfCompleteness)
    {
        LOG("VERBOSE_2") << "Degree of completeness: " << actualDegreeOfCompleteness << " < " << degreeOfCompleteness << std::endl;
        LOG("VERBOSE_2") << "Going to add transitions." << std::endl;
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
                stringstream ss;
ss << "meetDegreeOfCompleteness(): Could not create requested number of transitions. This shouldn't happen.";
std::cerr << ss.str();
throw ss.str();
            }

            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);

            LOG("VERBOSE_1") << "Created transition " << transition->str() << std::endl;
            actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
        }
        metRequirement = true;
    }
    else if (actualDegreeOfCompleteness > degreeOfCompleteness)
    {
        LOG("VERBOSE_2") << "Degree of completeness: " << actualDegreeOfCompleteness << " > " << degreeOfCompleteness << std::endl;
        LOG("VERBOSE_2") << "Going to remove transitions." << std::endl;
        while (actualDegreeOfCompleteness >= degreeOfCompleteness && !metRequirement)
        {

            if (doesMeetDegreeOfCompleteness(degreeOfCompleteness, nodePool))
            {
                LOG("VERBOSE_1") << "Won't remove any transition, as degree of completeness would be too small afterwards." << std::endl;
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
                LOG("VERBOSE_2") << "Selected node " << node->getName() << std::endl;
                vector<shared_ptr<FsmTransition>> detTrans = node->getDeterminisitcTransitions();
                LOG("VERBOSE_2") << "Found " << detTrans.size() << " deterministic transitions." << std::endl;
                if (!detTrans.empty())
                {
                    size_t transIndex = static_cast<size_t>(rand()) % detTrans.size();
                    const shared_ptr<FsmTransition>& tr = detTrans.at(transIndex);
                    LOG("VERBOSE_2") << "Removing transition " << tr->str() << std::endl;
                    if (!node->removeTransition(tr))
                    {
                        stringstream ss;
ss << "Could not remove transition " << tr->str() << " from node " << node->getName();
std::cerr << ss.str();
throw ss.str();
                    }
                    detTrans.erase(detTrans.begin()
                                   + static_cast<vector<shared_ptr<FsmTransition>>::difference_type>(transIndex));
                    removedTransition = true;
                }
                else
                {
                    // No deterministic transition found. Trying next node.
                    LOG("VERBOSE_2") << "No deterministic transitions found. Trying next node." << std::endl;
                    selectFrom.erase(selectFrom.begin()
                                     + static_cast<vector<shared_ptr<FsmTransition>>::difference_type>(nodeIdx));
                }

            }

            if (removedTransition)
            {
                actualDegreeOfCompleteness = getDegreeOfCompleteness(0, nodePool);
                continue;
            }

            LOG("VERBOSE_2") << "Could not find any node with deterministic transitions. "
                    << "Going to remove several transitions instead.";

            selectFrom = nodePool;
            removedTransition = false;
            while (!selectFrom.empty() && !removedTransition)
            {
                nodeIdx = static_cast<size_t>(rand()) % selectFrom.size();
                node = selectFrom.at(nodeIdx);
                LOG("VERBOSE_2") << "Selected node " << node->getName() << std::endl;
                vector<shared_ptr<FsmTransition>> transitions = node->getTransitions();
                LOG("VERBOSE_2") << "Found " << transitions.size() << " transitions." << std::endl;
                if (!transitions.empty())
                {
                    const shared_ptr<FsmTransition>& trans = transitions.at(static_cast<size_t>(rand()) % transitions.size());
                    int input = trans->getLabel()->getInput();
                    LOG("VERBOSE_2") << "Removing all transitions with input " << input << std::endl;

                    vector<shared_ptr<FsmTransition>> keepTrans;

                    for (const shared_ptr<FsmTransition>& t : transitions)
                    {
                        if (t->getLabel()->getInput() != input)
                        {
                            keepTrans.push_back(t);
                        }
                        else
                        {
                            LOG("VERBOSE_2") << "Removing transition " << t->str() << std::endl;
                        }
                    }
                    node->setTransitions(keepTrans);
                    removedTransition = true;
                }
                else
                {
                    // No deterministic transition found. Trying next node.
                    LOG("VERBOSE_2") << "No transitions found. Trying next node." << std::endl;
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
                LOG("ERROR") << "Could not comply with the required degree of completeness." << std::endl;
                break;
            }

        }
        if (!metRequirement)
        {
            LOG("ERROR") << "Could not comply with the required degree of completeness." << std::endl;
        }
    }
    LOG("VERBOSE_2") << "Finished with new degree of completeness: " << actualDegreeOfCompleteness << std::endl;
    return metRequirement;
}

bool Fsm::doesMeetDegreeOfCompleteness(const float& degreeOfCompleteness, vector<shared_ptr<FsmNode>> nodePool) const
{
    LOG("VERBOSE_2") << "doesMeetDegreeOfCompleteness()" << std::endl;

    if (degreeOfCompleteness <= 0)
    {
        return true;
    }

    if (nodePool.empty())
    {
        nodePool = nodes;
    }

    const float current = getDegreeOfCompleteness(0, nodePool);
    LOG("VERBOSE_2") << "Current degree of completeness: " << current << std::endl;

    bool yes;
    if (current >= degreeOfCompleteness)
    {
        float newDegreeOfCompleteness = getDegreeOfCompleteness(1, nodePool);
        yes = current >= degreeOfCompleteness && newDegreeOfCompleteness < degreeOfCompleteness;
        if (newDegreeOfCompleteness < degreeOfCompleteness)
        {
            LOG("VERBOSE_2") << "Fsm does meet degree of completeness, as removing one transition would "
                    << "reduce degree to " << newDegreeOfCompleteness << " < " << degreeOfCompleteness;
        }
    }
    else
    {
        LOG("VERBOSE_2") << "Fsm does not meet degree of completeness." << std::endl;
        yes = false;
    }
    return yes;
}

void Fsm::meetNumberOfStates(const int& maxState,
                             const float& maxDegreeOfNonDeterminism,
                             const bool& observable,
                             vector<shared_ptr<FsmNode>>& createdNodes)
{
    LOG("VERBOSE_1") << "**meetNumberOfStates()" << std::endl;
    LOG("VERBOSE_2") << "maxState: " << maxState << std::endl;

    int numIn = maxInput + 1;
    int numOut = maxOutput + 1;
    int numStates = maxState + 1;

    LOG("VERBOSE_2") << "numIn: " << numIn << std::endl;
    LOG("VERBOSE_2") << "numOut: " << numOut << std::endl;
    LOG("VERBOSE_2") << "numStates: " << numStates << std::endl;

    int currentNumberNodes = static_cast<int>(size());
    int missingStates = numStates - currentNumberNodes;
    LOG("VERBOSE_2") << "missingStates: " << missingStates << std::endl;
    LOG("VERBOSE_2") << "currentNumberNodes: " << currentNumberNodes << std::endl;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        LOG("VERBOSE_2") << "  Node at index " << i << " has ID " << nodes.at(i)->getId() << std::endl;
    }

    // Produce the nodes and put them into a vector.
    vector<shared_ptr<FsmNode>> unReachedNodes;
    unReachedNodes.reserve(static_cast<size_t>(missingStates));
    int lowestId = currentNumberNodes;
    int highestId = lowestId + missingStates - 1;
    for (int n = highestId; n >=lowestId; --n) {
        shared_ptr<FsmNode> node = make_shared<FsmNode>(n, name, presentationLayer);
        LOG("VERBOSE_2") << "Created node " << node->getName() << " with id " << n << " (" << node << ")" << std::endl;
        unReachedNodes.push_back(node);
        createdNodes.push_back(node);
    }

    // Connecting all nodes.
    while (unReachedNodes.size() > 0)
    {
        const shared_ptr<FsmNode>& targetNode = unReachedNodes.back();
        LOG("VERBOSE_2") << "targetNode: " << targetNode->getName() << std::endl;

        shared_ptr<FsmNode> srcNode;
        shared_ptr<FsmLabel> label;

        selectRandomNodeAndCreateLabel(nodes, maxDegreeOfNonDeterminism, false, observable, srcNode, label);

        // We could not find a source node or a valid label.
        if (!srcNode || !label)
        {
            if (srcNode)
            {
                LOG("VERBOSE_1") << "Could not create requested number of transitions." << std::endl;
                LOG("VERBOSE_1") << "Going to change the target of an existing one instead" << std::endl;
                const vector<shared_ptr<FsmTransition>>& transitions = srcNode->getTransitions();
                shared_ptr<FsmTransition> transition = transitions.at(static_cast<size_t>(rand()) % transitions.size());
                LOG("VERBOSE_2") << "Selected transition: " << transition->str() << std::endl;
                LOG("VERBOSE_2") << "Replacing target node " << transition->getTarget()->getName() << " with node " << targetNode->getName() << std::endl;
                transition->setTarget(targetNode);
                nodes.push_back(targetNode);
                LOG("VERBOSE_2") << "Modified transition: " << transition->str() << std::endl;
                unReachedNodes.pop_back();
            }
            else
            {
                stringstream ss;
ss << "meetNumberOfStates(): Could not create requested number of transitions. This shouldn't happen.";
std::cerr << ss.str();
throw ss.str();
            }
        }
        else
        {
            shared_ptr<FsmTransition> transition = make_shared<FsmTransition>(srcNode, targetNode, label);
            srcNode->addTransition(transition);
            nodes.push_back(targetNode);
            unReachedNodes.pop_back();
            LOG("VERBOSE_1") << "Created transition " << transition->str() << std::endl;
        }
    }

    LOG("VERBOSE_2") << "Checking node IDs:" << std::endl;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        LOG("VERBOSE_2") << "Node at index " << i << " has ID " << nodes.at(i)->getId() << std::endl;
        if (i != static_cast<size_t>(nodes.at(i)->getId()))
        {
            stringstream ss;
            ss << "Node at index " << i << " has ID " << nodes.at(i)->getId() << ". ";
            ss << "This is an invalid internal state and should not happen!.";
            LOG("FATAL") << ss.str() << std::endl;
            throw ss.str();
        }
    }

    LOG("VERBOSE_2") << "Connected all nodes." << std::endl;
}

shared_ptr<FsmLabel> Fsm::createRandomLabel(const shared_ptr<FsmNode>& srcNode,
                                            const float& maxDegreeOfNonDeterminism,
                                            const bool& onlyNonDeterministic,
                                            const bool& observable) const
{
    LOG("VERBOSE_2") << "createRandomLabel()" << std::endl;

    const int numIn = maxInput + 1;
    const int numOut = maxOutput + 1;

    const bool couldAddMoreNonDet = getDegreeOfNonDeterminism(1, nodes) <= maxDegreeOfNonDeterminism;
    LOG("VERBOSE_2") << "couldAddMoreNonDet: " << boolalpha << couldAddMoreNonDet << std::endl;

    shared_ptr<FsmLabel> label;

    if (numIn == 0 || numOut == 0)
    {
        // We can't create a label without an input or an output.
        LOG("VERBOSE_2") << "No input or output allowed." << std::endl;
        return label;
    }

    if (onlyNonDeterministic && maxDegreeOfNonDeterminism <= 0)
    {
        stringstream ss;
ss << "Invalid choice of parameters.";
std::cerr << ss.str();
throw ss.str();
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
                stringstream ss;
                ss << "Requested to create only non-determinisitc lacels, ";
                ss << "but the degree of non-determinism is already too high.";
                LOG("FATAL") << ss.str() << std::endl;
                throw ss.str();
            }
            // Allow only inputs that are not defined in the source node,
            // but allow every output.
            // Observability isn't an issue in this case, since we choose only
            // inputs that are not yet defined.
            LOG("VERBOSE_2") << "Use only inputs that are not yet defined in node " << srcNode->getName() << std::endl;
            allowedInputs = srcNode->getNotDefinedInputs(maxInput);

            if (allowedInputs.empty())
            {
                LOG("VERBOSE_2") << "No input allowed. Impossible to create label." << std::endl;
                // The source node has no input left under the given circumstances.
                impossible = true;
                break;
            }

            // Allow every output.
            LOG("VERBOSE_2") << "All outputs allowed:" << std::endl;
            for (int o = 0; o < numOut; ++o)
            {
                LOG("VERBOSE_2") << "  " << presentationLayer->getOutId(static_cast<unsigned int>(o)) << std::endl;
                allowedOutputs.push_back(o);
            }
            int input = allowedInputs.at(static_cast<size_t>(rand()) % allowedInputs.size());
            int output = allowedOutputs.at(static_cast<size_t>(rand()) % allowedOutputs.size());
            LOG("VERBOSE_2") << "Selected input: " << presentationLayer->getInId(static_cast<unsigned int>(input)) << std::endl;
            LOG("VERBOSE_2") << "Selected output: " << presentationLayer->getOutId(static_cast<unsigned int>(output)) << std::endl;
            label = make_shared<FsmLabel>(input, output, presentationLayer);
        }
        else
        {
            // We can still create non-deterministic transitions.
            LOG("VERBOSE_2") << "Non-determinsism allowed. Allowed inputs:" << std::endl;
            for (int i = 0; i < numIn; ++i)
            {
                // Check if we have to create non-deterministic transitions only.
                if (!onlyNonDeterministic || srcNode->hasTransition(i))
                {
                    LOG("VERBOSE_2") << "  " << presentationLayer->getInId(static_cast<unsigned int>(i)) << std::endl;
                    allowedInputs.push_back(i);
                }
            }
            if (observable)
            {
                // But we have to stay observable. Therefore we have to pick an
                // input and see, if there is any non-defined output left for
                // that input.
                LOG("VERBOSE_2") << "Fsm has to be observable." << std::endl;
                while (!allowedInputs.empty())
                {
                    LOG("VERBOSE_2") << "Still inputs left." << std::endl;
                    size_t inputIndex = static_cast<size_t>(rand()) % allowedInputs.size();
                    int input = allowedInputs.at(inputIndex);
                    LOG("VERBOSE_2") << "Getting allowed outputs for input "
                            << presentationLayer->getInId(static_cast<unsigned int>(input));
                    allowedOutputs = srcNode->getNotDefinedOutputs(input, maxOutput);

                    if (allowedOutputs.empty())
                    {
                        // There are no more outputs left. Remove selected input from
                        // allowed inputs and try again.
                        LOG("VERBOSE_2") << "No outputs allowed for the given input. Trying next input." << std::endl;
                        allowedInputs.erase(allowedInputs.begin()
                                            + static_cast<vector<shared_ptr<int>>::difference_type>(inputIndex));
                        continue;
                    }
                    else
                    {
                        int output = allowedOutputs.at(static_cast<size_t>(rand()) % allowedOutputs.size());
                        LOG("VERBOSE_2") << "Selected output: " << presentationLayer->getOutId(static_cast<unsigned int>(output)) << std::endl;
                        label = make_shared<FsmLabel>(input, output, presentationLayer);
                        break;
                    }
                }
                if (allowedInputs.empty())
                {
                    // The source node has no input left under the given circumstances.
                    LOG("VERBOSE_2") << "No input allowed. Impossible to create label." << std::endl;
                    impossible = true;
                    break;
                }
            }
            else
            {
                // We can have non-determinism and the FSM does not have to
                // be observable. Any output is allowed.
                LOG("VERBOSE_2") << "No need for observability. All outputs allowed:" << std::endl;
                for (int o = 0; o < numOut; ++o)
                {
                    LOG("VERBOSE_2") << "  " << presentationLayer->getOutId(static_cast<unsigned int>(o)) << std::endl;
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
        LOG("VERBOSE_2") << "Created label: " << *label << std::endl;
    }
    else
    {
        LOG("VERBOSE_2") << "Could not create a label." << std::endl;
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
    LOG("VERBOSE_2") << "selectRandomNodeAndCreateLabel()" << std::endl;

    node = nullptr;
    label = nullptr;

    // Select a reached node at random.
    LOG("VERBOSE_2") << "Trying to find a source node. Allowed:" << std::endl;
    vector<shared_ptr<FsmNode>> allowedSourceNodes;
    for (const shared_ptr<FsmNode>& n : srcNodePool)
    {
        if (!onlyNonDeterministic || !n->getTransitions().empty())
        {
            allowedSourceNodes.push_back(n);
            LOG("VERBOSE_2") << "  " << n->getName() << std::endl;
        }
    }
    while ((!node || !label) && !allowedSourceNodes.empty())
    {
        size_t srcNodeIndex = static_cast<size_t>(rand()) % (allowedSourceNodes.size());
        node = allowedSourceNodes.at(srcNodeIndex);
        LOG("VERBOSE_2") << "Trying node " << node->getName() << std::endl;
        label = createRandomLabel(node, maxDegreeOfNonDeterminism, onlyNonDeterministic, observable);
        if (!label)
        {
            // No label found. We have to try another source node
            LOG("VERBOSE_2") << "Could not find a label. Discarding node " << node->getName() << std::endl;
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
    LOG("VERBOSE_1") << "removeUnreachableNodes()" << std::endl;
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
            LOG("VERBOSE_1") << "Removing node " << oldNames.at(n->getId()) << " (" << n->getId() << ", " << n << ")." << std::endl;
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



bool Fsm::distinguishable(const FsmNode& s1, const FsmNode& s2) {
    
    if ( ofsmTableLst.empty() ) {
        calcOFSMTables();
    }
    
    shared_ptr<OFSMTable> p = ofsmTableLst.back();
    
    return ( p->getS2C().at(s1.getId()) != p->getS2C().at(s2.getId()) );
    
}




