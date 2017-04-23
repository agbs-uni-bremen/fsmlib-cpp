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
#include "fsm/OFSMTable.h"
#include "sets/HittingSet.h"
#include "trees/TreeNode.h"
#include "trees/OutputTree.h"
#include "trees/Tree.h"
#include "trees/IOListContainer.h"
#include "trees/TestSuite.h"


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
            InputTrace i = leftNode->calcDistinguishingTrace(rightNode, ofsmTableLst, maxInput, maxOutput);
            shared_ptr<vector<vector<int>>> lli = make_shared<vector<vector<int>>>();
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
    
    return tcl;
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
    
    calcStateIdentificationSets();
    
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
        out << i << "[label=\"" << nodeName << "(" << i << ")\"];" << endl;
        
        if (i == fsm.initStateIdx)
        {
            out << endl << "node [shape = circle]" << endl;
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
            if ( x == x0 and rand() % 2 ) continue;
            
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






