/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <chrono>

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

shared_ptr<FsmNode> Fsm::newNode(const int id, const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p,
                                 shared_ptr<FsmPresentationLayer> pl)
{
    shared_ptr<FsmNode> n = make_shared<FsmNode>(id, pl->getStateId(id,""), pl);
    n->setPair(p);
    return n;
}

bool Fsm::contains(const deque<shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>>>& lst,
                   const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p)
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

shared_ptr<FsmNode> Fsm::findp(const vector<shared_ptr<FsmNode>>& lst,
                               const shared_ptr<pair<shared_ptr<FsmNode>, shared_ptr<FsmNode>>> p)
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

int Fsm::getMaxState() const {
	return maxState;
}

Fsm Fsm::intersect(const Fsm & f)
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
    
    //newPl->dumpState(cout);  commented because of test outputs
    
    return Fsm(f.getName(), maxInput, maxOutput, fsmInterNodes, newPl);
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
    return Fsm(name + "_O", maxInput, maxOutput, nodeLst, obsPl);
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

Fsm Fsm::minimiseObservableFSM()
{
    calcOFSMTables();

    // The last OFSMTable defined has classes corresponding to
    // the minimised FSM to be constructed*/
    shared_ptr<OFSMTable> tbl = ofsmTableLst.back();
    
    // Create the minimised FSM from tbl and return it
    Fsm fsm = tbl->toFsm(name + "_MIN");
    fsm.minimal = True;
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
    
	return createRandomFsmRepeatable(fsmName, maxInput, maxOutput, maxState, pl);
    
}

std::shared_ptr<Fsm>
Fsm::createRandomFsmRepeatable(const std::string & fsmName,
	                           const int maxInput,
	                           const int maxOutput,
	                           const int maxState,
	                           const std::shared_ptr<FsmPresentationLayer> pl) {
	// Produce the nodes and put them into a vector.
    // All nodes are marked 'white' by the costructor - this is now
    // used to mark unreachable states which have to be made reachable
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= maxState; n++) {
		lst.push_back(make_shared<FsmNode>(n, fsmName, pl));
	}

	// At index 0 of the vector, the initial state is store, and
	// this is reachable
	lst[0]->setColor(FsmNode::black);

	// We create transitions by starting from black (reachable) nodes
	// and trying to reach at least one white node from there.
	deque< shared_ptr<FsmNode> > bfsq;
	bfsq.push_back(lst[0]);

	while (not bfsq.empty()) {

		shared_ptr<FsmNode> srcNode = bfsq.front();
		bfsq.pop_front();

		// Generation part 1.
		// Select an uncovered node at random
		int whiteNodeIndex = rand() % (maxState + 1);
		shared_ptr<FsmNode> whiteNode = nullptr;
		shared_ptr<FsmNode> startNode = lst[whiteNodeIndex];
		shared_ptr<FsmNode> thisNode = startNode;

		do {

			if (thisNode->getColor() == FsmNode::white) {
				whiteNode = thisNode;
			}
			else {
				whiteNodeIndex = (whiteNodeIndex + 1) % (maxState + 1);
				thisNode = lst[whiteNodeIndex];
			}

		} while (whiteNode == nullptr and thisNode != startNode);

		// Link srcNode by random transition to thisNode
		// and mark thisNode as black. Also insert into BFS queue
		int x0 = -1;
		int y0 = -1;

		if (whiteNode != nullptr) {
			x0 = rand() % (maxInput + 1);
			y0 = rand() % (maxOutput + 1);
			auto theTrans =
				make_shared<FsmTransition>(srcNode, whiteNode,
					make_shared<FsmLabel>(x0, y0, pl));
			// Add transition to adjacency list of the source node
			srcNode->addTransition(theTrans);
			thisNode->setColor(FsmNode::black);
			bfsq.push_back(thisNode);
		}

		// Generation part 2.
		// Random transition generation.
		// At least one transition for every input, with
		// arbitrary target nodes.
		for (int x = 0; x <= maxInput; x++) {
			// If x equals x0 produced already above,
			// we may skip it at random
			if (x == x0 and (rand() % 2)) continue;

			// How many transitions do we want for input x?
			// We construct at most 2 of these transitions
			int numTrans = rand() % 2;
			for (int t = 0; t <= numTrans; t++) {
				// Which output do we want?
				int y = rand() % (maxOutput + 1);
				// Which target node?
				int tgtNodeId = rand() % (maxState + 1);
				auto tgtNode = lst[tgtNodeId];
				if (tgtNode->getColor() == FsmNode::white) {
					tgtNode->setColor(FsmNode::black);
					bfsq.push_back(tgtNode);
				}
				auto theTrans =
					make_shared<FsmTransition>(srcNode, tgtNode,
						make_shared<FsmLabel>(x, y, pl));
				// Add transition to adjacency list of the source node
				srcNode->addTransition(theTrans);
			}

		}

	}

	return make_shared<Fsm>(fsmName, maxInput, maxOutput, lst, pl);
}


shared_ptr<Fsm> Fsm::createMutant(const std::string & fsmName,
                                  const size_t numOutputFaults,
                                  const size_t numTransitionFaults){
    
    srand(getRandomSeed());
	return createMutantRepeatable(fsmName, numOutputFaults, numTransitionFaults); 
}

/*
  Original
*/
//std::shared_ptr<Fsm> Fsm::createMutantRepeatable(const std::string & fsmName,
//	const size_t numOutputFaults,
//	const size_t numTransitionFaults) {
//	// Create new nodes for the mutant.
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= maxState; n++) {
//		lst.push_back(make_shared<FsmNode>(n, fsmName, presentationLayer));
//	}
//
//	// Now add transitions that correspond exactly to the transitions in
//	// this FSM
//	for (int n = 0; n <= maxState; n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = nodes[n];
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//	}
//
//	// Now add transition faults to the new machine
//	for (size_t tf = 0; tf < numTransitionFaults; tf++) {
//		int srcNodeId = rand() % (maxState + 1);
//		int newTgtNodeId = rand() % (maxState + 1);
//		int trNo = rand() % lst[srcNodeId]->getTransitions().size();
//		auto tr = lst[srcNodeId]->getTransitions()[trNo];
//		if (tr->getTarget()->getId() == newTgtNodeId) {
//			newTgtNodeId = (newTgtNodeId + 1) % (maxState + 1);
//		}
//		lst[srcNodeId]->getTransitions()[trNo]->setTarget(lst[newTgtNodeId]);
//	}
//
//	// Now add output faults to the new machine
//	for (size_t of = 0; of < numOutputFaults; of++) {
//
//		int srcNodeId = rand() % (maxState + 1);
//		int trNo = rand() % lst[srcNodeId]->getTransitions().size();
//		auto tr = lst[srcNodeId]->getTransitions()[trNo];
//		int theInput = tr->getLabel()->getInput();
//		int newOutVal = rand() % (maxOutput + 1);
//		int originalNewOutVal = rand() % (maxOutput + 1);
//		bool newOutValOk;
//
//		// We don't want to modify this transition in such a way
//		// that another one with the same label and the same
//		// source/target nodes already exists.
//		do {
//
//			newOutValOk = true;
//
//			for (auto trOther : lst[srcNodeId]->getTransitions()) {
//				if (tr == trOther) continue;
//				if (trOther->getTarget()->getId() != tr->getTarget()->getId())
//					continue;
//				if (trOther->getLabel()->getInput() != theInput) continue;
//				if (trOther->getLabel()->getOutput() == newOutVal) {
//					newOutValOk = false;
//				}
//			}
//
//			if (not newOutValOk) {
//				newOutVal = (newOutVal + 1) % (maxOutput + 1);
//			}
//
//		} while ((not newOutValOk) and (originalNewOutVal != newOutVal));
//
//		if (newOutValOk) {
//
//			auto newLbl = make_shared<FsmLabel>(tr->getLabel()->getInput(),
//
//				newOutVal,
//				presentationLayer);
//
//			tr->setLabel(newLbl);
//		}
//	}
//	return make_shared<Fsm>(fsmName, maxInput, maxOutput, lst, presentationLayer);
//}

/*
	Corrected
*/
std::shared_ptr<Fsm> Fsm::createMutantRepeatable(const std::string & fsmName,
	                                             const size_t numOutputFaults,
	                                             const size_t numTransitionFaults) {
	// Create new nodes for the mutant.
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= maxState; n++) {
		lst.push_back(make_shared<FsmNode>(n, fsmName, presentationLayer));
	}

	// Now add transitions that correspond exactly to the transitions in
	// this FSM
	for (int n = 0; n <= maxState; n++) {
		auto theNewFsmNodeSrc = lst[n];
		auto theOldFsmNodeSrc = nodes[n];
		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
			int tgtId = tr->getTarget()->getId();
			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
			shared_ptr<FsmTransition> newTr =
				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
			theNewFsmNodeSrc->addTransition(newTr);
		}
	}

	// Now add transition faults to the new machine
	for (size_t tf = 0; tf < numTransitionFaults; tf++) {
		int srcNodeId = rand() % (maxState + 1);
		int newTgtNodeId = rand() % (maxState + 1);
		if (lst[srcNodeId]->getTransitions().empty()) continue; // Fix possible divide by zero error
		int trNo = rand() % lst[srcNodeId]->getTransitions().size();
		auto tr = lst[srcNodeId]->getTransitions()[trNo];
		if (tr->getTarget()->getId() == newTgtNodeId) {
			newTgtNodeId = (newTgtNodeId + 1) % (maxState + 1);
		}
		lst[srcNodeId]->getTransitions()[trNo]->setTarget(lst[newTgtNodeId]);
	}

	// Now add output faults to the new machine
	for (size_t of = 0; of < numOutputFaults; of++) {

		int srcNodeId = rand() % (maxState + 1);
		if (lst[srcNodeId]->getTransitions().empty()) continue; // Fix possible divide by zero error
		int trNo = rand() % lst[srcNodeId]->getTransitions().size();
		auto tr = lst[srcNodeId]->getTransitions()[trNo];
		int theInput = tr->getLabel()->getInput();
		int newOutVal = rand() % (maxOutput + 1);
		int originalNewOutVal = rand() % (maxOutput + 1);
		bool newOutValOk;

		// We don't want to modify this transition in such a way
		// that another one with the same label and the same
		// source/target nodes already exists.
		do {

			newOutValOk = true;

			for (auto trOther : lst[srcNodeId]->getTransitions()) {
				if (tr == trOther) continue;
				if (trOther->getTarget()->getId() != tr->getTarget()->getId())
					continue;
				if (trOther->getLabel()->getInput() != theInput) continue;
				if (trOther->getLabel()->getOutput() == newOutVal) {
					newOutValOk = false;
				}
			}

			if (not newOutValOk) {
				newOutVal = (newOutVal + 1) % (maxOutput + 1);
			}

		} while ((not newOutValOk) and (originalNewOutVal != newOutVal));

		if (newOutValOk) {

			auto newLbl = make_shared<FsmLabel>(tr->getLabel()->getInput(),

				newOutVal,
				presentationLayer);

			tr->setLabel(newLbl);
		}
	}
	return make_shared<Fsm>(fsmName, maxInput, maxOutput, lst, presentationLayer);
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


// Original
//bool Fsm::removeUnreachableNodes(std::vector<shared_ptr<FsmNode>>& unreachableNodes) {
//    
//    vector<shared_ptr<FsmNode>> newNodes;
//    FsmVisitor v;
//    
//    // Mark all reachable nodes as 'visited'
//    accept(v);
//    
//    // When removing nodes from the FSM, the node ids of all remaining nodes
//    // have to be adapted, in order to match the index in the list of nodes.
//    // This is necessary, because during minimisation with OFSM tables or
//    // Pk-tables, the algorithms rely on the range of row numbers being
//    // identical to the range of node ids of the reachable nodes.
//    int subtractFromId = 0;
//    for ( auto n : nodes ) {
//        if ( not n->hasBeenVisited() ) {
//            unreachableNodes.push_back(n);
//            presentationLayer->removeState2String(n->getId() - subtractFromId);
//            ++subtractFromId;
//        }
//        else {
//            n->setId(n->getId() - subtractFromId);
//            newNodes.push_back(n);
//        }
//    }
//    
//    nodes = newNodes;
//    
//    return (unreachableNodes.size() > 0);
//}

// Corrected
bool Fsm::removeUnreachableNodes(std::vector<shared_ptr<FsmNode>>& unreachableNodes) {

	vector<shared_ptr<FsmNode>> newNodes;
	FsmVisitor v;

	// Mark all reachable nodes as 'visited'
	accept(v);

	// get the number of unreachable states in the nodes list with smaller
	// ids than initStateIdx
	size_t smallerIDs = 0;
	for (size_t id = 0; id < initStateIdx; ++id) {
		if (not nodes.at(id)->hasBeenVisited()) ++smallerIDs;
	}

	// When removing nodes from the FSM, the node ids of all remaining nodes
	// have to be adapted, in order to match the index in the list of nodes.
	// This is necessary, because during minimisation with OFSM tables or
	// Pk-tables, the algorithms rely on the range of row numbers being
	// identical to the range of node ids of the reachable nodes.
	int subtractFromId = 0;
	for (auto n : nodes) {
		if (not n->hasBeenVisited()) {
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

	maxState = nodes.size() - 1;
	initStateIdx = initStateIdx - smallerIDs;

	return (unreachableNodes.size() > 0);
}



bool Fsm::distinguishable(const FsmNode& s1, const FsmNode& s2) {
    
    if ( ofsmTableLst.empty() ) {
        calcOFSMTables();
    }
    
    shared_ptr<OFSMTable> p = ofsmTableLst.back();
    
    return ( p->getS2C().at(s1.getId()) != p->getS2C().at(s2.getId()) );
    
}


//------------------------------------------------------------------------------------

bool Fsm::checkNodeIds() const {
	for (size_t i = 0; i < nodes.size(); ++i) {
		if (nodes.at(i)->getId() != i) return false;
	}
	return true;
}

bool Fsm::contains(const shared_ptr<FsmNode> node) const {
	for (auto n : nodes) {
		if (n == node) return true;
	}
	return false;
}

bool Fsm::checkAllTransitions() const {
	for (auto n : nodes) {
		for (auto tr : n->getTransitions()) {
			if (tr == nullptr || tr->getLabel() == nullptr || tr->getLabel()->getInput() > maxInput
				|| tr->getLabel()->getOutput() > maxOutput || tr->getLabel()->getInput() < 0 || tr->getLabel()->getOutput() < 0
				|| tr->getSource() != n
				|| not contains(tr->getTarget())) {
				return false;
			}
		}
	}
	return true;
}

/**
 * Returns true iff the Fsm class invariant holds for this object.
 */
bool Fsm::checkInvariant() const {
	//cout << "Fsm_inv" << endl;
	if (maxInput < 0) return false;
	if (maxOutput < 0) return false;
	if (nodes.size() < 1) return false;
	if (not checkNodeIds()) return false;
	if (contains(nullptr)) return false;
	if (not checkAllTransitions()) return false;
	if (maxState != nodes.size() - 1) return false;
	if (not(0 <= initStateIdx and initStateIdx <= maxState)) return false;
	return true;
}

