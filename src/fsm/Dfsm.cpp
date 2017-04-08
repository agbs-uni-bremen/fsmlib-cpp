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

using namespace std;

void Dfsm::createAtRandom()
{
    srand((unsigned int) time(0));
    
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        nodes [i] = make_shared<FsmNode>(i, presentationLayer);//insertion
    }
    
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        shared_ptr<FsmNode> source = nodes.at(i);
        
        for (int input = 0; input <= maxInput; ++ input)
        {
            int nTarget = rand() % nodes.size();
            shared_ptr<FsmNode> target = nodes.at(nTarget);
            int output = rand() % (maxOutput + 1);
            shared_ptr<FsmTransition> transition =
            make_shared<FsmTransition>(source,
                                       target,
                                       make_shared<FsmLabel>(input,
                                                             output,
                                                             presentationLayer));
            source->addTransition(transition);
        }
    }
}

shared_ptr<DFSMTable> Dfsm::toDFSMTable() const
{
    shared_ptr<DFSMTable> tbl
    = make_shared<DFSMTable>(nodes.size(), maxInput, presentationLayer);
    
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        if (nodes.at(i) == nullptr)
        {
            continue;
        }
        
        shared_ptr<DFSMTableRow> r = nodes.at(i)->getDFSMTableRow(maxInput);
        
        if (r == nullptr)
        {
            return nullptr;
        }
        tbl->setRow(i, r);
        
        
        
    }
    
    return tbl;
}

Dfsm::Dfsm(const string & fname, const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer> presentationLayer)
: Fsm(fname, fsmName, maxNodes, maxInput, maxOutput, presentationLayer)
{
    
}

Dfsm::Dfsm(const string& fname,
           const shared_ptr<FsmPresentationLayer> presentationLayer,
           const string & fsmName)
: Fsm(fname,presentationLayer,fsmName)
{
    
}

Dfsm::Dfsm(const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer> presentationLayer)
: Fsm(presentationLayer)
{
    name = fsmName;
    nodes.insert(nodes.end(), maxNodes, nullptr);
    initStateIdx = 0;
    this->maxInput = maxInput;
    this->maxOutput = maxOutput;
    currentParsedNode = nullptr;
    createAtRandom();
    ofstream out(getName() + ".txt");
    dumpFsm(out);
    out.close();
}

Dfsm::Dfsm(const string & fsmName, const int maxInput, const int maxOutput, const vector<shared_ptr<FsmNode>> lst, const shared_ptr<FsmPresentationLayer> presentationLayer)
: Fsm(fsmName, maxInput, maxOutput, lst, presentationLayer)
{
    
}

Dfsm::Dfsm(const Fsm & fsm)
: Fsm (fsm.getName(), fsm.getMaxInput(), fsm.getMaxOutput(), fsm.getNodes(), fsm.getPresentationLayer())
{
    initStateIdx = fsm.getInitStateIdx();;
    minimal = isMinimal();
    /*shared_ptr<FsmNode> currentParsedNode;
     vector<shared_ptr<OFSMTable>> ofsmTableLst;
     shared_ptr<Tree> characterisationSet;
     vector<shared_ptr<Tree>> stateIdentificationSets;*/
    
}

Dfsm Dfsm::minimise()
{
    dfsmTable = toDFSMTable();
    
    pktblLst.clear();
    shared_ptr<PkTable> p1 = dfsmTable->getP1Table();
    pktblLst.push_back(p1);
    shared_ptr<PkTable> pMin = p1;
    
    for (shared_ptr<PkTable> pk = p1->getPkPlusOneTable();
         pk != nullptr;
         pk = pk->getPkPlusOneTable())
    {
        pMin = pk;
        pktblLst.push_back(pk);
    }
    
    return pMin->toFsm(name);
}

void Dfsm::printTables() const
{
    ofstream file("tables.tex");
    if (dfsmTable != nullptr)
    {
        file << *dfsmTable;
    }
    
    for (unsigned int i = 0; i < pktblLst.size(); ++ i)
    {
        file << *pktblLst.at(i) << endl << endl;
    }
    file.close();
}

IOListContainer Dfsm::getCharacterisationSet()
{
    /*Create Pk-tables for the minimised FSM*/
    dfsmTable = toDFSMTable();
    pktblLst.clear();
    shared_ptr<PkTable> p1 = dfsmTable->getP1Table();
    pktblLst.push_back(p1);
    
    for (shared_ptr<PkTable> pk = p1->getPkPlusOneTable();
         pk != nullptr;
         pk = pk->getPkPlusOneTable())
    {
        pktblLst.push_back(pk);
    }
    
    /*Create an empty characterisation set as an empty InputTree instance*/
    shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
    
    /*Loop over all non-equal pairs of states. If they are not already distinguished by
     the input sequences contained in w, create a new input traces that distinguishes them
     and add it to w.*/
    for (unsigned int left = 0; left < nodes.size(); ++ left)
    {
        shared_ptr<FsmNode> leftNode = nodes.at(left);
        for (unsigned int right = left + 1; right < nodes.size(); ++ right)
        {
            shared_ptr<FsmNode> rightNode = nodes.at(right);
            
            if (leftNode->distinguished(rightNode, w) != nullptr)
            {
                continue;
            }
            
            /*We have to create a new input trace and add it to w, because
             leftNode and rightNode are not distinguished by the current
             input traces contained in w. This step is performed
             according to Gill's algorithm.*/
            InputTrace i = leftNode->calcDistinguishingTrace(rightNode, pktblLst, maxInput);
            shared_ptr<vector<vector<int>>> lli = make_shared<vector<vector<int>>>();
            lli->push_back(i.get());
            IOListContainer tcli = IOListContainer(lli, presentationLayer);
            w->addToRoot(tcli);
        }
    }
    
    /* Wrap list of lists by an IOListContainer instance */
    IOListContainer tcl = w->getIOLists();
    return tcl;
}

IOTrace Dfsm::applyDet(const InputTrace & i)
{
    OutputTrace o = OutputTrace(presentationLayer);
    
    shared_ptr<FsmNode> currentNode = nodes.at(initStateIdx);
    
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
        
        return IOTrace(InputTrace(vector<int>(ifirst, ilast), presentationLayer),
                       OutputTrace(vector<int>(ofirst, olast), presentationLayer));
        
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



IOListContainer Dfsm::wMethod(const unsigned int m) {
    
    Dfsm dfsmMin = minimise();
    return dfsmMin.wMethodOnMinimisedDfsm(m);
    
}


IOListContainer Dfsm::wMethodOnMinimisedDfsm(const unsigned int m)
{
    
    size_t mMinusN = ( m > nodes.size() ) ? (m - nodes.size()) : 0;
    
    shared_ptr<Tree> iTree = getTransitionCover();
    
    if (mMinusN > 0)
    {
        IOListContainer inputEnum = IOListContainer(maxInput,
                                                    1,
                                                    (int)mMinusN,
                                                    presentationLayer);
        iTree->add(inputEnum);
    }
    
    IOListContainer w = getCharacterisationSet();
    iTree->add(w);
    return iTree->getIOLists();
}

IOListContainer Dfsm::wpMethod(const int m)
{
    Fsm fMin = minimiseObservableFSM();
    return fMin.wpMethod(m);
    
}


IOListContainer Dfsm::tMethod()
{
    
    shared_ptr<Tree> iTree = getTransitionCover();
    
    return iTree->getIOLists();
    
}


void Dfsm::toCsv(const std::string& fname) {
    ofstream out(fname + ".csv");
    
    // Table heading contains input identifiers
    for ( int x = 0; x <= maxInput; x++ ) {
        out << " ; ";
        out << presentationLayer->getInId(x);
    }
    
    for ( size_t n = 0; n < nodes.size(); n++ ) {
        out << endl << "\"" << nodes[n]->getName() << "\"";
        
        for ( int x = 0; x <= maxInput; x++ ) {
            
            out << " ; ";
            
            for ( auto tr : nodes[n]->getTransitions() ) {
                if ( tr->getLabel()->getInput() == x ) {
                    out << "\"" << tr->getTarget()->getName()
                    << " / "
                    << presentationLayer->getOutId(tr->getLabel()->getOutput())
                    << "\"";
                    break;
                }
            }
            
        }
        
    }
    
    out << endl;
    out.close();
}














