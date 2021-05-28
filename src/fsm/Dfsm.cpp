/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <iostream>
#include <set>
#include <fstream>


#include "fsm/Dfsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmLabel.h"
#include "fsm/FsmTransition.h"
#include "fsm/PkTable.h"
#include "fsm/PkTableRow.h"
#include "fsm/DFSMTableRow.h"
#include "fsm/DFSMTable.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "fsm/SegmentedTrace.h"
#include "trees/Tree.h"
#include "trees/IOListContainer.h"
#include "trees/TreeNode.h"
#include "trees/TreeEdge.h"
#include "trees/OutputTree.h"
#include "json/json.h"

using namespace std;


/** key comparator for usage in string sets */
struct setCmp
{
    bool operator()(const string& s1, const string& s2) const
    {
        return (s1.compare(s2) < 0);
    }
};


shared_ptr<FsmPresentationLayer> Dfsm::createPresentationLayerFromCsvFormat(const string & fname) {
    

    // Set containing all output events defined so far
    set<string,setCmp> outStrSet;
    
    
    vector<string> in2String;
    vector<string> out2String;
    vector<string> state2String;
    
    ifstream inputFile(fname);
    
    if ( not inputFile.is_open() ) {
        cout << "Unable to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Insert a "no operation" output as first member of the
    // set to be filled with the output alphabet members
    outStrSet.insert(string("_nop"));
    
    string line;
    
    // ------------------------------------------------------------
    // Read input names from first line
    // ------------------------------------------------------------
    getline(inputFile,line);
    
    // Skip first field
    size_t pos = line.find(";");
    if ( pos == string::npos ) return nullptr;
    pos++;
    size_t posEnd;
    do {
        
        posEnd = line.find(";", pos);
        
        string newInput;
        if ( posEnd == string::npos ) {
            newInput = line.substr(pos);
        }
        else {
            newInput = line.substr(pos, posEnd - pos);
        }
        
        // Trim
        newInput.erase(0,newInput.find_first_not_of(" \n\r\t\""));
        newInput.erase(newInput.find_last_not_of(" \n\r\t\"")+1);
        
        // Add to vector of strings denoting input names
        in2String.push_back(newInput);
        
        pos = posEnd + 1;
        
    } while ( posEnd != string::npos );
    
    // ------------------------------------------------------------
    // Read state names and outputs from consecutive lines
    // ------------------------------------------------------------
    while (getline(inputFile, line))
    {
        // Get the state name
        posEnd = line.find(";");
        if ( posEnd == string::npos ) continue;
        string newState = line.substr(0,posEnd);
        newState.erase(0,newState.find_first_not_of(" \n\r\t\""));
        newState.erase(newState.find_last_not_of(" \n\r\t\"")+1);
        state2String.push_back(newState);
        
        // Look for new output names
        do {
            
            pos = posEnd + 1;
            posEnd = line.find(";",pos);
            
            string newEntry;
            if ( posEnd == string::npos ) {
                newEntry = line.substr(pos);
            }
            else {
                newEntry = line.substr(pos, posEnd - pos);
            }
            size_t outPos = newEntry.find("/");
            if ( outPos != string::npos ) {
                string outStr = newEntry.substr(outPos+1);
                outStr.erase(0,outStr.find_first_not_of(" \n\r\t\""));
                outStr.erase(outStr.find_last_not_of(" \n\r\t\"")+1);
                outStrSet.insert(outStr);
            }
            
        } while ( posEnd != string::npos );
        
    }
    inputFile.close();
    
    for (const auto &s : outStrSet ) {
        out2String.push_back(s);
    }
    
    maxInput = (int)in2String.size() - 1;
    maxOutput = (int)out2String.size() - 1;
    maxState = (int)state2String.size() - 1;
    initStateIdx = 0;
    
    return make_shared<FsmPresentationLayer>(in2String,out2String,state2String);
    
}



shared_ptr<FsmPresentationLayer>
Dfsm::createPresentationLayerFromCsvFormat(const string & fname,
                                           const shared_ptr<FsmPresentationLayer>& pl) {
    
    
    // Set containing all output events defined so far
    set<string,setCmp> outStrSet;
    
    
    vector<string> in2String(pl->getIn2String());
    vector<string> out2String(pl->getOut2String());
    vector<string> state2String;
    
    ifstream inputFile(fname);
    
    if ( not inputFile.is_open() ) {
        cout << "Unable to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Insert a "no operation" output as first member of the
    // output alphabet, if it's not already contained in pl
    if ( pl->out2Num("_nop") < 0 ) {
        outStrSet.insert(string("_nop"));
    }
    
    string line;
    
    // ------------------------------------------------------------
    // Read input names from first line
    // ------------------------------------------------------------
    getline(inputFile,line);
    
    // Skip first field
    size_t pos = line.find(";");
    if ( pos == string::npos ) return nullptr;
    pos++;
    size_t posEnd;
    do {
        
        posEnd = line.find(";", pos);
        
        string newInput;
        if ( posEnd == string::npos ) {
            newInput = line.substr(pos);
        }
        else {
            newInput = line.substr(pos, posEnd - pos);
        }
        
        // Trim
        newInput.erase(0,newInput.find_first_not_of(" \n\r\t\""));
        newInput.erase(newInput.find_last_not_of(" \n\r\t\"")+1);
        
        // Add to vector of strings denoting input names,
        // if this string is not already containe in pl
        if ( pl->in2Num(newInput) < 0 ) {
            in2String.push_back(newInput);
        }
        
        pos = posEnd + 1;
        
    } while ( posEnd != string::npos );
    
    // ------------------------------------------------------------
    // Read state names and outputs from consecutive lines
    // ------------------------------------------------------------
    while (getline(inputFile, line))
    {
        // Get the state name
        posEnd = line.find(";");
        if ( posEnd == string::npos ) continue;
        string newState = line.substr(0,posEnd);
        newState.erase(0,newState.find_first_not_of(" \n\r\t\""));
        newState.erase(newState.find_last_not_of(" \n\r\t\"")+1);
        state2String.push_back(newState);
        
        // Look for new output names
        do {
            
            pos = posEnd + 1;
            posEnd = line.find(";",pos);
            
            string newEntry;
            if ( posEnd == string::npos ) {
                newEntry = line.substr(pos);
            }
            else {
                newEntry = line.substr(pos, posEnd - pos);
            }
            size_t outPos = newEntry.find("/");
            if ( outPos != string::npos ) {
                string outStr = newEntry.substr(outPos+1);
                outStr.erase(0,outStr.find_first_not_of(" \n\r\t\""));
                outStr.erase(outStr.find_last_not_of(" \n\r\t\"")+1);
                
                // Insert to set only if not already contained in pl
                if ( pl->out2Num(outStr) < 0 ) {
                    outStrSet.insert(outStr);
                }
            }
            
        } while ( posEnd != string::npos );
        
    }
    inputFile.close();
    
    for (const auto &s : outStrSet ) {
        out2String.push_back(s);
    }
    
    maxInput = (int)in2String.size() - 1;
    maxOutput = (int)out2String.size() - 1;
    maxState = (int)state2String.size() - 1;
    initStateIdx = 0;
    
    return make_shared<FsmPresentationLayer>(in2String,out2String,state2String);
    
}



void Dfsm::createDfsmTransitionGraph(const string& fname) {
    
    ifstream inputFile(fname);
    size_t pos;
    size_t posEnd;
    string line;
    shared_ptr<FsmNode> tgtNode;
    
    // skip first line
    getline(inputFile,line);
    
    // Initialise nodes-vector with null pointers, so that states
    // not yet defined can be identified.
    for ( size_t n = 0; n <= (size_t)maxState; n++ ) nodes.push_back(nullptr);
    
    
    int nodeId = 0;
    while (getline(inputFile, line)) {
        
        currentParsedNode = nodes[nodeId];
        if ( currentParsedNode == nullptr ) {
            currentParsedNode =
            make_shared<FsmNode>(nodeId,
                                 presentationLayer->getStateId(nodeId,""),
                                 presentationLayer);
            nodes[nodeId] = currentParsedNode;
        }
        
        // Skip the first column
        pos = line.find(";");
        if ( pos++ == string::npos ) continue;
        int x = 0;
        // Create transitions emanating from currentParsedNode
        do {
            
            string tableEntry;
            
            posEnd = line.find(";",pos);
            if ( posEnd == string::npos ) {
                tableEntry = line.substr(pos);
            }
            else {
                tableEntry = line.substr(pos,posEnd-pos);
            }
            
            // Empty table entries lead to an x/_nop self loop
            // _nop has integer code 0
            tableEntry.erase(0,tableEntry.find_first_not_of(" \n\r\t\""));
            tableEntry.erase(tableEntry.find_last_not_of(" \n\r\t\"")+1);
            if ( tableEntry.empty() ) {
                tgtNode = currentParsedNode;
                shared_ptr<FsmLabel> lbl =
                make_shared<FsmLabel>(x,0,presentationLayer);
                shared_ptr<FsmTransition> tr =
                make_shared<FsmTransition>(currentParsedNode,tgtNode,lbl);
                currentParsedNode->addTransition(tr);
            }
            else {
                
                // get the target state from the table entry
                size_t i0;
                string tgtStateName;
                i0 = tableEntry.find("/");
                if ( i0 == string::npos ) {
                    tgtStateName = tableEntry;
                }
                else {
                    tgtStateName = tableEntry.substr(0,i0);
                }
                tgtStateName.erase(0,tgtStateName.find_first_not_of(" \n\r\t\""));
                tgtStateName.erase(tgtStateName.find_last_not_of(" \n\r\t\"")+1);
                
                int tgtStateId = presentationLayer->state2Num(tgtStateName);
                if ( tgtStateId >= 0 ) {
                    tgtNode = nodes[tgtStateId];
                    if ( tgtNode == nullptr ) {
                        tgtNode = make_shared<FsmNode>(tgtStateId,
                                                       tgtStateName,
                                                       presentationLayer);
                        nodes[tgtStateId] = tgtNode;
                    }
                    // Get the output associated with the current state
                    // and input
                    string outStr =
                    ( i0 == string::npos ) ? "" : tableEntry.substr(i0+1);
                    outStr.erase(0,outStr.find_first_not_of(" \n\r\t\""));
                    outStr.erase(outStr.find_last_not_of(" \n\r\t\"")+1);
                    
                    int y =
                    ( outStr.empty() )
                    ? 0
                    : presentationLayer->out2Num(outStr);
                    
                    if ( y >= 0 ) {
                        shared_ptr<FsmLabel> lbl =
                        make_shared<FsmLabel>(x,y,presentationLayer);
                        shared_ptr<FsmTransition> tr =
                        make_shared<FsmTransition>(currentParsedNode,tgtNode,lbl);
                        currentParsedNode->addTransition(tr);
                    }
                    
                }
                else {
                    cout << endl << "ERROR: undefined target state "
                    << tgtStateName << endl;
                }
            }
            
            x++;
            pos = posEnd + 1;
            
        } while ( posEnd != string::npos );
        
        
        
        nodeId++;
        
    }
    
    
    inputFile.close();
    
}

void Dfsm::initDistTraces() {
    
    distTraces.clear();
    
    for ( size_t n = 0; n < size(); n++ ) {
        // Create empty vector of for row n, to be extended in
        // the inner loop
        vector< vector< shared_ptr< vector<int> > > > thisRow;
        for ( size_t m = 0; m < size(); m++ ) {
            // Create empty vector of pointers to traces
            vector< shared_ptr< vector<int> > > v;
            thisRow.push_back(v);
        }
        // Add thisRow to distTraces
        distTraces.push_back(thisRow);
    }
    
}

Dfsm::Dfsm(const std::string & fname,
           const std::string & fsmName) : Fsm(nullptr), dfsmTable(nullptr) {
    name = fsmName;
    presentationLayer = createPresentationLayerFromCsvFormat(fname);
    createDfsmTransitionGraph(fname);
}

Dfsm::Dfsm(const std::string & fname,
           const std::string & fsmName,
           const std::shared_ptr<FsmPresentationLayer>& pl) : Fsm(nullptr), dfsmTable(nullptr) {
    name = fsmName;
    presentationLayer = createPresentationLayerFromCsvFormat(fname,pl);
    createDfsmTransitionGraph(fname);
}

void Dfsm::createAtRandom()
{
    srand(getRandomSeed());
    
    for (unsigned int i = 0; i < nodes.size(); ++ i)
    {
        nodes [i] = make_shared<FsmNode>(i, presentationLayer);
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

vector<shared_ptr<PkTable> > Dfsm::getPktblLst() const
{
    return pktblLst;
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

Dfsm::Dfsm(const string & fname, const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer>& presentationLayer)
: Fsm(fname, fsmName, maxNodes, maxInput, maxOutput, presentationLayer)
{
    dfsmTable = nullptr;
}

Dfsm::Dfsm(const string& fname,
           const shared_ptr<FsmPresentationLayer>& presentationLayer,
           const string & fsmName)
: Fsm(fname,presentationLayer,fsmName)
{
    dfsmTable = nullptr;
}

Dfsm::Dfsm(const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer>& presentationLayer)
: Fsm(presentationLayer)
{
    dfsmTable = nullptr;
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

Dfsm::Dfsm(const string & fsmName, const int maxInput, const int maxOutput, const vector<shared_ptr<FsmNode>>& lst, const shared_ptr<FsmPresentationLayer>& presentationLayer)
: Fsm(fsmName, maxInput, maxOutput, lst, presentationLayer)
{
    dfsmTable = nullptr;
}

Dfsm::Dfsm(const Fsm & fsm)
: Fsm (fsm.getName(), fsm.getMaxInput(), fsm.getMaxOutput(), fsm.getNodes(), fsm.getPresentationLayer()), dfsmTable(nullptr)
{
    initStateIdx = fsm.getInitStateIdx();;
    minimal = isMinimal();
    /*shared_ptr<FsmNode> currentParsedNode;
     vector<shared_ptr<OFSMTable>> ofsmTableLst;
     shared_ptr<Tree> characterisationSet;
     vector<shared_ptr<Tree>> stateIdentificationSets;*/
}



Dfsm::Dfsm(const Json::Value& fsmExport) :
Fsm(), dfsmTable(nullptr)
{

    if (!fsmExport.isObject()) {
        cerr << endl << "File format is JSON but NOT FSM-lib file structure.";
        return;
    }
    
    bool valid = true;
    Json::Value inputs = fsmExport["inputs"];
    Json::Value outputs = fsmExport["outputs"];
    Json::Value states = fsmExport["states"];
    Json::Value transitions = fsmExport["transitions"];
    Json::Value requirements = fsmExport["requirements"];
    
    map< string,shared_ptr<FsmNode> > name2node;
    
    // check JSON value for a valid FSM export
    if (inputs.isNull() || (!inputs.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM inputs from JSON export file structure.";
    }
    if (outputs.isNull() || (!outputs.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM outputs from JSON export file structure.";
    }
    if (states.isNull() || (!states.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM states from JSON export file structure.";
    }
    if (transitions.isNull() || (!transitions.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM transitions from JSON export file structure.";
    }
    if (requirements.isNull() || (!requirements.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of requirements from JSON export file structure.";
    }
    if (!valid) {
        return;
    }
    
    // iterate over all inputs
    vector<string> in2String;
    for (const auto &input : inputs) {
        in2String.push_back(input.asString());
    }
    
    // iterate over all outputs
    vector<string> out2String;
    int theNopNo;
    bool haveNop = false;
    for (unsigned int index = 0; index < outputs.size(); ++index ) {
        string outStr(outputs[index].asString());
        if ( outStr == "_nop" ) {
            haveNop = true;
            theNopNo = index;
        }
        out2String.push_back(outStr);
    }
    
    // Add a NOP output for the case where the FSM is incomplete
    if ( not haveNop ) {
        out2String.push_back("_nop");
        theNopNo = (int)out2String.size() - 1;
    }
    
    // iterate over all states, insert initial state at index 0
    // of the state2String vector.
    vector<string> state2String;
    for (const auto &state : states) {
        if (state["initial"].asBool()) {
            state2String.push_back(state["name"].asString());
            break;
        }
    }
    for (const auto &state : states) {
        if (state["initial"].asBool()) {
            continue; // Initial state has already been inserted
        }
        state2String.push_back(state["name"].asString());
    }
    
    // Create the presentation layer
    presentationLayer =
    make_shared<FsmPresentationLayer>(in2String,out2String,state2String);
    
    // Define basic attributes
    name = "FSM";
    currentParsedNode = nullptr;
    maxInput = (int)in2String.size() - 1;
    maxOutput = (int)out2String.size() - 1;
    maxState = (int)state2String.size() - 1;
    initStateIdx = 0;
    minimal = Maybe;
    
    
    // Create all FSM states
    for ( size_t s = 0; s < state2String.size(); s++ ) {
        shared_ptr<FsmNode> theNode =
        make_shared<FsmNode>((int)s,state2String[s],presentationLayer);
        nodes.push_back(theNode);
        name2node[state2String[s]] = theNode;
    }
    
    
    // Create all transitions
    for (const auto &transition : transitions) {
        // Handle source and target nodes
        string srcName(transition["source"].asString());
        string tgtName(transition["target"].asString());
        
        shared_ptr<FsmNode> srcNode = name2node[srcName];
        shared_ptr<FsmNode> tgtNode = name2node[tgtName];
        
        if ( srcNode == nullptr ) {
            cerr << "Cannot associated valid FSM node with source node name"
            << srcName << endl;
            exit(1);
        }
        
        if ( tgtNode == nullptr ) {
            cerr << "Cannot associated valid FSM node with target node name"
            << tgtName << endl;
            exit(1);
        }
        
        // Get the output
        string yString(transition["output"].asString());
        
        // Trim
        yString.erase(0,yString.find_first_not_of(" \n\r\t\""));
        yString.erase(yString.find_last_not_of(" \n\r\t\"")+1);
        
        int y = presentationLayer->out2Num(yString);
        
        if ( y < 0 ) {
            cerr << "Unidentified output symbol `"
            <<  yString
            << "' in transition "
            << srcName << " --> " << tgtName
            << endl;
            exit(1);
        }
        
        // For each input, create a separate transition
        // and add it to the source node
        Json::Value inputlist = transition["input"];
        for (const auto &inidx : inputlist) {
            
            string xString(inidx.asString());
            // Trim
            xString.erase(0,xString.find_first_not_of(" \n\r\t\""));
            xString.erase(xString.find_last_not_of(" \n\r\t\"")+1);
            int x = presentationLayer->in2Num(xString);
            if ( x < 0 ) {
                cerr << "Unidentified input symbol `"
                <<  xString
                << "' in transition "
                << srcName << " --> " << tgtName
                << endl;
                exit(1);
            }
            shared_ptr<FsmLabel> theLabel =
            make_shared<FsmLabel>(x,y,presentationLayer);
            shared_ptr<FsmTransition> tr =
            make_shared<FsmTransition>(srcNode,tgtNode,theLabel);
            
            
            // Record the requirements satisfied by the transition
            Json::Value satisfies = transition["requirements"];
            for (const auto &satisfie : satisfies) {
                tr->addSatisfies(satisfie.asString());
            }
            
            srcNode->addTransition(tr);
        }
        
    }
    
    // Add requirements to nodes
    for (const auto &state : states) {
        string nodeName(state["name"].asString());
        
        for (const auto &n : nodes ) {
            if ( n->getName() == nodeName ) {
                Json::Value satisfies = state["requirements"];
                for (const auto &satisfie : satisfies) {
                    n->addSatisfies(satisfie.asString());
                }
            }
        }
        
    }
    
    
    // Force DFSM to be completely defined:
    // For each node and each input event x that has not been used
    // in any outgoing transition of this node,
    // create a self-loop transition with label x/NOP
    for ( auto n : nodes ) {
        
        vector<bool> inputs;
        for ( int x = 0; x <= maxInput; x++ ) {
            inputs.push_back(false);
        }
        for (const auto &tr : n->getTransitions() ) {
            inputs[tr->getLabel()->getInput()] = true;
        }
        for ( int x = 0; x <= maxInput; x++ ) {
            if ( not inputs[x] ) {
                shared_ptr<FsmLabel> theLabel =
                make_shared<FsmLabel>(x,theNopNo,presentationLayer);
                shared_ptr<FsmTransition> tr =
                make_shared<FsmTransition>(n,n,theLabel);
                n->addTransition(tr);
            }
        }
        
    }
    
}



Dfsm::Dfsm(const Json::Value& fsmExport,
           const std::shared_ptr<FsmPresentationLayer>& pl) :
Fsm(), dfsmTable(nullptr)
{

    if (!fsmExport.isObject()) {
        cerr << endl << "File format is JSON but NOT FSM-lib file structure.";
        return;
    }
    
    if ( pl == nullptr ) {
        cerr << endl << "Undefined presentation layer.";
        return;
    }
    
    bool valid = true;
    Json::Value inputs = fsmExport["inputs"];
    Json::Value outputs = fsmExport["outputs"];
    Json::Value states = fsmExport["states"];
    Json::Value transitions = fsmExport["transitions"];
    Json::Value requirements = fsmExport["requirements"];
    
    map< string,shared_ptr<FsmNode> > name2node;
    
    // check JSON value for a valid FSM export
    if (inputs.isNull() || (!inputs.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM inputs from JSON export file structure.";
    }
    if (outputs.isNull() || (!outputs.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM outputs from JSON export file structure.";
    }
    if (states.isNull() || (!states.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM states from JSON export file structure.";
    }
    if (transitions.isNull() || (!transitions.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of FSM transitions from JSON export file structure.";
    }
    if (requirements.isNull() || (!requirements.isArray())) {
        valid = false;
        cout << endl << "Unable to extract expected array of requirements from JSON export file structure.";
    }
    if (!valid) {
        return;
    }
    
    
    
    // iterate over all inputs; add all inputs not already contained
    // in pl to in2String.
    vector<string> in2String(pl->getIn2String());
    for (const auto &input : inputs) {
        string theInput(input.asString());
        if ( pl->in2Num(theInput) < 0 ) {
            in2String.push_back(theInput);
        }
    }
    
    // iterate over all outputs
    vector<string> out2String(pl->getOut2String());
    for (const auto &output : outputs) {
        string theOutput(output.asString());
        if ( pl->out2Num(theOutput) < 0 ) {
            out2String.push_back(theOutput);
        }
    }
    // Check whether the _nop output is already contained in pl,
    // otherwise add it to out2String
    int theNopNo = pl->out2Num("_nop");
    if ( theNopNo < 0 ) {
        out2String.push_back("_nop");
    }
    
    
    // iterate over all states, insert initial state at index 0
    // of the state2String vector.
    vector<string> state2String;
    for (const auto &state : states) {
        if (state["initial"].asBool()) {
            state2String.push_back(state["name"].asString());
            break;
        }
    }
    for (const auto &state : states) {
        if (state["initial"].asBool()) {
            continue; // Initial state has already been inserted
        }
        state2String.push_back(state["name"].asString());
    }
    
    // Create the presentation layer
    presentationLayer =
    make_shared<FsmPresentationLayer>(in2String,out2String,state2String);
    
    // Define basic attributes
    name = "FSM";
    currentParsedNode = nullptr;
    maxInput = (int)in2String.size() - 1;
    maxOutput = (int)out2String.size() - 1;
    maxState = (int)state2String.size() - 1;
    initStateIdx = 0;
    minimal = Maybe;
    
    
    // Create all FSM states
    for ( size_t s = 0; s < state2String.size(); s++ ) {
        shared_ptr<FsmNode> theNode =
        make_shared<FsmNode>((int)s,state2String[s],presentationLayer);
        nodes.push_back(theNode);
        name2node[state2String[s]] = theNode;
    }
    
    
    // Create all transitions
    for (const auto &transition : transitions) {
        // Handle source and target nodes
        string srcName(transition["source"].asString());
        string tgtName(transition["target"].asString());
        
        shared_ptr<FsmNode> srcNode = name2node[srcName];
        shared_ptr<FsmNode> tgtNode = name2node[tgtName];
        
        if ( srcNode == nullptr ) {
            cerr << "Cannot associated valid FSM node with source node name"
            << srcName << endl;
            exit(1);
        }
        
        if ( tgtNode == nullptr ) {
            cerr << "Cannot associated valid FSM node with target node name"
            << tgtName << endl;
            exit(1);
        }
        
        // Get the output
        string yString(transition["output"].asString());
        
        // Trim
        yString.erase(0,yString.find_first_not_of(" \n\r\t\""));
        yString.erase(yString.find_last_not_of(" \n\r\t\"")+1);
        
        int y = presentationLayer->out2Num(yString);
        
        if ( y < 0 ) {
            cerr << "Unidentified output symbol `"
            <<  yString
            << "' in transition "
            << srcName << " --> " << tgtName
            << endl;
            exit(1);
        }
        
        // For each input, create a separate transition
        // and add it to the source node
        Json::Value inputlist = transition["input"];
        for (const auto &inidx : inputlist) {
            
            string xString(inidx.asString());
            // Trim
            xString.erase(0,xString.find_first_not_of(" \n\r\t\""));
            xString.erase(xString.find_last_not_of(" \n\r\t\"")+1);
            int x = presentationLayer->in2Num(xString);
            if ( x < 0 ) {
                cerr << "Unidentified input symbol `"
                <<  xString
                << "' in transition "
                << srcName << " --> " << tgtName
                << endl;
                exit(1);
            }
            shared_ptr<FsmLabel> theLabel =
            make_shared<FsmLabel>(x,y,presentationLayer);
            shared_ptr<FsmTransition> tr =
            make_shared<FsmTransition>(srcNode,tgtNode,theLabel);
            
            
            // Record the requirements satisfied by the transition
            Json::Value satisfies = transition["requirements"];
            for (const auto &satisfie : satisfies) {
                tr->addSatisfies(satisfie.asString());
            }
            
            srcNode->addTransition(tr);
        }
        
    }
    
    // Add requirements to nodes
    for (const auto &state : states) {
        string nodeName(state["name"].asString());
        
        for (const auto &n : nodes ) {
            if ( n->getName() == nodeName ) {
                Json::Value satisfies = state["requirements"];
                for (const auto &satisfie : satisfies) {
                    n->addSatisfies(satisfie.asString());
                }
            }
        }
        
    }
    
    
    // Force DFSM to be completely defined:
    // For each node and each input event x that has not been used
    // in any outgoing transition of this node,
    // create a self-loop transition with label x/NOP
    for ( auto n : nodes ) {
        
        vector<bool> inputs;
        for ( int x = 0; x <= maxInput; x++ ) {
            inputs.push_back(false);
        }
        for (const auto &tr : n->getTransitions() ) {
            inputs[tr->getLabel()->getInput()] = true;
        }
        for ( int x = 0; x <= maxInput; x++ ) {
            if ( not inputs[x] ) {
                shared_ptr<FsmLabel> theLabel =
                make_shared<FsmLabel>(x,theNopNo,presentationLayer);
                shared_ptr<FsmTransition> tr =
                make_shared<FsmTransition>(n,n,theLabel);
                n->addTransition(tr);
            }
        }
        
    }
    
}


void Dfsm::calcPkTables() {
    
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
    
#if 0
    cout << "DFSM-Table" << endl;
    cout << *dfsmTable << endl;
    cout << "PK-Tables" << endl;
    PkTable::counter = 0;
    for (auto p : pktblLst) {
        PkTable::counter++;
        cout << *p << endl;
        
    }
#endif
    
}


Dfsm Dfsm::minimise()
{
    
    vector<shared_ptr<FsmNode>> uNodes;
    removeUnreachableNodes(uNodes);
    
    calcPkTables();
    shared_ptr<PkTable> pMin = pktblLst[pktblLst.size()-1];
    
    auto dfsm = pMin->toFsm(name, maxOutput);
    dfsm.calcPkTables();
    return dfsm;
}

void Dfsm::printTables() const
{
    ofstream file("tables.tex");
    if (dfsmTable != nullptr)
    {
        file << *dfsmTable;
    }
    
    for (const auto &i : pktblLst) {
        file << *i << endl << endl;
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
    
    
#if 0
    
    cout << "Dfsm::getCharacterisationSet()" << endl;
    cout << *p1 << endl;
    for (auto p : pktblLst) {
        
        cout << *p << endl;
    }
    
#endif
    
    /*Create an empty characterisation set as an empty InputTree instance*/
    characterisationSet = make_shared<Tree>(make_shared<TreeNode>(), presentationLayer);
    
    /*Loop over all non-equal pairs of states. If they are not already distinguished by
     the input sequences contained in w, create a new input traces that distinguishes them
     and add it to w.*/
    for (unsigned int left = 0; left < nodes.size(); ++ left)
    {
        shared_ptr<FsmNode> leftNode = nodes.at(left);
        for (unsigned int right = left + 1; right < nodes.size(); ++ right)
        {
            shared_ptr<FsmNode> rightNode = nodes.at(right);
            
            if (leftNode->distinguished(rightNode, characterisationSet) != nullptr)
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
            characterisationSet->addToRoot(tcli);
            
#if 0
            cout << "Distinguishing trace for " << nodes[left]->getName()
            << " and " << nodes[right]->getName() << ": " << endl;
            cout << i << endl;
#endif
            
        }
    }
    
#if 0
    /* CHECK */
    for (unsigned int left = 0; left < nodes.size(); ++ left)
    {
        shared_ptr<FsmNode> leftNode = nodes.at(left);
        for (unsigned int right = left + 1; right < nodes.size(); ++ right)
        {
            shared_ptr<FsmNode> rightNode = nodes.at(right);
            
            if (leftNode->distinguished(rightNode, characterisationSet) == nullptr)
            {
                cerr << "ERROR: nodes " << leftNode->getName() << " and " << rightNode->getName() << " cannot be distinguished by W" << endl;
            }
            
        }
    }
#endif
    
    // Minimise characterisation set
    // minimiseCharSet(characterisationSet);
    
    /* Wrap list of lists by an IOListContainer instance */
    IOListContainer tcl = characterisationSet->getIOLists();
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



IOListContainer Dfsm::wMethod(const unsigned int numAddStates) {
    
    Dfsm dfsmMin = minimise();
    return dfsmMin.wMethodOnMinimisedDfsm(numAddStates);
    
}


IOListContainer Dfsm::wMethodOnMinimisedDfsm(const unsigned int numAddStates)
{
    
    shared_ptr<Tree> iTree = getTransitionCover();
    
    if (numAddStates > 0)
    {
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

IOListContainer Dfsm::wpMethod(const unsigned int numAddStates)
{
    Dfsm dfsmMin = minimise();
    return dfsmMin.wpMethodOnMinimisedDfsm(numAddStates);
}

IOListContainer Dfsm::wpMethodOnMinimisedDfsm(const unsigned int numAddStates)
{
    shared_ptr<Tree> scov = getStateCover();

    shared_ptr<Tree> tcov = getTransitionCover();

    tcov->remove(scov);
    const shared_ptr<Tree> &r = tcov;

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

IOListContainer Dfsm::hsiMethod(const unsigned int numAddStates)
{
    Fsm fMin = minimiseObservableFSM();
    return fMin.hsiMethod(numAddStates);
}

IOListContainer Dfsm::tMethod()
{
    
    shared_ptr<Tree> iTree = getTransitionCover();
    
    return iTree->getIOLists();
    
}


void Dfsm::toCsv(const string& fname) {
    ofstream out(fname + ".csv");
    
    // Table heading contains input identifiers
    for ( int x = 0; x <= maxInput; x++ ) {
        out << " ; ";
        out << presentationLayer->getInId(x);
    }
    
    for (const auto &node : nodes) {
        out << endl << "\"" << node->getName() << "\"";
        
        for ( int x = 0; x <= maxInput; x++ ) {
            
            out << " ; ";
            
            for (const auto &tr : node->getTransitions() ) {
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

InputTrace Dfsm::calcDistinguishingTrace(
        const shared_ptr<InputTrace> iAlpha,
        const shared_ptr<InputTrace> iBeta,
        const shared_ptr<Tree> tree)
{
    shared_ptr<FsmNode> s0 = getInitialState();
    shared_ptr<FsmNode> s1 = *s0->after(*iAlpha).begin();
    shared_ptr<FsmNode> s2 = *s0->after(*iBeta).begin();

    InputTrace gamma = calcDistinguishingTraceInTree(s1, s2, tree);
    if (!gamma.get().empty())
        return gamma;

    InputTrace gamma2 = calcDistinguishingTraceAfterTree(s1, s2, tree);
    if (!gamma2.get().empty())
        return gamma2;

    return s1->calcDistinguishingTrace(s2,
                                       pktblLst,
                                       maxInput);
}


vector<int> Dfsm::calcDistinguishingTrace(shared_ptr<SegmentedTrace> alpha,
                                               shared_ptr<SegmentedTrace> beta, const shared_ptr<TreeNode> treeNode) {
    
    
    shared_ptr<FsmNode> s0 = getInitialState();
    shared_ptr<FsmNode> s1 = alpha->getTgtNode();
    shared_ptr<FsmNode> s2 = beta->getTgtNode();
    
    shared_ptr<Tree> tree = make_shared<Tree>(treeNode,presentationLayer);
    
    InputTrace gamma = calcDistinguishingTraceInTree(s1, s2, tree);
    if (!gamma.get().empty())
        return gamma.get();
    
    InputTrace gamma2 = calcDistinguishingTraceAfterTree(s1, s2, tree);
    if (!gamma2.get().empty())
        return gamma2.get();
    
    return s1->calcDistinguishingTrace(s2,
                                       pktblLst,
                                       maxInput).get();
    
}




InputTrace Dfsm::calcDistinguishingTraceInTree(
        const shared_ptr<FsmNode> s_i,
        const shared_ptr<FsmNode> s_j,
        const shared_ptr<Tree> tree)
{
    shared_ptr<TreeNode> root = tree->getRoot();
    shared_ptr<TreeNode> currentNode = root;
    deque<shared_ptr<InputTrace>> q1;

    /* initialize queue */
    for (const shared_ptr<TreeEdge> &e : *currentNode->getChildren())
    {
        shared_ptr<InputTrace> itrc = make_shared<InputTrace>(presentationLayer);
        itrc->add(e->getIO());
        q1.push_back(itrc);
    }

    /* Breadth-first search */
    while(!q1.empty())
    {
        shared_ptr<InputTrace> itrc = q1.front();
        q1.pop_front();

        if(s_i->distinguished(s_j, itrc->get()))
        {
            return *itrc;
        }

        currentNode = root->after(itrc->cbegin(), itrc->cend());

        for (const shared_ptr<TreeEdge> &ne : *currentNode->getChildren())
        {
            shared_ptr<InputTrace> itrcTmp = make_shared<InputTrace>(itrc->get(), presentationLayer);
            vector<int>nItrc;
            nItrc.push_back(ne->getIO());
            itrcTmp->append(nItrc);
            q1.push_back(itrcTmp);
        }
    }
    // Return empty trace: no distinguishing trace found in tree
    return InputTrace(presentationLayer);
}

InputTrace Dfsm::calcDistinguishingTraceInTree(const shared_ptr<InputTrace> alpha, const shared_ptr<InputTrace> beta, const shared_ptr<Tree> tree)
{
    // Only one element in the set, since FSM is deterministic
    shared_ptr<FsmNode> s_i = *(getInitialState()->after(*alpha)).begin();
    shared_ptr<FsmNode> s_j = *(getInitialState()->after(*beta)).begin();
    return calcDistinguishingTraceInTree(s_i, s_j, tree);
}

InputTrace Dfsm::calcDistinguishingTraceAfterTree(
        const shared_ptr<FsmNode> s_i,
        const shared_ptr<FsmNode> s_j,
        const shared_ptr<Tree> tree)
{
    vector<shared_ptr<TreeNode>> leaves = tree->getLeaves();
    for(const shared_ptr<TreeNode> &leaf : leaves)
    {
        shared_ptr<InputTrace> itrc = make_shared<InputTrace>(leaf->getPath(), presentationLayer);
        shared_ptr<FsmNode> s_i_after_input = *(s_i->after(*itrc)).begin();
        shared_ptr<FsmNode> s_j_after_input = *(s_j->after(*itrc)).begin();
        
        // Cannot find distinguishing trace if the states reached
        // by the path are identical
        if(s_i_after_input == s_j_after_input ) continue;

        // Since we are dealing with a minimised DFSM, a distinguishing
        // trace can ALWAYS be found, if s_i_after_input != s_j_after_input
        InputTrace gamma = s_i_after_input->calcDistinguishingTrace(s_j_after_input ,
                                                                        pktblLst,
                                                                        maxInput);

        itrc->append(gamma.get());

        return *itrc;
    }
    
    // Return empty trace: could not find a tree extension
    // distinguishing s_i and s_j
    return InputTrace(presentationLayer);
}

IOListContainer Dfsm::hMethodOnMinimisedDfsm(const unsigned int numAddStates) {
    
    // Our initial state
    shared_ptr<FsmNode> s0 = getInitialState();
    
    // We need a valid set of DFSM table and Pk-Tables for this method
    if ( dfsmTable == nullptr ) {
        calcPkTables();
    }
    
    // Auxiliary state cover set needed for further computations
    shared_ptr<Tree> V = getStateCover();
    
    // Test suite is initialised with the state cover
    shared_ptr<Tree> iTree = getStateCover();
    
    IOListContainer inputEnum = IOListContainer(maxInput,
                                                (int)numAddStates+1,
                                                (int)numAddStates+1,
                                                presentationLayer);
    
    // Initial test suite set is V.Sigma^{m-n+1}, m-n = numAddStates
    iTree->add(inputEnum);
    
    // Step 1.
    // Add all alpha.gamma, beta.gamma where alpha, beta in V
    // and gamma distinguishes s0-after-alpha, s0-after-beta
    // (if alpha.gamma or beta.gamma are already in iTree, addition
    // will not lead to a new test case)
    IOListContainer iolcV = V->getIOListsWithPrefixes();
    shared_ptr<vector<vector<int>>> iolV = iolcV.getIOLists();
    
    for ( size_t i = 0; i < iolV->size(); i++ ) {
        
        shared_ptr<InputTrace> alpha =
        make_shared<InputTrace>(iolV->at(i),presentationLayer);

        for ( size_t j = i+1; j < iolV->size(); j++ ) {
            
            shared_ptr<InputTrace> beta =
            make_shared<InputTrace>(iolV->at(j),presentationLayer);

            shared_ptr<Tree> alphaTree = iTree->getSubTree(alpha);
            shared_ptr<Tree> betaTree = iTree->getSubTree(beta);
            shared_ptr<Tree> prefixRelationTree = alphaTree->getPrefixRelationTree(betaTree);

            InputTrace gamma = calcDistinguishingTrace(alpha, beta, prefixRelationTree);

            shared_ptr<InputTrace> iAlphaGamma = make_shared<InputTrace>(alpha->get(), presentationLayer);
            iAlphaGamma->append(gamma.get());

            shared_ptr<InputTrace> iBetaGamma = make_shared<InputTrace>(beta->get(), presentationLayer);
            iBetaGamma->append(gamma.get());

            iTree->addToRoot(iAlphaGamma->get());
            iTree->addToRoot(iBetaGamma->get());
        }
        
    }
    
    // Step 2.
    // For each sequence α.β, α ∈ Q, |β| = m – n + 1, and each non-empty prefix
    // β1 of β that takes the DFSM from s0 to state s,
    // add sequences α.β1.γ and ω.γ, where ω ∈ V and s0-after-ω ≠ s,
    // and γ is a distinguishing sequence of states s0-after-α.β1
    // and s0-after-ω.
    IOListContainer allBeta = IOListContainer(maxInput,
                                              1,
                                              (int)numAddStates+1,
                                              presentationLayer);
    
    shared_ptr<vector<vector<int>>> iolAllBeta = allBeta.getIOLists();
    
    for (const auto &beta : *iolAllBeta ) {
        
        for (const auto &alpha : *iolV ) {
            
            shared_ptr<InputTrace> iAlphaBeta =
                make_shared<InputTrace>(alpha,presentationLayer);
            iAlphaBeta->append(beta);
            unordered_set<shared_ptr<FsmNode>>
                s_alpha_betaSet = s0->after(*iAlphaBeta);
            shared_ptr<FsmNode> s_alpha_beta = *s_alpha_betaSet.begin();
            
            for ( auto omega : *iolV ) {
                shared_ptr<InputTrace>
                    iOmega = make_shared<InputTrace>(omega,presentationLayer);
                unordered_set<shared_ptr<FsmNode>>
                    s_omegaSet = s0->after(*iOmega);
                shared_ptr<FsmNode> s_omega = *s_omegaSet.begin();

                if ( s_alpha_beta == s_omega ) continue;

                shared_ptr<Tree> alphaBetaTree = iTree->getSubTree(iAlphaBeta);
                shared_ptr<Tree> trAfterOmega = iTree->getSubTree(iOmega);
                shared_ptr<Tree> prefixRelationTree = alphaBetaTree->getPrefixRelationTree(trAfterOmega);

                InputTrace gamma = calcDistinguishingTrace(iAlphaBeta, iOmega, prefixRelationTree);

                shared_ptr<InputTrace> iAlphaBetaGamma = make_shared<InputTrace>(iAlphaBeta->get(), presentationLayer);
                iAlphaBetaGamma->append(gamma.get());

                shared_ptr<InputTrace> iOmegaGamma = make_shared<InputTrace>(iOmega->get(), presentationLayer);
                iOmegaGamma->append(gamma.get());

                iTree->addToRoot(iAlphaBetaGamma->get());
                iTree->addToRoot(iOmegaGamma->get());
            }
            
        }
        
    }
    
    // Step 3.
    // For each sequence α.β,α∈Q,|β|=m–n+1, and each two
    // non-empty prefixes β1 and β2 of β that take the
    // DFSM from state s0-after-alpha
    // to two different states add sequences α.β1.γ and α.β2.γ,
    // where γ is a distinguishing sequence of states
    // s0-after-alpha.beta1 and s0-after-alpha.beta2.
    
    for ( auto alpha : *iolV ) {
        
        shared_ptr<InputTrace> iAlpha =
            make_shared<InputTrace>(alpha,presentationLayer);
        
        for ( auto beta : *inputEnum.getIOLists() ) {
        
            for ( size_t i = 0; i < beta.size() - 1; i++ ) {
                
                shared_ptr<InputTrace> iBeta_1 = make_shared<InputTrace>(presentationLayer);
                for ( size_t k = 0; k <= i; k++ ) {
                    iBeta_1->add(beta[k]);
                }
                
                for ( size_t j = i+1; j < beta.size(); j++ ) {
                    
                    shared_ptr<InputTrace> iBeta_2 =
                        make_shared<InputTrace>(presentationLayer);
                    for ( size_t k = 0; k <= j; k++ ) {
                        iBeta_2->add(beta[k]);
                    }
                    
                    shared_ptr<InputTrace> iAlphaBeta_1 =
                        make_shared<InputTrace>(alpha,presentationLayer);
                    iAlphaBeta_1->append(iBeta_1->get());
                    
                    shared_ptr<InputTrace> iAlphaBeta_2 =
                    make_shared<InputTrace>(alpha,presentationLayer);
                    iAlphaBeta_2->append(iBeta_2->get());
                    
                    unordered_set<shared_ptr<FsmNode>> s1Set =
                    s0->after(*iAlphaBeta_1);
                    shared_ptr<FsmNode> s1 = *s1Set.begin();
                    
                    unordered_set<shared_ptr<FsmNode>> s2Set =
                    s0->after(*iAlphaBeta_2);
                    shared_ptr<FsmNode> s2 = *s2Set.begin();
                    
                    if ( s1 == s2 ) continue;

                    shared_ptr<Tree> afterAlphaBeta1Tree = iTree->getSubTree(iAlphaBeta_1);
                    shared_ptr<Tree> afterAlphaBeta2Tree = iTree->getSubTree(iAlphaBeta_2);
                    shared_ptr<Tree> prefixRelationTree = afterAlphaBeta1Tree->getPrefixRelationTree(afterAlphaBeta2Tree);

                    InputTrace gamma = calcDistinguishingTrace(iAlphaBeta_1, iAlphaBeta_2, prefixRelationTree);

                    shared_ptr<InputTrace> iAlphaBeta_1Gamma = make_shared<InputTrace>(iAlphaBeta_1->get(), presentationLayer);
                    iAlphaBeta_1Gamma->append(gamma.get());

                    shared_ptr<InputTrace> iAlphaBeta_2Gamma = make_shared<InputTrace>(iAlphaBeta_2->get(), presentationLayer);
                    iAlphaBeta_2Gamma->append(gamma.get());

                    iTree->addToRoot(iAlphaBeta_1Gamma->get());
                    iTree->addToRoot(iAlphaBeta_2Gamma->get());
                }
            }
        
        }
    }

    return iTree->getIOLists();

}


bool Dfsm::distinguishable(const FsmNode& s1, const FsmNode& s2) {
    
    if ( pktblLst.empty() ) {
        calcPkTables();
    }
    
    shared_ptr<PkTable> p = pktblLst.back();
    
    return ( p->getClass(s1.getId()) != p->getClass(s2.getId()) );
    
}

void Dfsm::calculateDistMatrix() {
    initDistTraces();
    calcPkTables();
    
    for ( size_t n = 0; n < size(); n++ ) {
        for ( size_t m = n+1; m < size(); m++ ) {
            // Skip indistinguishable nodes
            if ( not distinguishable(*nodes[n], *nodes[m]) ) continue;
            
            vector< shared_ptr< vector<int> > > u = calcDistTraces(*nodes[n],*nodes[m]);
            distTraces[n][m] = u;
            distTraces[m][n] = u;
            
        }
    }
    
    
}

vector< shared_ptr< vector<int> > >  Dfsm::calcDistTraces(shared_ptr< vector<int> > trc,
                                                          int id1,
                                                          int id2) {
    
    vector< shared_ptr< vector<int> > > v;
    
    for ( int x = 0; x <= maxInput; x++ ) {
        int y1 = dfsmTable->getRow(id1)->getioSection()[x];
        int y2 = dfsmTable->getRow(id2)->getioSection()[x];
        
        if ( y1 != y2 ) {
            shared_ptr< vector<int> > newTrc = make_shared< vector<int> >(*trc);
            newTrc->push_back(x);
            v.push_back(newTrc);
        }
    }
    
    return v;
    
}

vector< shared_ptr< vector<int> > > Dfsm::calcDistTraces(size_t l,
                                                         shared_ptr< vector<int> > trc,
                                                         int id1,
                                                         int id2) {
    if ( l == 0 ) return calcDistTraces(trc,id1,id2);
    
    vector< shared_ptr< vector<int> > > v;
    shared_ptr<PkTable> thisPkTbl = pktblLst[l];
    shared_ptr<PkTable> prevPkTbl = pktblLst[l-1];
    
    for ( int x = 0; x <= maxInput; x++ ) {
        int idNext1 = thisPkTbl->getRow(id1)->getI2PMap()[x];
        int idNext2 = thisPkTbl->getRow(id2)->getI2PMap()[x];

        if ( prevPkTbl->getClass(idNext1) != prevPkTbl->getClass(idNext2) ) {
            shared_ptr< vector<int> > newTrc = make_shared< vector<int> >(*trc);
            newTrc->push_back(x);
            vector< shared_ptr< vector<int> > > w = calcDistTraces(l-1,newTrc,idNext1,idNext2);
            v.insert(v.end(),w.begin(),w.end());
        }
        
    }
    
    return v;
    
}


vector< shared_ptr< vector<int> > > Dfsm::calcDistTraces(FsmNode& s1,
                                                         FsmNode& s2) {
    
    int id1 = s1.getId();
    int id2 = s2.getId();
    
    // Find first Pk-table distPkTbl,  where s1 and s2 are distinguished
    size_t l = 0;
    shared_ptr<PkTable> distPkTbl = nullptr;
    for ( ; l < pktblLst.size(); l++ ) {
        distPkTbl = pktblLst[l];
        if ( distPkTbl->getClass(id1) != distPkTbl->getClass(id2) )
            break;
    }
    
    return calcDistTraces(l,make_shared< vector<int> >(),id1,id2);
}


vector< shared_ptr< vector<int> > > Dfsm::getDistTraces(FsmNode& s1,
                                                                  FsmNode& s2) {
    
    return distTraces[s1.getId()][s2.getId()];
    
}






size_t Dfsm::spyhEstimateGrowthOfTestSuite(const std::shared_ptr<FsmNode> u, const std::shared_ptr<FsmNode> v, int x) {
    // return 1 if x already distinguishes u and v
    InputTrace xTrace(x,presentationLayer);
    if (u->apply(xTrace) != v->apply(xTrace))
        return 1;
    
    std::shared_ptr<FsmNode> ux = *u->after(x).cbegin();
    std::shared_ptr<FsmNode> vx = *v->after(x).cbegin();
    
    // if after application of x both states coincide or remain as before 
    // the application or "switch", then a value larger then produced by 
    // any shortest distinguishing trace is returned
    if (ux == vx
        || (ux == u && vx == v)
        || (ux == v && vx == u))
        return 2 * (getMaxNodes()+1);

    // otherwise compute and use the shortest distinguishing trace of u and v
    auto distTraces = getDistTraces(*u,*v);
    size_t minEst = distTraces[0]->size();
    for (auto trace : distTraces) {
        if (trace->size() < minEst) {
            minEst = trace->size();
        }
    }
    return minEst * 2 + 1;

    // TODO: avoid recomputation of minimal length dist-traces
}


std::pair<size_t,std::stack<int>> Dfsm::spyhGetPrefixOfSeparatingTrace(const std::vector<int>& trace1, const std::vector<int>& trace2, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph) {
    // get shortest sequences in each class
    std::vector<int> u = trace1;
    for (auto trace : graph.getConvergentTraces(trace1)) {
        if (trace->getPath().size() < u.size()) {
            u = trace->getPath();
        }
    }
    std::vector<int> v = trace2;
    for (auto trace : graph.getConvergentTraces(trace2)) {
        if (trace->getPath().size() < u.size()) {
            u = trace->getPath();
        }
    }

    std::shared_ptr<FsmNode> stateU = *getInitialState()->after(u).begin();
    std::shared_ptr<FsmNode> stateV = *getInitialState()->after(v).begin();

    // get shortest distinguishing trace and use double its size as minEst
    auto distTraces = getDistTraces(*stateU,*stateV);
    size_t minEst = distTraces[0]->size();
    for (auto trace : distTraces) {
        if (trace->size() < minEst) {
            minEst = trace->size();
        }
    }
    minEst = 2 * minEst;

    // start with the empty prefix
    std::stack<int> bestPrefix;

    if (!graph.hasLeaf(u)) 
        minEst = minEst + u.size();
    if (!graph.hasLeaf(v)) 
        minEst = minEst + v.size();

    for (int x = 0; x <= getMaxInput(); ++x) {
        bool foundExtensionForU = false;        

        // check if there exist u' in [u] and trace w such that u'xw are contained in the test suite
        for (auto convU : graph.getConvergentTraces(u)) {
            
            // u'xw is contained in the test suite if and only if x is defined in the testsuite after u'
            if (!convU->isDefined(x))
                continue;

            foundExtensionForU = true;    
            bool foundExtensionForV = false;  

            // check if there exist v' in [v] and trace w' such that v'xw' are contained in the test suite
            for (auto convV : graph.getConvergentTraces(v)) {

                // u'xw is contained in the test suite if and only if x is defined in the testsuite after u'
                if (!convV->isDefined(x))
                    continue;

                foundExtensionForV = true;   
            
                // if x produces different outputs to the states reached by u' and v', then the
                // test suite already distinguishes u and v
                InputTrace xTrace(x,presentationLayer);
                if (stateU->apply(xTrace) != stateV->apply(xTrace))
                    return std::make_pair<size_t,std::stack<int>>(0,{});

                
                // if the states reached by u'x and v'x coincide, and the same output is produced on x after both u' and v', then 
                // x cannot distinguish the states reached by u' and v'
                if (stateU->after(x) == stateV->after(x))
                    continue;

                std::vector<int> nextU(convU->getPath());
                nextU.push_back(x);
                std::vector<int> nextV(convV->getPath());
                nextV.push_back(x);

                // recursive call on u'x and v'x
                auto targetPair = spyhGetPrefixOfSeparatingTrace(nextU, nextV, testSuite, graph);

                // if u'x and v'x are already distinguished, then so are u' and v'
                if (targetPair.first == 0)
                    return std::make_pair<size_t,std::stack<int>>(0,{});

                // if the target pair is better than the current best, then use the target pair
                if (targetPair.first <= minEst) {
                    minEst = targetPair.first;
                    bestPrefix = targetPair.second;
                    // input x must be added to the front of the best result for the target pair
                    bestPrefix.push(x);
                }                
            }

            // if no such v' has been found, estimate growth on u and v
            if (!foundExtensionForV) {
                size_t e = spyhEstimateGrowthOfTestSuite(stateU,stateV,x);   
                if (e != 1) {
                    if (graph.hasLeaf(u)) {
                       e = e+1; 
                    } else {
                        std::vector<int> nextU(u);
                        nextU.push_back(x);
                        if (!graph.hasLeaf(nextU)) {
                            e = e + 1 + u.size();
                        }                        
                    }                    
                }
                if (!graph.hasLeaf(v)) {
                    e = e + v.size();
                }                    
                if (e < minEst) {
                    minEst = e;
                    bestPrefix = std::stack<int>();
                    bestPrefix.push(x);
                }    
            }
        }

        // if no such u' has been found, check if some applicable v' exists and estimate growth on u and v
        if (!foundExtensionForU) {
            bool foundExtensionForV = false;  

            // check if there exist v' in [v] and trace w' such that v'xw' are contained in the test suite
            for (auto convV : graph.getConvergentTraces(v)) {

                // u'xw is contained in the test suite if and only if x is defined in the testsuite after u'
                if (!convV->isDefined(x))
                    continue;

                foundExtensionForV = true;
                break;
            }

            if (foundExtensionForV) {
                size_t e = spyhEstimateGrowthOfTestSuite(stateU,stateV,x);   
                if (e != 1) {
                    if (graph.hasLeaf(v)) {
                       e = e+1; 
                    } else {
                        std::vector<int> nextV(v);
                        nextV.push_back(x);
                        if (!graph.hasLeaf(nextV)) {
                            e = e + 1 + v.size();
                        }                        
                    }                    
                }
                if (!graph.hasLeaf(u)) {
                    e = e + u.size();
                }                    
                if (e < minEst) {
                    minEst = e;
                    bestPrefix = std::stack<int>();
                    bestPrefix.push(x);
                }    
            }  
        }
    }
    return std::make_pair(minEst,bestPrefix); 
}


void Dfsm::spyhAppendSeparatingSequence(const std::vector<int>& traceToAppendTo, const std::vector<int>& traceToAppend, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph) 
{
    // get shortest u' in [traceToAppendTo]
    std::vector<int> uBest = traceToAppendTo;
    for (auto trace : graph.getConvergentTraces(traceToAppendTo)) {
        if (trace->getPath().size() < uBest.size()) {
            uBest = trace->getPath();
        }
    }

    // note: this is initialised as -1 in the original SPYH method;
    //       below we compare idx and maxLength using >= instead of >,
    //       resulting in the same effect but enabling the use of unsigned 
    //       types for maxLength
    unsigned int maxLength = 0;
    for (auto u : graph.getConvergentTraces(traceToAppendTo)) {
        // get length of the longest prefix w of traceToAppend such that traceToAppendTo.w is in the test suite
        unsigned int idx = 0;
        std::shared_ptr<TreeNode> node = u;
        for (; idx < traceToAppend.size(); ++idx) {
            node = node->after(traceToAppend[idx]);
            if (node == nullptr) break;
        }

        // if traceToAppend has already been fully applied, nothing further needs to be done
        if (idx == traceToAppend.size())
            return;

        // update the best length only if the current length is better and requires no branching
        if (idx >= maxLength && (node == nullptr || node->isLeaf())) {
            uBest = u->getPath();
            maxLength = idx;
        }
    }    

    // append traceToAppend after uBest and insert the trace into testsuite and graph
    uBest.insert(uBest.end(), traceToAppend.begin(), traceToAppend.end());
    testSuite->addToRoot(uBest);
    graph.add(uBest);
}


void Dfsm::spyhDistinguish(const std::vector<int>& trace, const std::unordered_set<std::shared_ptr<InputTrace>> traces, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph) {
    std::shared_ptr<FsmNode> u = *getInitialState()->after(trace).begin();
    
    for (auto other : traces) {
        std::vector<int> otherTrace = other->get();
        std::shared_ptr<FsmNode> v = *getInitialState()->after(otherTrace).begin();
        
        // trace needs to be distinguished only from traces reaching other states
        if ( u == v) continue;

        auto distPair = spyhGetPrefixOfSeparatingTrace(trace,otherTrace,testSuite,graph);
        std::vector<int> w;
        while (!distPair.second.empty()) {
            int x = distPair.second.top();
            distPair.second.pop();
            w.push_back(x);
        }

        // the test suite needs to be extended only if it does not already r-dist 
        if (distPair.first > 0) {
            
            std::vector<int> uw = trace;
            uw.insert(uw.end(), w.begin(), w.end());
            std::vector<int> vw = otherTrace;
            uw.insert(vw.end(), w.begin(), w.end());
            std::shared_ptr<FsmNode> uwState = *getInitialState()->after(uw).begin();
            std::shared_ptr<FsmNode> vwState = *getInitialState()->after(vw).begin();

            // always chooses the first shortest distinguishing trace
            auto distTrace = getDistTraces(*uwState,*vwState)[0];   

            auto traceToAppend = w;
            traceToAppend.insert(traceToAppend.end(), distTrace->begin(), distTrace->end());    

            // append separating sequences to 
            spyhAppendSeparatingSequence(trace,traceToAppend,testSuite,graph);
            spyhAppendSeparatingSequence(otherTrace,traceToAppend,testSuite,graph);
        }
    }    
}

void Dfsm::spyhDistinguishFromSet(const std::vector<int>& u, const std::vector<int>& v, const std::unordered_set<std::shared_ptr<InputTrace>> stateCover, std::unordered_set<std::shared_ptr<InputTrace>> tracesToDistFrom, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph, unsigned int depth) {
    
    // dist u 
    spyhDistinguish(u,tracesToDistFrom,testSuite,graph);

    // check whether [v] contains no traces from the state cover
    bool notReferenced = true;
    auto tracesConvergentToV = graph.getConvergentTraces(v);
    for (auto coverTrace : stateCover) {
        auto s = coverTrace->get();
        for (auto traceConvergentToV : tracesConvergentToV) {
            if (traceConvergentToV->getPath() == s) {
                notReferenced = false;
                break;
            }
        }
        if (!notReferenced) 
            break;
    }

    // if v has not been referenced already, distinguish it too
    if (notReferenced) {
        spyhDistinguish(v,tracesToDistFrom,testSuite,graph);
    }

    if (depth > 0) {
        std::shared_ptr<InputTrace> uTrace = std::make_shared<InputTrace>(u,presentationLayer);
        std::shared_ptr<InputTrace> vTrace = std::make_shared<InputTrace>(v,presentationLayer);
        
        tracesToDistFrom.insert(uTrace);
        if (notReferenced) 
            tracesToDistFrom.insert(vTrace);

        for (int x = 0; x <= getMaxInput(); ++x) {
            spyhAppendSeparatingSequence(u,{x},testSuite,graph);
            spyhAppendSeparatingSequence(v,{x},testSuite,graph);
            std::vector<int> ux = u;
            ux.push_back(x);
            std::vector<int> vx = v;
            vx.push_back(x);
            spyhDistinguishFromSet(ux,vx,stateCover,tracesToDistFrom,testSuite,graph,depth-1);
        }


        if (notReferenced) 
            tracesToDistFrom.erase(vTrace);
        tracesToDistFrom.erase(uTrace);
    }
}



IOListContainer Dfsm::spyhMethodOnMinimisedCompleteDfsm(const unsigned int numAddStates) {
    // Our initial state
    shared_ptr<FsmNode> s0 = getInitialState();
    
    // We need a valid set of DFSM table, Pk-Tables, and dist-traces for this method
    if ( dfsmTable == nullptr ) {
        calcPkTables();
        calculateDistMatrix();
    }
    
    // collect all transitions to be verified
    std::unordered_set<std::shared_ptr<FsmTransition>> transitions;
    for (auto node : nodes) {
        for (auto t : node->getTransitions()) {
            transitions.insert(t);
        }
    }

    // State Cover
    shared_ptr<Tree> V = getStateCover();
    IOListContainer iolcV = V->getIOListsWithPrefixes();
    shared_ptr<vector<vector<int>>> iolV = iolcV.getIOLists();
    std::unordered_set<std::shared_ptr<InputTrace>> stateCover;
    std::vector<std::shared_ptr<InputTrace>> stateCoverAssignment;
    stateCoverAssignment.reserve(getMaxNodes()+1);
    for ( size_t i = 0; i < iolV->size(); i++ ) {
        shared_ptr<InputTrace> alpha = make_shared<InputTrace>(iolV->at(i),presentationLayer);
        stateCover.insert(alpha);
        std::shared_ptr<FsmNode> reachedState = *getInitialState()->after(iolV->at(i)).begin();
        stateCoverAssignment[reachedState->getId()] = alpha;
    }
    
    // Test suite is initialised with the state cover
    shared_ptr<Tree> testSuite = getStateCover();
    ConvergenceGraph graph(*this, testSuite);
    

    // distinguish traces in state cover and note already verified nodes
    for ( auto trace : stateCover ) {
        spyhDistinguish(trace->get(),stateCover,testSuite,graph);
    }

    // filter all already verified transitions - i.e. all transitions
    // (s,x,y,s') such that the state cover reaches s by some u and
    // s' by u.x
    // As the state cover generated by getStateCover() is minimal,
    // this holds for all transitions along prefixes of traces in the
    // state cover.
    for (auto trace : stateCover) {
        auto u = trace->get();
        std::shared_ptr<FsmNode> state = getInitialState();
        for (auto x : u) {
            for (auto t : state->getTransitions()) {
                if (t->getLabel()->getInput() == x) {
                    transitions.erase(t);
                    state = t->getTarget();
                    break;
                }
            }
        }
    }

    // verify all remaining transitions
    for (auto transition : transitions) {
        std::shared_ptr<FsmNode> source = transition->getSource();
        std::shared_ptr<FsmNode> target = transition->getTarget();
        int x = transition->getLabel()->getInput();


        auto sourceTrace = stateCoverAssignment[source->getId()];
        auto targetTrace = stateCoverAssignment[target->getId()];

        auto ux = sourceTrace->get();
        ux.push_back(x);
        testSuite->addToRoot(ux);
        graph.add(ux);

        auto v = targetTrace->get();

        // create a copy of the state cover 
        std::unordered_set<std::shared_ptr<InputTrace>> stateCoverCopy;
        for (auto elem : stateCover) {
            stateCoverCopy.insert(elem);
        }

        spyhDistinguishFromSet(ux,v,stateCover,stateCoverCopy,testSuite,graph,numAddStates);

        graph.merge(sourceTrace->get(),x,v);
    }

    
    return testSuite->getIOLists();
}



































