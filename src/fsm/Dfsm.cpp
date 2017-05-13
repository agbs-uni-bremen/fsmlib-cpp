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



shared_ptr<FsmPresentationLayer> Dfsm::createPresentationLayerFromCsvFormat(const std::string & fname) {
    
    // key comparator for usage in string sets
    class setCmp
    {
    public:
        bool operator()(const string& s1, const string& s2)
        {
            return (s1.compare(s2) < 0);
        }
    };
    // Set containing all output events defined so far
    set<string,setCmp> outStrSet;

    
    vector<std::string> in2String;
    vector<std::string> out2String;
    vector<std::string> state2String;
    
    ifstream inputFile(fname);
    
    if ( not inputFile.is_open() ) {
        cout << "Unable to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Insert a "no operation" output as first member of the
    // output alphabet
    out2String.push_back("_nop");
    outStrSet.insert(out2String[0]);
    
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
    
    for ( auto s : outStrSet ) {
        out2String.push_back(s);
    }
    
    maxInput = in2String.size() - 1;
    maxOutput = out2String.size() - 1;
    maxState = state2String.size() - 1;
    initStateIdx = 0;
    
    return make_shared<FsmPresentationLayer>(in2String,out2String,state2String);
    
}


void Dfsm::createDfsmTransitionGraph(const std::string& fname) {
    
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



Dfsm::Dfsm(const std::string & fname,
           const std::string & fsmName) : Fsm(nullptr)   {
    
    name = fsmName;
    presentationLayer = createPresentationLayerFromCsvFormat(fname);
    createDfsmTransitionGraph(fname);
    
}

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



Dfsm::Dfsm(const Json::Value& fsmExport) :
Fsm()
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
    for (unsigned int index = 0; index < inputs.size(); ++index ) {
        in2String.push_back(inputs[index].asString());
    }
    
    // iterate over all outputs
    vector<string> out2String;
    for (unsigned int index = 0; index < outputs.size(); ++index ) {
        out2String.push_back(outputs[index].asString());
    }
    // Add a NOP output for the case where the FSM is incomplete
    out2String.push_back("_nop");
    int theNopNo = out2String.size() - 1;
    
    
    
    // iterate over all states, insert initial state at index 0
    // of the state2String vector.
    vector<string> state2String;
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
        if (state["initial"].asBool()) {
            state2String.push_back(state["name"].asString());
            break;
        }
    }
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
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
    maxInput = in2String.size() - 1;
    maxOutput = out2String.size() - 1;
    maxState = state2String.size() - 1;
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
    for (unsigned int index = 0; index < transitions.size(); ++index ) {
        Json::Value transition = transitions[index];
        
        
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
        for (unsigned int inidx = 0; inidx < inputlist.size(); ++inidx ) {
            
            string xString(inputlist[inidx].asString());
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
            for (unsigned int ridx = 0; ridx < satisfies.size(); ++ridx ) {
                tr->addSatisfies(satisfies[ridx].asString());
            }
            
            srcNode->addTransition(tr);
        }
        
    }
    
    // Add requirements to nodes
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
        string nodeName(state["name"].asString());
        
        for ( auto n : nodes ) {
            if ( n->getName() == nodeName ) {
                Json::Value satisfies = state["requirements"];
                for (unsigned int ridx = 0; ridx < satisfies.size(); ++ridx ) {
                    n->addSatisfies(satisfies[ridx].asString());
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
        for ( auto tr : n->getTransitions() ) {
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
           const std::shared_ptr<FsmPresentationLayer> pl) :
Fsm()
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
    for (unsigned int index = 0; index < inputs.size(); ++index ) {
        string theInput(inputs[index].asString());
        if ( pl->in2Num(theInput) < 0 ) {
            in2String.push_back(theInput);
        }
    }
    
    // iterate over all outputs
    vector<string> out2String(pl->getOut2String());
    for (unsigned int index = 0; index < outputs.size(); ++index ) {
        string theOutput(outputs[index].asString());
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
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
        if (state["initial"].asBool()) {
            state2String.push_back(state["name"].asString());
            break;
        }
    }
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
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
    maxInput = in2String.size() - 1;
    maxOutput = out2String.size() - 1;
    maxState = state2String.size() - 1;
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
    for (unsigned int index = 0; index < transitions.size(); ++index ) {
        Json::Value transition = transitions[index];
        
        
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
        for (unsigned int inidx = 0; inidx < inputlist.size(); ++inidx ) {
            
            string xString(inputlist[inidx].asString());
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
            for (unsigned int ridx = 0; ridx < satisfies.size(); ++ridx ) {
                tr->addSatisfies(satisfies[ridx].asString());
            }
            
            srcNode->addTransition(tr);
        }
        
    }
    
    // Add requirements to nodes
    for (unsigned int index = 0; index < states.size(); ++index ) {
        Json::Value state = states[index];
        string nodeName(state["name"].asString());
        
        for ( auto n : nodes ) {
            if ( n->getName() == nodeName ) {
                Json::Value satisfies = state["requirements"];
                for (unsigned int ridx = 0; ridx < satisfies.size(); ++ridx ) {
                    n->addSatisfies(satisfies[ridx].asString());
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
        for ( auto tr : n->getTransitions() ) {
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
    
    return pMin->toFsm(name, maxOutput);
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
        }
    }
    
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
    Fsm fMin = minimiseObservableFSM();
    return fMin.wpMethod(numAddStates);
    
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














