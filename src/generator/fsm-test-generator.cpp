/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <string.h>

#include "interface/FsmPresentationLayer.h"
#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/IOTrace.h"
#include "fsm/FsmPrintVisitor.h"
#include "fsm/FsmSimVisitor.h"
#include "fsm/FsmOraVisitor.h"
#include "trees/IOListContainer.h"
#include "trees/OutputTree.h"
#include "trees/TestSuite.h"
#include "json/json.h"


using namespace std;
using namespace Json;


/**
 *   Program execution parameters and associated types
 */
typedef enum {
    FSM_CSV,
    FSM_JSON,
    FSM_BASIC
} model_type_t;

typedef enum {
    WMETHOD,
    WPMETHOD
} generation_method_t;

static model_type_t modelType;
static string modelFile;
static string plStateFile;
static string plInputFile;
static string plOutputFile;
static string fsmName;
static string testSuiteFileName;
static generation_method_t genMethod;
static unsigned int numAddStates;

static shared_ptr<FsmPresentationLayer> pl = nullptr;
static shared_ptr<Dfsm> dfsm = nullptr;
static shared_ptr<Fsm> fsm = nullptr;

static bool isDeterministic = false;


/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    cerr << "usage: " << name << " [-w] [-n fsmname] [-p infile outfile statefile] [-a additionalstates] [-t testsuitename] modelfile" << endl;
}

/**
 *  Determine the model type of a model specified in an *.fsm file.
 *
 *  @param modelFile Name of the *.fsm file containing the model
 *
 *  @return FSM_BASIC, if the model file contains the low-level
 *                     encoding.
 *
 *  @return FSM_JSON, if the model contains the JSON encoding.
 */
static model_type_t getModelType(const string& modelFile) {
    
    model_type_t t = FSM_BASIC;
    
    ifstream inputFile(modelFile);
    string line;
    getline(inputFile,line);
    // Basic encoding does not contain any { or [
    if ( line.find("{") != string::npos or
        line.find("[") != string::npos) {
        t = FSM_JSON;
    }
    
    inputFile.close();
    
    return t;
    
}

/**
 * Parse parameters, stop execution if parameters are illegal.
 *
 * @param argc parameter 1 from main() invocation
 * @param argv parameter 2 from main() invocation
 */
static void parseParameters(int argc, char* argv[]) {
    
    // Set default parameters
    genMethod = WPMETHOD;
    fsmName = string("FSM");
    testSuiteFileName = string("testsuite.txt");
    numAddStates = 0;
    
    for ( int p = 1; p < argc; p++ ) {
        if ( strcmp(argv[p],"-w") == 0 ) {
            genMethod = WMETHOD;
        }
        else if ( strcmp(argv[p],"-n") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing FSM name" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                fsmName = string(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-t") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing test suite name" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                testSuiteFileName = string(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-a") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing number of additional states" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                numAddStates = atoi(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-p") == 0 ) {
            if ( argc < p+4 ) {
                cerr << argv[0] << ": missing presentation layer files" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                plInputFile = string(argv[++p]);
                plOutputFile = string(argv[++p]);
                plStateFile = string(argv[++p]);
            }
        }
        else if ( strstr(argv[p],".csv")  ) {
            modelFile = string(argv[p]);
            modelType = FSM_CSV;
        }
        else if ( strstr(argv[p],".fsm")  ) {
            modelFile = string(argv[p]);
            modelType = getModelType(modelFile);
        }
        else if ( strstr(argv[p],".@todo")  ) {
            modelFile = string(argv[p]);
            modelType = FSM_JSON;
        }
        else {
            cerr << argv[0] << ": illegal parameter `" << argv[p] << "'" << endl;
            printUsage(argv[0]);
            exit(1);
        }
    }
    
    if ( modelFile.empty() ) {
        cerr << argv[0] << ": missing model file" << endl;
        printUsage(argv[0]);
        exit(1);
    }
    
}


/**
 *   Instantiate DFSM or FSM from input file according to
 *   the different input formats which are supported.
 */
static void readModel() {
    
    switch ( modelType ) {
        case FSM_CSV:
            isDeterministic = true;
            dfsm = make_shared<Dfsm>(modelFile,fsmName);
            pl = dfsm->getPresentationLayer();
            break;
            
        case FSM_JSON:
        {
            Reader jReader;
            Value root;
            stringstream document;
            ifstream inputFile(modelFile);
            document << inputFile.rdbuf();
            
            if ( jReader.parse(document.str(),root) ) {
                dfsm = make_shared<Dfsm>(root);
            }
            else {
                cerr << "Could not parse JSON model - exit." << endl;
                exit(1);
            }
        }
            break;
            
        case FSM_BASIC:
            if ( plStateFile.empty() ) {
                pl = make_shared<FsmPresentationLayer>();
            }
            else {
                pl = make_shared<FsmPresentationLayer>(plInputFile,plOutputFile,plStateFile);
            }
            fsm = make_shared<Fsm>(modelFile,pl,fsmName);
            if ( fsm->isDeterministic() ) {
                isDeterministic = true;
                dfsm = make_shared<Dfsm>(modelFile,pl,fsmName);
                fsm = nullptr;
            }
            break;
    }
    
    if ( fsm != nullptr ) {
        fsm->toDot(fsmName);
    }
    else if ( dfsm != nullptr ) {
        dfsm->toDot(fsmName);
        dfsm->toCsv(fsmName);
    }
    
}

static void generateTestSuite() {
    
    shared_ptr<TestSuite> testSuite =
    make_shared<TestSuite>();
    
    switch ( genMethod ) {
        case WMETHOD:
            if ( dfsm != nullptr ) {
                IOListContainer iolc = dfsm->wMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
            }
            else {
                IOListContainer iolc = fsm->wMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(fsm->apply(*itrc));
                }
            }
            break;
        case WPMETHOD:
            if ( dfsm != nullptr ) {
                IOListContainer iolc = dfsm->wpMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
            }
            else {
                IOListContainer iolc = fsm->wpMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(fsm->apply(*itrc));
                }
            }
            break;
    }
    
    testSuite->save(testSuiteFileName);
    
}

int main(int argc, char* argv[])
{
    
    parseParameters(argc,argv);
    readModel();
    generateTestSuite();
    
    exit(0);
    
}



