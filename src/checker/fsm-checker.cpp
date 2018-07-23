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

static model_type_t sutModelType;
static string sutmodelFileName;
static string testSuiteFileName;


static string fsmSutName;


static shared_ptr<FsmPresentationLayer> pl = nullptr;
static shared_ptr<Dfsm> dfsmSut = nullptr;

static bool isDeterministic = true;


/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    cerr << "usage: " << name
    << "sutmodelfile testsuite"
    << endl;
}


/**
 * Parse parameters, stop execution if parameters are illegal.
 *
 * @param argc parameter 1 from main() invocation
 * @param argv parameter 2 from main() invocation
 */
static void parseParameters(int argc, char* argv[]) {
    
    if ( argc < 3 ) {
        printUsage(argv[0]);
        exit(1);
    }
    
    sutmodelFileName = string(argv[1]);
    testSuiteFileName = string(argv[2]);
    
    if ( strstr(sutmodelFileName.c_str(),".csv")  ) {
        sutModelType = FSM_CSV;
    }
    else {
        sutModelType = FSM_JSON;
    }

    fsmSutName = string("SUT");
    
}



static void readSUTModel() {
    
    switch ( sutModelType ) {
        case FSM_CSV:
            isDeterministic = true;
            dfsmSut = make_shared<Dfsm>(sutmodelFileName,fsmSutName);
            pl = dfsmSut->getPresentationLayer();
            break;
            
        case FSM_JSON:
        {
            Reader jReader;
            Value root;
            stringstream document;
            ifstream inputFile(sutmodelFileName);
            document << inputFile.rdbuf();
            inputFile.close();
            
            if ( jReader.parse(document.str(),root) ) {
                dfsmSut = make_shared<Dfsm>(root);
            }
            else {
                cerr << "Could not parse JSON model - exit." << endl;
                exit(1);
            }
        }
            break;
            
        default:
            cerr << "Could not parse this model type - exit." << endl;
            exit(1);
            break;
    }
    
    if ( dfsmSut != nullptr ) {
        dfsmSut->toDot(fsmSutName);
    }
     
}

static void getNextIO(char** p, char** x, char** y) {
    
    *x = NULL;
    *y = NULL;
    
    char* aux;
    
    while ( **p != 0 && **p != '(' ) (*p)++;
    if ( **p == 0 ) return;
    
    (*p)++;
    aux = strchr(*p,'/');
    if ( aux == NULL ) return;
    
    *x = *p;
    *aux = 0;
    
    *p = aux+1;
    aux = strchr(*p,')');
    if ( aux == NULL ) return;
    
    *y = *p;
    *aux = 0;
    
    *p = aux + 1;
    
}

static void executeTestCase(const char* tcId, char* line) {
    
    char* p = line;
    char* x = 0;
    char* y = 0;
    int xInt;
    int yInt;
    string theLine(line);
    vector<int> inVec;
    vector<int> outVec;
    
    
    
    printf("%s",tcId);
    
    
    
    while ( *p ) {
        
        getNextIO(&p,&x,&y);
        
        if ( x != NULL && y != NULL ) {
            xInt = pl->in2Num(x);
            yInt = pl->out2Num(y);
        }
        else {
            cerr << "Could not parse test case " << theLine << endl;
            return;
        }
        
        if ( xInt < 0 ) {
            cerr << "Unknown input " << x
            << " i test case " << theLine << endl;
        }
        else if ( yInt < 0 ) {
            cout << "FAIL: SUT does not produce expected output "
            << y << " occurring in test case " << theLine << endl;
            return;
        }
        
        inVec.push_back(xInt);
        outVec.push_back(yInt);
        
    }
    
    
    InputTrace inTrace(inVec,pl);
    OutputTrace outTrace(outVec,pl);
    
    IOTrace io(inTrace,outTrace);
    
    
    cout << "Check IO Trace " << io << ": ";
    
    if ( dfsmSut->pass(io) ) {
        printf(" PASS\n");
    }
    else {
        printf(" FAIL\n");
    }
    
}


static void executeTestSuite(const char* fname) {
    
    const int lineSize = 100000;
    char* line = (char*)calloc(lineSize,1);
    FILE* f = fopen(fname,"r");
    if ( f == NULL ) {
        fprintf(stderr,"Could not open file %s - exit.\n",fname);
        exit(1);
    }
    
    int tcNum = 0;
    while ( fgets(line,lineSize,f) ) {
        
        size_t len = strlen(line);
        
        // Replace newline by null character
        if ( len > 1 ) {
            line[len-1] = 0;
            char tcId[100];
            *tcId = 0;
            sprintf(tcId,"TC-%d: ",++tcNum);
            executeTestCase(tcId,line);
        }
        
    }
    fclose(f);
    
    
}

int main(int argc, char* argv[])
{
    
    parseParameters(argc,argv);
    readSUTModel();
    executeTestSuite(testSuiteFileName.c_str());
    
    exit(0);
    
}



