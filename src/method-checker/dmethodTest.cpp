#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"

using namespace std;

unsigned int getRandomSeed() {

    return static_cast<unsigned int>
    (std::chrono::high_resolution_clock::now().time_since_epoch().count());

}

shared_ptr<FsmPresentationLayer> createPresentationLayer(const size_t maxInput, const size_t refMin, const size_t maxOutput)
{
    vector<string> in2string = vector<string>();
    for (size_t i = 0; i <= maxInput; i++)
    {
        in2string.push_back(to_string(i));
    }
    vector<string> out2string = vector<string>();
    for (size_t i = 0; i <= maxOutput; i++)
    {
        out2string.push_back(to_string(i));
    }
    vector<string>state2string = vector<string>();
    for (size_t i = 0; i < refMin; i++)
    {
        state2string.push_back(to_string(i));
    }
    shared_ptr<FsmPresentationLayer> pl{
            new FsmPresentationLayer(
                    in2string,
                    out2string,
                    state2string)};

    return pl;
}


void checkDsLength()
{
    //shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(1,6,1);

    //shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("../../../resources/lee94_no_pds.fsm",pl,"lee94");
    shared_ptr<FsmPresentationLayer> pl;
    shared_ptr<Dfsm> dfsm;
    vector<int> distinguishingSequence;

    while(!distinguishingSequence.size()) {
        pl = createPresentationLayer(3,30,3);
        dfsm = make_shared<Dfsm>("Dfsm", 30, 3, 3, pl);

        dfsm->minimise();

        distinguishingSequence = dfsm->createDistinguishingSequence();
    }
    distinguishingSequence = dfsm->createDistinguishingSequence();
    dfsm->toDot("dfsm_min");
    InputTrace ds(distinguishingSequence,dfsm->getPresentationLayer());

    cout << "ds-length: " << ds.size() << endl;
    auto nodesA = dfsm->getNodes();
    auto& nodesB = nodesA;
    cout << "is ds:" << dfsm->distinguishesAllStates(nodesB,ds) << endl;
    /*for(auto& node: nodesB) {
        vector<shared_ptr<OutputTrace>> outputs;
        node->getPossibleOutputs(ds,outputs);
        cout << "output for node " << node->getId() << endl;
        for(auto& output:outputs) {
            cout << *output << endl;
        }
    }*/
}

void testRandomFaultCoverage(const int numStates,const int numInput,const int numOutput)
{

    srand(getRandomSeed());

    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput,numStates,numOutput);
    auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
    shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());
    dfsmMin->toDot("dfsm_min_fc");
    IOListContainer ts = dfsmMin->dMethodOnMinimisedDfsm(0);


    vector<IOTrace> dTestsuite;
    for(vector<int> v: *ts.getIOLists()) {
        InputTrace i(v,pl);
        IOTrace io = dfsmMin->applyDet(i);
        dTestsuite.push_back(io);
    }

    unsigned int count_equals=0,
            count_dMethod_pass=0,
            count_mutants=1000;
    if(dTestsuite.size() > 0) {
        for (int i = 0; i < count_mutants; i++) {
            unsigned int numOutputFaults = rand() % ((dfsmMin->getMaxInput() + 1) * dfsmMin->getMaxNodes()),
                    numTransitionFaults = rand() % ((dfsmMin->getMaxInput() + 1) * dfsmMin->getMaxNodes());

            shared_ptr<Dfsm> mutant = dfsmMin->createMutant("Mutant", numOutputFaults, numTransitionFaults);


            bool result = dfsmMin->equivalenceCheck(*mutant);
            if (result) count_equals++;

            result = true;
            for (IOTrace io: dTestsuite) {
                result &= mutant->pass(io);
            }
            if (result) count_dMethod_pass++;
        }
        /*cout << "Printing the d-testsuite:" << endl;
        for(IOTrace io:dTestsuite) {
            cout << io << endl;
        }
        cout << "the ds:" << endl;
        auto ds = dfsmMin->createDistinguishingSequence();
        for(auto i:ds) {
            cout << i << ".";
        }
        cout << endl;
        InputTrace dstrace(ds,dfsm->getPresentationLayer());
        for(auto& node: dfsmMin->getNodes()) {
            vector<shared_ptr<OutputTrace>> outputs;
            node->getPossibleOutputs(dstrace,outputs);
            cout << "output for node " << node->getId() << endl;
            for(auto& output:outputs) {
                cout << *output << endl;
            }
        }*/
    }
    cout << "|Fault coverage test|" << endl;
    cout << "Number of mutants: " << count_mutants << endl;
    cout << "Number of equal Mutants: " << count_equals << endl;
    cout << "Number of mutants passing D-Method Testsuite: " << count_dMethod_pass << endl;
    cout << "Testsuite size: " << dTestsuite.size() << endl;
}

void testRandomApplicability(const int numStates,const int numInput,const int numOutput)
{

    int numberOfTests = 100;
    int numberOfDfsmWithPds = 0;
    float ratio = (float) numOutput/ (float) numStates;
    cout << "test " << numberOfTests << " dfsm for the existence of an ads" << endl;
    cout << "number of states: " << numStates << endl;
    cout << "number of inputs: " << numInput << endl;
    cout << "number of outputs: " << numOutput << endl;
    cout << "ratio (|outputs|/|states|): " << ratio << endl;

    for(int i=0;i<100;++i) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput,numStates,numOutput);
        auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
        shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());

        auto distinguishingSequence = dfsmMin->createDistinguishingSequence();

        if(!distinguishingSequence.empty()) {
            numberOfDfsmWithPds++;
        }

    }
    cout << "number of dfsm with pds: " << numberOfDfsmWithPds << endl;
}

int main(int argc, char* argv[])
{

    //Logging
    //LogCoordinator& logger = LogCoordinator::getStandardLogger();
    //logger.setDefaultStream(cout);

    //testRandomFaultCoverage(30,6,6);
    testRandomApplicability(30,7,7);

}

