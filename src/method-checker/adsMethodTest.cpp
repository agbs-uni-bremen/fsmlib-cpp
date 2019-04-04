#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"
#include "trees/InputOutputTree.h"
#include "trees/SplittingTree.h"
#include "trees/IOTreeContainer.h"

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

void testCustomAds()
{

    cout << "affirm that there is an ads for lee94_no_pds.fsm" << endl;
    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(1,6,1);
    shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("../../../resources/lee94_no_pds.fsm",pl,"lee94_no_pds");

    auto adaptiveDistinguishingSequence = dfsm->createAdaptiveDistinguishingSequence();
    if(!adaptiveDistinguishingSequence) {
        cout << "ads does not exist" << endl;
    } else {
        cout << "algorithm delivers ads with following hsi:" << endl;
        //cout << *adaptiveDistinguishingSequence << endl;
        auto ioll = adaptiveDistinguishingSequence->getHsi();
        for(auto& ioLst:*ioll) {
            string iolst_str = "";
            for(int input:ioLst) {
                iolst_str += to_string(input) + ".";
            }
            iolst_str = iolst_str.substr(0,iolst_str.size()-1);
            cout << iolst_str << endl;
        }
    }

    auto adsList = make_shared<vector<shared_ptr<InputOutputTree>>>();
    adsList->push_back(adaptiveDistinguishingSequence);
    IOTreeContainer adaptiveTestCases(adsList,dfsm->getPresentationLayer());

    auto nodesA = dfsm->getNodes();
    auto& nodesB = nodesA;
    cout << "is ads? -> " << dfsm->distinguishesAllStates(nodesB,nodesB,adaptiveTestCases) << endl;

    shared_ptr<FsmPresentationLayer> plOther = createPresentationLayer(1,3,1);
    auto dfsmOther = make_shared<Dfsm>("../../../resources/lee94_no_ads.fsm",plOther,"lee94_no_ads");

    auto otherNodesA = dfsmOther->getNodes();
    auto& otherNodesB = otherNodesA;
    cout << "is not ads of lee94_no_ads? -> " << !dfsmOther->distinguishesAllStates(otherNodesB,otherNodesB,adaptiveTestCases)
    << endl << endl;

    cout << "affirm that there is no ads for lee_94_no_ads.fsm -> " << !dfsmOther->createAdaptiveDistinguishingSequence() << endl << endl;

    cout << "affirm that there is no pds for lee_94_no_pds.fsm -> " << (dfsmOther->createDistinguishingSequence().size() == 0) << endl << endl;
}

void testRandomAds(const int numStates,const int numInput,const int numOutput)
{

    int numberOfTests = 100;
    int numberOfCorrectAds = 0;
    float ratio = (float) numOutput/ (float) numStates;
    cout << "test " << numberOfTests << " dfsm with ads, if the created ads really distinguishes each state from every other state" << endl;
    cout << "number of states: " << numStates << endl;
    cout << "number of inputs: " << numInput << endl;
    cout << "number of outputs: " << numOutput << endl;
    cout << "ratio (|outputs|/|states|): " << ratio << endl;

    for(int i=0;i<numberOfTests;++i) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput,numStates,numOutput);
        auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
        shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());

        auto adaptiveDistinguishingSequence = dfsmMin->createAdaptiveDistinguishingSequence();

        while (!adaptiveDistinguishingSequence) {
            dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
            dfsmMin = make_shared<Dfsm>(dfsm->minimise());
            adaptiveDistinguishingSequence = dfsmMin->createAdaptiveDistinguishingSequence();
        }

        auto adsList = make_shared<vector<shared_ptr<InputOutputTree>>>();
        adsList->push_back(adaptiveDistinguishingSequence);
        IOTreeContainer adaptiveTestCases(adsList,dfsm->getPresentationLayer());

        auto nodesA = dfsm->getNodes();
        auto& nodesB = nodesA;
        bool isAds = dfsm->distinguishesAllStates(nodesB,nodesB,adaptiveTestCases);

        shared_ptr<FsmPresentationLayer> plOther = createPresentationLayer(1,3,1);
        auto dfsmOther = make_shared<Dfsm>("../../../resources/lee94_no_ads.fsm",plOther,"lee94_no_ads");

        auto otherNodesA = dfsmOther->getNodes();
        auto& otherNodesB = otherNodesA;
        bool isNotAds = !dfsm->distinguishesAllStates(otherNodesB,otherNodesB,adaptiveTestCases);

        if(isAds && isNotAds) {
            numberOfCorrectAds++;
        }
    }
    cout << "number of correct ads: " << numberOfCorrectAds << endl;
    float percent = ((float) numberOfCorrectAds / (float) numberOfTests) * 100;
    cout << "percentage of correct ads: " << percent << endl << endl;


}

void testRandomApplicability(const int numStates,const int numInput,const int numOutput)
{

    int numberOfTests = 100;
    int numberOfDfsmWithAds = 0;
    int numberOfDfsmWithPds = 0;

    unsigned long totalMinLengthAds = 0,
        totalMaxLengthAds = 0;
    long double totalAvgLengthPds = 0,
        totalAvgLengthAds = 0;

    float ratio = (float) numOutput/ (float) numStates;
    cout << "test " << numberOfTests << " dfsm for the existence of an ads" << endl;
    cout << "number of states: " << numStates << endl;
    cout << "number of inputs: " << numInput << endl;
    cout << "number of outputs: " << numOutput << endl;
    cout << "ratio (|outputs|/|states|): " << ratio << endl;

    bool adsLongerThanPds = false;
    for(int i=0;i<numberOfTests;++i) {
        unsigned long lengthPds = 0,
            minLengthAds = 0,
            maxLengthAds = 0;
        long double avgLengthAds = 0;

        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput,numStates,numOutput);
        auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
        shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());

        auto adaptiveDistinguishingSequence = dfsmMin->createAdaptiveDistinguishingSequence();
        auto presetDistinguishingSequence = dfsmMin->createDistinguishingSequence();

        if(adaptiveDistinguishingSequence) {
            assert(!presetDistinguishingSequence.empty());
            numberOfDfsmWithAds++;
            numberOfDfsmWithPds++;

            lengthPds = presetDistinguishingSequence.size();

            auto inputLists = adaptiveDistinguishingSequence->getInputLists();
            auto ioLsts = inputLists.getIOLists();
            avgLengthAds = 0;
            maxLengthAds = 0;
            minLengthAds = 0;
            for(auto& ioLst:*ioLsts) {
                assert(ioLst.size() >0);
                if(minLengthAds > 0 )
                    minLengthAds = ioLst.size() < minLengthAds?ioLst.size():minLengthAds;
                else
                    minLengthAds = ioLst.size();
                maxLengthAds = ioLst.size() > maxLengthAds?ioLst.size():maxLengthAds;
                avgLengthAds += ioLst.size();
                assert(minLengthAds>0);
            }
            avgLengthAds = avgLengthAds / (long double) ioLsts->size();

            //if(avgLengthAds > lengthPds) {
             //   cout << "pds with size " << lengthPds << " smaller than ads with avg length " << avgLengthAds << endl;
           // }
            adsLongerThanPds |= avgLengthAds > lengthPds;

            totalAvgLengthPds += lengthPds;
            totalAvgLengthAds += avgLengthAds;

            totalMaxLengthAds = maxLengthAds > totalMaxLengthAds?maxLengthAds:totalMaxLengthAds;
            if(totalMinLengthAds > 0)
                totalMinLengthAds = minLengthAds < totalMinLengthAds?minLengthAds:totalMinLengthAds;
            else
                totalMinLengthAds = minLengthAds;
        }

    }
    totalAvgLengthPds = totalAvgLengthPds / (long double) numberOfDfsmWithPds;
    totalAvgLengthAds = totalAvgLengthAds / (long double) numberOfDfsmWithAds;

    cout << "is there an ads with avg length greater than a pds? -> " << adsLongerThanPds << endl;
    cout << "average length pds : " << totalAvgLengthPds << endl;
    cout << "average length ads : " << totalAvgLengthAds << endl;
    cout << "min length ads: " << totalMinLengthAds << endl;
    cout << "max length ads: " << totalMaxLengthAds << endl;

    cout << "number of dfsm with ads: " << numberOfDfsmWithAds << endl;
    cout << "number of dfsm with pds: " << numberOfDfsmWithPds << endl << endl;
}

void testRandomFaultCoverage(const int numStates,const int numInput,const int numOutput)
{

    srand(getRandomSeed());

    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput,numStates,numOutput);
    auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
    shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());
    dfsmMin->toDot("dfsm_min_fc");
    IOListContainer ts = dfsmMin->dMethodOnMinimisedDfsm(0, true);


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
    cout << "Number of mutants passing D-Method(ADS) Testsuite: " << count_dMethod_pass << endl;
    cout << "Testsuite size: " << ts.getFlatSize() << endl;
}

int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    //testCustomAds();
    //testRandomAds(15, 3, 3);
    testRandomFaultCoverage(30, 6, 6);

    /*
    testRandomApplicability(30,2,2);
    testRandomApplicability(30,4,4);
    testRandomApplicability(30,6,6);
    testRandomApplicability(30,10,10);

    testRandomApplicability(50,3,3);
    testRandomApplicability(50,6,6);
    testRandomApplicability(50,10,10);
    testRandomApplicability(50,17,17);
*/

}