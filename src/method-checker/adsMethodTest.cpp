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
        cout << "algorithm delivers ads" << endl;
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
    cout << "is not ads of lee94_no_ads? -> " << !dfsmOther->distinguishesAllStates(otherNodesB,otherNodesB,adaptiveTestCases) << endl;

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

    for(int i=0;i<100;++i) {
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
    cout << "percentage of correct ads: " << percent << endl;


}

void testRandomApplicability(const int numStates,const int numInput,const int numOutput)
{

    int numberOfTests = 100;
    int numberOfDfsmWithAds = 0;
    int numberOfDfsmWithPds = 0;
    int totalAvgLengthPds = 0;
    int totalMinLengthAds = 0;
    int totalMaxLengthAds = 0;
    int totalAvgLengthAds = 0;

    float ratio = (float) numOutput/ (float) numStates;
    cout << "test " << numberOfTests << " dfsm for the existence of an ads" << endl;
    cout << "number of states: " << numStates << endl;
    cout << "number of inputs: " << numInput << endl;
    cout << "number of outputs: " << numOutput << endl;
    cout << "ratio (|outputs|/|states|): " << ratio << endl;

    for(int i=0;i<100;++i) {
        int lengthPds = 0,
            minLengthAds = 0,
            maxLengthAds = 0,
            avgLengthAds = 0;

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
            for(auto ioLst:*ioLsts) {

            }
        }

    }
    cout << "number of dfsm with ads: " << numberOfDfsmWithAds << endl;
    cout << "number of dfsm with pds: " << numberOfDfsmWithPds << endl;
}

int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    testCustomAds();
    //testRandomAds(30, 4, 4);
    testRandomApplicability(30,2,2);
    testRandomApplicability(30,4,4);
    testRandomApplicability(30,6,6);
    testRandomApplicability(30,10,10);

    testRandomApplicability(50,3,3);
    testRandomApplicability(50,6,6);
    testRandomApplicability(50,10,10);
    testRandomApplicability(50,17,17);


}