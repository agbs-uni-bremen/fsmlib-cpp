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


void test() {
    srand(getRandomSeed());

    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(1,6,1);
    //shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("../../../resources/lee94_no_pds.fsm",pl,"lee94_no_pds");
    shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("../../../resources/hierons_multicomp.fsm",pl,"hierons_multicomp");

    dfsm->hieronsDMethodOnMinimisedDfsm(true);
}
void testRandom(const int numStates,const int numInput,const int numOutput)
{
    int numberOfTests = 100;

    for(int i=0;i<numberOfTests;++i) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput, numStates, numOutput);
        auto dfsm = make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl);
        shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());

        dfsm->hieronsDMethodOnMinimisedDfsm(true);
        cout << "next one " << endl;
    }
}

int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    test();
    //testRandom(6,1,1);
    //testCustomAds();
    //testRandomAds(15, 3, 3);
    //testRandomFaultCoverage(30, 6, 6);

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