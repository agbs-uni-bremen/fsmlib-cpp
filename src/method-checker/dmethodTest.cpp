#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"

using namespace std;

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


int main(int argc, char* argv[])
{

    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(3,5,3);

    shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("Dfsm",20,4,7,pl);

    dfsm->minimise();

    vector<int> distinguishingSequence = dfsm->createDistinguishingSequence();

    cout << "ds-length: " << distinguishingSequence.size() << endl;
}

