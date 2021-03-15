#ifndef FSM_ENUMERATOR_H_
#define FSM_ENUMERATOR_H_

// TODO: refine includes
#include <memory>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <deque>
#include <stdexcept>

#include "fsm/Fsm.h"

/**
 * A simple enumerator that enumerates all FSMs using a given number of states, inputs and outputs.
 * Returns connected and observable but possibly non-deterministic and/or partial FSMs.
 * 
 * A flag in the constructor defines whether all states of the generated FSMs must be reachable.
 * The generated FSMs may use fewer inputs than the maximum numbers given.
 */
class FsmEnumerator
{
private:
    const int maxInput;
    const int maxOutput;
    const int maxState;
    const bool generateSmallerFsms; 
    
    int candidateNum = 0;
    Fsm fsm;
    
    std::vector<std::vector<std::vector<int>>> currentTable;

    bool updateTable();
    Fsm generateFsmFromTable();

public:
    /**
     * @param maxInput The largest input to be used.
     * @param maxOutput The largest output to be used.
     * @param maxState The largest state index to be used.
     * @param presentationLayer The presentation layer to be used in returned FSMs.
     * @param generateSmallerFsms If true, then all FSMs with up to (maxState+1) reachable states are generated.
     *                            Otherwise, all FSMs with exactly (maxState+1) reachable states are generated.
     */ 
    FsmEnumerator(int maxInput, int maxOutput, int maxState, const std::shared_ptr<FsmPresentationLayer>& presentationLayer, bool generateSmallerFsms = false);

    bool hasNext();
    Fsm getNext();
};


#endif // FSM_ENUMERATOR_H_

