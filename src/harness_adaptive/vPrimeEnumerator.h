#ifndef VPRIMEENUMERATOR_H
#define VPRIMEENUMERATOR_H

#include <memory>
#include "../fsm/Fsm.h"
#include "../fsm/IOTrace.h"
#include "../fsm/IOTraceContainer.h"

using namespace std;

/**
 * Enumerator for calculating all V'' from V' used in adaptive state counting.
 *
 *  Concept: http://phrogz.net/lazy-cartesian-product
 */
class VPrimeEnumerator
{
private:
    vector<vector<shared_ptr<const IOTrace>>> vPrime;
    size_t current;
    size_t length;
    vector<pair<size_t,size_t>> pairs;

public:
    /**
     * Constructor for lazy calculation of V'.
     * @param detStateCover The deterministic state cover for the specification fsm.
     * @param iut The implementation under test.
     */
    VPrimeEnumerator(const vector<IOTraceContainer>& detStateCoverResponses);

    /**
     * Returns the `n`-th V'' for a given index `n`.
     * @param n The given index.
     * @return The `n`-th V''.
     */
    IOTraceContainer getVDoublePrime(size_t n) const;

    /**
     * Returns the next V''. If there is no next V'', an exception will be thrown.
     * @return The next V''.
     */
    IOTraceContainer getNext();

    /**
     * Determines, wether there is some next V''.
     * @return `true`, if there is some next V'', `false`, otherwise.
     */
    bool hasNext() const;

    void reset();

    /**
     * Returns the size of V' (number of elements, aka V'').
     * @return Number of elements in V'.
     */
    size_t size() const { return length; }
};

#endif // VPRIMEENUMERATOR_H
