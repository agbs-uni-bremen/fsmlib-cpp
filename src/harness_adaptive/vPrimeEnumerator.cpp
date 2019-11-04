#include "vPrimeEnumerator.h"
#include "../utils/Logger.hpp"
#include "../fsm/FsmNode.h"

VPrimeEnumerator::VPrimeEnumerator(const vector<IOTraceContainer>& detStateCoverResponses):
    vPrime(vector<vector<shared_ptr<const IOTrace>>>())
{

    
    current = 0;
    LOG("VERBOSE_1") << "detStateCoverResponses:" << std::endl;
    
    // TODO: verify
    for (unsigned int i = 0; i < detStateCoverResponses.size(); ++i) {
        vector<shared_ptr<const IOTrace>> temp;
        for (auto it = detStateCoverResponses.at(i).cbegin(); it != detStateCoverResponses.at(i).cend(); ++it)
        {
            temp.push_back(*it);
        }
        
    }

    size_t i = vPrime.size();
    pairs.resize(i);
    size_t f = 1;
    size_t l = 0;
    while (i > 0)
    {
        --i;
        l = vPrime.at(i).size();
        pairs.at(i) = (make_pair(f, l));
        f *= l;
    }
    length = f;
}


IOTraceContainer VPrimeEnumerator::getVDoublePrime(size_t n) const
{
    size_t i = vPrime.size();
    IOTraceContainer c;
    while (i > 0)
    {
        --i;
        size_t idx = n / pairs.at(i).first % pairs.at(i).second;
        c.add(vPrime.at(i).at(idx));
    }
    return c;
}

IOTraceContainer VPrimeEnumerator::getNext()
{
    return getVDoublePrime(++current);
}

bool VPrimeEnumerator::hasNext() const
{
    return current < length;
}

void VPrimeEnumerator::reset()
{
    current = 0;
}
