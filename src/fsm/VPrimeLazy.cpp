#include "VPrimeLazy.h"
#include "utils/Logger.hpp"
#include "fsm/FsmNode.h"

VPrimeLazy::VPrimeLazy(const InputTraceSet& detStateCover, const Fsm& iut):
    allPossibleTraces(vector<vector<shared_ptr<const IOTrace>>>())
{
    current = 0;
    LOG("VERBOSE_1") << "allPossibleTraces:" << std::endl;
    for (const shared_ptr<InputTrace>& input : detStateCover)
    {
        vector<shared_ptr<OutputTrace>> producedOutputs;
        vector<shared_ptr<FsmNode>> reached;
        iut.getInitialState()->getPossibleOutputs(*input, producedOutputs, reached);
        vector<shared_ptr<const IOTrace>> producedIOTraces;
        for (size_t j = 0; j < producedOutputs.size(); ++j)
        {
            LOG("VERBOSE_1") << *producedOutputs.at(j) << std::endl;
            producedIOTraces.push_back(make_shared<const IOTrace>(*input, *producedOutputs.at(j), reached.at(j)));
        }
        allPossibleTraces.push_back(producedIOTraces);
        LOG("VERBOSE_1") << "--------------" << std::endl;
    }

    size_t i = allPossibleTraces.size();
    pairs.resize(i);
    size_t f = 1;
    size_t l = 0;
    while (i > 0)
    {
        --i;
        l = allPossibleTraces.at(i).size();
        pairs.at(i) = (make_pair(f, l));
        f *= l;
    }
    length = f;
}


IOTraceContainer VPrimeLazy::getVDoublePrime(size_t n) const
{
    size_t i = allPossibleTraces.size();
    IOTraceContainer c;
    while (i > 0)
    {
        --i;
        size_t idx = n / pairs.at(i).first % pairs.at(i).second;
        c.add(allPossibleTraces.at(i).at(idx));
    }
    return c;
}

IOTraceContainer VPrimeLazy::getNext()
{
    return getVDoublePrime(++current);
}

bool VPrimeLazy::hasNext() const
{
    return current < length;
}

void VPrimeLazy::reset()
{
    current = 0;
}
