#include "VPrimeLazy.h"
#include "logging/easylogging++.h"
#include "fsm/FsmNode.h"

VPrimeLazy::VPrimeLazy(const vector<shared_ptr<InputTrace>>& detStateCover, const Fsm& iut):
    allPossibleTraces(vector<vector<shared_ptr<IOTrace>>>())
{
    current = 0;
    VLOG(1) << "allPossibleTraces:";
    for (size_t i = 0; i < detStateCover.size(); ++i)
    {
        const shared_ptr<InputTrace>& input = detStateCover.at(i);
        vector<shared_ptr<OutputTrace>> producedOutputs;
        vector<shared_ptr<FsmNode>> reached;
        iut.getInitialState()->getPossibleOutputs(*input, producedOutputs, reached);
        vector<shared_ptr<IOTrace>> producedIOTraces;
        for (size_t j = 0; j < producedOutputs.size(); ++j)
        {
            VLOG(1) << *producedOutputs.at(j);
            producedIOTraces.push_back(make_shared<IOTrace>(*input, *producedOutputs.at(j), reached.at(j)));
        }
        allPossibleTraces.push_back(producedIOTraces);
        VLOG(1) << "--------------";
    }

    size_t i = allPossibleTraces.size();
    dm.resize(i);
    size_t f = 1;
    size_t l = 0;
    while (i > 0)
    {
        --i;
        l = allPossibleTraces.at(i).size();
        dm.at(i) = (make_pair(f, l));
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
        size_t idx = n / dm.at(i).first % dm.at(i).second;
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
