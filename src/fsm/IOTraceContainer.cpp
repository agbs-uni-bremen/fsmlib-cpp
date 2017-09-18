#include "IOTraceContainer.h"

using namespace std;

IOTraceContainer::IOTraceContainer(const shared_ptr<FsmPresentationLayer> presentationLayer):
    list(make_shared<std::vector<IOTrace>>()), presentationLayer(presentationLayer)
{

}

IOTraceContainer::IOTraceContainer(shared_ptr<vector<IOTrace>>& list, const shared_ptr<FsmPresentationLayer> presentationLayer):
    list(list), presentationLayer(presentationLayer)
{

}

IOTraceContainer::IOTraceContainer(shared_ptr<IOTrace> trace, const std::shared_ptr<FsmPresentationLayer> presentationLayer):
    list(make_shared<std::vector<IOTrace>>()), presentationLayer(presentationLayer)
{
    list->push_back(*trace);
}


shared_ptr<vector<IOTrace>> IOTraceContainer::getList() const
{
    return list;
}

void IOTraceContainer::addUnique(IOTrace& trc)
{
    for(auto inLst : *list)
    {
        if (trc == inLst){
            return;
        }
    }
    list->push_back(trc);
}


void IOTraceContainer::add(IOTrace& trc)
{
    list->push_back(trc);
}

void IOTraceContainer::addUnique(IOTraceContainer& container)
{
    for (IOTrace& trace : *container.getList())
    {
        addUnique(trace);
    }
}

void IOTraceContainer::add(IOTraceContainer& container)
{
    for (IOTrace& trace : *container.getList())
    {
        add(trace);
    }
}

bool IOTraceContainer::contains(IOTrace& trace) const
{
    for (IOTrace& t : *list)
    {
        if (t == trace)
        {
            return true;
        }
    }
    return false;
}

void IOTraceContainer::concatenate(IOTrace& trace)
{
    for (IOTrace& t : *list)
    {
        t.append(trace);
    }
}

void IOTraceContainer::concatenate(IOTraceContainer& container)
{
    for (IOTrace& thisTrace : *list)
    {
        for (IOTrace& otherTrace : *container.getList())
        {
            thisTrace.append(otherTrace);
        }
    }
}

void IOTraceContainer::remove (IOTrace& trace)
{
    for (auto it = list->begin(); it != list->end();)
    {
        if (*it == trace)
        {
            //TODO fix;
            it = list->erase(it);
        }
        else
        {
            ++it;
        }
    }
}

void IOTraceContainer::remove (IOTraceContainer& container)
{
    for (IOTrace& trace : *container.getList())
    {
        remove(trace);
    }
}

std::ostream & operator<<(std::ostream & out, const IOTraceContainer & iot)
{
    out << "{ ";

    bool isFirst = true;
    for (IOTrace& trace : *iot.list)
    {
        if (!isFirst)
        {
            out << "," << std::endl << "  ";
        }

        out << trace;
        isFirst = false;
    }
    out << " }";
    return out;
}
