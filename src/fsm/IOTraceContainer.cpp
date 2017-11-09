#include "IOTraceContainer.h"

using namespace std;

IOTraceContainer::IOTraceContainer():
    list(make_shared<std::vector<IOTrace>>())
{

}

IOTraceContainer::IOTraceContainer(shared_ptr<vector<IOTrace>>& list):
    list(list)
{

}

IOTraceContainer::IOTraceContainer(shared_ptr<IOTrace> trace):
    list(make_shared<std::vector<IOTrace>>())
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

void IOTraceContainer::addUniqueRemovePrefixes(const IOTrace& trc)
{
    removeRealPrefixes(trc);
    for (IOTrace& trace : *list)
    {
        if (trc == trace || trace.isPrefix(trc))
        {
            // The list contains either the trace itself or the trace as prefix of another trace.
            return;
        }
    }
    // No other trace has been found that contains the trace as prefix.
    list->push_back(trc);
}

void IOTraceContainer::addUniqueRemovePrefixes(const IOTraceContainer& cont)
{
    for (IOTrace& trc : *cont.getList())
    {
        addUniqueRemovePrefixes(trc);
    }
}

void IOTraceContainer::removeRealPrefixes(const IOTrace & trc)
{
    auto it = list->begin();
    while (it  != list->end())
    {
        if (it->isPrefixOf(trc))
        {
            it = list->erase(it);
        }
        else
        {
            ++it;
        }
    }

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

void IOTraceContainer::addUnique(OutputTree& tree)
{
    std::vector<IOTrace> iOTraces;
    tree.toIOTrace(iOTraces);
    for (IOTrace& trace : iOTraces)
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

void IOTraceContainer::add(OutputTree& tree)
{
    std::vector<IOTrace> iOTraces;
    tree.toIOTrace(iOTraces);
    for (IOTrace& trace : iOTraces)
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

void IOTraceContainer::concatenateToFront(InputTrace& inputTrace, OutputTrace outputTrace)
{
    IOTrace newIOTrace(inputTrace, outputTrace);
    concatenateToFront(newIOTrace);
}

void IOTraceContainer::concatenateToFront(IOTrace& iOTrace)
{
    for (IOTrace& iOT : *list)
    {
        iOT.prepend(iOTrace);
    }
}

void IOTraceContainer::clear()
{
    list = make_shared<std::vector<IOTrace>>();
}

void IOTraceContainer::remove (IOTrace& trace)
{
    for (auto it = list->begin(); it != list->end();)
    {
        if (*it == trace)
        {
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

vector<OutputTrace> IOTraceContainer::getOutputTraces() const
{
    vector<OutputTrace> result;
    for (IOTrace& iOTrace : *list)
    {
        result.push_back(iOTrace.getOutputTrace());
    }
    return result;
}

void IOTraceContainer::addUnique(std::vector<IOTraceContainer>& container, const IOTraceContainer& elem)
{
    for (const IOTraceContainer& cont : container)
    {
        if (cont == elem)
        {
            return;
        }
    }
    container.push_back(elem);
}

void IOTraceContainer::remove(std::vector<IOTraceContainer>& container, const IOTraceContainer& elem)
{
    for (auto it = container.begin(); it != container.end();)
    {
        if (*it == elem)
        {
            it = container.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

std::ostream & operator<<(std::ostream & out, const IOTraceContainer & iot)
{
    out << "{\n";

    bool isFirst = true;
    for (IOTrace& trace : *iot.list)
    {
        if (!isFirst)
        {
            out << ",\n";
        }

        out << "  " << trace;
        isFirst = false;
    }
    out << "\n}";
    return out;
}

bool operator==(IOTraceContainer const & cont1, IOTraceContainer const & cont2)
{
    if (cont1.size() != cont2.size())
    {
        return false;
    }
    for (size_t i = 0; i < cont1.size(); ++i)
    {
        shared_ptr<vector<IOTrace>> list = cont2.getList();
        if (find(list->begin(), list->end(), cont1.getList()->at(i)) == list->end())
        {
            return false;
        }
    }
    return true;
}
