#include <ostream>

#include "IOTraceContainer.h"
#include "trees/OutputTree.h"
#include "fsm/IOTrace.h"

using namespace std;

IOTraceContainer::IOTraceContainer()
{

}

IOTraceContainer::IOTraceContainer(const shared_ptr<IOTraceCont>& list):
    list(*list)
{

}

IOTraceContainer::IOTraceContainer(const shared_ptr<IOTrace>& trace)
{
    list.insert(trace);
}

void IOTraceContainer::add(const shared_ptr<const IOTrace>& trc)
{
    list.insert(trc);
}

void IOTraceContainer::addRemovePrefixes(const shared_ptr<const IOTrace>& trc)
{
    removePrefixes(trc);
    add(trc);
}

void IOTraceContainer::addRemovePrefixes(const IOTraceContainer& cont)
{
    for (auto it = cont.cbegin(); it != cont.cend(); ++it)
    {
        addRemovePrefixes(*it);
    }
}

void IOTraceContainer::removePrefixes(const shared_ptr<const IOTrace>& trc)
{
    auto it = list.begin();
    while (it != list.end())
    {
        if ((*it)->isPrefixOf(*trc))
        {
            it = list.erase(it);
        }
        else
        {
            ++it;
        }
    }

}

void IOTraceContainer::add(const IOTraceContainer& container)
{
    list.insert(container.cbegin(), container.cend());
}

void IOTraceContainer::add(OutputTree& tree)
{
    std::vector<shared_ptr<IOTrace>> iOTraces;
    tree.toIOTrace(iOTraces);
    list.insert(iOTraces.begin(), iOTraces.end());
}

bool IOTraceContainer::contains(const shared_ptr<const IOTrace>& trace) const
{
    return list.find(trace) != list.end();
}

IOTraceCont::const_iterator IOTraceContainer::get(const InputTrace& inputTrace) const
{
    for (auto it = list.cbegin(); it != list.cend(); ++it)
    {
        if ((*it)->getInputTrace() == inputTrace)
        {
            return it;
        }
    }
    return list.cend();
}

void IOTraceContainer::concatenate(IOTrace& trace)
{
    IOTraceCont newList;
    for (const shared_ptr<const IOTrace>& t : list)
    {
        shared_ptr<const IOTrace> newTrace = make_shared<const IOTrace>(*t, trace);
        newList.insert(newTrace);
    }
    list = newList;
}

void IOTraceContainer::concatenate(IOTraceContainer& container)
{
    IOTraceCont newList;
    for (const shared_ptr<const IOTrace>& thisTrace : list)
    {
        for (auto it = container.cbegin(); it != container.cend(); ++it)
        {
            shared_ptr<const IOTrace> newTrace = make_shared<const IOTrace>(*thisTrace, **it);
            newList.insert(newTrace);
        }
    }
    list = newList;
}

void IOTraceContainer::concatenateToFront(const shared_ptr<InputTrace>& inputTrace, const shared_ptr<OutputTrace>& outputTrace)
{
    shared_ptr<IOTrace> newIOTrace = make_shared<IOTrace>(*inputTrace, *outputTrace);
    concatenateToFront(newIOTrace);
}

void IOTraceContainer::concatenateToFront(const shared_ptr<const IOTrace>& iOTrace)
{
    IOTraceCont newList;
    for (const shared_ptr<const IOTrace>& iOT : list)
    {
        shared_ptr<const IOTrace> newTrace = make_shared<const IOTrace>(*iOT, *iOTrace, true);
        newList.insert(newTrace);
    }
    list = newList;
}

void IOTraceContainer::clear()
{
    list = IOTraceCont();
}

bool IOTraceContainer::remove(const shared_ptr<const IOTrace>& trace)
{
    return list.erase(trace) > 0;
}

vector<OutputTrace> IOTraceContainer::getOutputTraces() const
{
    vector<OutputTrace> result;
    for (const shared_ptr<const IOTrace>& iOTrace : list)
    {
        result.push_back(iOTrace->getOutputTrace());
    }
    return result;
}

shared_ptr<IOTrace> IOTraceContainer::getLongestTrace() const
{
    shared_ptr<const IOTrace> longest;
    for (const shared_ptr<const IOTrace>& trace : list)
    {
        if (!longest || longest->size() < trace -> size())
        {
            longest = trace;
        }
    }
    if (longest)
    {
        return make_shared<IOTrace>(*longest);
    }
    else
    {
        return nullptr;
    }

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

void IOTraceContainer::remove(std::unordered_set<IOTraceContainer>& container, const IOTraceContainer& elem)
{
    container.erase(elem);
}

std::ostream & operator<<(std::ostream & out, const IOTraceContainer & iot)
{
    out << "{\n";

    bool isFirst = true;
    for (const shared_ptr<const IOTrace>& trace : iot.list)
    {
        if (!isFirst)
        {
            out << ",\n";
        }

        out << "  " << *trace;
        isFirst = false;
    }
    out << "\n}";
    return out;
}

bool operator==(IOTraceContainer const & cont1, IOTraceContainer const & cont2)
{
    return cont1.list == cont2.list;
}

bool std::operator==(IOTraceCont const& x, IOTraceCont const& y)
{
    if (x.size() != y.size())
    {
        return false;
    }
    for (auto e : x)
    {
        if (y.find(e) == y.end())
        {
            return false;
        }
    }
    return true;
}

bool operator!=(IOTraceContainer const & cont1, IOTraceContainer const & cont2)
{
    return !(cont1 == cont2);
}
