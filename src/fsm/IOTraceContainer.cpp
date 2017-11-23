#include "IOTraceContainer.h"
#include "logging/easylogging++.h"

using namespace std;

IOTraceContainer::IOTraceContainer()
{

}

IOTraceContainer::IOTraceContainer(unordered_set<shared_ptr<IOTrace>>& list):
    list(list)
{

}

IOTraceContainer::IOTraceContainer(shared_ptr<IOTrace> trace)
{
    list.insert(trace);
}

void IOTraceContainer::add(const shared_ptr<IOTrace>& trc)
{
    list.insert(trc);
}

void IOTraceContainer::addRemovePrefixes(const shared_ptr<IOTrace>& trc)
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

void IOTraceContainer::removePrefixes(const shared_ptr<IOTrace>& trc)
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

void IOTraceContainer::add(IOTraceContainer& container)
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(1));
    for (auto it = container.cbegin(); it != container.cend(); ++it)
    {
        add(*it);
    }
}

void IOTraceContainer::add(OutputTree& tree)
{
    std::vector<shared_ptr<IOTrace>> iOTraces;
    tree.toIOTrace(iOTraces);
    for (const shared_ptr<IOTrace>& trace : iOTraces)
    {
        add(trace);
    }
}

bool IOTraceContainer::contains(const shared_ptr<IOTrace>& trace) const
{
    for (const shared_ptr<IOTrace>& t : list)
    {
        if (*t == *trace)
        {
            return true;
        }
    }
    return false;
}

unordered_set<shared_ptr<IOTrace>>::const_iterator IOTraceContainer::get(const InputTrace& inputTrace) const
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
    for (const shared_ptr<IOTrace>& t : list)
    {
        t->append(trace);
    }
}

void IOTraceContainer::concatenate(IOTraceContainer& container)
{
    for (const shared_ptr<IOTrace>& thisTrace : list)
    {
        for (auto it = container.cbegin(); it != container.cend(); ++it)
        {
            thisTrace->append(**it);
        }
    }
}

void IOTraceContainer::concatenateToFront(const shared_ptr<InputTrace>& inputTrace, const shared_ptr<OutputTrace>& outputTrace)
{
    TIMED_FUNC_IF(timerObj, VLOG_IS_ON(1));
    shared_ptr<IOTrace> newIOTrace = make_shared<IOTrace>(*inputTrace, *outputTrace);
    concatenateToFront(newIOTrace);
}

void IOTraceContainer::concatenateToFront(const shared_ptr<IOTrace>& iOTrace)
{
    for (const shared_ptr<IOTrace>& iOT : list)
    {
        iOT->prepend(*iOTrace);
    }
}

void IOTraceContainer::clear()
{
    list = std::unordered_set<shared_ptr<IOTrace>>();
}

void IOTraceContainer::remove (IOTrace& trace)
{
    for (auto it = list.cbegin(); it != list.cend();)
    {
        if (**it == trace)
        {
            it = list.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

void IOTraceContainer::remove (IOTraceContainer& container)
{
    for (auto it = container.cbegin(); it != container.cend(); ++it)
    {
        remove(**it);
    }
}

vector<OutputTrace> IOTraceContainer::getOutputTraces() const
{
    vector<OutputTrace> result;
    for (const shared_ptr<IOTrace>& iOTrace : list)
    {
        result.push_back(iOTrace->getOutputTrace());
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
    for (const shared_ptr<IOTrace>& trace : iot.list)
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

bool operator!=(IOTraceContainer const & cont1, IOTraceContainer const & cont2)
{
    return !(cont1 == cont2);
}
