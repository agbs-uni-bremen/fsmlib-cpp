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
