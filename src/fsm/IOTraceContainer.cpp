#include "IOTraceContainer.h"

using namespace std;

IOTraceContainer::IOTraceContainer(const shared_ptr<FsmPresentationLayer> presentationLayer):
    presentationLayer(presentationLayer)
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
