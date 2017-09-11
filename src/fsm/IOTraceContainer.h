#ifndef IOTRACECONTAINER_H
#define IOTRACECONTAINER_H

#include "fsm/IOTrace.h"

class IOTraceContainer
{
private:
    std::shared_ptr<std::vector<IOTrace>> list;
    const std::shared_ptr<FsmPresentationLayer> presentationLayer;
public:
    IOTraceContainer(const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    IOTraceContainer(std::shared_ptr<std::vector<IOTrace>>& list, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    std::shared_ptr<std::vector<IOTrace>> getList() const;
    void addUnique(IOTrace& trc);
};

#endif // IOTRACECONTAINER_H
