#ifndef RDISTINGUISHABILITY_H
#define RDISTINGUISHABILITY_H

#include <vector>
#include <memory>
#include <map>

class FsmNode;
class InputOutputTree;
class FsmPresentationLayer;


class RDistinguishability
{
private:
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
    bool hBeenCalculated;
public:
    RDistinguishability(const std::shared_ptr<FsmPresentationLayer> presentationLayer);

protected:
    std::map<size_t, std::vector<std::shared_ptr<FsmNode>>> rDistinguishableWith;
    std::map<size_t, std::vector<std::shared_ptr<FsmNode>>> notRDistinguishableWith;
    std::map<std::shared_ptr<FsmNode>, std::shared_ptr<InputOutputTree>> adaptiveIOSequences;
public:
    std::vector<std::shared_ptr<FsmNode>>::iterator removeNotDistinguishable(size_t i, std::shared_ptr<FsmNode> node);
    void addDistinguishable(size_t i, std::shared_ptr<FsmNode> node);
    void addNotDistinguishable(size_t i, std::shared_ptr<FsmNode> node);
    void addNotDistinguishable(size_t i);
    void addAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode, std::shared_ptr<InputOutputTree> tree);
    std::vector<std::shared_ptr<FsmNode>> getRDistinguishableWith(size_t i);
    std::vector<std::shared_ptr<FsmNode>> getNotRDistinguishableWith(size_t i);
    bool isRDistinguishableWith(size_t i, std::shared_ptr<FsmNode> node);
    bool isRDistinguishableWith(std::shared_ptr<FsmNode> node);
    bool isRDistinguishableWith(std::vector<std::shared_ptr<FsmNode>> nodes);
    bool isNotRDistinguishable();
    std::shared_ptr<InputOutputTree> getAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode);
    void inheritDistinguishability(size_t i);
    bool hasBeenCalculated() const;
    void hasBeenCalculated(bool hasBeen);

};

#endif // RDISTINGUISHABILITY_H
