#include <utility>

#include "RDistinguishability.h"
#include "trees/OutputTree.h"

using namespace std;
RDistinguishability::RDistinguishability()
{

}

vector<shared_ptr<FsmNode>>::iterator RDistinguishability::removeNotDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    vector<shared_ptr<FsmNode>>& notDist = notRDistinguishableWith.at(i);
    for (vector<shared_ptr<FsmNode>>::iterator it = notDist.begin(); it != notDist.end(); ++it)
    {
        if (*it == node) {
            return notDist.erase(it);
        }
    }
    return notDist.end();
}

void RDistinguishability::addDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    auto it = rDistinguishableWith.find(i);
    if (it == rDistinguishableWith.end())
    {
        rDistinguishableWith.insert(pair<size_t, std::vector<std::shared_ptr<FsmNode>>>(i, {node}));
    }
    else
    {
        it->second.push_back(node);
    }
}

void RDistinguishability::addNotDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    auto it = notRDistinguishableWith.find(i);
    if (it == notRDistinguishableWith.end())
    {
        notRDistinguishableWith.insert(pair<size_t, std::vector<std::shared_ptr<FsmNode>>>(i, {node}));
    }
    else
    {
        it->second.push_back(node);
    }
}

void RDistinguishability::addNotDistinguishable(size_t i)
{
    notRDistinguishableWith.insert(pair<size_t, std::vector<std::shared_ptr<FsmNode>>>(i, {}));
}

void RDistinguishability::addAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode, std::shared_ptr<OutputTree> tree)
{
    auto it = adaptiveIOSequences.find(otherNode);
    if (it == adaptiveIOSequences.end())
    {
        adaptiveIOSequences.insert(pair<shared_ptr<FsmNode>, std::shared_ptr<OutputTree>>(otherNode, {tree}));
    }
    else
    {

    }
}

vector<shared_ptr<FsmNode>> RDistinguishability::getRDistinguishableWith(size_t i)
{
    return rDistinguishableWith.at(i);
}

vector<shared_ptr<FsmNode>> RDistinguishability::getNotRDistinguishableWith(size_t i)
{
    return notRDistinguishableWith.at(i);
}

OutputTree RDistinguishability::getAdaptiveIOSequence(shared_ptr<FsmNode> otherNode)
{
    return *adaptiveIOSequences.at(otherNode);
}
