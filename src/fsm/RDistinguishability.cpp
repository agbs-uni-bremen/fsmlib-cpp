#include <utility>

#include "RDistinguishability.h"
#include "trees/InputOutputTree.h"
#include "trees/AdaptiveTreeNode.h"

using namespace std;
RDistinguishability::RDistinguishability(shared_ptr<FsmPresentationLayer> presentationLayer) : presentationLayer(presentationLayer)
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

void RDistinguishability::addAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode, std::shared_ptr<InputOutputTree> tree)
{
    adaptiveIOSequences.insert(pair<shared_ptr<FsmNode>, std::shared_ptr<InputOutputTree>>(otherNode, tree));
}

vector<shared_ptr<FsmNode>> RDistinguishability::getRDistinguishableWith(size_t i)
{
    return rDistinguishableWith.at(i);
}

vector<shared_ptr<FsmNode>> RDistinguishability::getNotRDistinguishableWith(size_t i)
{
    return notRDistinguishableWith.at(i);
}

bool RDistinguishability::isRDistinguishableWith(size_t i, std::shared_ptr<FsmNode> node)
{
    for (size_t j = i; i > 0; --i)
    {
        auto dist = rDistinguishableWith.at(j);
        if ( std::find(dist.begin(), dist.end(), node) != dist.end() )
        {
            return true;
        }
    }
    return false;
}

bool RDistinguishability::isNotRDistinguishableWith(size_t i, std::shared_ptr<FsmNode> node)
{
    for (size_t j = i; i > 0; --i)
    {
        auto notDist = notRDistinguishableWith.at(j);
        if ( std::find(notDist.begin(), notDist.end(), node) != notDist.end() )
        {
            return true;
        }
    }
    return false;
}

InputOutputTree RDistinguishability::getAdaptiveIOSequence(shared_ptr<FsmNode> otherNode)
{
    auto it = adaptiveIOSequences.find(otherNode);
    if (it == adaptiveIOSequences.end())
    {
        return InputOutputTree(make_shared<AdaptiveTreeNode>(), presentationLayer);
    }
    return *it->second;
}

void RDistinguishability::inheritDistinguishability(size_t i)
{
    auto it = rDistinguishableWith.find(i - 1);
    if (it != rDistinguishableWith.end())
    {
        rDistinguishableWith.insert(make_pair(i, it->second));
    }
    it = notRDistinguishableWith.find(i - 1);
    if (it != notRDistinguishableWith.end())
    {
        notRDistinguishableWith.insert(make_pair(i, it->second));
    }
}
