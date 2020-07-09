#include <utility>
#include <algorithm>

#include "RDistinguishability.h"
#include "trees/InputOutputTree.h"
#include "trees/AdaptiveTreeNode.h"
#include "fsm/FsmNode.h"
#include "utils/Logger.hpp"

using namespace std;
RDistinguishability::RDistinguishability(const shared_ptr<FsmPresentationLayer>& presentationLayer) : presentationLayer(presentationLayer)
{
    hBeenCalculated = false;
}

vector<int>::iterator RDistinguishability::removeNotRDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    vector<int>& notDist = notRDistinguishableWith.at(i);
    for (vector<int>::iterator it = notDist.begin(); it != notDist.end(); ++it)
    {
        if (*it == node->getId()) {
            return notDist.erase(it);
        }
    }
    return notDist.end();
}

void RDistinguishability::initRDistinguishable(size_t i)
{
    auto it = rDistinguishableWith.find(i);
    if (it == rDistinguishableWith.end())
    {
        rDistinguishableWith.insert(pair<size_t, std::vector<int>>(i, {}));
    }
    else
    {
        LOG("WARNING") << "Overwriting r-distinguishability." << std::endl;
        rDistinguishableWith[i] = {};
    }
}

void RDistinguishability::addRDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    auto it = rDistinguishableWith.find(i);
    if (it == rDistinguishableWith.end())
    {
        rDistinguishableWith.insert(pair<size_t, std::vector<int>>(i, {node->getId()}));
    }
    else
    {
        it->second.push_back(node->getId());
    }
}

void RDistinguishability::addNotRDistinguishable(size_t i, std::shared_ptr<FsmNode> node)
{
    auto it = notRDistinguishableWith.find(i);
    if (it == notRDistinguishableWith.end())
    {
        notRDistinguishableWith.insert(pair<size_t, std::vector<int>>(i, {node->getId()}));
    }
    else
    {
        it->second.push_back(node->getId());
    }
}

void RDistinguishability::addNotRDistinguishable(size_t i)
{
    notRDistinguishableWith[i] = {};
}

void RDistinguishability::addAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode, std::shared_ptr<InputOutputTree> tree)
{
    adaptiveIOSequences.insert(pair<int, shared_ptr<InputOutputTree>>(otherNode->getId(), tree));
}

vector<int> RDistinguishability::getRDistinguishableWith(size_t i)
{
    return rDistinguishableWith.at(i);
}

vector<int> RDistinguishability::getRDistinguishableWith()
{
    if (rDistinguishableWith.rbegin() != rDistinguishableWith.rend())
    {
        return rDistinguishableWith.rbegin()->second;
    }
    else
    {
        return vector<int>();
    }

}

vector<int> RDistinguishability::getNotRDistinguishableWith(size_t i)
{
    return notRDistinguishableWith.at(i);
}

bool RDistinguishability::isNotRDistinguishable()
{
    return notRDistinguishableWith.rbegin()->second.size() != 0;
}

bool RDistinguishability::isRDistinguishableWith(size_t i, std::shared_ptr<FsmNode> node)
{
    for (; i > 0; --i)
    {
        try {
            auto dist = rDistinguishableWith.at(i);
            if ( std::find(dist.begin(), dist.end(), node->getId()) != dist.end() )
            {
                return true;
            }
        } catch (out_of_range &e) {
            continue;
        }

    }
    return false;
}

bool RDistinguishability::isRDistinguishableWith(std::shared_ptr<FsmNode> node)
{
    for (auto it = rDistinguishableWith.rbegin(); it != rDistinguishableWith.rend(); ++it)
    {
        vector<int> dist = it->second;
        if ( std::find(dist.begin(), dist.end(), node->getId()) != dist.end() )
        {
            return true;
        }
    }
    return false;
}

bool RDistinguishability::isRDistinguishableWith(vector<shared_ptr<FsmNode>> nodes)
{
    for (shared_ptr<FsmNode> node : nodes)
    {
        if (!isRDistinguishableWith(node))
        {
            return false;
        }
    }
    return true;
}

shared_ptr<InputOutputTree> RDistinguishability::getAdaptiveIOSequence(shared_ptr<FsmNode> otherNode)
{
    auto it = adaptiveIOSequences.find(otherNode->getId());
    if (it == adaptiveIOSequences.end())
    {
        return make_shared<InputOutputTree>(make_shared<AdaptiveTreeNode>(), presentationLayer);
    }
    return it->second->Clone();
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

bool RDistinguishability::hasBeenCalculated() const
{
    return hBeenCalculated;
}

void RDistinguishability::hasBeenCalculated(bool hasBeen)
{
    hBeenCalculated = hasBeen;
}
