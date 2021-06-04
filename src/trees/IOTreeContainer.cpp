#include "IOTreeContainer.h"
#include "IOListContainer.h"
#include "utils/Logger.hpp"
#include "trees/InputOutputTree.h"

using namespace std;

IOTreeContainer::IOTreeContainer(const std::shared_ptr<FsmPresentationLayer>& presentationLayer):
    list(make_shared<vector<shared_ptr<InputOutputTree>>>()), presentationLayer(presentationLayer)
{
}

IOTreeContainer::IOTreeContainer(const std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>>& list, const std::shared_ptr<FsmPresentationLayer>& presentationLayer):
    list(list), presentationLayer(presentationLayer)
{

}

shared_ptr<vector<shared_ptr<InputOutputTree>>> IOTreeContainer::getList() const
{
    return list;
}


IOListContainer IOTreeContainer::toIOList() const
{
    LOG("VERBOSE_8") << "toIOList()" << std::endl;
    IOListContainer result = IOListContainer(presentationLayer);
    bool containsEmpty = false;
    for (shared_ptr<InputOutputTree> tree : *list)
    {
        LOG("VERBOSE_8") << "tree: " << tree->str() << std::endl;
        if (tree->isEmpty())
        {
            LOG("VERBOSE_8") << "  Tree is empty." << std::endl;
            if (!containsEmpty)
            {
                LOG("VERBOSE_8") << "  Adding empty tree." << std::endl;
                result.addUnique((Trace({})));
            }
            continue;
        }
        IOListContainer container = tree->getInputLists();
        LOG("VERBOSE_8") << "  Tree as input list: " << container << std::endl;
        shared_ptr<vector<vector<int>>> set = container.getIOLists();
        LOG("VERBOSE_8") << "  Tree as IO set: " << std::endl;
        for (vector<int> e : *set)
        {
            stringstream ss;
            ss << "[";
            for (int i : e)
            {
                ss << i << ",";
            }
            ss << "]";
            LOG("VERBOSE_8") << ss.str() << std::endl;
        }
        for (vector<int> trace : *set)
        {
            stringstream ss;
            ss << "Adding trace: [";
            for (int i : trace)
            {
                ss << i << ",";
            }
            ss << "]";
            LOG("VERBOSE_8") << ss.str() << std::endl;
            result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
            LOG("VERBOSE_8") << "  result: " << result << std::endl;
        }
    }
    return result;
}

void IOTreeContainer::addUnique(std::shared_ptr<InputOutputTree> tree)
{
    for(auto inLst : *list)
    {
        if (*tree == *inLst){
            return;
        }
    }
    list->push_back(tree);
}

void IOTreeContainer::removePrefixes(std::shared_ptr<InputOutputTree> tree)
{
    auto it = list->begin();
    while (it != list->end())
    {
        if (tree->contains(**it))
        {
            it = list->erase(it);
            continue;
        }
        ++it;
    }
}

void IOTreeContainer::addUniqueRemovePrefixes(std::shared_ptr<InputOutputTree> tree)
{
    removePrefixes(tree);
    for (shared_ptr<InputOutputTree> t : *list)
    {
        if (t->contains(*tree))
        {
            return;
        }
    }
    list->push_back(tree);
}

size_t IOTreeContainer::size() const
{
    return list->size();
}

ostream & operator<<(ostream & out, const IOTreeContainer & ot)
{
    out << "{";
    bool first = true;
    for (shared_ptr<InputOutputTree> tree : *ot.getList())
    {
        if (!first)
        {
            out << ", ";
        }
        out << "[";
        out << *tree;
        out << "]";
        first = false;
    }
    out << "}";
    return out;
}
