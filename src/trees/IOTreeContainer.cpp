#include "IOTreeContainer.h"
#include "logging/easylogging++.h"

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
    VLOG(8) << "toIOList()";
    IOListContainer result = IOListContainer(presentationLayer);
    bool containsEmpty = false;
    for (shared_ptr<InputOutputTree> tree : *list)
    {
        VLOG(8) << "tree: " << tree->str();
        if (tree->isEmpty())
        {
            VLOG(8) << "  Tree is empty.";
            if (!containsEmpty)
            {
                VLOG(8) << "  Adding empty tree.";
                result.addUnique((Trace({})));
            }
            continue;
        }
        IOListContainer container = tree->getInputLists();
        VLOG(8) << "  Tree as input list: " << container;
        shared_ptr<vector<vector<int>>> set = container.getIOLists();
        VLOG(8) << "  Tree as IO set: ";
        for (vector<int> e : *set)
        {
            stringstream ss;
            ss << "[";
            for (int i : e)
            {
                ss << i << ",";
            }
            ss << "]";
            VLOG(8) << ss.str();
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
            VLOG(8) << ss.str();
            result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
            VLOG(8) << "  result: " << result;
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
