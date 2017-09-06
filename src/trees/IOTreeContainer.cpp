#include "IOTreeContainer.h"

using namespace std;

IOTreeContainer::IOTreeContainer(const std::shared_ptr<FsmPresentationLayer> presentationLayer):
    list(make_shared<vector<shared_ptr<InputOutputTree>>>()), presentationLayer(presentationLayer)
{
}

IOTreeContainer::IOTreeContainer(const std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> list, const std::shared_ptr<FsmPresentationLayer> presentationLayer):
    list(list), presentationLayer(presentationLayer)
{

}

shared_ptr<vector<shared_ptr<InputOutputTree>>> IOTreeContainer::getList() const
{
    return list;
}


IOListContainer IOTreeContainer::toIOList() const
{
    IOListContainer result = IOListContainer(presentationLayer);
    for (shared_ptr<InputOutputTree> tree : *list)
    {
        if (tree->isEmpty())
        {
            result.addUnique((Trace({})));
            continue;
        }
        IOListContainer container = tree->getInputLists();
        auto set = container.getIOLists();
        for (auto trace : *set)
        {
            result.addUniqueRemovePrefixes(Trace(trace, presentationLayer));
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
        cout << "[";
        out << *tree;
        cout << "]";
        first = false;
    }
    out << "}\n";
    return out;
}
