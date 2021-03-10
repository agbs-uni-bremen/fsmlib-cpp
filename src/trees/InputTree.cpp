/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/InputTree.h"
#include "fsm/IOTrace.h"
#include "trees/TreeNode.h"
#include "trees/TreeEdge.h"
#include "trees/IOListContainer.h"
#include "interface/FsmPresentationLayer.h"

#include <algorithm>
#include <fstream>

using namespace std;

void InputTree::printChildrenInput(ostream& out, const shared_ptr<TreeNode>& top, const shared_ptr<int>& idNode, const int idInput) const
{
	int idNodeBase = *idNode;
	for (shared_ptr<TreeEdge> edge : *top->getChildren())
	{
		out << idNodeBase << " -> " << ++ *idNode << "[label = \"" << edge->getIO() << "\" ];" << endl;
		printChildrenInput(out, edge->getTarget(), idNode, idInput + 1);
	}
}

InputTree::InputTree(const shared_ptr<TreeNode>& root,
                       const shared_ptr<FsmPresentationLayer>& presentationLayer)
	: Tree(root, presentationLayer)
{

}

InputTree::InputTree(const InputTree* other):
    Tree(other)
{

}

InputTree::InputTree(const shared_ptr<FsmPresentationLayer>& presentationLayer)
	: Tree(make_shared<TreeNode>(), presentationLayer)
{

}


bool InputTree::contains(const InputTree& ot) const
{
    //Get all input traces
    vector<InputTrace> myInputs = getInputTraces();
    vector<InputTrace> otherInputs = ot.getInputTraces();
    
    //Sort both sequences of input traces to obtain sorted sets of input
    //traces ( O(nlog(n)) )
    sort(myInputs.begin(), myInputs.end());
    sort(otherInputs.begin(), otherInputs.end());
    vector<InputTrace>::const_iterator myLast = unique(myInputs.begin(), myInputs.end());
    vector<InputTrace>::const_iterator otherLast = unique(otherInputs.begin(), otherInputs.end());
    //Return whether the set of input traces of this tree contains the set of
    //input traces of the other ( O(n) )
    return includes(myInputs.cbegin(), myLast, otherInputs.cbegin(), otherLast);
}

vector<InputTrace> InputTree::getInputTraces() const {
    //Get all traces
    auto lli = getIOLists().getIOLists();
    vector<InputTrace> traces;
    //Reserve enough space in the result vector
    traces.reserve(lli->size());
    //Move each trace into the result vector and construct an InputTrace for
    //each
    for(auto &trace : *lli) {
        traces.emplace_back(std::move(trace), presentationLayer);
    }
    return traces;
}



void InputTree::toDot(ostream& out) const
{
	out << "digraph InputTree {" << endl;
	out << "\trankdir=TB;" << endl;//Top -> Bottom, to create a vertical graph
	out << "\tnode [shape = circle];" << endl;
	shared_ptr<int> id = make_shared<int>(0);
	printChildrenInput(out, root, id, 0);
	out << "}";
}

void InputTree::store(ofstream& file)
{
	vector<vector<int>> lli = *getIOLists().getIOLists();
	for (vector<int> lst : lli)
	{
		for (unsigned int i = 0; i < lst.size(); ++ i)
		{
			if (i != 0)
			{
				file << ".";
			}

			file << lst.at(i);
		}
	}
}

InputTree* InputTree::_clone() const
{
    return new InputTree( this );
}

std::shared_ptr<InputTree> InputTree::Clone() const
{
    return std::shared_ptr<InputTree>(_clone());
}

ostream& operator<<(ostream& out, InputTree& ot)
{
	vector<vector<int>> lli = *ot.getIOLists().getIOLists();
	for (vector<int> lst : lli)
	{
		for (unsigned int i = 0; i < lst.size(); ++ i)
		{
            
            if ( i > 0 ) out << ".";

			out << ot.presentationLayer->getInId(lst.at(i));
		}
		out << endl;
	}
	return out;
}

bool operator==(InputTree const &inputTree1, InputTree const &inputTree2)
{
    
    return ( inputTree1.contains(inputTree2) and inputTree2.contains(inputTree1) );    
}

bool operator!=(InputTree const &inputTree1, InputTree const &inputTree2)
{
    return not (inputTree1 == inputTree2);
}

shared_ptr<InputTree> InputTree::getSubtreeForInput(int input){
    shared_ptr<TreeNode> targetRoot = root->after(input);
    shared_ptr<TreeNode> cpyNode = targetRoot->clone();
    return make_shared<InputTree>(cpyNode,presentationLayer);
}

std::vector<int> InputTree::getInputsAtRoot() const{
    std::vector<int> inputs;
    for (auto edge : *(root->getChildren())) {
        inputs.push_back(edge->getIO());
    }
    return inputs;
}