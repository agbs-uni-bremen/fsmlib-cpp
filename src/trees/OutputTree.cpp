/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/OutputTree.h"

void OutputTree::printChildrenOutput(std::ostream & out, const std::shared_ptr<TreeNode> top, const std::shared_ptr<int> idNode, const int idInput) const
{
	int idNodeBase = *idNode;
	for (std::shared_ptr<TreeEdge> edge : *top->getChildren())
	{
		out << idNodeBase << " -> " << ++ *idNode << "[label = \"" << inputTrace.get().at(idInput) << "/" << edge->getIO() << "\" ];" << std::endl;
		printChildrenOutput(out, edge->getTarget(), idNode, idInput + 1);
	}
}

OutputTree::OutputTree(const std::shared_ptr<TreeNode> root, const InputTrace & inputTrace, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Tree(root, presentationLayer), inputTrace(inputTrace)
{

}

bool OutputTree::contains(const OutputTree & ot) const
{
	/* If the associated input traces differ,
	the output trees can never be equal*/
	if (!(inputTrace == ot.inputTrace))
	{
		return false;
	}

	return getRoot()->superTreeOf(ot.getRoot());
}

void OutputTree::toDot(std::ostream & out) const
{
	out << "digraph OutputTree {" << std::endl;
	out << "\trankdir=TB;" << std::endl;//Top -> Bottom, to create a vertical graph
	out << "\tnode [shape = circle];" << std::endl;
	std::shared_ptr<int> id = std::make_shared<int>(0);
	printChildrenOutput(out, root, id, 0);
	out << "}";
}

void OutputTree::store(std::ofstream & file)
{
	std::vector<std::vector<int>> lli = *getIOLists().getIOLists();
	for (std::vector<int> lst : lli)
	{
		for (unsigned int i = 0; i < lst.size(); ++ i)
		{
			if (i != 0)
			{
				file << ".";
			}

			file << "(" << inputTrace.get().at(i) << "," << lst.at(i) <<")";
		}
	}
}

std::ostream & operator<<(std::ostream & out, OutputTree & ot)
{
	std::vector<std::vector<int>> lli = *ot.getIOLists().getIOLists();
	for (std::vector<int> lst : lli)
	{
		bool isFirst = true;

		for (unsigned int i = 0; i < lst.size(); ++ i)
		{
			if (isFirst)
			{
				isFirst = false;
			}
			else
			{
				out << ".";
			}

			out << "(" << ot.presentationLayer->getInId(ot.inputTrace.get().at(i)) << "/" << ot.presentationLayer->getOutId(lst.at(i)) << ")";
		}
		out << std::endl;
	}
	return out;
}

bool operator==(OutputTree const & outputTree1, OutputTree const & outputTree2)
{
	/*If the associated input traces differ,
	the output trees can never be equal*/
	if (!(outputTree1.inputTrace == outputTree2.inputTrace))
	{
		return false;
	}

	/*Since outputTree1 and outputTree2 are two output trees over the same
	input trace, the trees have the same height. We only have
	to check whether each corresponding node has the same number
	of children, and whether corresponding edges carry the same labels
	(output values).*/
	return *outputTree1.getRoot() == *outputTree2.getRoot();
}
