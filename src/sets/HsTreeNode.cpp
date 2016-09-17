/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "sets/HsTreeNode.h"

std::unordered_set<int> HsTreeNode::hSmallest;

int HsTreeNode::maxNodeNum = 0;

HsTreeNode::HsTreeNode(const std::unordered_set<int>& x, const std::vector<std::unordered_set<int>>& s)
	: x(x), s(s), nodeNum(maxNodeNum ++)
{

}

size_t HsTreeNode::size() const
{
	return x.size();
}

bool HsTreeNode::isHittingSet() const
{
	for (std::unordered_set<int> z : s)
	{
		std::unordered_set<int> z1;

		/*Not present in java, it is because the set in C++ doesn't have any retainAll method*/
		for (int i : z)
		{
			if (std::find(x.cbegin(), x.cend(), i) != x.cend())
			{
				z1.insert(i);
			}
		}

		if (z1.empty())
		{
			return false;
		}
	}

	/*This node is a valid hitting set
	Check whether x is smaller than hSmallest,
	or whether hSmallest is undefined yet.*/
	if (x.size() < hSmallest.size())
	{
		hSmallest = std::unordered_set<int>(x);
	}
	return true;
}

void HsTreeNode::add(const HsTreeNode & node)
{
	children.push_back(node);
}

void HsTreeNode::expandNode()
{
	for (int a : x)
	{
		std::unordered_set<int> xNew(x);
		xNew.erase(a);

		if (!xNew.empty())
		{
			HsTreeNode newNode = HsTreeNode(xNew, s);
			if (newNode.isHittingSet())
			{
				add(newNode);
				newNode.expandNode();
			}
		}
	}
}

std::string HsTreeNode::toDot()
{
	std::stringstream ss;
	ss << std::endl << "\t" << nodeNum << "[label=\"[";

	/*Not present in java, it is because the vector in C++ doesn't have any toString method*/
	for (auto it = x.cbegin(); it != x.cend(); ++ it)
	{
		if (it != x.cbegin())
		{
			ss << ", ";
		}
		ss << *it;
	}

	ss << "]\"];";

	for (HsTreeNode & n : children)
	{
		ss << std::endl << "\t" << nodeNum << " -> " << n.nodeNum << ";";
		ss << n.toDot();
	}
	std::string str;
	ss >> str;
	return str;
}
