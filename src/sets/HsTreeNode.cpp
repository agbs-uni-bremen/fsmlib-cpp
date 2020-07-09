/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "sets/HsTreeNode.h"

#include <sstream>

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
        
        bool emptyIntersection = true;

        // Check whether x has a non-empty intersection with z
		for (int i : z)
		{
			if ( x.find(i) != x.cend())
			{
                emptyIntersection = false;
                break;
			}
		}

        if ( emptyIntersection ) return false;
		 
	}

    // If x is smaller than hSmallest, x becomes the new
    // smallest hitting set candidate.
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
