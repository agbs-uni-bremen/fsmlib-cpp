/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/OutputTrace.h"
#include "fsm/OFSMTable.h"
#include "fsm/DFSMTableRow.h"
#include "fsm/PkTable.h"
#include "trees/TreeEdge.h"
#include "trees/TreeNode.h"
#include "trees/OutputTree.h"
#include "trees/Tree.h"
#include "trees/TreeNode.h"
#include "trees/IOListContainer.h"
#include "interface/FsmPresentationLayer.h"

FsmNode::FsmNode(const int id, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: id(id), visited(false), color(white), presentationLayer(presentationLayer), derivedFromPair(nullptr)
{

}

FsmNode::FsmNode(const int id, const std::string & name, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: FsmNode(id, presentationLayer)
{
	this->name = name;
}

void FsmNode::addTransition(const FsmTransition & transition)
{
	transitions.push_back(transition);
}

std::vector<FsmTransition> FsmNode::getTransitions() const
{
	return transitions;
}

int FsmNode::getId() const
{
	return id;
}

std::string FsmNode::getName() const
{
	return presentationLayer->getStateId(id, name);
}

bool FsmNode::hasBeenVisited() const
{
	return visited;
}

void FsmNode::setVisited()
{
	visited = true;
}

void FsmNode::setPair(const std::shared_ptr<FsmNode> l, const std::shared_ptr<FsmNode> r)
{
	derivedFromPair = std::make_shared<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>(l, r);
}

void FsmNode::setPair(const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p)
{
	derivedFromPair = p;
}

bool FsmNode::isDerivedFrom(const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p) const
{
	return derivedFromPair != nullptr && *derivedFromPair == *p;
}

std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> FsmNode::getPair() const
{
	return derivedFromPair;
}

std::shared_ptr<FsmNode> FsmNode::apply(const int e, OutputTrace & o)
{
	for (FsmTransition & tr : transitions)
	{
		if (tr.getLabel().getInput() == e)
		{
			o.add(tr.getLabel().getOutput());
			return tr.getTarget();
		}
	}
	return nullptr;
}

OutputTree FsmNode::apply(const InputTrace & itrc)
{
	std::vector<std::shared_ptr<TreeNode>> tnl;
	std::unordered_map<std::shared_ptr<TreeNode>, std::shared_ptr<FsmNode>> t2f;

	std::shared_ptr<TreeNode> root = std::make_shared<TreeNode>();
	OutputTree ot = OutputTree(root, itrc, presentationLayer);

	if (itrc.get().size() == 0)
	{
		return ot;
	}

	t2f [root] = shared_from_this();//insertion

	for (auto it = itrc.cbegin(); it != itrc.cend(); ++ it)
	{
		int x = *it;
		tnl = ot.getLeaves();

		while (!tnl.empty())
		{
			std::shared_ptr<TreeNode> thisTreeNode = tnl.front();
			tnl.erase(tnl.begin());

			std::shared_ptr<FsmNode> thisState = t2f.at(thisTreeNode);

			for (FsmTransition & tr : thisState->getTransitions())
			{
				if (tr.getLabel().getInput() == x)
				{
					int y = tr.getLabel().getOutput();
					std::shared_ptr<FsmNode> tgtState = tr.getTarget();
					std::shared_ptr<TreeNode> tgtNode = std::make_shared<TreeNode>();
					std::shared_ptr<TreeEdge> te = std::make_shared<TreeEdge>(y, tgtNode);
					thisTreeNode->add(te);
					t2f [tgtNode] = tgtState;//insertion
				}
			}
		}
	}
	return ot;
}

std::unordered_set<std::shared_ptr<FsmNode>> FsmNode::after(const InputTrace & itrc)
{
	std::unordered_set<std::shared_ptr<FsmNode>> nodeSet;
	nodeSet.insert(shared_from_this());

	for (auto it = itrc.cbegin(); it != itrc.cend(); ++ it)
	{
		int x = *it;
		std::unordered_set<std::shared_ptr<FsmNode>> newNodeSet;

		for (std::shared_ptr<FsmNode> n : nodeSet)
		{
			std::unordered_set<std::shared_ptr<FsmNode>> ns = n->afterAsSet(x);
			newNodeSet.insert(ns.begin(), ns.end());
		}
		nodeSet = newNodeSet;
	}
	return nodeSet;
}

std::vector<std::shared_ptr<FsmNode>> FsmNode::after(const int x)
{
	std::vector<std::shared_ptr<FsmNode>> lst;
	for (FsmTransition & tr : transitions)
	{
		if (tr.getLabel().getInput() == x)
		{
			lst.push_back(tr.getTarget());
		}
	}
	return lst;
}

std::unordered_set<std::shared_ptr<FsmNode>> FsmNode::afterAsSet(const int x)
{
	std::unordered_set<std::shared_ptr<FsmNode>> lst;

	for (FsmTransition & tr : transitions)
	{
		if (tr.getLabel().getInput() == x)
		{
			lst.insert(tr.getTarget());
		}
	}
	return lst;
}

void FsmNode::setColor(const int pcolor)
{
	color = pcolor;
}

int FsmNode::getColor()
{
	return color;
}

std::shared_ptr<DFSMTableRow> FsmNode::getDFSMTableRow(const int maxInput)
{
	std::shared_ptr<DFSMTableRow> r = std::make_shared<DFSMTableRow>(id, maxInput);

	IOMap& io = r->getioSection();
	I2PMap& i2p = r->geti2postSection();

	for (FsmTransition tr : transitions)
	{
		int x = tr.getLabel().getInput();

		/*Check whether transitions from this state are nondeterministic.
		This is detected when detecting a second transition triggered
		by the same input. In this case we cannot calculate a  DFSMTableRow.*/
		if (io.at(x) >= 0)
		{
			std::cout << "Cannot calculated DFSM table for nondeterministic FSM." << std::endl;
			return nullptr;
		}

		io[x] = tr.getLabel().getOutput();
		i2p[x] = tr.getTarget()->getId();
	}
	return r;
}

bool FsmNode::distinguished(const std::shared_ptr<FsmNode> otherNode, const std::vector<int>& iLst)
{
	InputTrace itr = InputTrace(iLst, presentationLayer);
	OutputTree ot1 = apply(itr);
	OutputTree ot2 = otherNode->apply(itr);
	return !(ot1 == ot2);
}

std::shared_ptr<InputTrace> FsmNode::distinguished(const std::shared_ptr<FsmNode> otherNode, std::shared_ptr<Tree> w)
{
	IOListContainer iolc = w->getIOLists();
	std::shared_ptr<std::vector<std::vector<int>>> inputLists = iolc.getIOLists();

	for (std::vector<int>& iLst : *inputLists)
	{
		if (distinguished(otherNode, iLst))
		{
			return std::make_shared<InputTrace>(iLst, presentationLayer);
		}
	}
	return nullptr;
}

InputTrace FsmNode::calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<PkTable>>& pktblLst, const int maxInput)
{
	/*Determine the smallest l >= 1, such that this and otherNode are
	distinguished by P_l, but not by P_(l-1).
	Note that table P_n is found at pktblLst.get(n-1)*/
	unsigned int l;
	for (l = 1; l < pktblLst.size(); ++ l)
	{
		/*Two nodes are distinguished by a Pk-table, if they
		reside in different Pk-table classes.*/
		std::shared_ptr<PkTable> pk = pktblLst.at(l - 1);
		if (pk->getClass(this->getId()) != pk->getClass(otherNode->getId()))
		{
			break;
		}
	}

	std::shared_ptr<FsmNode> qi = shared_from_this();
	std::shared_ptr<FsmNode> qj = otherNode;

	InputTrace itrc = InputTrace(presentationLayer);

	for (int k = 1; l - k > 0; ++ k)
	{
		std::shared_ptr<PkTable> plMinK = pktblLst.at(l - k - 1);
		/*Determine input x such that qi.after(x) is distinguished
		from qj.after(x) in plMinK*/

		for (int x = 0; x <= maxInput; ++ x)
		{
			/*We are dealing with completely defined DFSMs,
			so after() returns an ArrayList containing exactly
			one element.*/
			std::shared_ptr<FsmNode> qiNext = qi->after(x).front();
			std::shared_ptr<FsmNode> qjNext = qj->after(x).front();

			if (plMinK->getClass(qiNext->getId() != plMinK->getClass(qjNext->getId())))
			{
				qi = qiNext;
				qj = qjNext;
				itrc.add(x);
				break;
			}
		}
	}

	/*Now the case l == k. qi and qj must be distinguishable by at least
	one input*/
	for (int x = 0; x <= maxInput; ++ x)
	{
		OutputTrace oti = OutputTrace(presentationLayer);
		OutputTrace otj = OutputTrace(presentationLayer);
		qi->apply(x, oti);
		qj->apply(x, otj);
		if (oti.get().front() != otj.get().front())
		{
			itrc.add(x);
			break;
		}
	}
	return itrc;
}

InputTrace FsmNode::calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput)
{
	InputTrace itrc = InputTrace(presentationLayer);
	int q1 = this->getId();
	int q2 = otherNode->getId();

	/*Now we know that this and otherNode are NOT distinguished by OFSM-Table-0.
	Determine the smallest l >= 1, such that this and otherNode are
	distinguished by OFSM-Table l, but not by OFSM-table (l-1).
	Note that table OFSM-table n is found at ofsmTblLst.get(n).*/
	unsigned int l;
	for (l = 1; l < ofsmTblLst.size(); ++ l)
	{
		/*Two nodes are distinguished by a OFSM-table, if they
		reside in different OFSM-table classes.*/
		std::shared_ptr<OFSMTable> ot = ofsmTblLst.at(l);
		if (ot->getS2C().at(q1) != ot->getS2C().at(q2))
		{
			break;
		}
	}

	for (int k = 1; l - k > 0; ++ k)
	{
		std::shared_ptr<OFSMTable> ot = ofsmTblLst.at(l - k);

		/*Determine IO x/y such that qi.after(x/y) is distinguished
		from qj.after(x/y) in ot*/
		for (int x = 0; x <= maxInput; ++ x)
		{
			for (int y = 0; y <= maxOutput; ++ y)
			{
				int q1Post = ot->get(q1, x, y);
				int q2Post = ot->get(q2, x, y);

				if (q1Post < 0 || q2Post < 0)
				{
					continue;
				}

				if (ot->getS2C().at(q1Post) != ot->getS2C().at(q2Post))
				{
					itrc.add(x);

					/*Set q1,q2 to their post-states under x/y*/
					q1 = q1Post;
					q2 = q2Post;
					break;
				}
			}
		}
	}

	/*Now the case l == k. q1 and q2 must be distinguishable by at least
	one IO in OFSM-Table-0*/
	std::shared_ptr<OFSMTable> ot0 = ofsmTblLst.front();
	for (int x = 0; x <= maxInput; ++ x)
	{
		for (int y = 0; y <= maxOutput; ++ y)
		{
			if ((ot0->get(q1, x, y) < 0 && ot0->get(q2, x, y) >= 0) || (ot0->get(q1, x, y) >= 0 && ot0->get(q2, x, y) < 0))
			{
				itrc.add(x);
				return itrc;
			}
		}
	}
	return itrc;
}

bool FsmNode::isObservable() const
{
	/*If a label already contained in lblSet is also
	attached to at least one other transition, the
	node violates the observability condition*/
	std::unordered_set<FsmLabel> lblSet;
	for (FsmTransition tr : transitions)
	{
		FsmLabel lbl = tr.getLabel();
		if (!lblSet.insert(lbl).second)
		{
			return false;
		}
	}
	return true;
}

bool FsmNode::isDeterministic() const
{
	std::unordered_set<int> inputSet;

	/*Check if more than one outgoing transition
	is labelled with the same input value*/
	for (FsmTransition tr : transitions)
	{
		int inp = tr.getLabel().getInput();
		if (!inputSet.insert(inp).second)
		{
			return false;
		}
	}
	return true;
}

std::ostream & operator<<(std::ostream & out, const FsmNode & node)
{
	for (FsmTransition tr : node.transitions)
	{
		out << tr << std::endl;
	}
	return out;
}

bool operator==(FsmNode const & node1, FsmNode const & node2)
{
	if (node1.id == node2.id)
	{
		return true;
	}
	return false;
}
