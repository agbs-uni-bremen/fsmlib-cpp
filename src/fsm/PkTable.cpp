/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/PkTable.h"
#include "fsm/FsmNode.h"
#include "fsm/Dfsm.h"
#include "fsm/FsmLabel.h"
#include "fsm/FsmTransition.h"

using namespace std;

PkTable::PkTable(const int numStates, const int maxInput, const shared_ptr<FsmPresentationLayer> presentationLayer)
	: s2c(numStates), maxInput(maxInput), presentationLayer(presentationLayer)
{
	rows.insert(rows.end(), numStates, nullptr);
}

PkTable::PkTable(const int numStates, const int maxInput, const vector<shared_ptr<PkTableRow>> rows, const shared_ptr<FsmPresentationLayer> presentationLayer)
	: rows(rows), s2c(numStates), maxInput(maxInput), presentationLayer(presentationLayer)
{

}

void PkTable::setRow(const int s, const shared_ptr<PkTableRow> row)
{
	rows[s] = row;
}

void PkTable::setClass(const int n, const int c)
{
	s2c[n] = c;
}

int PkTable::getClass(const int n) const
{
	return s2c.at(n);
}

int PkTable::maxClassId() const
{
	int id = 0;
	for (unsigned int i = 0; i < s2c.size(); ++i)
	{
		if (s2c.at(i) > id)
		{
			id = s2c.at(i);
		}
	}
	return id;
}

shared_ptr<PkTable> PkTable::getPkPlusOneTable() const
{
	shared_ptr<PkTable> pkp1 = make_shared<PkTable>(rows.size(), maxInput, rows, presentationLayer);

	int thisClass = 0;
	int thisNewClassId = maxClassId() + 1;
	shared_ptr<PkTableRow> refRow;
	shared_ptr<PkTableRow> newClassRefRow;
	bool haveNewClasses = false;

	do
	{
		refRow = nullptr;
		for (unsigned int i = 0; i < rows.size(); ++i)
		{
			if (s2c.at(i) != thisClass)
			{
				continue;
			}

			if (pkp1->getClass(i) >= 0)
			{
				continue;
			}

			if (refRow == nullptr)
			{
				refRow = rows.at(i);
				pkp1->setClass(i, thisClass);
				continue;
			}

			if ( refRow->isEquivalent(*rows.at(i), s2c) )
			{
				pkp1->setClass(i, thisClass);
			}
			else
			{
				haveNewClasses = true;
				newClassRefRow = rows.at(i);
				pkp1->setClass(i, thisNewClassId);

				for (unsigned int j = i + 1; j < rows.size(); ++j)
				{
					if ( s2c.at(j) == thisClass and
                         newClassRefRow->isEquivalent(*rows.at(j), s2c) )
					{
						pkp1->setClass(j, thisNewClassId);
					}
				}

				newClassRefRow = nullptr;
				++thisNewClassId;
			}
		}
		++thisClass;
	} while (refRow != nullptr);

	return haveNewClasses ? pkp1 : nullptr;
}

Dfsm PkTable::toFsm(string name, const int maxOutput)
{
    string minFsmName("");
	vector<shared_ptr<FsmNode>> nodeLst;
    
    /* We need a new presentation layer.
     * Input and output names are the same as for the original FSM,
     * but states should have new names including the set of 
     *  original nodes that are equivalent.
     */
    vector<string> minState2String;
    for (int i = 0; i <= maxClassId(); ++i) {
        string newName(minFsmName + getMembers(i));
        minState2String.push_back(newName);
    }
    
    shared_ptr<FsmPresentationLayer> minPl =
    make_shared<FsmPresentationLayer>(presentationLayer->getIn2String(),
                                      presentationLayer->getOut2String(),
                                      minState2String);
    

	/*Create the FSM states, one for each class*/
	for (int i = 0; i <= maxClassId(); ++i)
	{
		shared_ptr<FsmNode> newNode = make_shared<FsmNode>(i, "", minPl);
		nodeLst.push_back(newNode);
	}

	/*For each FSM state, add outgoing transitions*/
	for (shared_ptr<FsmNode> srcNode : nodeLst)
	{
		int classId = srcNode->getId();
		shared_ptr<PkTableRow> row = nullptr;

		for (unsigned int i = 0; i < rows.size() && row == nullptr; ++i)
		{
			if (classId == s2c.at(i))
			{
				row = rows.at(i);
			}
		}

		for (int x = 0; x <= maxInput; ++ x)
		{
			int y = row->getIOMap().at(x);
			int cAux = row->getI2PMap().at(x);
            
            // If the Fsm is not completely specified,
            // it may be the case that no transition for this
            // input exists.
            if ( cAux < 0 ) continue;
            
			int cTarget = s2c.at(cAux);
			shared_ptr<FsmNode> tgtNode = nullptr;
            
            // Find the new FsmNode in the minimised FSM
            // which has cTarget as node id
			for (shared_ptr<FsmNode> node : nodeLst)
			{
				if (node->getId() == cTarget)
				{
					tgtNode = node;
					break;
				}
			}
			shared_ptr<FsmLabel> lbl = make_shared<FsmLabel>(x, y, minPl);
			srcNode->addTransition(make_shared<FsmTransition>(srcNode,
                                                              tgtNode,
                                                              lbl));
		}
	}

	return Dfsm(minFsmName, maxInput, maxOutput, nodeLst, minPl);
}

string PkTable::getMembers(const int c) const
{
	string memSet = "{";
	bool first = true;
	for (unsigned int i = 0; i < rows.size(); ++i)
	{
		if (s2c.at(i) != c)
		{
			continue;
		}

		if (!first)
		{
			memSet += ",";
		}
		first = false;
        memSet += presentationLayer->getStateId(i,"");
	}
	memSet += "}";
	return memSet;
}



ostream & operator<<(ostream & out, const PkTable & pkTable)
{
    
    // Create the table header
    out << endl << "\\begin{center}" << endl << "\\begin{tabular}{|c|c||";
    for (int i = 0; i <= pkTable.maxInput; ++i)
    {
        out << "c|";
    }
    out << "|";
    
    for (int i = 0; i <= pkTable.maxInput; ++i)
    {
        out << "c|";
    }
    
    out << "}\\hline\\hline" << endl;
    out << " & & \\multicolumn{" << pkTable.maxInput + 1;
    out << "}{|c||}{\\bf I2O} & \\multicolumn{" << pkTable.maxInput + 1;
    out << "}{|c|}{\\bf I2P}" << endl << "\\\\\\hline" << endl;
    out << "{\\bf [q]} & {\\bf q} ";
    
    for (int i = 0; i <= pkTable.maxInput; ++i)
    {
        out << " & " << i;
    }
    for (int i = 0; i <= pkTable.maxInput; ++i)
    {
        out << " & " << i;
    }
    out << "\\\\\\hline\\hline" << endl;
    
    
    // Output each table row
    for (unsigned int i = 0; i < pkTable.rows.size(); ++i)
    {
        if (pkTable.rows.at(i) == nullptr)
        {
            continue;
        }
        out << pkTable.s2c.at(i)  << " & " << i << " " << *pkTable.rows.at(i);
    }
    
    // Create the table footer
    out << "\\hline" << endl << "\\end{tabular}" << endl << "\\end {center}" << endl << endl;
    return out;
}
