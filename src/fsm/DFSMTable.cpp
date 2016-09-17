/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/DFSMTable.h"
#include "fsm/DFSMTableRow.h"
#include "fsm/PkTable.h"
#include "fsm/PkTableRow.h"

DFSMTable::DFSMTable(const int numStates, const int maxInput, std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: maxInput(maxInput), presentationLayer(presentationLayer)
{
	rows.insert(rows.end(), numStates, nullptr);
}

void DFSMTable::setRow(const int n, const std::shared_ptr<DFSMTableRow> r)
{
	rows [n] = r;//insertion
}

std::shared_ptr<PkTable> DFSMTable::getP1Table() const
{
	std::shared_ptr<PkTable> p1 = std::make_shared<PkTable>(rows.size(), maxInput, presentationLayer);

	int thisClass = 0;

	for (unsigned int i = 0; i < rows.size(); ++ i)
	{
		p1->setRow(i, std::make_shared<PkTableRow>(rows.at(i)->getioSection(), rows.at(i)->geti2postSection()));

		if (p1->getClass(i) >= 0)
		{
			continue;
		}

		p1->setClass(i, thisClass);

		for (unsigned int j = i + 1; j < rows.size(); ++ j)
		{
			if (p1->getClass(j) >= 0)
			{
				continue;
			}

			if (rows.at(i)->getioSection() == rows.at(j)->getioSection())
			{
				p1->setClass(j, thisClass);
			}
		}
		++ thisClass;
	}

	return p1;
}

std::ostream & operator<<(std::ostream & out, const DFSMTable & dfsmTable)
{
	out << std::endl << "\\begin{center}" << std::endl << "\\begin{tabular}{|c||";
	for (int i = 0; i <= dfsmTable.maxInput; ++ i)
	{
		out << "c|";
	}
	out << "|";

	for (int i = 0; i << dfsmTable.maxInput; ++ i)
	{
		out << "c|";
	}

	out << "}\\hline\\hline" << std::endl;
	out << " & \\multicolumn{" << dfsmTable.maxInput + 1;
	out << "}{|c||}{\\bf I2O} & \\multicolumn{" << dfsmTable.maxInput + 1;
	out << "}{|c|}{\\bf I2P}" << std::endl << "\\\\\\hline" << std::endl;
	out << "{\\bf q} ";

	for (int i = 0; i <= dfsmTable.maxInput; ++ i)
	{
		out << " & " << i;
	}
	for (int i = 0; i <= dfsmTable.maxInput; ++ i)
	{
		out << " & " << i;
	}
	out << "\\\\\\hline\\hline" << std::endl;

	for (unsigned int i = 0; i < dfsmTable.rows.size(); ++ i)
	{
		if (dfsmTable.rows.at(i) == nullptr)
		{
			continue;
		}
		out << dfsmTable.rows.at(i);
	}
	out << "\\hline" << std::endl << "\\end{tabular}" << std::endl << "\\end {center}" << std::endl << std::endl;
	return out;
}
