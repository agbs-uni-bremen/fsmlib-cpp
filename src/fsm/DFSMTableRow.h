/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_DFSMTABLEROW_H_
#define FSM_FSM_DFSMTABLEROW_H_

#include <iostream>

#include "fsm/Int2IntMap.h"
#include "fsm/typedef.inc"

/**
Class representing one row of a DFSM table
*/
class DFSMTableRow
{
private:
	/**
	The state associated with this row
	*/
	int state;

	/**
	The input-to-output section, indexed over the input alphabet
	*/
	IOMap ioSection;

	/**
	The input-to-post-state section, indexed over the input alphabet
	*/
	I2PMap i2postSection;
public:
	/**
	Create a DFSMTableRow
	\param q The number of states
	\param maxInput The number of input
	*/
	DFSMTableRow(const int q, const int maxInput);

	/**
	Getter for the input/output section
	\return The input/output section
	*/
	IOMap getioSection() const;

	/**
	Getter for the input/post section
	\return The input/post section
	*/
	I2PMap geti2postSection() const;

	/**
	Output the DFSMTableRow to a standard output stream
	\param out The standard output stream to use
	\param dfsmTableRow The DFSMTableRow to print
	\return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const DFSMTableRow & dfsmTableRow);
};
#endif //FSM_FSM_DFSMTABLEROW_H_