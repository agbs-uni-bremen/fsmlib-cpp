/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_OFSMTABLEROW_H_
#define FSM_FSM_OFSMTABLEROW_H_

#include <memory>
#include <vector>

#include "fsm/Int2IntMap.h"
#include "fsm/typedef.inc"

/**
Class representing one table row of an OFSMTable
*/
class OFSMTableRow
{
private:
	/**
	The maximal input
	*/
	int maxInput;

	/**
	The maximal output
	*/
	int maxOutput;

	/**
	The OFSMTableRow itself
	*/
	std::vector<std::vector<int>> matrix;
public:
	/**
	Create a OFSMTableRow
	@param maxInput The maximal input
	@param maxOutput The maximal output
	*/
	OFSMTableRow(const int maxInput, const int maxOutput);

	/**
	Set for the element at the position i / j as postState
	@param i The line number
	@param j The column number
	@param postState The value
	@return The element
	*/
	void set(const int i, const int j, const int postState);

	/**
	Getter for the element at the position i / j
	@param i The line number
	@param j The column number
	@return The element
	*/
	int get(const int i, const int j) const;

	/**
	Return false if and only if this row represents a state that
	is mapped to a post state by input/output i/o, while the other state
	represented by r has no emanating transition labelled with i/o or
	vice versa.
	*/
    bool ioEquals(const std::shared_ptr<OFSMTableRow>& r) const;

	/**
	Return false if and only if this row represents a state that
	is mapped to a post state by input/output i/o which is associated
	with another equivalence class than the post state r.get(i,o).
	*/
    bool classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow>& r);
};
#endif //FSM_FSM_OFSMTABLEROW_H_
