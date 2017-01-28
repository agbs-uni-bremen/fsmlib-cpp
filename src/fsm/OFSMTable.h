/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_OFSMTABLE_H_
#define FSM_FSM_OFSMTABLE_H_

#include <memory>
#include <vector>

#include "fsm/typedef.inc"

class OFSMTableRow;
class FsmPresentationLayer;
class FsmNode;
class Fsm;

/**
Table representation for observable FSMs

The table structure is as follows.
1. One row for each FSM state, the FSM state number in 0..nodes.length-1
   is the table index.
2. Each row r has the following matrix structure.
   (a) The matrix is indexed over inputs i in 0..maxInput and outputs o in 0..maxOutput
   (b) If FSM state r has no transition labelled with i/o, then the matrix contains
       value -1 at position [i][o]
   (c) If FSM state r transits with i/o to post-state r', then the matrix contains
       value r' at position [i][o]

\note This representation is well-defined if and only if the FSM is observable.

Additionally, each OFSMTable contains an S2CMap which maps FSM states to their equivalence
class associated with the current OFSMTable
*/
class OFSMTable
{
private:
	/**
	The number of states
	*/
	int numStates;

	/**
	The maximal input
	*/
	int maxInput;

	/**
	The maximal output
	*/
	int maxOutput;

    /** 
     *  The number of this OFSM table: The first one to be created
     *  gets tblId 0, the last one created during the mnimisation 
     *  process gets the highest tblId and carries the final mapping
     *  of states to equivalence classes in s2c.
     */
	int tblId;

	/** Mapping from a given state to its current class */
	S2CMap s2c;

	/** Rows of the OFSM table */
	std::vector<std::shared_ptr<OFSMTableRow>> rows;

	/**
	The presentation layer used by the OFSMTable
	*/
	const std::shared_ptr<FsmPresentationLayer> presentationLayer;

	/**
	Create new OFSMTable from the initial one.
	@return OFSMTable where two FSM nodes are associated with the same
	class if and only if they have outgoing transitions for
	exactly the same set of input/output labels
	*/
	std::shared_ptr<OFSMTable> nextAfterZero();
public:
	/**
	This constructor creates the initial OFSMTable for an observable FSM.
	All nodes are associated with equivalence class 0.
	@param nodes    array of FSM states
	@param maxInput input alphabet is in range 0..maxInput
	@param maxOutput output alphabet is in range 0..maxOutput
	@param presentationLayer The presentation layer used by the OFSMTable
	*/
	OFSMTable(const std::vector<std::shared_ptr<FsmNode>>& nodes, const int maxInput, const int maxOutput, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
	
	//TODO
	OFSMTable(const int numStates, const int maxInput, const int maxOutput, const std::vector<std::shared_ptr<OFSMTableRow>>& rows, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
	
	//TODO
	int getId();

	//TODO
	S2CMap getS2C();

	//TODO
	void setS2C(const S2CMap & ps2c);

	/**
	Get table entry for given state and I/O
	@param id state id in range 0..(nodes.length-1)
	@param x  input value in range 0..maxInput
	@param y  output value in range 0..maxOutput
	@return -1 if FSM state id does not have an outgoing transition labelled by x/y
	n >= 0 if FSM state id has an outgoing transition labelled by x/y,
	which ends at FSM state n
	*/
	int get(const int id, const int x, const int y);

	//TODO
	int maxClassId() const;

	/**
	Create the next OFSMTable on the basis of the current table.
	@return The next OFSMTable if it exists,
	null otherwise
	*/
	std::shared_ptr<OFSMTable> next();

	/**
	Return members of an equivalence class c as set string
	*/
	std::string getMembers(const int c) const;
    
    
    /** 
     *  Compare two columns of the OFSM table. Recall that each
     *  column of an OFSM table is uniquely identified by a pair
     *  of input value x and output value y.
     *
     *  @param x1 Input value associated with the first column
     *  @param y1 Output value associated with the first column
     *  @param x2 Input value associated with the second column
     *  @param y2 Output value associated with the second column
     *  @return true iff the first column has the same entries as the second one
     */
    bool compareColumns(int x1, int y1, int x2, int y2);

	/**
	Create minimised FSM from OFSM-Table
	@param name the name for the FSM to be created
	@return the FSM instance
	\note This operation expects that this OFSMTable is already the last
	one in a sequence of OFSMTable transformations, so that it really represents
	the minimised FSM.
	*/
	Fsm toFsm(const std::string & name) const;

	/**
	Produce a stream given a tabular LaTeX representation of
	the OFSMTable, to be included in a LaTeX document.
	*/
	friend std::ostream & operator<<(std::ostream & out, const OFSMTable & ofsmTable);
};
#endif //FSM_FSM_OFSMTABLE_H_
