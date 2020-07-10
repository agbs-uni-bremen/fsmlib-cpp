/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_DFSMTABLE_H_
#define FSM_FSM_DFSMTABLE_H_

#include <memory>
#include <vector>

#include "fsm/typedef.inc"
#include "fsm/DFSMTableRow.h"
#include "fsm/PkTable.h"
#include "interface/FsmPresentationLayer.h"

class DFSMTableRow;
class PkTable;
class FsmPresentationLayer;

/**
 Class representing DFSM tables
 */
class DFSMTable
{
private:
    /**
     * For each state, one DFSM table row is created
     */
    std::vector< std::shared_ptr<DFSMTableRow> > rows;
    
    /**
     * Maximal value of the input alphabet in range 0..maxInput
     */
    int maxInput;
    
    /**
     * The presentation layer used by the FSMTable
     */
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
public:
    /**
     * Create a DFSMTable
     * @param numStates Size of the DFSMTable
     * @param maxInput Maximum Input of the FSM
     * @param presentationLayer The presentation layer used by the DFSM
     */
    DFSMTable(const int numStates, const int maxInput, std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     * Modify a row of the DFSMTable
     * @param n The id of the row to change
     * @param r The row to insert
     */
    void setRow(const int n, const std::shared_ptr<DFSMTableRow>& r);
    
    std::shared_ptr<DFSMTableRow> getRow(int n) { return rows[n]; }
    
    /**
     * Modify a row of the DFSMTable
     * @return The PkTablep1
     */
    std::shared_ptr<PkTable> getP1Table() const;
    
    /**
     * Output the DFSMTable to a standard output stream in LaTeX
     *  tabular format.
     * @param out The standard output stream to use
     * @param dfsmTable The DFSMTable to print
     * @return The standard output stream used, to allow user to cascade <<
     */
    friend std::ostream & operator<<(std::ostream & out, const DFSMTable & dfsmTable);
};
#endif //FSM_FSM_DFSMTABLE_H_
