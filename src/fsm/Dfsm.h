/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_DFSM_H_
#define FSM_FSM_DFSM_H_

#include <memory>
#include <string>
#include <vector>
#include <stack>

#include "fsm/Fsm.h"
#include "fsm/ConvergenceGraph.h"

class PkTable;
class IOTrace;
class TreeNode;
class DFSMTable;
class SegmentedTrace;


namespace Json {
    class Value;
}

class Dfsm : public Fsm
{
private:
	//TODO
	std::shared_ptr<DFSMTable> dfsmTable;

	//TODO
	std::vector<std::shared_ptr<PkTable>> pktblLst;

	/**
	Create a DFSMTable from the DFSM
	@return The DFSMTable created
	*/
	std::shared_ptr<DFSMTable> toDFSMTable() const;
    
    std::shared_ptr<FsmPresentationLayer> createPresentationLayerFromCsvFormat(const std::string& fname);
    
    std::shared_ptr<FsmPresentationLayer> createPresentationLayerFromCsvFormat(const std::string& fname,const std::shared_ptr<FsmPresentationLayer>& pl);
    
    void createDfsmTransitionGraph(const std::string& fname);
    
    /**
     *   distTraces[n][m] contains a vector of pointers to
     *   traces distinguishing this from FsmNode number n
     */
    std::vector< std::vector< std::vector< std::shared_ptr< std::vector<int> > > > > distTraces;
    
    void initDistTraces();
    
    std::vector< std::shared_ptr< std::vector<int> > > calcDistTraces(FsmNode& s1,
                                                                      FsmNode& s2);
    
    std::vector< std::shared_ptr< std::vector<int> > > calcDistTraces(size_t l,
                                                                      std::shared_ptr< std::vector<int> > trc,
                                                                      int id1,
                                                                      int id2);
    
    std::vector< std::shared_ptr< std::vector<int> > > calcDistTraces(std::shared_ptr< std::vector<int> > trc,
                                                                      int id1,
                                                                      int id2);



    // TODO: doc
    void spyhDistinguish(const std::vector<int>& trace, std::unordered_set<std::shared_ptr<InputTrace>> traces, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph);
    std::pair<size_t,std::stack<int>> spyhGetPrefixOfSeparatingTrace(const std::vector<int>& trace1, const std::vector<int>& trace2, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph);
    size_t spyhEstimateGrowthOfTestSuite(const std::shared_ptr<FsmNode> u, const std::shared_ptr<FsmNode> v, int input);
    void spyhAppendSeparatingSequence(const std::vector<int>& traceToAppendTo, const std::vector<int>& traceToAppend, std::shared_ptr<Tree> testSuite, ConvergenceGraph& graph);

public:
	/**
	Create a DFSM from a file description
	@param fname The name of the file containing the FSM informations
	@param fsmName The name of the DFSM
	@param maxNodes Number of DFSM states
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param presentationLayer The presentation layer used by the DFSM
     
     \note This constructor is deprecated and will be removed in a future version
	*/
	Dfsm(const std::string & fname,
         const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    /**
     Create a DFSM from a file description
     @param fname The name of the file containing the FSM informations
     @param presentationLayer The presentation layer used by the DFSM
     @param fsmName The name of the DFSM
     
     The parameters maxNodes, maxInput, and maxOutput used in the constructor above
     are determined from the input file specifying the FSM.
     */
    Dfsm(const std::string & fname,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer,
         const std::string & fsmName);

    
    /**
     Create a random mutant of this DFSM
     */
    void createAtRandom();
    
	/**
	Random creation of a completely defined deterministic FSM
	@param fsmName The name of the DFSM
	@param maxNodes Number of DFSM states
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	Create a DFSM from a list of nodes
	@param fsmName The name of the DFSM
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param lst The list of nodes of the DFSM
	@param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxInput,
         const int maxOutput,
         const std::vector<std::shared_ptr<FsmNode>>& lst,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	Create a DFSM from the equivalent deterministic FSM
	@param fsm The deterministic FSM
	*/
	Dfsm(const Fsm & fsm);
    
    
    /**
     * Create a DFSM from a file written in csv format.
     * @item fields are separated by ;
     * @item First row starts with empty field, followed by
     *       input event names to be used by presentation layer
     * @item Consecutive rows are ordered by DFSM states,
     *       starting with the initial state
     * @item Each of these rows contain the state name in the first 
     *       column, as it should be used by the presentation layer.
     * @item In row indexed by state s, column indexed by input x
     *       contains entry 
     *                  s' / y
     *       where s' is the post-state reached when applying x in state s,
     *       and y is the associated output.
     *
     * @param fname Filename with DFSM specification in csv format.
     * @param fsmName Name of the DFSM to be used in output files etc.
     */
    Dfsm(const std::string& fname,
         const std::string& fsmName);
    
    /**
     *  This constructor acts like the previous, but uses a
     *  a given presentation layer. Typically, this is required
     *  if an implementation DFSM is read from a file where
     *  inputs and outputs occur in a different order, and
     *  where some inputs/outputs occurring in the reference model
     *  do not occur in the implementation. Also, the implementation
     *  might use additional inputs and/or outputs that do not
     *  occur in the reference model. These are appended to the respective
     *  lists of the given presentation layer.
     *
     *  The constructor ensures that both implementation model and
     *  reference model use the same internal numbers for the
     *  same inputs and outputs.
     */
    Dfsm(const std::string& fname,
         const std::string& fsmName,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    
    /**
     *  Construct an DFSM from a json model file. The DFSM is completely
     *  specified: in every state of the json model, inputs x that do
     *  not occur in any outgoing transition, are associated with
     *  self-loop transitions that have an additional NO OPERATION
     *  action as output. The NOP action is added as last action
     *  to the output alphabet.
     */
    Dfsm(const Json::Value& jsonModel);
    
    /**
     *  This constructor acts like the previous, but uses a 
     *  a given presentation layer. Typically, this is required
     *  if an implementation DFSM is read from a file where
     *  inputs and outputs occur in a different order, and 
     *  where some inputs/outputs occurring in the reference model
     *  do not occur in the implementation. Also, the implementation
     *  might use additional inputs and/or outputs that do not
     *  occur in the reference model. These are appended to the respective 
     *  lists of the given presentation layer.
     *
     *  The constructor ensures that both implementation model and
     *  reference model use the same internal numbers for the
     *  same inputs and outputs.
     */
    Dfsm(const Json::Value& jsonModel,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);


	/**
	Minimise this DFSM
	As a side effect, create the DFSM table and all Pk tables needed.
	@return the new minimised DFSM
	*/
	Dfsm minimise();

	/**
	Output the DFSM table and the Pk-tables in LaTeX format.
	The output file name is hard-coded as tables.tex.
	\note This requires that the operation minimise() has been called before.
	*/
	void printTables() const;

    void calcPkTables();

	/**
	Calculate characterisation set according to Gill's algorithm
	based on Pk-tables
	\pre the DFSM has already been minimised, and it is completely defined
	@return Characterisation set
	*/
	IOListContainer getCharacterisationSet();

	/**
	Apply input trace to the initial state of the FSM.
	\note This operation is only applicable to deterministic FSMs
	@param i input trace
	@return IOTrace containing the maximal sub-trace of i that could be processed as inputs
	(if FSM is completely defined this is always the complete trace i) and the associated
	outputs of the same length (the outputs are uniquely determined by the inputs, since the FSM
	is deterministic.
	*/
	IOTrace applyDet(const InputTrace & i);

	/**
	Check whether IOTrace is in the language of the DFSM
	@param io IOTrace to be checked against the DFSM
	@return True if io is in the language of the DFSM,
	False otherwise
	*/
	bool pass(const IOTrace & io);

   /**
	* Perform test generation by means of the W-Method.
    * The DFSM this method is applied to is regarded as the reference
    * model. If an implementation DFSM passes this test suite, this proves
    * language equivalence (I/O-equivalence) between reference DFSM and
    * implementation DFSM, provided that the implementation DFSM in minimised
    * form does not have more than numAddStates additional states, when compared
    * to the minimised reference DFSM. If this assumption does not hold,
    * the test suite may not uncover certain errors in the implementation DFSM.
    *
    * The reference DFSM will first be minimised, and then the
    * proper W-Method is applied to the minimised DFSM.
	* @param numAddStates The maximal number of additional states,
    *                     which the implementation DFSM in minimised 
    *                     for may have, when compared to the reference
    *                     model in minimised form.
	* @return A test suite
    *
    * @note The size of the test suite to be produced grows exponentially
    *        with numAddStates
    *
    * @note If it is already known that the reference DFSM is minimal,
    *       then method wMethodOnMinimisedDfsm() should rather be used,
    *       since it avoids unnecessary minisation steps.
	*/
	IOListContainer wMethod(const unsigned int numAddStates);

    /**
     *  Apply the W-Method on a DFSM that is already minimised
     */
    IOListContainer wMethodOnMinimisedDfsm(const unsigned int numAddStates);


	/**
	* Perform test generation by means of the Wp Method. The algorithm
	* we have implemented is applicable to both nondeterministic and
	* deterministic FSMs. It relies, however, on the existence of
	* OFSM tables which we only calculate for nondeterministic FSMs.
	* Therefore we need a wrapper method for DFSMs which first
	* calculates OFSM tables (by means of a call to minimiseObervableFSM())
	* and then calls the wpMethod() operation of the super class Fsm.
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     form may have, when compared to the reference
     *                     model in minimised form.
	* @return A test suite
	*/
	IOListContainer wpMethod(const unsigned int numAddStates);

    /**
     * Apply the Wp Method on a DFSM that is already minimised
     */
    IOListContainer wpMethodOnMinimisedDfsm(const unsigned int numAddStates);

    /**
     * WORK IN PROGRESS
     * Perform test generation by means of the HSI-Method. The algorithm
     * we have implemented is applicable to both nondeterministic and
     * deterministic FSMs.
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     form may have, when compared to the reference
     *                     model in minimised form.
     * @return a test suite
     */
    IOListContainer hsiMethod(const unsigned int numAddStates);
    /**
     * Perform test generation by means of the T-Method. The algorithm 
     * applies to deterministic FSMs which are completely defined.
     * It produces test cases from the transition cover. The resulting 
     * test suite is complete for I/O-equivalence in the fault domain
     * of all completely defined DFSMs M' differing from the reference
     * model M only by zero or more output faults. The number of in M and M'
     * must be identical, and the transition arrows and input labels
     * must be identical as well to ensure completeness.
     */
    IOListContainer tMethod();
    
    
    /**
     *  Perform test generation by means of the H-Method, as 
     *  described in Theorem 1 of
     *    Rita Dorofeeva1, Khaled El-Fakih, and Nina Yevtushenko:
     *    An Improved Conformance Testing Method. 
     *    F. Wang (Ed.): FORTE 2005, LNCS 3731, pp. 204 – 218, 2005. 
     *    IFIP International Federation for Information Processing 2005
     *
     *  This implementation requires the DFSM to be already minimised and
     *  completely specified.
     *
     *  @note This implementation is still under construction.
     *  @note Further operations implementing the variants of the
     *        H-Method for nondeterministic FSMs which are not
     *        completely specified will be added to the library 
     *        in the future.
     */
    IOListContainer hMethodOnMinimisedDfsm(const unsigned int numAddStates);
    
    
    /**
     *  Output DFSM in tabular format as *.csv file
     *  The table consists of one row per state and one column per input.
     *  The table entry for state s and input x consists of the pair
     *  s'/y, where s' is the post-state of the transition triggered
     *  by x in state s, and y is the corresponding output.
     *  @param fname Name of the output file without extension; 
     *         extension .csv is added internally to the filename.
     */
    void toCsv(const std::string& fname);

    /**
     * This method searches in a given Tree for a trace
     * distinguishing s_i=s_0-after(alpha) and s_j=s_0-after(beta).
     * If no such trace could be found, this method tries to lengthen
     * an existing trace in tree so that it distinguishes s_i and s_j.
     * If this is not possible, a new distinguishing trace is calculated.
     * @param iAlpha InputTrace that leads to state s_i
     * @param iBeta InputTrace that leads to state s_j
     * @param Search in this Tree for a distinguishing Trace for s_i and s_j
     *        or search in the leaves of this tree for a trace that
     *        can be lenghtened to distinguish s_i and s_j.
     * @return a distinguishing trace for s_i and s_j
     */
    InputTrace calcDistinguishingTrace(const std::shared_ptr<InputTrace> iAlpha, const std::shared_ptr<InputTrace> iBeta, const std::shared_ptr<Tree> tree);
    
    std::vector<int> calcDistinguishingTrace(std::shared_ptr<SegmentedTrace> alpha,
                                             std::shared_ptr<SegmentedTrace> beta, const std::shared_ptr<TreeNode> treeNode);

    /**
     *  Breadth-first search in a given Tree for a Trace that
     *  distinguishes the states s_i=s0->after(alpha) and s_j=s0->after(beta)
     *  of this DFSM. The distinguishing Trace is returned.
     *  This method is used to avoid appending distinguishing traces
     *  that widen the tree by creating new tree branches.
     *
     *  @param alpha InputTrace that leads to state s_i
     *  @param beta InputTrace that leads to state s_j
     *  @param tree Search in this Tree for a distinguishing Trace
     *  @return distinguishing InputTrace for s_i and s_j.
     *   returns an empty InputTrace if no distinguishing trace could be found
     */
    InputTrace calcDistinguishingTraceInTree(const std::shared_ptr<FsmNode> s_i, const std::shared_ptr<FsmNode> s_j, const std::shared_ptr<Tree> tree);
    InputTrace calcDistinguishingTraceInTree(const std::shared_ptr<InputTrace> alpha, const std::shared_ptr<InputTrace> beta, const std::shared_ptr<Tree> tree);

    /**
     *  Calculate trace that distinguishes the states s_i=s0->after(alpha) and
     *  s_j=s0->after(beta), starting with a path of tree.
     *  The distinguishing trace is returned.
     *  This method tries to lengthen an InputTrace in tree so that it
     *  distinguishes s_i and s_j.
     *  @param s_i state reached by some input trace alpha
     *  @param s_j state reached by some input trace beta
     *  @param tree Search in the paths of this tree for a trace
     *         that can be extended to distinguish s_i and s_j.
     *  @return distinguishing InputTrace for s_i and s_j.
     *   returns an empty InputTrace if no distinguishing trace could be found
     */
    InputTrace calcDistinguishingTraceAfterTree(const std::shared_ptr<FsmNode> s_i, const std::shared_ptr<FsmNode> s_j, const std::shared_ptr<Tree> tree);

    std::vector<std::shared_ptr<PkTable> > getPktblLst() const;
    std::shared_ptr<DFSMTable> getDFSMTable() const { return dfsmTable; }
    
    
    /**
     *  Return true if and only if the two FSM states are distinguishable
     */
    virtual bool distinguishable(const FsmNode& s1, const FsmNode& s2);
    
    /**
     *   Calculate the distinguishability matrix
     */
    void calculateDistMatrix();
    
    /**
     * Return the vector of shortest traces distinguishing s1 and s2
     */
    std::vector< std::shared_ptr< std::vector<int> > > getDistTraces(FsmNode& s1,
                                                                     FsmNode& s2);


                                                                    
    /**
     *  Computes a test suite using the SPYH-method.
     *
     *  This implementation requires the DFSM to be already minimised and
     *  completely specified.
     * 
     * @param numAddStates The number of additional states such that the
     *                     generated test suite is complete for testing
     *                     against SUTs containing up to 
     *                     (numAddStates + size of this fsm) states.
     * 
     * @return A complete test suite.
     */
    IOListContainer spyhMethodOnMinimisedCompleteDfsm(const unsigned int numAddStates);
};
#endif //FSM_FSM_DFSM_H_
