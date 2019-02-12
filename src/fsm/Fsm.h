/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_FSM_H_
#define FSM_FSM_FSM_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <deque>

#include "fsm/FsmVisitor.h"


class Dfsm;
class FsmNode;
class Tree;
class OutputTree;
class InputTrace;
class FsmPresentationLayer;
class OFSMTable;
class IOListContainer;
class TestSuite;

enum Minimal
{
    True, False, Maybe
};

class Fsm
{
protected:
    
    /**
     *  Default constructors without effect - needed by sub-classes
     */
    Fsm(const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    Fsm();
    
    /** Name of the FSM -- appears in nodes when printing the FSM as a dot graph */
    std::string name;
    
    /** FSM states */
    std::vector<std::shared_ptr<FsmNode>> nodes;
    
    std::shared_ptr<FsmNode> currentParsedNode;
    
    /** Maximal value of the input alphabet in range 0..maxInput */
    int maxInput;
    
    /** Maximal value of the output alphabet in range 0..maxOutput */
    int maxOutput;
    
    /** Maximal value of the state id in range 0..maxState */
    int maxState;
    
    /** Integer id of the initial state */
    int initStateIdx;
    
    std::shared_ptr<Tree> characterisationSet;
    Minimal minimal;

    
    std::vector<std::shared_ptr<OFSMTable>> ofsmTableLst;
    std::vector<std::shared_ptr<Tree>> stateIdentificationSets;
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
    std::shared_ptr<FsmNode> newNode(const int id, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p,
                                     std::shared_ptr<FsmPresentationLayer> pl);
    bool contains(const std::deque<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p);
    bool contains(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<FsmNode> n);
    std::shared_ptr<FsmNode> findp(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p);
    void parseLine(const std::string & line);
    void readFsm(const std::string & fname);
    
    void parseLineInitial (const std::string & line);
    void readFsmInitial (const std::string & fname);
    
    
    std::string labelString(std::unordered_set<std::shared_ptr<FsmNode>>& lbl) const;
    
    /**
     *  Return a random seed to be used for random generation
     * of FSMs by public methods createRandomFsm() and
     * createMutant().
     */
    static unsigned int getRandomSeed();
    
    /** 
     * calculate eqivalent inputs for an FSM which already is a prime 
     * machine 
     */
    std::vector< std::unordered_set<int> > getEquivalentInputsFromPrimeMachine();

    /**
     *   Calculate OFSM tables and store them in ofsmTableLst
     */
    void calcOFSMTables();
    
public:
    
    
    /** Copy constructor
     *  A deep copy is created, so that the new Fsm does not
     *  point to any nodes, transitions, or labels of the old FSM
     */
    Fsm(const Fsm& other);
    
    /**
     *  Constructor creating an FSM from file - used only internally
     */
    Fsm(const std::string & fname,
        const std::string & fsmName,
        const int maxNodes,
        const int maxInput,
        const int maxOutput,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     *  Constructor creating an FSM specified in a file.
     *  @param fname Filename of a text file (typically with extension *.fsm),
     *  where each row is formatted as
     *       pre-state input output post-state
     * \item pre-state is a number in range 0..(number of states -1) specifying
     *       the FSM pre-state
     * \item input is a number in range 0..maxInput specifying the input applied to
     *       the pre-state
     * \item output is a number in range 0..maxOutput specifying the output
     *       produced by the FSM when transiting to post-state on 'input'
     * \item post-state is a number in range 0..(number of states -1)
     * The initial state is expected to be the pre-state of the first row.
     * The parameters 'number of states', maxInput, and maxOutput are
     * determined when parsing the input file fname.
     *
     * @param presentationLayer Pointer to instance of a presentation layer
     *        which associates each input number in range 0..maxInput with
     *        an input name, each output number in range 0..maxOutput with
     *        an output name, and each state number in range 0..(number of states -1)
     *        with a state name.
     *
     * @param fsmName name of the FSM
     *
     */
    Fsm(const std::string& fname,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer,
        const std::string& fsmName);
    
    
    /**
     Constructor for creating an FSM from a list of FsmNodes that have
     been created beforehand, together with the outgoing transitions of each FsmNode/.
     @param fsmName  name of the FSM (appears in every node)
     @param maxInput maximal value of the (integer) input alphabet - admissible
     values are 0..maxInput
     @param maxOutput maximal value of (integer) output alphabet - admissible
     values are 0..maxOutput
     */
    Fsm(const std::string & fsmName,
        const int maxInput,
        const int maxOutput,
        const std::vector<std::shared_ptr<FsmNode>> lst,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    
    /**
     *  Create a completely specified FSM at random. Every state in the FSM 
     *  will be reachable, but the FSM may be nondeterministic, non-observable,
     *  and not minimal.
     *   @param fsmName Name of the FSM to be created
     *   @param maxInput Maximal value of the input alphabet, ranging from
     *                   0 to maxInput
     *   @param maxOutput Maximal value of the output alphabet in range 0..maxOutput
     *   @param seed If 0, a "real" random seed will be internally calculated 
     *               using getRandom(). Otherwise the seed value > 0 will
     *               be taken to initialise the random number generation 
     *               with srand(). The first variant is useful for
     *               creating different FSMs every time the constructor
     *               is used. The second variant is useful for producing
     *               repeatable generation sequences of pseudo random FSMs.
     *               Also note for the second variant,
     *               that for every new random instance, another seed
     *               value needs to be supplied by the user, because
     *               otherwise the constructor will always produce the same FSM.
     *   @param maxState  Maximal value of the states in range 0..maxState
     *   @return an FSM created at random according to these specifications.
     */
    static std::shared_ptr<Fsm>
    createRandomFsm(const std::string & fsmName,
                    const int maxInput,
                    const int maxOutput,
                    const int maxState,
                    const std::shared_ptr<FsmPresentationLayer>
                    presentationLayer,
                    const unsigned seed = 0);
    
    
    /**
     *  Create a mutant of the FSM, producing output faults
     *  and/or transition faults only.
     *
     *  The number of states remains the same. If FSM is completely
     *  specified, the same will hold for the mutant.
     */
    std::shared_ptr<Fsm> createMutant(const std::string & fsmName,
                                      const size_t numOutputFaults,
                                      const size_t numTransitionFaults);
    
    
    /**
     * Write FSM to text file using the standard format in each line
     * which is also used for instantiating an FSM from file:
     *    pre-state input output post-state
     * \item pre-state is a number in range 0..(number of states -1) specifying
     *       the FSM pre-state
     * \item input is a number in range 0..maxInput specifying the input applied to
     *       the pre-state
     * \item output is a number in range 0..maxOutput specifying the output
     *       produced by the FSM when transiting to post-state on 'input'
     * \item post-state is a number in range 0..(number of states -1)
     */
    void dumpFsm(std::ofstream & outputFile) const;
    std::shared_ptr<FsmNode> getInitialState() const;
    
    /**
     This is the getter for the name of the FSM
     @return Name of the FSM
     */
    std::string getName() const;
    virtual int getMaxNodes() const;
    int getMaxInput() const;
    int getMaxOutput() const;
    std::vector<std::shared_ptr<FsmNode>> getNodes() const;
    std::shared_ptr<FsmPresentationLayer> getPresentationLayer() const;
    int getInitStateIdx() const;
    void resetColor();
    void toDot(const std::string & fname);

	int getMaxState() const;
	void setMaxState(int maxState);
    
    /**
     Create a new FSM that represents the intersection of this and the other FSM
     @param f the other FSM
     @return a new FSM which equals the intersection of this and f
     */
    Fsm intersect(const Fsm & f);
    
    /**
     * Generate the state cover of an arbitrary FSM
     * (deterministic or nondeterministic, completely specified 
     *  or not, observable or not, minimised or not)
     */
    std::shared_ptr<Tree> getStateCover();
    
    /**
     * Generate the transition cover of an arbitrary FSM
     */
    std::shared_ptr<Tree> getTransitionCover();
    
    /**
     *  Apply an input trace to an FSM and return its
     *  resulting output tree.
     *
     *  @param itrc Input trace to be process on the FSM, starting in
     *              the FSM's initial state
     *  @param  markAsVisited If true, every FSM node visited 
     *              while executing the input trace is marked 
     *              as visited by setting atrubute 'visited' to true
     *
     *  @return The set of outputs created by input trace itrc;
     *          the set is encoded as an OutputTree.
     *
     */
    OutputTree apply(const InputTrace & itrc, bool markAsVisited = false);
    
    /**
     *  Transform an FSM to its observable equivalent.
     */
    Fsm transformToObservableFSM() const;
    
    /**
     Check this FSM with respect to observability
     @return true if and only if the FSM is observable
     */
    bool isObservable() const;
    Minimal isMinimal() const;
    
    /**
     *   Check for unreachable states and remove them from the 
     *   FSM
     *
     *   @param unreachableNodes On termination, this vector
     *          contains pointers to the FSM
     *          nodes that have been removed, due to unreachability.
     *
     *   @return true if and only if at least one unreachable node
     *                has been found and removed.
     */
    bool removeUnreachableNodes(std::vector<std::shared_ptr<FsmNode>>& unreachableNodes);
    
    /**
     Create the minimal observable FSM which is equivalent to this FSM.
     \pre This method can only be applied to an observable OFSM
     @return minimal observable FSM which is equivalent to this FSM
     */
    Fsm minimiseObservableFSM();
    
    /**
     Create the minimal observable FSM which is equivalent to this FSM.
     If this FSM is not observable, an observable equivalent is created
     first. The observable FSM is then minimised.
     @return minimal observable FSM which is equivalent to this FSM
     */
    Fsm minimise();
    bool isCharSet(const std::shared_ptr<Tree> w) const;
    void minimiseCharSet(const std::shared_ptr<Tree> w);
    
    /**
     Calculate the characterisation set W of a (possibly nondeterministic) FSM.
     In addition, calculate the state identification
     sets W_i containing prefixes of input traces
     of W, such that W_i distinguishes FSM state number i from every other
     FSM state.
     \pre The FSM must be observable and minimised
     @return The characterisation set W, represented by an instance of
     IOListContainer
     \note As a side effect, the state identification sets W_i
     are stored in the list stateIdentificationSets.
     */
    IOListContainer getCharacterisationSet();
    
    /**
     * Calculate the state identification sets. The sets are stored
     * in list stateIdentificationSets, ordered by the FSM state numbers.
     * \pre The FSM must be observable and minimal, and the characterisation
     * set must have been previously calculated using operation
     * getCharacterisationSet().
     */
    void calcStateIdentificationSets();
    void calcStateIdentificationSetsFast();

    void appendStateIdentificationSets(const std::shared_ptr<Tree> Wp2) const;
    
    /**
     * Perform test generation by means of the W Method, as applicable
     * to nondeterministc FSMs that do not need to be completely specified.
     *
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     for may have, when compared to the reference
     *                     model in minimised form.
     *
     * The reference machine is
     * first transformed into an observable minimised one. Then the
     * W-Method is applied on the minimised machine, using wMethodOnMinimisedFsm().
     *
     * @return A test suite
     *
     */
    IOListContainer wMethod(const unsigned int numAddStates);
    
    
    /**
     * Perform test generation by means of the W Method, as applicable
     * to nondeterministc FSMs that do not need to be completely specified.
     * @param m Maximum number of states in the observable, minimised FSM
     * reflecting the implementation behaviour. It is assumed that the
     * method is called on an observable, minimised FSM
     *
     * @return A test suite
     *
     */
    IOListContainer wMethodOnMinimisedFsm(const unsigned int m);
    
    
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
     * WORK IN PROGRESS
     * Perform test generation by means of the HSI-Method. The algorithm
     * we have implemented is applicable to both nondeterministic and
     * deterministic FSMs.
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     form may have, when compared to the reference
     *                     model in minimised form.
     * @return A test suite
     */
    IOListContainer hsiMethod(const unsigned int numAddStates);

    
    /**
     *  Return a test suite from a collection of input test cases.
     *  Typically, the IOListContainer has been created by some test
     *  generation method such as the W-Method or the Wp-Method.
     *
     *  Method createTestSuite() executes these input sequences
     *  agains the FSM the method is called on. A TestSuite
     *  instance consists of a vector of OutputTrees. Each OutputTree
     *  instance in the TestSuite contains a single input trace
     *  and a tree whose edges specify the possible outputs resulting
     *  from the input trace, when applied to this FSM.
     *
     *  Recall that the OutputTree is reduced to just an output list for
     *  the given input trace, if the FSM is deterministic.
     *
     *  Recall further that two TestSuite instances can be compared
     *  with respect to I/O-equivalence, using TestSuite method isEquivalentTo().
     *  They can also be checked with respect to a reduction relationship
     *  by using the isReductionOf() method defined for TestSuite instances.
     */
    TestSuite createTestSuite(const IOListContainer & testCases);
    
    /**
     *  Identify equivalent inputs for completely specified FSMs.
     *  Two members x1, x2 of the FSM input alphabet are equivalent,
     *  if and only if from every state s in its prime machine
     *  and for every output y, labels x1/y and x2/y lead to the same
     *  target state or are both undefined in s.
     *
     *   @return A vector of length less or equal the alphabet size.
     *           Each vector element is a set of equivalent inputs.
     *           Of course, the sets are disjoint.
     */
    std::vector< std::unordered_set<int> > getEquivalentInputs();
    
    /**
     *  Return true if the FSM is completely specified.
     *  This means that in every state, for every value
     *  of the input alphabet, at least one outgoing transition
     *  labelled by this input is specified.
     */
    bool isCompletelyDefined() const;
    
    /**
     * Check if FSM is deterministic
     * @return true if FSM is deterministic
     */
    bool isDeterministic() const;
    
    
    void setPresentationLayer(const std::shared_ptr<FsmPresentationLayer> ppresentationLayer);
    
    
    /** Return the number of states in this FSM */
    size_t size() const { return nodes.size(); }
    
    
    
    friend std::ostream & operator<<(std::ostream & out, const Fsm & fsm);
    
    
    /**
     *  Accept an FsmVisitor and initiate a BFS traversal of the Fsm.
     *  All Fsm nodes will be set to "unvisited", before the traversal
     *  starts. Nodes that already have been visited will ignore the 
     *  the accept command.
     */
    void accept(FsmVisitor& v);
    
    /**
     *  Return true if and only if the two FSM states are distinhuishable
     */
    virtual bool distinguishable(const FsmNode& s1, const FsmNode& s2);


	//--------------------------------------------
	const std::vector<std::shared_ptr<OFSMTable>>& calcAndGetOfsmTbls();

	const std::vector<std::shared_ptr<Tree>>& getStateIdentificationSets();

	/**
	 * Test function for InputTrace FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
                                             const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
                                             const int maxInput,
                                             const int maxOutput).
	   Needs access to ofsmTables of this class.
	*/
	friend void testCalcDistinguishingTrace2(Fsm &m);
	//--------------------------------------------

};
#endif //FSM_FSM_FSM_H_
