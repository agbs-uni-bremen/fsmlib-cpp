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
#include <map>

#include "fsm/FsmVisitor.h"
#include "fsm/FsmLabel.h"


class Dfsm;
class FsmNode;
class Tree;
class OutputTree;
class InputTrace;
class FsmPresentationLayer;
class OFSMTable;
class IOListContainer;
class IOTreeContainer;
class TestSuite;
class OutputTrace;
class InputOutputTree;
class IOTrace;
class IOTraceContainer;

enum Minimal
{
    True, False, Maybe
};

/**
 * Modes that can be applied to completely specify a FSM.
 */
enum CompleteMode
{
    /**
      For every missing transition, add a transition with the missing input and output "error"
      to an error state that can never be left.
      */
    ErrorState,
    /**
      For every missing transition, add a transition with the missing input and output "ε" to
      the state that is missing the transition.
      */
    SelfLoop
};

class Fsm
{
protected:
    
    /**
     *  Default constructors without effect - needed by sub-classes
     */
    Fsm(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
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
    
    int failOutput = FsmLabel::UNDEFINED_OUTPUT;
    std::shared_ptr<Tree> characterisationSet;
    Minimal minimal;
    bool complete;

    std::vector<std::shared_ptr<FsmTransition>> transitions;
    std::vector<std::shared_ptr<OFSMTable>> ofsmTableLst;
    std::vector<std::shared_ptr<Tree>> stateIdentificationSets;
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
    std::shared_ptr<FsmNode> newNode(const int id, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);
    bool contains(const std::vector<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);
    bool contains(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<FsmNode>& n);
    std::shared_ptr<FsmNode> findp(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);
    void parseLine(const std::string & line);
    void readFsm(const std::string & fname);
    
    void parseLineInitial (const std::string & line);
    void readFsmInitial (const std::string & fname);
    /**
     * Reads a dot file and creates an FSM according to the dot definitions.
     * @param fname The path to the dot file.
     */
    void readFsmFromDot (const std::string & fname);
    
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
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
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
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer,
        const std::string& fsmName);

    /**
     * Constructs an FSM based on the specification from a dot file.
     * @param dotFileName The path to the dot file.
     * @param fsmName Name of the FSM.
     */
    Fsm(const std::string& dotFileName,
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
        const std::vector<std::shared_ptr<FsmNode>>& lst,
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

    Fsm(const std::string & fsmName,
        const int maxInput,
        const int maxOutput,
        const std::vector<std::shared_ptr<FsmNode>>& lst,
        const int initStateIdx,
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    
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
     *   @param observable  When `true`, an observable FSM will be created.
     *   @return an FSM created at random according to these specifications.
     */
    static std::shared_ptr<Fsm>
    createRandomFsm(const std::string & fsmName,
                    const int maxInput,
                    const int maxOutput,
                    const int maxState,
                    const std::shared_ptr<FsmPresentationLayer>&
                    presentationLayer,
                    const bool observable = false,
                    const unsigned seed = 0);

    /**
     * Creates the product machine for two given FSM.
     * The product machine behaves like `iut` where it is consistent with
     * `reference`. Otherwise it moves to the state `Fail` and stays there.
     * @param reference The reference model
     * @param iut The implementation under test
     * @param fsmName Name of the product machine to be created
     * @return The product machine
     */
    static std::shared_ptr<Fsm> createProductMachine(std::shared_ptr<Fsm> reference, std::shared_ptr<Fsm> iut, const std::string & fsmName = "_prod");

    /**
     * Creates the product machine for two given FSM.
     * The product machine behaves like `iut` where it is consistent with
     * `reference`. Otherwise it moves to the state `Fail` and stays there.
     * @param reference The reference model
     * @param iut The implementation under test
     * @param fsmName Name of the product machine to be created
     * @return The product machine
     */
    static std::shared_ptr<Fsm> createProductMachine(Fsm reference, Fsm iut, const std::string & fsmName = "_prod");
    
    
    /**
     *  Create a mutant of the FSM, producing output faults
     *  and/or transition faults only.
     *
     *  The number of states remains the same. If FSM is completely
     *  specified, the same will hold for the mutant.
     */
    std::shared_ptr<Fsm> createMutant(const std::string & fsmName,
                                      const size_t numOutputFaults,
                                      const size_t numTransitionFaults,
                                      const unsigned seed = 0);
    
    
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
    std::vector<std::shared_ptr<FsmNode>> getDReachableStates(std::vector<std::shared_ptr<InputTrace>>& detStateCover);
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
    
    /**
     Create a new FSM that represents the intersection of this and the other FSM
     @param f the other FSM
     @return a new FSM which equals the intersection of this and f
     */
    Fsm intersect(const Fsm & f);
    
    /**
     * Generate the deterministic state cover.
     * @return A tree representing the deterministic state cover of this FSM.
     */
    std::shared_ptr<Tree> getDeterministicStateCover();

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
     * Calculates each output that can be generated by a given input trace and the corresponding target nodes.
     * @param input The given input trace.
     * @param producedOutputs The calculated outputs that can be produced by the given input trace.
     * @param reachedNodes The calculated nodes that can be reached by the given input trace.
     */
    void apply(const InputTrace& input, std::vector<std::shared_ptr<OutputTrace>>& producedOutputs, std::vector<std::shared_ptr<FsmNode>>& reachedNodes) const;
    
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
    bool isComplete() const;
    
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

    /**
     * Returns a new, completely specified FSM, based on this FSM.
     * @param mode The way this FSM will be completely specified.
     * @return The new completely specified FSM
     */
    Fsm makeComplete(CompleteMode mode);

    /**
     * Returns the output that is associated with "fail" (used for adaptive
     * state counting).
     * @return The output that is associated with "fail".
     */
    int getFailOutput() const;

    bool isCharSet(const std::shared_ptr<Tree>& w) const;
    void minimiseCharSet(const std::shared_ptr<Tree>& w);
    
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
     * Returns all outputs that may occur on <b>both</b> given states when applying the
     * given input.
     * @param q1 First state.
     * @param q2 Second state.
     * @param x The input that will be applied to both states.
     * @return All outputs that may occur on <b>both</b> states with the given input.
     */
    std::vector<std::shared_ptr<OutputTrace>> getOutputIntersection(std::shared_ptr<FsmNode> q1, std::shared_ptr<FsmNode> q2, int x) const;

    /**
     * Calculates for every state the r(1)-distinguishable states.
     */
    void calcROneDistinguishableStates();

    /**
     * Calculates for every state the r-distinguishable states.
     */
    void calcRDistinguishableStates();

    /**
     * Calculates the state characterisation set for a given state, based on the
     * previsously calculated r-distinguishability.
     *
     * The state characterisation set distinguishes all r-distinguishable states
     * from the given state, if the given state is r-distinguishable.
     *
     * @param node The given state
     * @return The state characterisation set for the given state in the form of
     * a list of input sequences.
     */
    IOListContainer getRStateCharacterisationSet(std::shared_ptr<FsmNode> node) const;

    /**
     * Calculates the adaptive state characterisation set for a given state, based on the
     * previsously calculated r-distinguishability.
     *
     * The adaptive state characterisation set distinguishes all r-distinguishable states
     * from the given state, if the given state is r-distinguishable.
     *
     * @param node The given state
     * @return The adaptive state characterisation set for the given state in the form of
     * a list of input sequences.
     */
    IOTreeContainer getAdaptiveRStateCharacterisationSet(std::shared_ptr<FsmNode> node) const;

    /**
     * Calculates the state characterisation set that r-distinguishes all r-distinguishable
     * states.
     * @return The state characterisation set
     */
    IOListContainer getRCharacterisationSet() const;

    /**
     * Calculates the adaptive state characterisation set that r-distinguishes all
     * r-distinguishable states.
     * @return The adaptive state characterisation set
     */
    IOTreeContainer getAdaptiveRCharacterisationSet() const;

    /**
     * Returns all possible input/output sequences that can be produced and all
     * nodes that can be reached when applying a given adaptive test case tree
     * to a given state.
     * @param node The given state
     * @param tree The given adaptive test case
     * @return All input/output traces that can be produced
     */
    IOTraceContainer getPossibleIOTraces(std::shared_ptr<FsmNode> node,
                                         std::shared_ptr<InputOutputTree> tree,
                                         const bool cleanTrailingEmptyTraces = true) const;

    /**
     * Returns all possible input/output sequences that can be produced and all nodes
     * that can be reached when applying each element of a given adaptive test case
     * tree container to a given state.
     * @param node The given state
     * @param treeContainer The given container with adaptive test cases
     * @return All input/output traces that can be produced
     */
    IOTraceContainer getPossibleIOTraces(std::shared_ptr<FsmNode> node,
                                         const IOTreeContainer& treeContainer,
                                         const bool cleanTrailingEmptyTraces = true) const;

    /**
     * Returns a set of input/output sequences that can be produced by this
     * FSM when applying each element of the given adaptive test cases to
     * the state that gets reached by applying the given trace `trace`
     * to this FSM.
     * @param adaptiveTestCases The given adaptive test cases
     * @param trace The given trace
     * @return A set of all input/output traces that can be produced
     */
    IOTraceContainer bOmega(const IOTreeContainer& adaptiveTestCases, const IOTrace& trace) const;

    /**
     * Returns a set of input/output sequences that can be produced by this
     * FSM when applying each element of the given adaptive test cases to
     * each state that can be reached by applying each input trace from the given
     * input traces `inputTraces` to this FSM.
     * @param adaptiveTestCases The given adaptive test cases
     * @param inputTraces The given input traces
     * @return A set of all input/output traces that can be produced
     */
    IOTraceContainer bOmega(const IOTreeContainer& adaptiveTestCases, const std::vector<std::shared_ptr<InputTrace>>& inputTraces) const;

    /**
     * Calculates all possible output traces for the deterministic state cover of the
     * fsm and returns all possible combinations of the corresponding input traces
     * and the related output traces.
     * @return List of possible combinations of input traces and output traces, based
     * on the deterministic state cover.
     */
    std::vector<IOTraceContainer> getVPrime(const std::vector<std::shared_ptr<InputTrace>>& detStateCover);

    /**
     * Calculates all prefixes of `base.suffix`, that reach the given node and that
     * extend `base`.
     * @param node The given node
     * @param base Given base input output trace
     * @param suffix Given suffix input output trace
     * @return List of input output traces
     */
    IOTraceContainer r(std::shared_ptr<FsmNode> node,
                       const IOTrace& base,
                       const IOTrace& suffix) const;

    /**
     * Calculates all prefixes of `base.suffix`, that reach the given node and that
     * extend `base` plus the input output sequence from `vDoublePrime` that reaches
     * the given node (if there is such a sequence).
     * @param node The given node
     * @param base Given base input output trace
     * @param suffix Given suffix input output trace
     * @return List of input output traces
     */
    IOTraceContainer rPlus(std::shared_ptr<FsmNode> node,
                       const IOTrace& prefix,
                       const IOTrace& suffix,
                       const IOTraceContainer& vDoublePrime) const;

    /**
     * Calculates the lower bound that may be placed on the number of states
     * of the FSM if there has been no repetition in the states of the product
     * machine for a given input/output sequence.
     *
     * Let x/y be the given input/output sequence, being observed when applying
     * a determinisitc state cover input sequence.
     * @param base Input/output sequence with the input sequence being an element
     * of the deterministic state cover and the output sequence being a possible
     * output to said input sequence. `base` is a prefix of x/y.
     * @param suffix Suffix x/y with base.suffix = x/y
     * @param takenInputs Set of input sequences that have been followed during
     * adaptive state counting.
     * @param states A maximum set of r-distinguishable states, that has been used
     * during the calculation of `base` and `suffix`.
     * @param adaptiveTestCases The corresponding adaptive test cases.
     * @param vDoublePrime A set of input/output sequences that represents one
     * possible combination of all input traces from the deterministic state cover
     * and theirs corresponding output sequences.
     * @param dReachableStates All d-reachable states.
     * @return The lower bound that may be placed on the number of states
     * of the FSM if there has been no repetition in the states of the product
     * machine for a given input/output sequence.
     */
    static size_t lowerBound(const IOTrace& base,
                             const IOTrace& suffix,
                             const std::vector<std::shared_ptr<InputTrace>>& takenInputs,
                             const std::vector<std::shared_ptr<FsmNode>>& states,
                             const IOTreeContainer& adaptiveTestCases,
                             const IOTraceContainer& vDoublePrime,
                             const std::vector<std::shared_ptr<FsmNode>>& dReachableStates,
                             const Fsm& spec,
                             const Fsm& iut);

    /**
     * Calculates a test suite that determines if a given IUT is a reduction of the given
     * specification..
     *
     * It is assumed that the IUT behaves like some unknown element of a fault domain
     * of completely specified observable FSMs with the same input and output alphabets
     * as the specification.
     * @param spec The given specification
     * @param iut The given IUT
     * @param observedTraces Return parameter for the observed traces during the test
     * suite creation.
     * @return `true`, if no failure has been observed during the creation of the test suite,
     * `false`, otherwise.
     */
    static bool adaptiveStateCounting(Fsm& spec, Fsm& iut, const size_t m, IOTraceContainer& observedTraces);

    /**
     * Determines if the given adaptive test cases distinguish all states from
     * {@code nodesA} from all states from {@code nodeB}.
     * @param nodesA First set of states
     * @param nodesB Second set of states
     * @param adaptiveTestCases The adaptive test cases that will be used
     * @return {@code} true, if {@code adaptiveTestCases} distuinguishes all states
     * from {@code nodesA} from all states from {@code nodesB}; {@code false}, otherwise
     */
    bool rDistinguishesAllStates(std::vector<std::shared_ptr<FsmNode>>& nodesA,
                                std::vector<std::shared_ptr<FsmNode>>& nodesB,
                                const IOTreeContainer& adaptiveTestCases) const;

    /**
     * Determines if the given adaptive test cases distinguish state {@code nodesA}
     * from state {@code nodeB}.
     * @param nodeA First state
     * @param nodesB Second state
     * @param adaptiveTestCases The adaptive test cases that will be used
     * @return {@code} true, if {@code adaptiveTestCases} distuinguishes state
     * {@code nodeA} from state {@code nodeB}; {@code false}, otherwise
     */
    bool rDistinguishes(std::shared_ptr<FsmNode> nodeA,
                      std::shared_ptr<FsmNode> nodeB,
                      const IOTreeContainer& adaptiveTestCases) const;

    /**
     * Determines if the given adaptive test case distinguish state {@code nodesA}
     * from state {@code nodeB}.
     * @param nodeA First state
     * @param nodesB Second state
     * @param adaptiveTestCase The adaptive test case that will be used
     * @return {@code} true, if {@code adaptiveTestCase} distuinguishes state
     * {@code nodeA} from state {@code nodeB}; {@code false}, otherwise
     */
    bool rDistinguishes(std::shared_ptr<FsmNode> nodeA,
                      std::shared_ptr<FsmNode> nodeB,
                      std::shared_ptr<InputOutputTree> adaptiveTestCase) const;

    /**
     * Calculates a set of maximal sets of r-distinguishable states.
     * @return A set of maximal sets of r-distinguishable states
     */
    std::vector<std::vector<std::shared_ptr<FsmNode>>> getMaximalSetsOfRDistinguishableStates() const;
    
    /**
     * Calculate the state identification sets. The sets are stored
     * in list stateIdentificationSets, ordered by the FSM state numbers.
     * \pre The FSM must be observable and minimal, and the characterisation
     * set must have been previously calculated using operation
     * getCharacterisationSet().
     */
    void calcStateIdentificationSets();
    void calcStateIdentificationSetsFast();

    void appendStateIdentificationSets(const std::shared_ptr<Tree>& Wp2) const;
    
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
    
    
    void setPresentationLayer(const std::shared_ptr<FsmPresentationLayer>& ppresentationLayer);
    
    
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

};
#endif //FSM_FSM_FSM_H_
