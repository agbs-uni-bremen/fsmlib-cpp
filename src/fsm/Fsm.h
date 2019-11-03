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
#include <map>

#include "fsm/FsmVisitor.h"
#include "fsm/FsmLabel.h"
#include "fsm/InputTrace.h"


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


class too_many_transition_faults : public std::runtime_error
{
public:
    too_many_transition_faults(const std::string& msg);
};

class too_many_output_faults : public std::runtime_error
{
public:
    too_many_output_faults(const std::string& msg);
};

class unexpected_reduction : public std::runtime_error
{
public:
    unexpected_reduction(const std::string& msg);
};

class reduction_not_possible : public std::runtime_error
{
public:
    reduction_not_possible(const std::string& msg);
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
    
    std::shared_ptr<Tree> characterisationSet;
    std::vector<std::shared_ptr<FsmNode>> dReachableStates;
    Minimal minimal;
    bool complete;

    std::vector<std::shared_ptr<OFSMTable>> ofsmTableLst;
    std::vector<std::shared_ptr<Tree>> stateIdentificationSets;
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
    std::shared_ptr<FsmNode> newNode(const int id, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p,
                                     const std::shared_ptr<FsmPresentationLayer>& pl);
    bool contains(const std::deque<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);
    bool contains(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<FsmNode>& n);
    bool contains(const std::vector<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);

    std::shared_ptr<FsmNode> findp(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>& p);

    void parseLine(const std::string & line);
    void readFsm(const std::string & fname);
    
    /**
     *  Parse a line from from a *.fsm file supposed to be in raw format
     *  <source-state-number> <input-number> <output-number> <target-state-number>
     *  The line is checked w.r.t. consisting of exactly 4 numbers. If the line has another format,
     *  false is returned.
     */
    bool checkRawFormat(const std::string& line);
    void parseLineInitial (const std::string & line);
    void readFsmInitial (const std::string & fname);
    /**
     * Reads a dot file and creates an FSM according to the dot definitions.
     * @param fname The path to the dot file.
     */
    void readFsmFromDot (const std::string & fname, const std::string name = "");
    
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
    
    void addRandomTransitions(const float& maxDegreeOfNonDeterminism,
                              const bool& onlyNonDeterministic,
                              const bool& observable,
                              const float& factor,
                              std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>());
    bool meetDegreeOfCompleteness(const float& degreeOfCompleteness,
                                  const float& maxDegreeOfNonDeterminism,
                                  const bool& observable,
                                  std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>());
    bool doesMeetDegreeOfCompleteness(const float& degreeOfCompleteness, std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    void meetNumberOfStates(const int& maxState, const float& maxDegreeOfNonDeterminism, const bool& observable,
                            std::vector<std::shared_ptr<FsmNode>>& createdNodes);
    std::shared_ptr<FsmLabel> createRandomLabel(
            const std::shared_ptr<FsmNode>& srcNode,
            const float& maxDegreeOfNonDeterminism,
            const bool& onlyNonDeterministic,
            const bool& observable) const;

    void selectRandomNodeAndCreateLabel(
            const std::vector<std::shared_ptr<FsmNode>> srcNodePool,
            const float& maxDegreeOfNonDeterminism,
            const bool& onlyNonDeterministic,
            const bool& observable,
            std::shared_ptr<FsmNode>& node,
            std::shared_ptr<FsmLabel>& label) const;

    bool moreTransitionsPossible(const float& maxDegreeOfNonDeterminism,
                                 const bool& onlyNonDeterministic,
                                 std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;

    int getNumberOfPossibleTransitions(std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    int getNumberOfNotDefinedDeterministicTransitions(std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    int getNumberOfNonDeterministicTransitions(std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    int getNumberOfTotalTransitions(std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;

    std::vector<std::shared_ptr<FsmTransition>> getNonDeterministicTransitions() const;

public:
    
    
    /** Copy constructor
     *  A deep copy is created, so that the new Fsm does not
     *  point to any nodes, transitions, or labels of the old FSM
     */
    Fsm(const Fsm& other);

    Fsm(const Fsm& other,
        const std::string& fsmName,
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    /**
     *  Constructor creating an FSM from file - used only internally
     */
    Fsm(const std::string & fname,
        const std::string & fsmName,
        const int maxNodes,
        const int maxInput,
        const int maxOutput,
        const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    /** Destructor */
    virtual ~Fsm();
    
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
                    const std::shared_ptr<FsmPresentationLayer>& presentationLayer,
                    const bool observable = false,
                    const unsigned seed = 0);


    static std::shared_ptr<Fsm>
    createRandomFsm(const std::string & fsmName,
                    const int& maxInput,
                    const int& maxOutput,
                    const int& maxState,
                    const std::shared_ptr<FsmPresentationLayer>& presentationLayer,
                    const float& degreeOfCompleteness,
                    const float& maxDegreeOfNonDeterminism,
                    const bool& forceNonDeterminism,
                    const bool& minimal,
                    const bool& observable,
                    const unsigned& seed = 0);

    /**
     *  Create a mutant of the FSM, producing output faults
     *  and/or transition faults only.
     *
     *  The number of states remains the same. If FSM is completely
     *  specified, the same will hold for the mutant.
     */
    std::shared_ptr<Fsm> createMutant(const std::string & fsmName,
                                      const int numOutputFaults,
                                      const int numTransitionFaults,
                                      const bool keepObservability = false,
                                      const unsigned seed = 0,
                                      const std::shared_ptr<FsmPresentationLayer>& pLayer = nullptr);

    std::shared_ptr<Fsm> createReduction(const std::string& fsmName,
                                         const bool& force,
                                         int& removedTransitions,
                                         const unsigned seed = 0,
                                         const std::shared_ptr<FsmPresentationLayer>& pLayer = nullptr) const;
    
    
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
    std::vector<std::shared_ptr<FsmNode>> getDReachableStates() { return dReachableStates; }
    std::vector<std::shared_ptr<FsmNode>> calcDReachableStates(InputTraceSet& detStateCover);
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
    std::shared_ptr<FsmNode> getNode(int id) const;
    std::shared_ptr<FsmPresentationLayer> getPresentationLayer() const;
    int getInitStateIdx() const;
    void resetColor();
    
    /**
     * Create a dot (GraphViz)-File from this FSM and store it in the
     *  working directory with file name fname.
     */
    void toDot(const std::string & fname);
    
    /**
     * Store the FSM in internal file format consisting of 4-tuples
     *          <source state> <input> <output> <target state>
     * All states, inputs, outputs ar encoded as integers in range 0,1,2,...
     * @param fname Basename  of the file to be created. Extension
     * .fsm will be added if not already used as extension in fname.
     */
    void toInternalFsmFormat(const std::string & fname);
    
    float getDegreeOfCompleteness(const int& minus = 0, std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    float getDegreeOfNonDeterminism(const int& diff = 0, std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;
    int getNumberOfDifferentInputTransitions(std::vector<std::shared_ptr<FsmNode>> nodePool = std::vector<std::shared_ptr<FsmNode>>()) const;

    /**
     Create a new FSM that represents the intersection of this and the other FSM
     @param f the other FSM
     @return a new FSM which equals the intersection of this and f
     */
    Fsm intersect(const Fsm & f, std::string name = "");

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
    Fsm transformToObservableFSM(const std::string& nameSuffix = "_O") const;
    
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
    Fsm minimiseObservableFSM(const std::string& nameSuffix = "_MIN", bool prependFsmName = true);
    
    /**
     Create the minimal observable FSM which is equivalent to this FSM.
     If this FSM is not observable, an observable equivalent is created
     first. The observable FSM is then minimised.
     @return minimal observable FSM which is equivalent to this FSM
     */
    Fsm minimise(const std::string& nameSuffixMin = "_MIN",
                 const std::string& nameSuffixObs = "_0",
                 bool prependFsmName = true);

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


    void addPossibleIOTraces(std::shared_ptr<FsmNode> node,
                             std::shared_ptr<InputOutputTree> tree,
                             IOTraceContainer& iOTraceContainer,
                             const bool cleanTrailingEmptyTraces = true) const;


    void addPossibleIOTraces(std::shared_ptr<FsmNode> node,
                             const IOTreeContainer& treeContainer,
                             IOTraceContainer& iOTraceContainer,
                             const bool cleanTrailingEmptyTraces = true) const;

    bool hasFailure() const;

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
    void bOmega(const IOTreeContainer& adaptiveTestCases,
                const InputTraceSet& inputTraces,
                std::unordered_set<IOTraceContainer>& result) const;

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
                           const IOTraceContainer& vDoublePrime,
                           const bool onlyPlusPortion = false) const;

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
                             const std::vector<std::shared_ptr<FsmNode>>& states,
                             const IOTreeContainer& adaptiveTestCases,
                             std::unordered_set<IOTraceContainer> bOmegaT,
                             const IOTraceContainer& vDoublePrime,
                             const std::vector<std::shared_ptr<FsmNode>>& dReachableStates,
                             const Fsm& spec,
                             const Fsm& iut);

    static bool exceedsBound(const size_t m,
                             const IOTrace& base,
                             const IOTrace& suffix,
                             const std::vector<std::shared_ptr<FsmNode>>& states,
                             const IOTreeContainer& adaptiveTestCases,
                             std::unordered_set<IOTraceContainer> bOmegaT,
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
     * @param observedTraces Return parameter for the trace that caused a failure during the test
     * suite creation.
     * @return `true`, if no failure has been observed during the creation of the test suite,
     * `false`, otherwise.
     */
    static bool adaptiveStateCounting(Fsm& spec, Fsm& iut, const size_t m,
                                      IOTraceContainer& observedTraces,
                                      std::shared_ptr<IOTrace>& failTrace,
                                      int& iterations);

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
     * Determines if the given adaptive test cases distinguish all states from
     * {@code nodesA} from all states from `nodeB`.
     * @param nodesA First set of states
     * @param nodesB Second set of states
     * @param adaptiveTestCases The adaptive test cases that will be used
     * @return `true`, if 'adaptiveTestCases` distuinguishes all states
     * from `nodesA` from all states from `nodesB`; `false`, otherwise
     */
    bool distinguishesAllStates(std::vector<std::shared_ptr<FsmNode>>& nodesA,
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


    bool rDistinguishes(std::shared_ptr<FsmNode> nodeA,
                      std::shared_ptr<FsmNode> nodeB,
                      const IOListContainer& testCases) const;


    /**
     * Determines if the given adaptive test cases distinguish state `nodesA`
     * from state `nodeB`.
     * @param nodeA First state
     * @param nodesB Second state
     * @param adaptiveTestCases The adaptive test cases that will be used
     * @return `true`, if `adaptiveTestCases` distuinguishes state
     * `nodeA` from state `nodeB`; `false`, otherwise
     */
    bool distinguishes(std::shared_ptr<FsmNode> nodeA,
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

    bool rDistinguishes(std::shared_ptr<FsmNode> nodeA,
                        std::shared_ptr<FsmNode> nodeB,
                        const std::vector<int>& list) const;

    /**
     * Determines if the given adaptive test case distinguishes state `nodesA`
     * from state `nodeB`.
     * @param nodeA First state
     * @param nodesB Second state
     * @param adaptiveTestCase The adaptive test case that will be used
     * @return `true`, if `adaptiveTestCase` distuinguishes state
     * `nodeA` from state `nodeB`; `false`, otherwise
     */
    bool distinguishes(std::shared_ptr<FsmNode> nodeA,
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
    
    /**
     *  Return true if and only if the two FSM states are distinhuishable
     */
    virtual bool distinguishable(const FsmNode& s1, const FsmNode& s2);

};
#endif //FSM_FSM_FSM_H_
