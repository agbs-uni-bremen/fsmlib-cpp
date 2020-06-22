#ifndef __GENERALIZED_H_METHOD_HPP__
#define __GENERALIZED_H_METHOD_HPP__

#include "fsm/Dfsm.h"
#include "fsm/InputTrace.h"
#include "fsm/FsmNode.h"
#include "trees/OutputTree.h"
#include <functional>
#include <cassert>

//Base FSM type helper template
//Does not provide any types, making template specialization for the FSM type necessary
//TODO: Rename to adaptor
template<typename FSM>
struct FSM_t {
};

template<>
struct FSM_t<Fsm> {
    struct tag {};
    typedef int InputType;
    typedef std::shared_ptr<FsmNode> StateType;
};

template<>
struct FSM_t<Dfsm> {
    struct tag {};
    typedef int InputType;
    typedef std::shared_ptr<FsmNode> StateType;
};

//template<typename FSM, typename InputTrace, typename Tag>
//bool inputTraceDefinedOnFSM(FSM &&fsm, InputTrace &&trace, Tag);
//template<typename FSM, typename StateType, typename InputTrace, typename Tag>
//bool inputTraceDefinedOnState(FSM &&fsm, StateType &&state, InputTrace &&trace, Tag);
//template<typename FSM, typename StateType, typename InputTrace, typename Tag>
//bool isDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, InputTrace &&trace, Tag);
//template<typename FSM, typename InputTrace, typename StateType, typename Tag>
//InputTrace getDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, Tag);
//template<typename FSM, typename InputTrace, typename StateType, typename Tag>
//StateType getStateAfter(FSM &&fsm, InputTrace &&trace, Tag);
//template<typename FSM, typename TestSuiteType, typename Tag>
//TestSuiteType getStateCover(FSM &&fsm, Tag);
//template<typename FSM, typename IterableInputAlphabet, typename Tag>
//IterableInputAlphabet getInputAlphabet(FSM &&fsm, Tag);
template<typename FSM, typename InputTrace>
bool inputTraceDefinedOnFSM(FSM &&fsm, InputTrace &&trace, struct FSM_t<Fsm>::tag) {
    return inputTraceDefinedOnState(std::forward<FSM>(fsm), fsm.getInitialState(), std::forward<InputTrace>(trace));
}
template<typename FSM, typename StateType, typename InputTrace = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType> >
bool inputTraceDefinedOnState(FSM &&fsm, StateType &&state, InputTrace &&trace, struct FSM_t<Fsm>::tag) {
    //NOTE: Due to the assumption of harmonized traces, all output traces should have the same length.
    return (state->apply(::InputTrace(trace, nullptr), false).getOutputTraces().front().size() == trace.size());
}
template<typename FSM, typename StateType, typename InputTrace = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType> >
bool isDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, InputTrace &&trace, struct FSM_t<Fsm>::tag) {
    return (state1->apply(::InputTrace(trace, nullptr), false) != state2->apply(::InputTrace(trace, nullptr), false));
}
template<typename FSM, typename InputTrace, typename StateType>
InputTrace getDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, struct FSM_t<Fsm>::tag tag) {
    //Assumption: There is a distinguishing sequence
    assert(fsm.distinguishable(*state1, *state2));
    return fsm.calculateDistinguishingTrace(state1, state2).get();
}
template<typename FSM, typename InputTrace, typename StateType>
std::vector<StateType> getStatesAfter(FSM &&fsm, InputTrace &&trace, struct FSM_t<Fsm>::tag) {
    std::vector<std::shared_ptr<OutputTrace>> outputs;
    std::vector<StateType> states;
    fsm.apply(::InputTrace(std::forward<InputTrace>(trace), nullptr), outputs, states);
    return states;
}

//TODO: When having template parameters setting the type for a return parameter, make sure it is a decayed type and a container type
template<typename FSM, typename TestSuiteType>
TestSuiteType getStateCover(FSM &&fsm, struct FSM_t<Fsm>::tag) {
    auto stateCover = fsm.getStateCover();
    auto iolists = stateCover->getIOListsWithPrefixes();
    auto traces = iolists.getIOLists();
    TestSuiteType result;
    std::copy(traces->begin(), traces->end(), std::inserter(result, result.end()));
    return result;
}
template<typename FSM, typename IterableInputAlphabet>
IterableInputAlphabet getInputAlphabet(FSM &&fsm, struct FSM_t<Fsm>::tag) {
    IterableInputAlphabet result;
    for(int symbol = 0; symbol <= fsm.getMaxInput(); ++symbol) {
        result.emplace(result.end(), symbol);
    }
    return result;
}
template<typename FSM>
bool isMinimal(FSM &&fsm, struct FSM_t<Fsm>::tag) {
    return fsm.isMinimal() == True;
}
template<typename FSM>
bool isHarmonized(FSM &&fsm, struct FSM_t<Fsm>::tag) {
    return fsm.isHarmonized();
}



template<typename FSM, typename InputTrace>
bool inputTraceDefinedOnFSM(FSM &&fsm, InputTrace &&trace, struct FSM_t<Dfsm>::tag) {
    return inputTraceDefinedOnState(std::forward<FSM>(fsm), fsm.getInitialState(), std::forward<InputTrace>(trace));
}
template<typename FSM, typename StateType, typename InputTrace = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType> >
bool inputTraceDefinedOnState(FSM &&fsm, StateType &&state, InputTrace &&trace, struct FSM_t<Dfsm>::tag) {
    return (state->apply(::InputTrace(trace, nullptr), false).getOutputTraces().front().size() == trace.size());
}
template<typename FSM, typename StateType, typename InputTrace = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType> >
bool isDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, InputTrace &&trace, struct FSM_t<Dfsm>::tag) {
    return (state1->apply(::InputTrace(trace, nullptr), false) != state2->apply(::InputTrace(trace, nullptr), false));
}
template<typename FSM, typename InputTrace, typename StateType>
InputTrace getDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, struct FSM_t<Dfsm>::tag) {
    //Assumption: There is a distinguishing sequence
    fsm.calculateDistMatrix();
    return *(fsm.getDistTraces(*state1, *state2).front());
}
template<typename FSM, typename InputTrace, typename StateType>
std::vector<StateType> getStatesAfter(FSM &&fsm, InputTrace &&trace, struct FSM_t<Dfsm>::tag) {
    std::vector<std::shared_ptr<OutputTrace>> outputs;
    std::vector<StateType> states;
    fsm.apply(::InputTrace(std::forward<InputTrace>(trace), nullptr), outputs, states);
    return states;
}

//TODO: When having template parameters setting the type for a return parameter, make sure it is a decayed type and a container type
template<typename FSM, typename TestSuiteType>
TestSuiteType getStateCover(FSM &&fsm, struct FSM_t<Dfsm>::tag) {
    auto stateCover = fsm.getStateCover();
    auto iolists = stateCover->getIOListsWithPrefixes();
    auto traces = iolists.getIOLists();
    TestSuiteType result;
    std::copy(traces->begin(), traces->end(), std::inserter(result, result.end()));
    return result;
}
template<typename FSM, typename IterableInputAlphabet>
IterableInputAlphabet getInputAlphabet(FSM &&fsm, struct FSM_t<Dfsm>::tag) {
    IterableInputAlphabet result;
    for(int symbol = 0; symbol <= fsm.getMaxInput(); ++symbol) {
        result.emplace(result.end(), symbol);
    }
    return result;
}
template<typename FSM>
bool isMinimal(FSM &&fsm, struct FSM_t<Dfsm>::tag) {
    return fsm.isMinimal() == True;
}
template<typename FSM>
bool isHarmonized(FSM &&fsm, struct FSM_t<Dfsm>::tag) {
    return true;
}

template<typename FSM, typename InputTrace>
bool inputTraceDefinedOnFSM(FSM &&fsm, InputTrace &&trace) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return inputTraceDefinedOnFSM(std::forward<FSM>(fsm), std::forward<InputTrace>(trace), tag);
}
template<typename FSM, typename StateType, typename InputTrace>
bool inputTraceDefinedOnState(FSM &&fsm, StateType &&state, InputTrace &&trace) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return inputTraceDefinedOnState(std::forward<FSM>(fsm), std::forward<StateType>(state), std::forward<InputTrace>(trace), tag);
}
template<typename FSM, typename StateType, typename InputTrace>
bool isDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2, InputTrace &&trace) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return isDistinguishingSequence(std::forward<FSM>(fsm), std::forward<StateType>(state1), std::forward<StateType>(state2), std::forward<InputTrace>(trace), tag);
}
template<typename FSM, typename StateType = typename FSM_t<typename std::decay<FSM>::type>::StateType, typename InputTrace = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType>>
InputTrace getDistinguishingSequence(FSM &&fsm, StateType &&state1, StateType &&state2) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return getDistinguishingSequence<FSM, InputTrace, StateType>(std::forward<FSM>(fsm), std::forward<StateType>(state1), std::forward<StateType>(state2), tag);
}
template<typename FSM, typename InputTrace, typename StateType = typename FSM_t<typename std::decay<FSM>::type>::StateType>
std::vector<StateType> getStatesAfter(FSM &&fsm, InputTrace &&trace) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return getStatesAfter<FSM, InputTrace, StateType>(std::forward<FSM>(fsm), std::forward<InputTrace>(trace), tag);
}
template<typename FSM, typename TestSuiteType = std::vector<std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType>>>
TestSuiteType getStateCover(FSM &&fsm) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return getStateCover<FSM, TestSuiteType>(std::forward<FSM>(fsm), tag);
}
template<typename FSM, typename IterableInputAlphabet = std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType>>
IterableInputAlphabet getInputAlphabet(FSM &&fsm) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return getInputAlphabet<FSM, IterableInputAlphabet>(std::forward<FSM>(fsm), tag);
}
template<typename FSM>
bool isApplicable(FSM &&fsm) {
    struct FSM_t<typename std::decay<FSM>::type>::tag tag;
    return isMinimal(std::forward<FSM>(fsm), tag)
       and isHarmonized(std::forward<FSM>(fsm), tag);
}



template<typename IIter1, typename IIter2, typename Value = typename std::decay<decltype(*std::declval<IIter1>())>::type>
std::vector<Value> concatenateTraces(IIter1 begin1, IIter1 end1, IIter2 begin2, IIter2 end2) {
    std::vector<Value> result(begin1, end1);
    std::copy(begin2, end2, std::back_inserter(result));
    return result;
}

template<typename Container1,
         typename Container2,
         typename TraceType = typename std::decay<decltype(concatenateTraces(std::declval<Container1>().begin()->begin(), std::declval<Container1>().begin()->end(),
                                                                std::declval<Container2>().begin()->begin(), std::declval<Container2>().begin()->end()))>::type>
std::vector<TraceType> getTraceProduct(Container1 container1, Container2 container2) {
    static_assert(std::is_same<typename std::decay<decltype(*container1.begin())>::type, typename std::decay<decltype(*container2.begin())>::type>::value);
    std::vector<TraceType> result;
    result.reserve(container1.size() * container2.size());
    for(auto &trace1 : container1) {
        for(auto &trace2 : container2) {
            result.push_back(concatenateTraces(trace1.begin(), trace1.end(), trace2.begin(), trace2.end()));
        }
    }
    static_assert(std::is_same<typename std::decay<decltype(*result.begin())>::type, typename std::decay<decltype(*container1.begin())>::type>::value);
    static_assert(std::is_same<typename std::decay<decltype(*result.begin())>::type, typename std::decay<decltype(*container2.begin())>::type>::value);
    return result;
}

template<typename Container,
         typename TraceContainer = typename std::decay<decltype(getTraceProduct(std::declval<Container>(), std::declval<Container>()))>::type>
TraceContainer getTracePower(Container container, unsigned int power) {
    TraceContainer result;
    static_assert(std::is_same<typename std::decay<TraceContainer>::type, typename std::decay<Container>::type>::value);
    result.emplace_back();
    for(unsigned int i = 0; i < power; ++i) {
        result = getTraceProduct(result, container);
        static_assert(std::is_same<typename std::decay<decltype(result)>::type, typename std::decay<Container>::type>::value);
    }
    return result;
}

template<typename IIter,
         typename Value = typename std::decay<decltype(*std::declval<IIter>())>::type>
std::vector<std::vector<Value>> getPrefixes(IIter begin, IIter end) {
    std::vector<std::vector<Value>> result;
    for(IIter intermediate = begin; intermediate != end; ++intermediate) {
        result.emplace_back(begin, intermediate);
    }
    return result;
}

template<typename IterableInputAlphabet, typename SingletonTraceContainer = std::vector<std::vector<typename std::decay<decltype(*std::declval<IterableInputAlphabet>().begin())>::type>>>
SingletonTraceContainer getSingletonTraces(IterableInputAlphabet &&inputAlphabet) {
    SingletonTraceContainer result;
    for(auto const &symbol : inputAlphabet) {
        typename std::decay<decltype(*result.begin())>::type singletonElement;
        auto inserter = std::inserter(singletonElement, singletonElement.end());
        *inserter = symbol;
        result.emplace_back(std::move(singletonElement));
    }
    return result;
}

template<typename IIter1, typename IIter2>
bool isPrefix(IIter1 prefixBegin, IIter1 prefixEnd, IIter2 begin, IIter2 end) {
    return std::search(begin, end, prefixBegin, prefixEnd) == begin;
}
template<typename IIter1, typename IIter2>
bool isSuffix(IIter1 suffixBegin, IIter1 suffixEnd, IIter2 begin, IIter2 end) {
    auto occurencePosition = std::search(begin, end, suffixBegin, suffixEnd);
    return (occurencePosition + std::distance(suffixBegin, suffixEnd)) == end;
}
template<typename IIter1, typename IIter2, typename InsertIterator>
void copyPrefixSequences(IIter1 begin, IIter1 end, IIter2 prefixBegin, IIter2 prefixEnd, InsertIterator &&inserter) {
    std::copy_if(begin, end, std::forward<InsertIterator>(inserter), [&prefixBegin, &prefixEnd](decltype(*begin) const &sequence) -> bool {
        return isPrefix(prefixBegin, prefixEnd, sequence.begin(), sequence.end());
    });
}
template<typename IIter1, typename IIter2, typename InsertIterator>
void copySuffixOfPrefixSequences(IIter1 begin, IIter1 end, IIter2 prefixBegin, IIter2 prefixEnd, InsertIterator &&inserter) {
    std::transform(begin, end, std::forward<InsertIterator>(inserter), [&prefixBegin, &prefixEnd](decltype(*begin) const &sequence) {
        auto occurencePosition = std::search(sequence.begin(), sequence.end(), prefixBegin, prefixEnd);
        if(occurencePosition != sequence.end()) {
            auto copyStartPosition = occurencePosition;
            std::advance(copyStartPosition, std::distance(prefixBegin, prefixEnd));
            typename std::decay<decltype(*std::declval<typename std::decay<decltype(inserter)>::type::container_type>().begin())>::type result;
            std::copy(copyStartPosition, sequence.end(), std::inserter(result, result.end()));
            return result;
        }
        return sequence;
    });
}

template<typename IIter1, typename IIter2, typename IIter3, typename InsertIterator>
void copyCommonSuffixesOfSequencesInTestSuitePrefixedBySequences1And2(IIter1 traceCollectionBegin, IIter1 traceCollectionEnd, IIter2 sequence1Begin, IIter2 sequence1End, IIter3 sequence2Begin, IIter3 sequence2End, InsertIterator &&inserter) {
    std::vector<typename std::decay<decltype(*traceCollectionBegin)>::type> prefixedBySequence1;
    std::vector<typename std::decay<decltype(*traceCollectionBegin)>::type> prefixedBySequence2;
    copyPrefixSequences(traceCollectionBegin, traceCollectionEnd, sequence1Begin, sequence1End, std::back_inserter(prefixedBySequence1));
    copyPrefixSequences(traceCollectionBegin, traceCollectionEnd, sequence2Begin, sequence2End, std::back_inserter(prefixedBySequence2));
    std::vector<typename std::decay<decltype(*traceCollectionBegin)>::type> suffixesOfSequencesPrefixedBySequence1;
    std::vector<typename std::decay<decltype(*traceCollectionBegin)>::type> suffixesOfSequencesPrefixedBySequence2;
    copySuffixOfPrefixSequences(prefixedBySequence1.begin(), prefixedBySequence1.end(), sequence1Begin, sequence1End, std::back_inserter(suffixesOfSequencesPrefixedBySequence1));
    copySuffixOfPrefixSequences(prefixedBySequence2.begin(), prefixedBySequence2.end(), sequence2Begin, sequence2End, std::back_inserter(suffixesOfSequencesPrefixedBySequence2));
    std::sort(suffixesOfSequencesPrefixedBySequence1.begin(), suffixesOfSequencesPrefixedBySequence1.end());
    std::sort(suffixesOfSequencesPrefixedBySequence2.begin(), suffixesOfSequencesPrefixedBySequence2.end());
    std::set_intersection(suffixesOfSequencesPrefixedBySequence1.begin(), suffixesOfSequencesPrefixedBySequence1.end(),
                          suffixesOfSequencesPrefixedBySequence2.begin(), suffixesOfSequencesPrefixedBySequence2.end(),
                          std::forward<InsertIterator>(inserter));
}

template<typename FSM, typename StateType, typename SequenceIter1, typename SequenceIter2, typename TestSuiteIter, typename InsertIterator>
void addDistinguishingSequencesIfNotAlreadyContained(FSM &&fsm, 
                                                     StateType &&state1, StateType &&state2,
                                                     SequenceIter1 prefix1Begin, SequenceIter1 prefix1End,
                                                     SequenceIter2 prefix2Begin, SequenceIter2 prefix2End,
                                                     TestSuiteIter testSuiteBegin, TestSuiteIter testSuiteEnd,
                                                     InsertIterator &&inserter) {
    auto distinguishingSequence = getDistinguishingSequence(std::forward<FSM>(fsm), state1, state2);
    auto trace1 = concatenateTraces(prefix1Begin, prefix1End, distinguishingSequence.begin(), distinguishingSequence.end());
    auto trace2 = concatenateTraces(prefix2Begin, prefix2End, distinguishingSequence.begin(), distinguishingSequence.end());
    bool insertTrace1 = (std::find(testSuiteBegin, testSuiteEnd, trace1) == testSuiteEnd);
    bool insertTrace2 = (std::find(testSuiteBegin, testSuiteEnd, trace2) == testSuiteEnd);
    if(insertTrace1) {
        *inserter = trace1;
    }
    if(insertTrace2) {
        *inserter = trace2;
    }
}



template<typename FSM, typename TestSuiteType = std::vector<std::vector<typename FSM_t<typename std::decay<FSM>::type>::InputType>>>
TestSuiteType generateHMethodTestSuite(FSM &&specification, unsigned int additionalStates) {
    auto stateCover = getStateCover(specification);
    auto inputAlphabet = getSingletonTraces(getInputAlphabet(specification));
    //TODO: Check that each trace in stateCover covers a different state

    std::vector<typename std::decay<decltype(getTracePower(inputAlphabet, 0))>::type> inputAlphabetPowers;
    static_assert(std::is_same<typename std::decay<decltype(inputAlphabet)>::type, typename std::decay<typename decltype(inputAlphabetPowers)::value_type>::type>::value);
    for(unsigned int power = 0; power <= additionalStates+1; ++power) {
        inputAlphabetPowers.push_back(getTracePower(inputAlphabet, power));
    }

    static_assert(std::is_same<typename std::decay<typename decltype(stateCover)::value_type>::type,
                               typename std::decay<typename std::decay<decltype(inputAlphabetPowers.back())>::type::value_type>::type
                              >::value);
    auto testSuiteComplete = getTraceProduct(stateCover, inputAlphabetPowers.back());
    static_assert(std::is_same<typename std::decay<decltype(*stateCover.begin())>::type, typename std::decay<decltype(*testSuiteComplete.begin())>::type>::value);
    static_assert(std::is_same<typename std::decay<decltype(testSuiteComplete)>::type, typename std::decay<decltype(stateCover)>::type>::value);
    decltype(testSuiteComplete) testSuiteDefined;
    static_assert(std::is_same<typename std::decay<decltype(testSuiteDefined)>::type, typename std::decay<decltype(stateCover)>::type>::value);
    std::copy_if(testSuiteComplete.begin(), testSuiteComplete.end(), std::back_inserter(testSuiteDefined),
                    [&specification](decltype(*testSuiteComplete.begin()) const &testCase) {
                        return inputTraceDefinedOnFSM(std::forward<FSM>(specification), testCase);
    });

    for(auto const &stateCoverSequence1 : stateCover) {
        for(auto const &stateCoverSequence2 : stateCover) {
            if(stateCoverSequence1 != stateCoverSequence2) {
                auto const &sequence1 = stateCoverSequence1;
                auto const &sequence2 = stateCoverSequence2;
                //Implicit assumption: state cover finally reaches states exactly once
                auto states1 = getStatesAfter(std::forward<FSM>(specification), stateCoverSequence1);
                auto states2 = getStatesAfter(std::forward<FSM>(specification), stateCoverSequence2);
                for(auto const &state1 : states1) {
                    for(auto const &state2 : states2) {
                        if(state1 != state2) {
                            decltype(stateCover) commonSuffixes;
                            copyCommonSuffixesOfSequencesInTestSuitePrefixedBySequences1And2(testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                             sequence1.begin(), sequence1.end(),
                                                                                             sequence2.begin(), sequence2.end(),
                                                                                             std::inserter(commonSuffixes, commonSuffixes.end()));
                            if(std::none_of(commonSuffixes.begin(), commonSuffixes.end(), [&specification, &state1, &state2](decltype(*commonSuffixes.begin()) const &suffix) -> bool {
                                return isDistinguishingSequence(std::forward<FSM>(specification), state1, state2, suffix);
                            })) {
                                addDistinguishingSequencesIfNotAlreadyContained(std::forward<FSM>(specification), state1, state2,
                                                                                sequence1.begin(), sequence1.end(),
                                                                                sequence2.begin(), sequence2.end(),
                                                                                testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                std::inserter(testSuiteDefined, testSuiteDefined.end()));
                            }
                        }
                    }
                }
            }
        }
    }

    for(auto const &stateCoverSequence1 : stateCover) {
        //NOTE: Implicit (usually justified) assumption, that stateCover sequences are all defined
        //TODO: Encode all implicit assumptions appropriately as assert statements. Preferably using the GSL
        auto baseStates = getStatesAfter(std::forward<FSM>(specification), stateCoverSequence1);
        for(auto const &stateCoverSequence2 : stateCover) {
            for(auto const &power : inputAlphabetPowers) {
                for(auto const &traceInAlphabetPower : power) {
                    if(std::all_of(baseStates.begin(), baseStates.end(), [&specification, &traceInAlphabetPower](decltype(*baseStates.begin()) const &baseState) {
                            return inputTraceDefinedOnState(std::forward<FSM>(specification), baseState, traceInAlphabetPower);
                        })) {
                        auto sequence1 = concatenateTraces(stateCoverSequence1.begin(), stateCoverSequence1.end(), traceInAlphabetPower.begin(), traceInAlphabetPower.end());
                        auto const &sequence2 = stateCoverSequence2;
                        auto states1 = getStatesAfter(std::forward<FSM>(specification), sequence1);
                        auto states2 = getStatesAfter(std::forward<FSM>(specification), sequence2);
                        for(auto const &state1 : states1) {
                            for(auto const &state2 : states2) {
                                if(state1 != state2) {
                                    std::vector<decltype(sequence1)> commonSuffixes;
                                    copyCommonSuffixesOfSequencesInTestSuitePrefixedBySequences1And2(testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                                     sequence1.begin(), sequence1.end(),
                                                                                                     sequence2.begin(), sequence2.end(),
                                                                                                     std::inserter(commonSuffixes, commonSuffixes.end()));
                                    if(std::none_of(commonSuffixes.begin(), commonSuffixes.end(), [&specification, &state1, &state2](typename std::decay<decltype(*commonSuffixes.begin())>::type const &suffix) {
                                        return isDistinguishingSequence(std::forward<FSM>(specification), state1, state2, suffix);
                                    })) {
                                        addDistinguishingSequencesIfNotAlreadyContained(std::forward<FSM>(specification), state1, state2,
                                                                                        sequence1.begin(), sequence1.end(),
                                                                                        sequence2.begin(), sequence2.end(),
                                                                                        testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                        std::inserter(testSuiteDefined, testSuiteDefined.end()));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for(auto const &stateCoverSequence : stateCover) {
        //NOTE: Implicit (usually justified) assumption, that stateCover sequences are all defined
        //TODO: Encode all implicit assumptions appropriately as assert statements. Preferably using the GSL
        auto baseStates = getStatesAfter(std::forward<FSM>(specification), stateCoverSequence);
        for(auto const &power : inputAlphabetPowers) {
            for(auto const &traceInAlphabetPower : power) {
                if(traceInAlphabetPower.size() > 0
                   and std::all_of(baseStates.begin(), baseStates.end(), [&specification, &traceInAlphabetPower](decltype(*baseStates.begin()) const &baseState) {
                            return inputTraceDefinedOnState(std::forward<FSM>(specification), baseState, traceInAlphabetPower);
                        })) {
                    for(auto iter = traceInAlphabetPower.begin() + 1; iter != traceInAlphabetPower.end(); ++iter) {
                        auto sequence1 = concatenateTraces(stateCoverSequence.begin(), stateCoverSequence.end(), traceInAlphabetPower.begin(), traceInAlphabetPower.end());
                        auto sequence2 = concatenateTraces(stateCoverSequence.begin(), stateCoverSequence.end(), traceInAlphabetPower.begin(), iter);
                        auto states1 = getStatesAfter(std::forward<FSM>(specification), sequence1);
                        auto states2 = getStatesAfter(std::forward<FSM>(specification), sequence2);
                        for(auto const &state1 : states1) {
                            for(auto const &state2 : states2) {
                                if(state1 != state2) {
                                    std::vector<decltype(sequence1)> commonSuffixes;
                                    copyCommonSuffixesOfSequencesInTestSuitePrefixedBySequences1And2(testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                                     sequence1.begin(), sequence1.end(),
                                                                                                     sequence2.begin(), sequence2.end(),
                                                                                                     std::inserter(commonSuffixes, commonSuffixes.end()));
                                    if(std::none_of(commonSuffixes.begin(), commonSuffixes.end(), [&specification, &state1, &state2](typename std::decay<decltype(*commonSuffixes.begin())>::type const &suffix) {
                                        return isDistinguishingSequence(std::forward<FSM>(specification), state1, state2, suffix);
                                    })) {
                                        addDistinguishingSequencesIfNotAlreadyContained(std::forward<FSM>(specification), state1, state2,
                                                                                        sequence1.begin(), sequence1.end(),
                                                                                        sequence2.begin(), sequence2.end(),
                                                                                        testSuiteDefined.begin(), testSuiteDefined.end(),
                                                                                        std::inserter(testSuiteDefined, testSuiteDefined.end()));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    decltype(testSuiteDefined) prefixFreeDefinedTestSuite;
    std::copy_if(testSuiteDefined.begin(), testSuiteDefined.end(), std::inserter(prefixFreeDefinedTestSuite, prefixFreeDefinedTestSuite.end()),
        [&testSuiteDefined](decltype(*testSuiteDefined.begin()) const &potentialTestCase) {
        return 1 >= std::count_if(testSuiteDefined.begin(), testSuiteDefined.end(), [&potentialTestCase](decltype(*testSuiteDefined.begin()) const &potentialSupersequence) {
            return isPrefix(potentialTestCase.begin(), potentialTestCase.end(), potentialSupersequence.begin(), potentialSupersequence.end());
        });
    });

    return prefixFreeDefinedTestSuite;
}

#endif //__GENERALIZED_H_METHOD_HPP__

