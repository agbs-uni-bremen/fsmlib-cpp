#ifndef __FSM_LIB_GENERIC_EQUIVALENCE_CLASS_CALCULATION_HPP__
#define __FSM_LIB_GENERIC_EQUIVALENCE_CLASS_CALCULATION_HPP__

template<typename StateEquivalencePartitioningType, typename StateType>
typename std::decay<decltype(std::declval<StateEquivalencePartitioningType>().begin())>::type
mapStateToEquivalenceClass(StateEquivalencePartitioningType &&partitioning, StateType &&state) {
    return std::find_if(partitioning.begin(), partitioning.end(), [&state](typename std::decay<decltype(*partitioning.begin())>::type const &partition) {
        return std::find(partition.begin(), partition.end(), state) != partition.end();
    });
}

template<typename StateEquivalencePartitioningType, typename SymbolContainerType, typename TransitionFunctionType>
bool areStatesEquivalent(StateEquivalencePartitioningType &&prePartitioning,
                         SymbolContainerType &&symbolSet,
                         TransitionFunctionType &&transitionFunction,
                         decltype(*prePartitioning.begin()->begin()) state1,
                         decltype(*prePartitioning.begin()->begin()) state2) {
    typedef typename std::decay<decltype(*symbolSet.begin())>::type SymbolType;
    typedef typename std::decay<decltype(prePartitioning.begin())>::type IteratorType;

    return std::all_of(symbolSet.begin(), symbolSet.end(), [&prePartitioning, &transitionFunction, &state1, &state2](SymbolType const &symbol) {
        auto transitionResultState1 = transitionFunction(state1, symbol);
        auto transitionResultState2 = transitionFunction(state2, symbol);
        std::set<IteratorType> equivalenceClassesState1;
        std::set<IteratorType> equivalenceClassesState2;
        for(auto const &postState : transitionResultState1) {
            equivalenceClassesState1.insert(mapStateToEquivalenceClass(prePartitioning, postState));
        }
        for(auto const &postState : transitionResultState2) {
            equivalenceClassesState2.insert(mapStateToEquivalenceClass(prePartitioning, postState));
        }
        return equivalenceClassesState1 == equivalenceClassesState2;
    });
}

template<typename StateEquivalencePartitioningType, typename SymbolContainerType, typename TransitionFunctionType>
typename std::decay<StateEquivalencePartitioningType>::type
refineStatePartitioning(StateEquivalencePartitioningType &&prePartitioning,
                        SymbolContainerType &&symbolSet,
                        TransitionFunctionType &&transitionFunction) {
    typedef typename std::decay<StateEquivalencePartitioningType>::type PartitioningType;
    typedef typename std::decay<decltype(*std::declval<StateEquivalencePartitioningType>().begin())>::type PartitionType;
    PartitioningType subpartitioning;
    for(auto const &partition : prePartitioning) {
        auto inserter = std::inserter(subpartitioning, subpartitioning.end());
        for(auto const &state : partition) {
            auto matchingPartition = std::find_if(subpartitioning.begin(), subpartitioning.end(), [&prePartitioning, &symbolSet, &transitionFunction, &state](PartitionType const &subpartition) {
                auto const &firstElement = *subpartition.begin();
                return areStatesEquivalent(prePartitioning, symbolSet, transitionFunction, firstElement, state);
            });
            if(matchingPartition == subpartitioning.end()) {
                PartitionType newSubpartition;
                auto subPartitionInserter = std::inserter(newSubpartition, newSubpartition.end());
                *subPartitionInserter = state;
                *inserter = newSubpartition;
            } else {
                auto subPartitionInserter = std::inserter(*matchingPartition, matchingPartition->end());
                *subPartitionInserter = state;
            }
        }
    }

    return subpartitioning;
}

template<typename StateCollectionType, typename PredicateFunction>
std::vector<typename std::decay<StateCollectionType>::type> createInitialPartitioning(StateCollectionType &&stateSet, PredicateFunction &&predicateFunction) {
    typedef typename std::decay<StateCollectionType>::type DecayedStateCollectionType;
    std::map<typename std::decay<decltype(predicateFunction(*stateSet.begin()))>::type, DecayedStateCollectionType> mapEquivalentStatesToSameCollection;
    for(auto const &state : stateSet) {
        auto &collection = mapEquivalentStatesToSameCollection[predicateFunction(state)];
        auto inserter = std::inserter(collection, collection.end());
        *inserter = state;
    }

    std::vector<DecayedStateCollectionType> result;
    for(auto const &kvp : mapEquivalentStatesToSameCollection) {
        result.emplace_back(kvp.second);
    }
    return result;
}


template<typename Data, typename Fn, typename... Args>
typename std::decay<Data>::type fixpoint(Fn &&function, Data &&initialData, Args&& ...additionalArgs) {
    typedef typename std::decay<Data>::type DataType;
    DataType currentData;
    DataType nextData = initialData;
    do {
        currentData = nextData;
        nextData = function(currentData, std::forward<Args>(additionalArgs)...);
    } while(nextData != currentData);
    return currentData;
}

#endif //__FSM_LIB_GENERIC_EQUIVALENCE_CLASS_CALCULATION_HPP__

