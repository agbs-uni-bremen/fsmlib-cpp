#include "fsm/IOTraceHash.h"

#include "fsm/IOTrace.h"

namespace std {
    size_t hash<IOTrace>::operator()(const IOTrace& trace) const
    {

        size_t seed = 0;
        std::hash<size_t> hasher;

        size_t iH = hash<InputTrace>()(trace.getInputTrace());
        size_t oH = hash<OutputTrace>()(trace.getOutputTrace());

        seed ^= hasher(iH) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(oH) + 0x9e3779b9 + (seed<<6) + (seed>>2);

        return seed;
    }
}

