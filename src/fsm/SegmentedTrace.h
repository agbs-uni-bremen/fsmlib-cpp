/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_SEGMENTEDTRACE_H_
#define FSM_FSM_SEGMENTEDTRACE_H_

#include <iostream>
#include <vector>

class TraceSegment {
    
private:
    std::shared_ptr< std::vector<int> > segment;
    size_t prefix;
    
public:
    
    TraceSegment();
    TraceSegment(std::shared_ptr< std::vector<int> > segment,
                 size_t prefix = std::string::npos);
    
    /** Shallow copy */
    TraceSegment(const TraceSegment& other);
    
    void setPrefix(size_t pref);
    
    size_t getPrefix() const { return prefix; }
    
    std::shared_ptr< std::vector<int> > get() { return segment; }
    
    std::vector<int> getCopy();
    
    size_t size();
    
    int at(size_t n);
    
};

class SegmentedTrace
{
private:
    
    std::vector< std::shared_ptr<TraceSegment> > segments;
    
	 
public:
	 
};
#endif  
