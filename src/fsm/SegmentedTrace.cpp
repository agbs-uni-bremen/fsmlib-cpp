/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/SegmentedTrace.h"

using namespace std;

TraceSegment::TraceSegment() {
    segment = make_shared< vector<int> >();
    prefix = string::npos;
    tgtNode = nullptr;
}

TraceSegment::TraceSegment(std::shared_ptr< std::vector<int> > segment,
                           size_t prefix,
                           std::shared_ptr<FsmNode> tgtNode)
{
    this->segment = segment;
    this->prefix = prefix;
    this->tgtNode = tgtNode;
}

TraceSegment::TraceSegment(const TraceSegment& other) {
    segment = other.segment;
    prefix = other.prefix;
    tgtNode = other.tgtNode;
}


void TraceSegment::setPrefix(size_t pref) {
    prefix = pref;
}


 
vector<int> TraceSegment::getCopy() {
    
    if ( prefix == string::npos or
         prefix >= segment->size() ) {
        vector<int> v(segment->begin(),segment->end());
        return v;
    }

    if ( prefix == 0 ) {
        vector<int> v;
        return v;
    }
    
    vector<int> w(segment->begin(),segment->begin() + prefix);
    
    return w;
    
}

size_t TraceSegment::size() const {
    return (prefix == string::npos or prefix > segment->size()) ?
             segment->size() : prefix;
}

int TraceSegment::at(size_t n) {
    
    if ( n >= segment->size() ) return -1;
    if ( prefix != string::npos and prefix <= n ) return -1;
    return segment->at(n);
    
}

ostream & operator<<(ostream & out, const TraceSegment& seg)
{
    
    if ( seg.segment->size() == 0 ) {
        out << "eps";
        return out;
    }
    
    out << seg.segment->at(0);
    
    for (size_t i = 1; i < seg.size(); i++ ) {
        out << "." << seg.segment->at(i);
    }
        
    return out;
    
}



// **************************************************************************

SegmentedTrace::SegmentedTrace(std::deque< std::shared_ptr<TraceSegment> > segments) {
    this->segments = segments;
}

SegmentedTrace::SegmentedTrace(const SegmentedTrace& other) {
    segments = other.segments;
}

void SegmentedTrace::add(std::shared_ptr<TraceSegment> seg) {
    segments.push_back(seg);
}

vector<int> SegmentedTrace::getCopy() {
    vector<int> v;
    for ( auto s : segments ) {
        vector<int> svec = s->getCopy();
        v.insert(v.end(),svec.begin(),svec.end());
    }
    return v;
}

shared_ptr<FsmNode> SegmentedTrace::getTgtNode() {
    
    if ( segments.size() == 0 ) return nullptr;
    
    return segments.back()->getTgtNode();
    
}

ostream & operator<<(ostream & out, const SegmentedTrace& trc)
{
    
    if ( trc.segments.size() == 0 ) {
        out << "eps";
        return out;
    }
    
    out << *trc.segments.at(0);
    
    for ( size_t i = 1; i < trc.segments.size(); i++ ) {
        cout << "." << *trc.segments[i];
    }
    return out;
    
}












