/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FsmVisitor_h
#define FsmVisitor_h


class Fsm;
class FsmNode;
class FsmTransition;
class FsmLabel;

class FsmVisitor {
    
protected:
    
    bool isNew;

    
public:
    
    FsmVisitor() { };
    ~FsmVisitor() { } ;
    
    
    /**
     *  @note These visitor methods are not declared
     *        abstract, because the empty methods can 
     *        already be used to mark reachable nodes
     *        (this is done by the accept() methods).
     */
    virtual void visit(Fsm& f) { }
    virtual void visit(FsmNode& n) { }
    virtual void visit(FsmTransition& t) { }
    virtual void visit(FsmLabel& t) { }
    
    
    void setNew(bool b) { isNew = b; }
    bool getNew() const { return isNew; }
    
};


#endif /* FsmVisitor_h */
