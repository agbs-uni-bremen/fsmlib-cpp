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
    
    FsmVisitor() { };
    bool isNew;

    
public:
    
    ~FsmVisitor() { } ;
    
    virtual void visit(Fsm& f) = 0;
    virtual void visit(FsmNode& n) = 0;
    virtual void visit(FsmTransition& t) = 0;
    virtual void visit(FsmLabel& t) = 0;
    
    
    void setNew(bool b) { isNew = b; }
    bool getNew() const { return isNew; }
    
};


#endif /* FsmVisitor_h */
