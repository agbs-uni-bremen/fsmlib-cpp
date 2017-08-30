#ifndef ICLONEABLE_H
#define ICLONEABLE_H

#include<memory>

struct ICloneable {
private:
    virtual ICloneable* clone() const = 0;
public:
    std::shared_ptr<ICloneable> Clone() const;
    virtual ~ICloneable();



};

#endif // ICLONEABLE_H
