#ifndef ICLONEABLE_H
#define ICLONEABLE_H

struct ICloneable {
    virtual ~ICloneable();
    virtual ICloneable* clone() const = 0;
};

#endif // ICLONEABLE_H
