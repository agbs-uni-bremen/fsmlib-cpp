#ifndef ICLONEABLE_H
#define ICLONEABLE_H

#include<memory>

/**
 * Interface for cloneable objects/classes.
 */
struct ICloneable {
private:
    /**
     * Helper method that returns a raw pointer to the cloned object.
     * @return Raw pointer to the cloned object.
     */
    virtual ICloneable* _clone() const = 0;
public:
    /**
     * Creates a clone and returns it through a shared pointer.
     * @return Shared pointer to the cloned object.
     */
    std::shared_ptr<ICloneable> Clone() const;
    virtual ~ICloneable();

};

#endif // ICLONEABLE_H
