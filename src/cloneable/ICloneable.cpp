#include "cloneable/ICloneable.h"

std::shared_ptr<ICloneable> ICloneable::Clone() const
{
    return std::shared_ptr<ICloneable>(clone());
}

ICloneable::~ICloneable()
{

}
