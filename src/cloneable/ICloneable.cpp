#include "cloneable/ICloneable.h"

std::shared_ptr<ICloneable> ICloneable::Clone() const
{
    return std::shared_ptr<ICloneable>(_clone());
}

ICloneable::~ICloneable()
{

}
