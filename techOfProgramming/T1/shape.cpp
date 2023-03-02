#include <limits>
#include <stdexcept>

#include "shape.hpp"

void safr::Shape::scale(double k)
{
    if (k <= std::numeric_limits< double >::epsilon())
    {
        throw std::range_error("bad scale");
    }
    scaling(k);
}
