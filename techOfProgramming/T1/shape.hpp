#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "base-types.hpp"

namespace safr
{
    class Shape
    {
    public:
        virtual ~Shape() = default;
        virtual double getArea() const = 0;
        virtual rectangle_t getFrameRect() const = 0;
        virtual void move(const point_t& newPoint) = 0;
        virtual void move(double newX, double newY) = 0;
        void scale(double scaleCoef);

    protected:
        virtual void scaling(double k) = 0;
    };
}

#endif
