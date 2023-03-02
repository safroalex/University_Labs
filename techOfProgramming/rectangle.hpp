#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <iomanip>
#include <memory>

#include "shape.hpp"
namespace safr
{
    class Rectangle : public Shape
    {
    public:
        Rectangle(std::shared_ptr< double[] > points);
        ~Rectangle();
        double getArea() const;
        rectangle_t getFrameRect() const;
        void move(const point_t& newPoint) override;
        void move(double dX, double dY) override;

    private:
        void scaling(double scaleCoef) override;
        rectangle_t rectangle;
    };
}

#endif
