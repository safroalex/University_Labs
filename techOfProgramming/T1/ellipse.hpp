#ifndef ELLIPSE_HPP
#define ELLIPSE_HPP

#include <memory>

#include "shape.hpp"

namespace safr
{
    class Ellipse : public Shape
    {
    public:
        Ellipse(std::shared_ptr< double[] > points, size_t ppts);
        ~Ellipse();
        double getArea() const;
        rectangle_t getFrameRect() const;
        void move(const point_t& newPoint) override;
        void move(double dX, double dY);

    private:
        void scaling(double scaleCoef) override;
        rectangle_t rectangle;
    };
}

#endif
