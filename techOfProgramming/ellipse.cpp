#include <memory>
#include <stdexcept>

#include "ellipse.hpp"

const double m_pi = 3.1415926535;

namespace safr
{
    Ellipse::Ellipse(std::shared_ptr< double[] > points, size_t ppts)
    {
        if (ppts != 4)
        {
            throw std::range_error("not enough points");
        }

        if (points[3] <= 0 || points[2] <= 0)
        {
            throw std::logic_error("incorrect radius\n");
        }
        double x = points[0];
        double y = points[1];
        rectangle.pos = point_t{ x, y };
        rectangle.width = 2.0 * points[3];
        rectangle.height = 2.0 * points[2];
    }
    Ellipse::~Ellipse() {}
    double Ellipse::getArea() const
    {
        return (rectangle.height / 2.0) * (rectangle.width / 2.0) * m_pi;
    }
    void Ellipse::move(double moveX, double moveY)
    {
        rectangle.pos.x += moveX;
        rectangle.pos.y = moveY;
    }
    void Ellipse::move(const point_t& point)
    {
        rectangle.pos = point;
    }
    void Ellipse::scaling(double coef)
    {
        rectangle.height *= coef;
        rectangle.width *= coef;
    }
    rectangle_t Ellipse::getFrameRect() const
    {
        return rectangle;
    }
}
