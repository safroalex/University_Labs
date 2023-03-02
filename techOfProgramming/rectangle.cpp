#include <cfloat>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#include "rectangle.hpp"

safr::Rectangle::Rectangle(std::shared_ptr< double[] > points)
{
    rectangle.height = points[3] - points[1];
    rectangle.width = points[2] - points[0];
    double x = points[0] + (rectangle.width) / 2;
    double y = points[1] + (rectangle.height) / 2;
    rectangle.pos = point_t{ x, y };
    if (rectangle.width <= DBL_EPSILON || rectangle.height <= DBL_EPSILON)


    {
        throw std::range_error("incorrect input");


    }
}

double safr::Rectangle::getArea() const
{
    return rectangle.height * rectangle.width;
}

safr::rectangle_t safr::Rectangle::getFrameRect() const
{
    return rectangle;
}

void safr::Rectangle::move(const point_t& newPoint)
{
    rectangle.pos = newPoint;
}

void safr::Rectangle::move(double dX, double dY)
{
    rectangle.pos.x += dX;
    rectangle.pos.y += dY;
}

void safr::Rectangle::scaling(double scaleCoef)
{
    rectangle.width *= scaleCoef;
    rectangle.height *= scaleCoef;
}

safr::Rectangle::~Rectangle() {}
