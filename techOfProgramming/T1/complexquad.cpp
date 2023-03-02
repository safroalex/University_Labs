#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "complexquad.hpp"

safr::Complexquad::Complexquad(std::shared_ptr< double[] > points, size_t ppts) :
    one({ points[0], points[1] }),
    two({ points[2], points[3] }),
    three({ points[4], points[5] }),
    four({ points[6], points[7] })
{
    if (ppts != 8)
    {
        throw std::range_error("not enough points");
    }

    if ((one.x == two.x && one.y == two.y) || (one.x == three.x && one.y == three.y) ||
        (two.x == three.x && two.y == three.y) || (one.x == four.x && one.y == four.y) ||
        (two.x == four.x && two.y == four.y) || (three.x == four.x && three.y == four.y))
    {
        throw std::logic_error("Identical points");
    }

    if ((three.x - one.x) / (two.x - one.x) == (three.y - one.y) / (two.y - one.y) ||
        (four.x - two.x) / (three.x - two.x) == (four.y - two.y) / (three.y - two.y) ||
        (four.x - one.x) / (two.x - one.x) == (four.y - one.y) / (two.y - one.y) ||
        (four.x - one.x) / (three.x - one.x) == (four.y - one.y) / (three.y - one.y))
    {
        throw std::logic_error("At least 3 points lying on the same line\n");
    }
}

safr::Complexquad::~Complexquad() {}

safr::point_t safr::Complexquad::getCenterPos() const
{
    double x = (((one.x * two.y - one.y * two.x) * (three.x - four.x) -
        (one.x - two.x) * (three.x * four.y - three.y * four.x)) /
        ((one.x - two.x) * (three.y - four.y) - (one.y - two.y) * (three.x - four.x)));
    double y = (((one.x * two.y - one.y * two.x) * (three.y - four.y) -
        (one.y - two.y) * (three.x * four.y - three.y * four.x)) /
        ((one.x - two.x) * (three.y - four.y) - (one.y - two.y) * (three.x - four.x)));
    return point_t{ x, y };
}

double safr::Complexquad::getArea() const
{
    double centerSide1 = getDistance(one, getCenterPos());
    double centerSide2 = getDistance(two, getCenterPos());
    double centerSide3 = getDistance(three, getCenterPos());
    double centerSide4 = getDistance(four, getCenterPos());
    double side14 = getDistance(one, four);
    double side23 = getDistance(two, three);
    double p1 = 0.5 * (side14 + centerSide1 + centerSide4);
    double p2 = 0.5 * (side23 + centerSide2 + centerSide3);
    return std::sqrt(p1 * (p1 - side14) * (p1 - centerSide1) * (p1 - centerSide4)) +
        std::sqrt(p2 * (p2 - side23) * (p2 - centerSide2) * (p2 - centerSide3));
}

safr::rectangle_t safr::Complexquad::getFrameRect() const
{
    double minX = one.x;
    double maxX = one.x;
    double minY = one.y;
    double maxY = one.y;
    if (two.x < minX)
    {
        minX = two.x;
    }
    if (three.x < minX)
    {
        minX = three.x;
    }
    if (four.x < minX)
    {
        minX = four.x;
    }
    if (two.y < minY)
    {
        minY = two.y;
    }
    if (three.y < minY)
    {
        minY = three.y;
    }
    if (four.y < minY)
    {
        minY = four.y;
    }
    if (two.x > maxX)
    {
        maxX = two.x;
    }
    if (three.x > maxX)
    {
        maxX = three.x;
    }
    if (four.x > maxX)
    {
        maxX = four.x;
    }
    if (two.y > maxY)
    {
        maxY = two.y;
    }
    if (three.y > maxY)
    {
        maxY = three.y;
    }
    if (four.y > maxY)
    {
        maxY = four.y;
    }
    double width = maxX - minX;
    double height = maxY - minY;
    double x = 0.5 * (minX + maxX);
    double y = 0.5 * (minY + maxY);
    return { width, height, point_t{ x, y } };
}

void safr::Complexquad::move(const point_t &destination)
{
    double moveX = destination.x - getFrameRect().pos.x;
    double moveY = destination.y - getFrameRect().pos.y;
    one.x += moveX;
    one.y += moveY;
    two.x += moveX;
    two.y += moveY;
    three.x += moveX;
    three.y += moveY;
    four.x += moveX;
    four.y += moveY;
}

void safr::Complexquad::move(double moveX, double moveY)
{
    one.x += moveX;
    one.y += moveY;
    two.x += moveX;
    two.y += moveY;
    three.x += moveX;
    three.y += moveY;
    four.x += moveX;
    four.y += moveY;
}

void safr::Complexquad::scaling(double scaleCoef)
{
    one.x *= scaleCoef;
    one.y *= scaleCoef;
    two.x *= scaleCoef;
    two.y *= scaleCoef;
    three.x *= scaleCoef;
    three.y *= scaleCoef;
    four.x *= scaleCoef;
    four.y *= scaleCoef;
}

double safr::Complexquad::getDistance(safr::point_t first, safr::point_t second)
{
    return std::sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}
