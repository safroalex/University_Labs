#ifndef COMPLEXQUAD_HPP
#define COMPLEXQUAD_HPP

#include <memory>

#include "shape.hpp"

namespace safr
{
    class Complexquad : public Shape
    {
    public:
        Complexquad(std::shared_ptr< double[] > points, size_t ppts);
        ~Complexquad();
        double getArea() const;
        rectangle_t getFrameRect() const;
        void move(const point_t& destination) override;
        void move(double moveX, double moveY) override;

    private:
        void scaling(double scaleCoef) override;
        point_t one, two, three, four;
        point_t getCenterPos() const;
        static double getDistance(point_t first, point_t second);
    };
}
#endif
