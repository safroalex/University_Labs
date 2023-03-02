#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "complexquad.hpp"
#include "ellipse.hpp"
#include "rectangle.hpp"

size_t figureInfo(std::shared_ptr< double[] >& pointsArray)
{
    std::string s = "";
    size_t cntPoint = 0;
    std::getline(std::cin, s, '\n');
    for (size_t i = 1; i < s.size() - 1; i++)
    {
        if (s[i - 1] >= '0' && s[i - 1] <= '9' && s[i] == '.' && s[i + 1] >= 0 && s[i + 1] <= '9')
        {
            cntPoint++;
        }
    }
    pointsArray = std::shared_ptr< double[] >(new double[cntPoint]);

    size_t firstSpace = 1;
    size_t secondSpace = s.find(' ', firstSpace);
    for (size_t i = 0; i < cntPoint; i++)
    {
        pointsArray[i] = strtod(s.substr(firstSpace, secondSpace).c_str(), nullptr);
        firstSpace = secondSpace + 1;
        secondSpace = s.find(' ', firstSpace);
    }
    return cntPoint;
}

std::shared_ptr< std::shared_ptr< safr::Shape >[] >
addFigure(std::shared_ptr< std::shared_ptr< safr::Shape >[] > current, size_t size) {
    std::shared_ptr< std::shared_ptr< safr::Shape >[] > temp(new std::shared_ptr< safr::Shape >[size]);
    for (size_t i = 0; i < size - 1; i++)
    {
        temp[i] = current[i];
    }
    return temp;
}

void write(safr::rectangle_t frame)
{
    double pntX1 = frame.pos.x + frame.width / 2;
    double pntY1 = frame.pos.y + frame.height / 2;
    double pntX2 = frame.pos.x - frame.width / 2;
    double pntY2 = frame.pos.y - frame.height / 2;
    std::cout << pntX2 << ' ' << pntY2 << ' ' << pntX1 << ' ' << pntY1;
}

int main()
{
    bool scaled = false;
    bool error = false;
    std::shared_ptr< std::shared_ptr< safr::Shape >[] > figures;
    std::string s;
    size_t cntFigure = 0;
    while (std::cin >> s)
    {
        try
        {
            std::shared_ptr< double[] > points(new double[0]);
            size_t ppts = 0;
            if (s == "RECTANGLE")
            {
                ppts = figureInfo(points);
                figures = addFigure(figures, cntFigure + 1);
                figures[cntFigure] = std::shared_ptr< safr::Rectangle >(new safr::Rectangle(points));
                cntFigure++;
            }
            else if (s == "ELLIPSE")
            {
                ppts = figureInfo(points);
                figures = addFigure(figures, cntFigure + 1);
                figures[cntFigure] = std::shared_ptr< safr::Ellipse >(new safr::Ellipse(points, ppts));
                cntFigure++;
            }
            else if (s == "COMPLEXQUAD")
            {
                ppts = figureInfo(points);
                figures = addFigure(figures, cntFigure + 1);
                figures[cntFigure] = std::shared_ptr< safr::Complexquad >(new safr::Complexquad(points, ppts));
                cntFigure++;
            }
            else if (s == "SCALE")
            {
                double x, y, coef;
                std::cin >> x >> y >> coef;
                if (!std::cin || coef < 0 || !cntFigure)
                {
                    std::cerr << "invalid scale\n";
                    return 1;
                }
                else
                {
                    scaled = true;
                    std::cout << std::fixed << std::setprecision(1);
                    double area = 0;
                    for (size_t i = 0; i < cntFigure; i++)
                    {
                        area += figures[i]->getArea();
                    }
                    std::cout << area << ' ';
                    for (size_t i = 0; i < cntFigure; i++)
                    {
                        write(figures[i]->getFrameRect());
                        if (i + 1 != cntFigure)
                        {
                            std::cout << ' ';
                        }
                    }
                    std::cout << std::endl;
                    for (size_t i = 0; i < cntFigure; i++)
                    {
                        safr::point_t lastPos = figures[i]->getFrameRect().pos;
                        figures[i]->move(safr::point_t{ x, y });
                        figures[i]->scale(coef);
                        double dx = lastPos.x - x;
                        double dy = lastPos.y - y;
                        figures[i]->move(dx * coef, dy * coef);
                    }
                    double areaTwo = 0;
                    for (size_t i = 0; i < cntFigure; i++)
                    {
                        areaTwo = areaTwo + figures[i]->getArea();
                    }
                    std::cout << areaTwo << ' ';
                    for (size_t i = 0; i < cntFigure; i++)
                    {
                        write(figures[i]->getFrameRect());
                        if (i + 1 != cntFigure)
                        {
                            std::cout << ' ';
                        }
                    }
                    std::cout << std::endl;
                    if (error)
                    {
                        std::cerr << "error in shapes\n";
                    }
                }
                break;
            }
        }
        catch (const std::exception& e)
        {
            error = true;
        }
    }
    if (!scaled)
    {
        std::cerr << "scale?\n";
        return 2;
    }
    return 0;
}
