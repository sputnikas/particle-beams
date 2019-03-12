#include <iostream>
#include <vector>
#include <cmath>

#include "vec2.h"
#include "expression.h"

////////////////////////////////////////////////////////////////////////////////
// s2d - Space 2d - содержит классы для работы с линиями с 
////////////////////////////////////////////////////////////////////////////////

namespace s2d {

////////////////////////////////////////////////////////////////////////////////
// s2d::Point - точка двухмерного пространства 
// Curve - класс линии, от которого всё наследуется 
//      : простые линии
//      Line - отрезок
//      Arc - дуга окружности
//      Parametric - параметрическая кривая
//      Bezier - кривая Безье
//      Spline - простой сплайн 3-го порядка
//      : сложные линии
//      CompositeCurve - сложная кривая, состоящая из нескольких простых
////////////////////////////////////////////////////////////////////////////////

typedef Vec2<double> Point;

struct Curve {
    int id;
    int style;
    virtual Point point(const double &t) = 0; // return point in t (0 <= t <= 1);
};

struct Line : public Curve {
    int start_point;
    int end_point;

    Line(const int &start_point, const int &end_point);
    Line(const Point &start, const Point &end);
    ~Line();

    Point point(const double &t);
    double measure();
};

struct Arc : public Curve {
    int center;
    double radius;
    double start_alpha;
    double end_alpha;
    int start_point;
    int end_point;
    int geq;

    Arc(const int &center, const double &radius = 1, const double &start_alpha = 0, const double &end_alpha = M_PI);
    Arc(const Point &center, const double &radius = 1, const double &start_alpha = 0, const double &end_alpha = M_PI);
    Arc(const int &center, const int &start, const int &end, const double &radius = 1, const int geq = 1);
    Arc(const Point &center, const Point &start, const Point &end, const double &radius = 1, const int geq = 1);
    ~Arc();

    Point point(const double &t);
    double measure();
};

struct Parametric : public Curve {
    Expression *x_t, *y_t;

    Point point(const double &t);
    double measure();
};

struct CompositeCurve {
    std::vector<Curve*> curves;

    CompositeCurve(Curve *a);
    CompositeCurve(const int &n, Curve *a, ...);
};

}  // namespace s2d

////////////////////////////////////////////////////////////////////////////////
// Sketch - класс эскиза. На нём можно рисовать векторную графику
////////////////////////////////////////////////////////////////////////////////

class Sketch {
    std::vector<s2d::Point> points; // points
    std::vector<int> points_aux_or_not; // every point remember its type - aux or not
    std::vector<std::vector<int>> points_on_curves; // every point know its curve
    std::vector<int> id_curve;  // first - last id, next - free id
    std::vector<int> id_sketch; // first - last id, next - free id
    std::vector<s2d::Curve *> curves;
    int id;
public:
    void AddCurve(s2d::Curve *);
    
};