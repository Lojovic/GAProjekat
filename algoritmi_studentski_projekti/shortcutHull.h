#ifndef SHORTCUTHULL_H
#define SHORTCUTHULL_H

#include "algoritambaza.h"
#include <csetjmp>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <QPoint>
#include <QLine>

struct Point;
struct Edge;
class GeometricGraph;

struct Point : QPoint{
    int index;

    Point(double x, double y) : QPoint(x, y), index(-1) {}

    bool operator<=(const struct Point& other) const {
        return index <= other.index;
    }

    bool operator>=(const Point& other) const {
        return index >= other.index;
    }

};

struct Edge : QLine{
    Point *start, *end;
    double cost;
    std::vector<Point*> pocketVertices;
    bool isShortcut;

    Edge(Point* start, Point* end) : QLine(*start, *end), start(start), end(end), cost(0.0), isShortcut(false) {}

    ~Edge();

    void calculateCost(const double lambda);
    double calculateTriangleArea( Point* a,  Point* b,  Point* c);
    double calculateAreaOfPocket(const std::vector<Point*>& pocketVertices);

};

class GeometricGraph {


public:

    std::unordered_map<Point*, std::vector<Edge*>> edges;
    std::vector<Point*> vertices;

    ~GeometricGraph();

    void addVertex( Point* vertex);

    void removeTopmostVertex();

    void addEdge( Point* start,  Point* end, bool isShortcut);

    void calculateAllCosts(double lambda);

    void findShortestPath(std::vector<Point*> *path);
    bool edgeIntersectsAnyEdge(Point* start, Point* end);
    bool isPointInsidePolygon(Point* p);
};

double distanceBetweenPoints(const Point* a, const Point* b);

class ShortcutHull : public AlgoritamBaza {
public:
    ShortcutHull(QWidget *, int,
                        const bool & = false,
                        std::string = "",
                        int = BROJ_SLUCAJNIH_OBJEKATA);
    virtual ~ShortcutHull() override;

    std::vector<Point*> getGlavni() const;

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *) const final;


    GeometricGraph _graph;

    std::vector<Point*> points;

    /*Ovde podesiti parametar lambda (0, 1)*/

    double lambda = 0.5;

    std::vector<Point*> shortestPath;

private:

    void initVertices();
    void initEdges();
    void initShortcuts();
    void generateAllShortcuts();


#ifndef GA06_BENCHMARK
    /* Staticki bafer za daleki skok */
    static jmp_buf _buf;
#endif
};

#endif // SHORTCUTHULL_H
