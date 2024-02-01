#include "shortcutHull.h"

#include <QDebug>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

ShortcutHull::ShortcutHull(QWidget *pCrtanje,
                                         int pauzaKoraka,
                                         const bool &naivni,
                                         std::string imeDatoteke,
                                         int brojTemena)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{
    std::vector<QPoint> temp_points;
    if (imeDatoteke == "")
        temp_points = generisiNasumicneTacke(brojTemena);
    else
        temp_points = ucitajPodatkeIzDatoteke(imeDatoteke);

    for(auto &p : temp_points)
        points.push_back(new Point(p.x(), p.y()));

}

ShortcutHull::~ShortcutHull()
{
    for (auto &p : points) {
        delete p;
    }
}

GeometricGraph::~GeometricGraph() {
    for (auto& pair : edges) {
        for (auto& edge : pair.second) {
            delete edge;
        }
    }
}

Edge::~Edge(){
    delete start;
    delete end;
}

/* Pomocne funkcije */

double calculateAngle(const Point* P1, const Point* Pi) {
    return std::atan2(Pi->y() - P1->y(), Pi->x() - P1->x());
}
double calculateDistance(const Point* a, const Point* b) {
    return sqrt(pow(a->x() - b->x(), 2) + pow(a->y() - b->y(), 2));
}

double Edge::calculateTriangleArea( Point* a,  Point* b,  Point* c) {
    return 0.5 * abs(a->x() * (b->y() - c->y()) + b->x() * (c->y() - a->y()) + c->x() * (a->y() - b->y()));
}

int orientation(Point* p, Point* q, Point* r) {
    int val = (q->y() - p->y()) * (r->x() - q->x()) -
              (q->x() - p->x()) * (r->y() - q->y());

    if (val == 0) return 0;
    return (val > 0) ? 1 : 2;
}


bool onSegment(Point* p, Point* q, Point* r) {
    if (q->x() <= std::max(p->x(), r->x()) && q->x() >= std::min(p->x(), r->x()) &&
        q->y() <= std::max(p->y(), r->y()) && q->y() >= std::min(p->y(), r->y()))
        return true;

    return false;
}

bool edgesIntersect(Point* p1, Point* q1, Point* p2, Point* q2) {

    // Trazimo 4 orijentacije koje su nam potrebne
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // Opsti slucaj
    if (o1 != o2 && o3 != o4)
        return true;

    // Specijalni slucajevi
    // p1, q1 i p2 su kolinearne i p2 lezi na duzi p1q1, slicno za ostale
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

bool GeometricGraph::isPointInsidePolygon(Point* p) {

    if (vertices.size() < 3) return false;

    // Tacka koja je sigurno van dimenzija ekrana
    Point extreme(3000, 3000);

    /* Brojimo preseke duzi p-ekstreme sa stranicama poligona. Ako je
     * broj preseka neparan, tacka je unutra, u suprotnom je spolja.
     * Imamo specijalan slucaj kada je neka od stranica kolinearna
     * sa duzi p-extreme i tada proveravamo da li je p bas na toj stranici */
    int count = 0;
    for (int i = 0; i < vertices.size(); ++i) {
        int next = (i + 1) % vertices.size();

        if (edgesIntersect(vertices[i], vertices[next], p, &extreme)) {
            if (orientation(vertices[i], p, vertices[next]) == 0){
                if(onSegment(vertices[i], p, vertices[next])){
                    return true;
                }
            }
            else
                count++;
        }
    }

    if(count % 2 == 1)
        return true;
    else
        return false;
}


/* Funkcija za racunanje povrsine dzepa. Prolazi se kroz cvorove u dzepu, pri cemu
 * su oni uredjeni kao u polaznom poligonu. Ukupna povrsina je upravo zbir trouglova dobijenih
 * od pocetnog temena i dva susedna tekuca temena */

double Edge::calculateAreaOfPocket(const std::vector<Point*>& pocketVertices) {
    double area = 0.0;
    for (int i = 1; i + 1 < pocketVertices.size(); ++i) {
        area += calculateTriangleArea(pocketVertices[0], pocketVertices[i], pocketVertices[i + 1]);
    }
    return area;
}

/*Funkcija cene za odredjenu precicu. Deljenje kod povrsine postoji zbog
 * vrste poligona na kojima pokrecemo nasumicne primere. Naime, oni su takvi da
 * su povrsine u dzepovima precica cesto velike, pa je povrsina uvek favorizovana.
 * kako bismo videli efekat promene lambda, smanjujemo uticaj povrsine. Na realnim mapama
 * povrsine u dzepovima precica su dosta manje pa ne bi bilo ovog problema */

void Edge::calculateCost(const double lambda) {
    double length = calculateDistance(start, end);
    double area = calculateAreaOfPocket(pocketVertices);

    cost = lambda * length + (1 - lambda) * area/70;
}


void GeometricGraph::addVertex(Point* vertex) {
    vertices.push_back(vertex);
}

/* Funkcija za dodavanje nove ivice. Svakoj ivici dodajemo i listu
 * temena u njenom dzepu */

void GeometricGraph::addEdge( Point* start,  Point* end, bool isShortcut) {
    edges[start].push_back(new Edge(start, end));
    edges[start].back()->isShortcut = isShortcut;
    for (auto& vertex : vertices) {
        if (*vertex >= *start && *vertex <= *end) {
            edges[start].back()->pocketVertices.push_back(vertex);
        }
    }

}

void GeometricGraph::calculateAllCosts(double lambda) {
    for (auto& vertex : vertices) {
        for(auto& edge : edges[vertex])
            edge->calculateCost(lambda);
    }

}

/* Glavna funkcija u algoritmu koja trazi najkraci put u DAG-u. Algoritam je modifikacija
 * Dijkstrinog algoritma koja prolazi kroz sve cvorove i za svakog od njihovih suseda
 * azurira najkraci put ukoliko je moguce (SSSP problem za DAG) */

void GeometricGraph::findShortestPath(std::vector<Point *> *path) {
    std::vector<double> dist(vertices.size(), std::numeric_limits<double>::max());
    std::vector<int> prev(vertices.size(), -1);
    dist[0] = 0;
    for (int i = 0; i < vertices.size(); ++i) {
        for (const auto& edge : edges[vertices[i]]) {
            int endIndex = edge->end->index;
            if (dist[endIndex] > dist[i] + edge->cost) {
                dist[endIndex] = dist[i] + edge->cost;
                prev[endIndex] = i;
            }

        }
    }

    for (int at = vertices.size() - 1; at != -1; at = prev[at]) {
        path->push_back(vertices[at]);

    }
    path->push_back(*path->begin());
    std::reverse(path->begin(), path->end());
}

/* Algoritam za kreiranje prostog poligona. Ideja algoritma je preuzeta iz
 * knjige profesora Janicica, sa malom modifikacijom kako bi poligoni lepse izgledali
 * (Tacka u odnosu na koju se sortira se uzima iz sredine ekrana). */

void ShortcutHull::initVertices() {

    Point *P1 = new Point(_pCrtanje->width()/2, _pCrtanje->height()/2);

    std::sort(points.begin(), points.end(),
              [P1](const Point* a, const Point* b) {
                  double angleA = calculateAngle(P1, a);
                  double angleB = calculateAngle(P1, b);
                  if (angleA == angleB) {
                      return calculateDistance(P1, a) < calculateDistance(P1, b);
                  }
                  return angleA < angleB;
              });


    for (int i = 0; i < points.size(); i++) {
        points[i]->index = i;
        _graph.addVertex(points[i]);
    }

    delete P1;
}


void ShortcutHull::pokreniAlgoritam(){

    initVertices();
    AlgoritamBaza_updateCanvasAndBlock();
    initEdges();
    AlgoritamBaza_updateCanvasAndBlock();
    initShortcuts();
    AlgoritamBaza_updateCanvasAndBlock();
    _graph.calculateAllCosts(lambda);
    AlgoritamBaza_updateCanvasAndBlock();
    _graph.findShortestPath(&shortestPath);
    AlgoritamBaza_updateCanvasAndBlock();

    for (const auto& point : shortestPath) {
        std::cout << "Vertex: (" << point->x() << ", " << point->y() << ")" << std::endl;
   }
}

void ShortcutHull::crtajAlgoritam(QPainter *painter) const
{
    QPen pen(Qt::black, 2);
    painter->setPen(pen);

    for(int i = 0; i < points.size(); i++){
        painter->drawPoint(*points[i]);
    }


    for(auto& vertex : _graph.vertices){
        auto it = _graph.edges.find(vertex);
        if(it != _graph.edges.end()){
            for(const auto& edge : it->second){
                if(!edge->isShortcut)
                    painter->drawLine(*edge);
            }
        }
    }

    QPen pen1(Qt::blue, 2);
    painter->setPen(pen1);

    for(auto& vertex : _graph.vertices){
        auto it = _graph.edges.find(vertex);
        if(it != _graph.edges.end()){
            for(const auto& edge : it->second){
                if(edge->isShortcut)
                    painter->drawLine(*edge);
            }
        }
    }


    QPen pen2(Qt::red, 5);
    painter->setPen(pen2);

    for(int i = 0; i + 1 < shortestPath.size() ; i++){
       painter->drawLine(*shortestPath[i], *shortestPath[i+1]);
    }

}

void ShortcutHull::pokreniNaivniAlgoritam(){}
void ShortcutHull::crtajNaivniAlgoritam(QPainter *) const {}



void ShortcutHull::initEdges(){
    for (size_t i = 0; i < points.size(); ++i) {
        Point *start = points[i];
        Point *end = points[(i + 1) % points.size()];
        _graph.addEdge(start, end, false);

    }
}

void ShortcutHull::initShortcuts(){
    generateAllShortcuts();
}


bool GeometricGraph::edgeIntersectsAnyEdge(Point* start, Point* end) {
    for(auto& vertex : vertices){
        for(auto& edge : edges[vertex]){
            if (edge->start == start || edge->end == end || edge->start == end || edge->end == start)
                continue;
            if (edgesIntersect(edge->start, edge->end, start, end)) {
                return true;
            } else {
                Point *P = new Point((start->x() + end->x())/2, (start->y() + end->y())/2);
                if(isPointInsidePolygon(P)){
                    delete P;
                    return true;
                }
                delete P;
            }
        }
    }
    return false;
}

/* Funkcija koja generise sve moguce precice zarad simulacije.
 * Ona vuce najveci deo slozenosti jer se za svaku mogucu precicu
 * proverava da li sece neku od postojecih ivica. */

void ShortcutHull::generateAllShortcuts() {
    for (int i = 0; i < points.size(); ++i) {
        for (int j = i + 2; j < points.size(); ++j) {
            if (i == 0 && j == points.size() - 1 || i == j) {
                continue;
            }

            Point* start = points[i];
            Point* end = points[j];

            if (!_graph.edgeIntersectsAnyEdge(start, end)) {
                _graph.addEdge(start, end, true);
            }
        }
    }
}



