#ifndef GA07_DATASTRUCTURES_H
#define GA07_DATASTRUCTURES_H

#include "ga06_dcel.h"
#include "pomocnefunkcije.h"

///
/// \brief The EventQueueComp struct
/// Brisuca prava se krece od jedne ka drugoj EventTacki
///
struct EventQueueCompTriangulation
{
    bool operator()(const Vertex* levi, const Vertex* desni) const
    {
        return ((levi->coordinates().y() > desni->coordinates().y())
       || ((fabs(levi->coordinates().y()-desni->coordinates().y())<EPS)
            && (levi->coordinates().x()<desni->coordinates().x())));
    }
};

///
/// \brief The StatusQueueCompTriangulation struct
///
/// Status za Monotono Particionisanje
///
struct StatusQueueCompTriangulation
{
private:
    const double *const yBrisucaPrava;

public:
    StatusQueueCompTriangulation(double *y)
        : yBrisucaPrava(y)
    {}

    bool operator()(const HalfEdge* line1, const HalfEdge* line2) const
    {
       QPointF presek1;
       QPointF presek2;

       QLineF line11(line1->origin()->coordinates(), line1->twin()->origin()->coordinates());
       QLineF line22(line2->origin()->coordinates(), line2->twin()->origin()->coordinates());

       QLineF brisucaPrava(0, *yBrisucaPrava, 100, *yBrisucaPrava);

       pomocneFunkcije::presekDuzi(line11, brisucaPrava, presek1);
       pomocneFunkcije::presekDuzi(line22, brisucaPrava, presek2);

       return presek1.x() < presek2.x();
    }
};

/// TODO
struct DiagonalsAddDECELComp
{
private:
    Vertex* _vertex;

public:
    DiagonalsAddDECELComp(Vertex* vertex = nullptr)
        :_vertex(vertex)
    {}

    bool operator()(const HalfEdge* line1, const HalfEdge* line2) const
    {
        return true;
    }
};

///
/// \brief The EventQueueCompTriangulation2 struct
/// U Event za Triangulaciju pamtimo HalfEdge, a ne vertex
/// (kao sto je dato u knjizi)
/// zbog problema u implementaciji koji postoje ako se koristi vertex
///
/// TODO
struct EventQueueCompTriangulation2
{
    bool operator()(const HalfEdge* lhs, const HalfEdge* rhs) const
    {
        return true;
    }
};


#endif // GA07_DATASTRUCTURES_H
