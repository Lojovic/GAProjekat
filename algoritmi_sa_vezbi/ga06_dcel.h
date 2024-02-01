#ifndef DCEL_H
#define DCEL_H

#include <QPoint>
#include <QWidget>

class Vertex;
class HalfEdge;
class Field;

///
/// \brief The DCEL class
/// Dvostruko povezna lista (Doubly-Connected Edge List)
/// Vise procitati u knjizi Computational Geometry
///
class DCEL
{
public:
    DCEL(std::string imeDatoteke, int h, int w);

    ///
    /// \brief DCEL
    /// Poligon zadatak tackama u pozitivnom matematickom smeru
    /// bez rupa
    ///
    DCEL(const std::vector<QPointF> &tacke);
    DCEL() = default;
    virtual ~DCEL();

    Vertex *vertex(size_t i) const;
    const std::vector<Vertex *> &vertices() const;
    HalfEdge *edge(size_t i) const;
    const std::vector<HalfEdge *> &edges() const;
    Field *field(size_t i) const;
    const std::vector<Field *> &fields() const;
    void setFields(const std::vector<Field *> &fileds);
    void loadData(const std::vector<QPointF> &tacke);
    void insertEdge(HalfEdge *e);
    void insertFiled(Field *f);
    HalfEdge *findEdge(Vertex *start, Vertex *end);

private:
    std::vector<Vertex *> _vertices;
    std::vector<HalfEdge *> _edges;
    std::vector<Field *> _fields;
};

///
/// \brief The VertexType enum
/// Potrebno za Triangulaciju
/// Sluzi za zadavanje tipa za teme
///
enum class VertexType {START, SPLIT, END, MERGE, REGULAR};

///
/// \brief The Vertex class
/// Teme poligona
/// Sadrzi koordinate
/// i jednu incident (proizvoljnu) HalfEdge
class Vertex{
public:
    Vertex();
    Vertex(float x, float y);
    Vertex(const QPointF &coordinates, HalfEdge *incidentEdge = nullptr);

    qreal x() const;
    qreal y() const;
    const QPointF &coordinates() const;
    void setCoordinates(const QPointF &coordinates);
    HalfEdge* incidentEdge() const;
    void setIncidentEdge(HalfEdge *incidentEdge);

    VertexType type() const { return _type; }
    void setType(VertexType type) { _type = type; }

private:
    QPointF _coordinates;
    HalfEdge *_incidentEdge;

private:
    /* za triangulaciju ga07_triangulation */
    VertexType _type;
};


///
/// \brief The HalfEdge class
/// Ovo je zapravo jedna ivica poligona, ali ivica moze imati dva smera, pa je
/// sa HalfEdge oznacen jedan smer, a ona na koju pokazuje (twin) je ista ivica, ali drugi smer
/// IncidentFace je face (strana) koja je sa leve strane posmatrane HalfEdge
/// Next i Prev su prethodna i sledeca HalfEdge
/// Origin je Vrtex koji je pocetak HalfEdge
/// Vise o ovoj strukturi procitati u Computational Geoemtry
///
class HalfEdge{
public:
    HalfEdge();
    HalfEdge(Vertex *origin, HalfEdge *twin = nullptr,
             HalfEdge *next = nullptr, HalfEdge *prev = nullptr,
             Field *incidentFace = nullptr);

    Vertex *origin() const;
    void setOrigin(Vertex *origin);
    HalfEdge *twin() const;
    void setTwin(HalfEdge *twin);
    HalfEdge *next() const;
    void setNext(HalfEdge *next);
    HalfEdge *prev() const;
    void setPrev(HalfEdge *prev);
    Field *incidentFace() const;
    void setIncidentFace(Field *incidentFace);

private:
    Vertex *_origin;
    HalfEdge *_twin;
    HalfEdge *_next;
    HalfEdge *_prev;
    Field *_incidentFace;
};

///
/// \brief The Field class
/// Sadrzi dva niza HalfEdge
/// Prvi je OuterComponent i pokazuje na one HalfEdge koje su na granici (boundary) posmatrane stranice
/// Za one stranice koje nemaju granicu, tj. koje su nepovezane taj niz je nil (tj. prazan)
/// Drugi niz je InnerComponent i pokazuje na one HalfEdge koje cine "rupu" u okviru posmatrane stranice
/// Moze biti prazan
/// Vise o ovoj strukturi procitati u Computational Geoemtry
///
class Field{
public:
    Field();
    Field(HalfEdge *outerComponent, const std::vector<HalfEdge *> &innerComponents);

    HalfEdge *outerComponent() const;
    void setOuterComponent(HalfEdge *outerComponent);
    const std::vector<HalfEdge *> &innerComponents() const;
    HalfEdge *innerComponent() const;
    void setInnerComponents(const std::vector<HalfEdge *> &innerComponents);
    void setInnerComponent(HalfEdge *innerComponent);

private:
    HalfEdge *_outerComponent;
    std::vector<HalfEdge *> _innerComponents;
};

#endif // DCEL_H
