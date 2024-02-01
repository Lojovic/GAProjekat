#include "ga07_triangulation.h"
#include "pomocnefunkcije.h"
#include <fstream>

Triangulation::Triangulation(QWidget *pCrtanje,
                         int pauzaKoraka,
                         const bool &naivni,
                         std::string imeDatoteke,
                         int brojTacaka)
   : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni),
     _brisucaPravaY(_pCrtanje ? _pCrtanje->height()-3 : 0),
     _statusQueue(StatusQueueCompTriangulation(&_brisucaPravaY))
{
    std::vector<QPointF> tacke;

    /* Ucitavanje tacka */
    if (imeDatoteke != "")
        tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        tacke = ucitajNasumicneTacke(brojTacaka);

    _monotone = true;
    _polygon.loadData(tacke);
}


void Triangulation::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    /* Podesavanje razlicitih olovki za crtanje */
    QPen regular = painter->pen();
    regular.setWidth(10);

    QPen magneta = painter->pen();
    magneta.setColor(Qt::magenta);
    magneta.setWidth(2);

    QPen blue = painter->pen();
    blue.setColor(Qt::blue);
    magneta.setWidth(2);

    if (_monotone) {
        /* Crtanje pravougaonika */
        for (auto i=0ul; i < _polygon.vertices().size() - 1; i++)
            painter->drawLine(_polygon.vertex(i)->coordinates(),
                              _polygon.vertex(i+1)->coordinates());
        painter->drawLine(_polygon.vertex(_polygon.vertices().size()-1)->coordinates(),
                          _polygon.vertex(0)->coordinates());


        /* Crtanje brisuce prave */
        painter->setPen(blue);
        painter->drawLine(QPointF(0, _brisucaPravaY),
                          QPointF(_pCrtanje->width(), _brisucaPravaY));


        /* Crtanje tacaka u redu dogadaja i njihovih tipova. */
        painter->save();
        painter->resetTransform();

        painter->setPen(regular);
        painter->drawText(600, 500, "regular");

        QPen split = regular;
        split.setColor(Qt::darkRed);
        painter->setPen(split);
        painter->drawText(600, 450, "split");

        QPen merge = regular;
        merge.setColor(Qt::darkCyan);
        painter->setPen(merge);
        painter->drawText(600, 400, "merge");

        QPen start = regular;
        start.setColor(Qt::darkYellow);
        painter->setPen(start);
        painter->drawText(600, 350, "start");

        QPen end = regular;
        end.setColor(Qt::darkGreen);
        painter->setPen(end);
        painter->drawText(600, 300, "end");

        painter->restore();

        for(auto v : _eventQueue) {
            switch(v->type())
            {
                case VertexType::START: painter->setPen(start); break;
                case VertexType::REGULAR: painter->setPen(regular); break;
                case VertexType::SPLIT: painter->setPen(split); break;
                case VertexType::END: painter->setPen(end); break;
                case VertexType::MERGE: painter->setPen(merge); break;
            }

            painter->drawPoint(v->coordinates());
        }

        /* Crtanje dijagonala. */
        painter->setPen(magneta);
        for(auto d : _allDiagonals)
            painter->drawLine(d.first->coordinates(), d.second->coordinates());
    } else {
        regular.setWidth(3);
        painter->setPen(regular);
        for (auto f : _polygon.fields()) {
            HalfEdge* v_start = f->outerComponent();
            HalfEdge* v_next = v_start;

            if (v_start == nullptr) continue;

            do {
                painter->drawLine(v_next->origin()->coordinates(), v_next->twin()->origin()->coordinates());
                v_next = v_next->next();
            } while (v_next != v_start);
        }

        /* Crtanje brisuce prave */
        painter->setPen(blue);
        painter->drawLine(QPointF(0, _brisucaPravaY),
                          QPointF(_pCrtanje->width(), _brisucaPravaY));

        /* Crtanje dijagonala */
        painter->setPen(magneta);
        for(auto d : _allDiagonals)
            painter->drawLine(d.first->coordinates(), d.second->coordinates());
    }
}

void Triangulation::pokreniAlgoritam()
{
    initialiseEventQueue();
    monotonePartition();
//    _monotone = false;
//    connectDiagonalsDCEL();
//    _allDiagonals.clear();
//    AlgoritamBaza_updateCanvasAndBlock();
//    for (auto f : _polygon.fields()) {
//        if (f->outerComponent() == nullptr) continue;
//        triangulacija(f);
//    }
    emit animacijaZavrsila();
}

/**********************************************************************************/
/*                                  Y-MONOTON                                     */
/**********************************************************************************/

void Triangulation::initialiseEventQueue()
{
     for(unsigned i=0;i<_polygon.vertices().size(); i++){
            Vertex * v = _polygon.vertex(i);
            Vertex * v_sledeci = v->incidentEdge()->next()->origin();
            Vertex * v_prethodni = v->incidentEdge()->prev()->origin();

            if (pomocneFunkcije::ispod(v->coordinates(), v_sledeci->coordinates())
                    && pomocneFunkcije::ispod(v->coordinates(), v_prethodni->coordinates())){

                    if (pomocneFunkcije::konveksan(v->coordinates(), v_prethodni->coordinates(), v_sledeci->coordinates())){
                           v->setType(VertexType::START);
                    } else
                           v->setType(VertexType::SPLIT);
            } else if (pomocneFunkcije::ispod(v_sledeci->coordinates(),v->coordinates())
                       && pomocneFunkcije::ispod(v_prethodni->coordinates(), v->coordinates())){
                    if (pomocneFunkcije::konveksan(v->coordinates(), v_prethodni->coordinates(), v_sledeci->coordinates())){
                            v->setType(VertexType::END);
                    } else
                            v->setType(VertexType::MERGE);
            } else v->setType(VertexType::REGULAR);

            _eventQueue.insert(v);
            AlgoritamBaza_updateCanvasAndBlock();
     }
}

void Triangulation::monotonePartition()
{
    while(!_eventQueue.empty()){
        Vertex* novi_dogadjaj = *_eventQueue.begin();
        _eventQueue.erase(_eventQueue.begin());

        _brisucaPravaY = novi_dogadjaj->coordinates().y();
        AlgoritamBaza_updateCanvasAndBlock();

        auto tip = novi_dogadjaj->type();
        if (tip == VertexType::START){
            handleStartVertex(novi_dogadjaj);
        } else if (tip==VertexType::SPLIT){
            handleSplitVertex(novi_dogadjaj);
        } else if (tip==VertexType::END){
            handleEndVertex(novi_dogadjaj);
        } else if (tip==VertexType::MERGE){
            handleMergeVertex(novi_dogadjaj);
        } else if (tip==VertexType::REGULAR){
            handleRegularVertex(novi_dogadjaj);
        } else return;
    }
}

HalfEdge *Triangulation::levaStranica(Vertex *v)
{
    HalfEdge* e=new HalfEdge(v);
    HalfEdge* e_twin = new HalfEdge(v,e);
    e->setTwin(e_twin);

    const auto it = _statusQueue.lower_bound(e);
    return it==_statusQueue.end() ? nullptr : *it;
}

void Triangulation::handleStartVertex(Vertex *v)
{
    _helpers[v->incidentEdge()] = v;
    _statusQueue.insert(v->incidentEdge());
}

void Triangulation::handleEndVertex(Vertex *v)
{
    if (_helpers[v->incidentEdge()->prev()]->type()==VertexType::MERGE){
        _allDiagonals.emplace_back(v, _helpers[v->incidentEdge()->prev()]);
        AlgoritamBaza_updateCanvasAndBlock();
    }
    _statusQueue.erase(v->incidentEdge()->prev());
}

void Triangulation::handleSplitVertex(Vertex *v)
{
    HalfEdge * prva_levo = levaStranica(v);
    if (prva_levo == nullptr) return;

    _allDiagonals.emplace_back(_helpers[prva_levo], v);
    AlgoritamBaza_updateCanvasAndBlock();

    _helpers[prva_levo]=v;
    _statusQueue.emplace(v->incidentEdge());
    _helpers[v->incidentEdge()]=v;
}

void Triangulation::handleMergeVertex(Vertex *v)
{
    auto e = v->incidentEdge()->prev();
    if (_helpers[e]->type()==VertexType::MERGE){
        _allDiagonals.emplace_back(_helpers[e], v);
        AlgoritamBaza_updateCanvasAndBlock();
    }

    _statusQueue.erase(e);

    auto leva_prva = levaStranica(v);
    if (leva_prva == nullptr) return;
    if (_helpers[leva_prva]->type()==VertexType::MERGE){
        _allDiagonals.emplace_back(_helpers[leva_prva],v);
        AlgoritamBaza_updateCanvasAndBlock();
    }

    _helpers[leva_prva]=v;


}

void Triangulation::handleRegularVertex(Vertex *v)
{

}

/**********************************************************************************/
/*                             TRIANGULACIJA                                      */
/**********************************************************************************/

void Triangulation::triangulacija(Field *f)
{
    /* radi vizuelizacije */
    _brisucaPravaY = _pCrtanje->height();
    AlgoritamBaza_updateCanvasAndBlock()


}

bool Triangulation::istiLanac(HalfEdge* e1, HalfEdge* e2)
{
    return true;
}

bool Triangulation::leviLanac(HalfEdge *e1, HalfEdge *e2)
{
    QPointF a = e1->origin()->coordinates();
    QPointF a_next = e1->twin()->origin()->coordinates();

    QPointF b = e2->origin()->coordinates();
    QPointF b_next = e2->twin()->origin()->coordinates();

    if (pomocneFunkcije::ispod(a_next, a) &&
        pomocneFunkcije::ispod(b_next, b))
        return true;
    else
        return false;
}

bool Triangulation::desniLanac(HalfEdge *e1, HalfEdge *e2)
{
    return true;
}


/**********************************************************************************/
/*                          POMOCNE FUNKCIJE                                      */
/**********************************************************************************/

void Triangulation::connectDiagonalsDCEL()
{
    std::map<Vertex*, std::set<HalfEdge*, DiagonalsAddDECELComp>> allDiagonals;
    std::map<Vertex*, std::set<HalfEdge*, DiagonalsAddDECELComp>>::iterator it;


    for(auto& pair : _allDiagonals) {
        HalfEdge* ei1 = new HalfEdge(pair.first);
        HalfEdge* ei2 = new HalfEdge(pair.second, ei1);
        ei1->setTwin(ei2);

        _polygon.insertEdge(ei1);
        _polygon.insertEdge(ei2);

        it = allDiagonals.find(pair.first);
        if (it != allDiagonals.end()) {
            allDiagonals[pair.first].emplace(ei1);
        } else {
            std::set<HalfEdge*, DiagonalsAddDECELComp> s1(DiagonalsAddDECELComp(pair.first));
            s1.emplace(ei1);
            allDiagonals[pair.first] = s1;
        }

        it = allDiagonals.find(pair.second);
        if (it != allDiagonals.end()) {
            allDiagonals[pair.second].emplace(ei2);
        } else {
            std::set<HalfEdge*, DiagonalsAddDECELComp> s2(DiagonalsAddDECELComp(pair.second));
            s2.emplace(ei2);
            allDiagonals[pair.second] = s2;
        }
    }

    Field* f_old = _polygon.field(1);
    for(auto& v : allDiagonals) {
        HalfEdge* v_prev = v.first->incidentEdge()->prev();
        HalfEdge* v_next = v.first->incidentEdge();

        for(auto& d : v.second) {
            v_prev->setNext(d);
            d->setPrev(v_prev);

            v_prev = d->twin();

            /* Dodavanje polja (face) ne radi kako treba.
             * Ponekad pravi greske i za outerComponent za face se zadaje pogresne
             * ivice sto dovodi do greske. */
            if (d->incidentFace() == nullptr) {
                Field* f = new Field();
                _polygon.insertFiled(f);

                if (sameDirectionVectors(f_old->outerComponent(), d)) {
                    f->setOuterComponent(d);
                    f_old->setOuterComponent(d->twin());
                    d->setIncidentFace(f);
                    d->twin()->setIncidentFace(f_old);
                } else {
                    f->setOuterComponent(d->twin());
                    f_old->setOuterComponent(d);
                    d->twin()->setIncidentFace(f);
                    d->setIncidentFace(f_old);
                }

                f_old = f;
            }
        }
        v_prev->setNext(v_next);
        v_next->setPrev(v_prev);
    }
}

bool Triangulation::sameDirectionVectors(HalfEdge *e1, HalfEdge *e2)
{
    Vertex* a = e1->origin();
    Vertex* b = e1->twin()->origin();

    Vertex* c = e2->origin();
    Vertex* d = e2->twin()->origin();

    return  (b->x() - a->x()) *
            (d->x() - c->x()) +
            (b->y() - a->y()) *
            (d->y() - c->y()) > 0;
}

std::vector<QPointF> Triangulation::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QPointF> points;
    int x, y;
    while(inputFile >> x >> y)
        points.emplace_back(x, y);
    return points;
}

std::vector<QPointF> Triangulation::generisiNasumicneTacke(int brojTacaka) const
{
    static int constexpr DRAWING_BORDER = 10;

    srand(static_cast<unsigned>(time(nullptr)));
    int xMax;
    int yMax;

    if (_pCrtanje)
    {
        xMax = _pCrtanje->width() - DRAWING_BORDER;
        yMax = _pCrtanje->height() - DRAWING_BORDER;
    }
    else
    {
        xMax = CANVAS_WIDTH;
        yMax = CANVAS_HEIGHT;
    }

    int xMin = DRAWING_BORDER;
    int yMin = DRAWING_BORDER;

    std::vector<QPointF> randomPoints;

    int xDiff = xMax-xMin;
    int yDiff = yMax-yMin;
    for(int i=0; i < brojTacaka; i++)
        randomPoints.emplace_back(xMin + rand()%xDiff, yMin + rand()%yDiff);

    return randomPoints;
}

std::vector<QPointF> Triangulation::ucitajNasumicneTacke(int brojTacaka) const
{
    std::vector<QPointF> tacke = generisiNasumicneTacke(brojTacaka);

    QPointF maxTacka = tacke[0];

    for (auto i = 1ul; i < tacke.size(); i++) {
        if (tacke[i].x() > maxTacka.x() ||
           (fabs(tacke[i].x() - maxTacka.x()) < EPS && tacke[i].y() < maxTacka.y()))
            maxTacka = tacke[i];
    }

    std::sort(tacke.begin(), tacke.end(), [&](const auto& lhs, const auto& rhs) {
        return pomocneFunkcije::konveksan(maxTacka, lhs, rhs);
    });

    return tacke;
}

/**********************************************************************************/
/*                         NAIVNI ALGORITAM                                       */
/**********************************************************************************/

void Triangulation::pokreniNaivniAlgoritam()
{

}

void Triangulation::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;
}
