#include "ga05_preseciduzi.h"

#include <fstream>

PreseciDuzi::PreseciDuzi(QWidget *pCrtanje,
                         int pauzaKoraka,
                         const bool &naivni,
                         std::string imeDatoteke,
                         int brojDuzi)
   : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni),
     _brisucaPravaY(_pCrtanje ? _pCrtanje->height()-3 : 0),
     _redDuzi(poredjenjeDuzi(&_brisucaPravaY))
{
    if (imeDatoteke != "")
        _duzi = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        _duzi = generisiNasumicneDuzi(brojDuzi);
    _i = _j = _duzi.size();
}

void PreseciDuzi::pokreniAlgoritam()
{
    /* Slozenost ovakvog algoritma: O(nlogn + klogn).
     * Izlazni parametar k je broj preseka. */
    for (auto &duz : _duzi) {
        _redDogadjaja.emplace(duz.p1(), tipDogadjaja::POCETAK_DUZI, &duz, nullptr);
        _redDogadjaja.emplace(duz.p2(), tipDogadjaja::KRAJ_DUZI, &duz, nullptr);
    }

    while (!_redDogadjaja.empty()) {
        auto td = *_redDogadjaja.begin();
        _redDogadjaja.erase(_redDogadjaja.begin());

        _brisucaPravaY = td.tacka.y();
        if (td.tip == tipDogadjaja::POCETAK_DUZI) {
            AlgoritamBaza_updateCanvasAndBlock()

            auto trenutna = _redDuzi.emplace(td.duz1).first;
            if (trenutna != _redDuzi.begin()) {
                auto prethodna= std::prev(trenutna);

                QPointF presek;
                if (pomocneFunkcije::presekDuzi(**trenutna, **prethodna, presek))
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK, *trenutna, *prethodna);
            }

            auto sledeca = std::next(trenutna);
            if (sledeca != _redDuzi.end()) {
                QPointF presek;
                if (pomocneFunkcije::presekDuzi(**trenutna, **sledeca, presek))
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK, *trenutna, *sledeca);
            }
        }
        else if (td.tip == tipDogadjaja::KRAJ_DUZI) {
            AlgoritamBaza_updateCanvasAndBlock()

            auto tr_duz = _redDuzi.find(td.duz1);
            if (tr_duz == _redDuzi.end()) continue;

            auto sledeca = std::next(tr_duz);
            if (tr_duz != _redDuzi.begin() && sledeca != _redDuzi.end()) {
                auto prethodna = std::prev(tr_duz);
                QPointF presek;
                if (pomocneFunkcije::presekDuzi(**prethodna, **sledeca, presek)
                        && presek.y() <= _brisucaPravaY)
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK, *prethodna, *sledeca);
            }

            _redDuzi.erase(tr_duz);
         }
         else /*if (td.tip == tipDogadjaja::PRESEK)*/ {
            _preseci.push_back(td.tacka);
            AlgoritamBaza_updateCanvasAndBlock()

            _redDuzi.erase(td.duz1);
            _redDuzi.erase(td.duz2);

            auto duz1 = _redDuzi.insert(td.duz1).first;
            auto duz2 = _redDuzi.insert(td.duz2).first;

            QPointF presek;

            if (duz1 != _redDuzi.begin()) {
                auto prethodna = std::prev(duz1);
                if (pomocneFunkcije::presekDuzi(**duz1, **prethodna, presek)
                        && presek != td.tacka && presek.y() <= _brisucaPravaY)
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK,*duz1,*prethodna);
            }

            auto sledeca = std::next(duz1);
            if (sledeca != _redDuzi.end()) {
                if(pomocneFunkcije::presekDuzi(**duz1, **sledeca, presek)
                        && presek != td.tacka && presek.y() <= _brisucaPravaY)
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK, *duz1, *sledeca);
            }

            if (duz2 != _redDuzi.begin()) {
                auto prethodna= std::prev(duz2);
                if (pomocneFunkcije::presekDuzi(**duz2, **prethodna, presek)
                        && presek != td.tacka && presek.y() <= _brisucaPravaY)
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK,*duz2,*prethodna);
            }

            sledeca = std::next(duz2);
            if (sledeca != _redDuzi.end()) {
                if (pomocneFunkcije::presekDuzi(**duz2, **sledeca, presek)
                        && presek != td.tacka && presek.y() <= _brisucaPravaY)
                    _redDogadjaja.emplace(presek, tipDogadjaja::PRESEK, *duz2, *sledeca);
            }
        }
    }

    _redDuzi.clear();
    _brisucaPravaY = 3;
    AlgoritamBaza_updateCanvasAndBlock()

    emit animacijaZavrsila();
}

void PreseciDuzi::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen magneta = painter->pen();
    magneta.setColor(Qt::magenta);
    magneta.setWidth(3);

    QPen blue = painter->pen();
    blue.setColor(Qt::blue);
    magneta.setWidth(3);

    QPen yellow = painter->pen();
    yellow.setColor(Qt::yellow);
    yellow.setWidth(10);

    QPen green = painter->pen();
    green.setColor(Qt::darkGreen);
    green.setWidth(10);

    QPen default_pen = painter->pen();

    //Iscrtavamo sve duzi
    for(const QLineF& d : _duzi) {
        painter->setPen(default_pen);
        painter->drawLine(d);
        painter->setPen(green);
        painter->drawPoint(d.p1());
        painter->drawPoint(d.p2());
    }

    //Iscrtavamo trenutnu poziciju brisuce prave
    painter->setPen(blue);
    painter->drawLine(QPointF(0, _brisucaPravaY),
                      QPointF(_pCrtanje->width(), _brisucaPravaY));

    //Iscrtavamo trenutno stanje statusa
    painter->setPen(magneta);
    //int i = 1; //brojac za ispis redosleda duzi u statusu
    for(const QLineF* l : _redDuzi)
    {
        painter->drawLine(*l);
        //painter.drawText(l->center(), QString::number(i++));
    }

    //Isrtavamo sve presecne tacke koje smo nasli
    painter->setPen(yellow);
    for(const QPointF& t : _preseci)
        painter->drawPoint(t);
}

void PreseciDuzi::pokreniNaivniAlgoritam()
{
    /* Slozenost naivnog algoritma: O(n^2). Ona je
     * asimptotski optimalna za najgori slucaj. */
    QPointF presek;
    for (_i = 0; _i < _duzi.size(); _i++) {
        for (_j = _i+1; _j < _duzi.size(); _j++) {
            if (pomocneFunkcije::presekDuzi(_duzi[_i],
                                            _duzi[_j],
                                            presek))
                _naivniPreseci.push_back(presek);
            AlgoritamBaza_updateCanvasAndBlock()
        }
    }
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

void PreseciDuzi::naglasiTrenutnu(QPainter *painter, unsigned long i) const
{
    /* Ako je trenutna duz u redu */
    if (i < _duzi.size()) {
        /* Transformacija cetkice */
        painter->save();
        painter->scale(1, -1);
        painter->translate(0, -2*_duzi[i].center().y()+5);

        /* Podesavanje stila olovke */
        auto olovka = painter->pen();
        olovka.setColor(Qt::darkGreen);
        painter->setPen(olovka);

        /* Oznacavanje trenutne duzi */
        painter->drawText(_duzi[i].center(), "_duzi[" + QString::number(i) + "]");

        /* Ponistavanje transformacija */
        painter->restore();

        /* Podesavanje stila olovke */
        olovka.setColor(Qt::red);
        painter->setPen(olovka);

        /* Iscrtavanje duzi */
        painter->drawLine(_duzi[i]);
    }
}

void PreseciDuzi::crtajNaivniAlgoritam(QPainter *painter) const
{
    /* Odustajanje u slucaju greske */
    if (!painter) return;

    /* Iscrtavanje svih duzi */
    for (const auto &duz : _duzi) {
        painter->drawLine(duz);
    }

    /* Podesavanje stila fonta */
    auto font = painter->font();
    font.setWeight(font.Bold);
    font.setPointSizeF(1.3*font.pointSizeF());
    painter->setFont(font);

    /* Naglasavanje trenutnih duzi */
    naglasiTrenutnu(painter, _i);
    naglasiTrenutnu(painter, _j);

    /* Podesavanje stila olovke */
    auto olovka = painter->pen();
    olovka.setColor(Qt::yellow);
    olovka.setWidth(2*olovka.width());
    painter->setPen(olovka);

    /* Iscrtavanje svih preseka */
    for (const auto &presek: _naivniPreseci) {
        painter->drawPoint(presek);
    }
}

std::vector<QLineF> PreseciDuzi::generisiNasumicneDuzi(int brojDuzi) const
{
    srand(static_cast<unsigned>(time(nullptr)));

    std::vector<QLineF> randomDuzi;

    std::vector<QPoint> tacke_za_duzi = generisiNasumicneTacke(2*brojDuzi);

    for(auto i=0ul; i < tacke_za_duzi.size()-1; i+=2){
        auto x1 = tacke_za_duzi[i].x();
        auto x2 = tacke_za_duzi[i+1].x();
        auto y1 = tacke_za_duzi[i].y();
        auto y2 = tacke_za_duzi[i+1].y();

        if (y1 < y2 || (y1 == y2 && x2 < x1)) {
            auto tmp_y = y1;
            auto tmp_x = x1;
            y1 = y2;
            x1 = x2;
            y2 = tmp_y;
            x2 = tmp_x;
        }
        randomDuzi.emplace_back(x1, y1, x2, y2);
    }

    return randomDuzi;
}

std::vector<QLineF> PreseciDuzi::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QLineF> duzi;

    double x1, y1, x2, y2;

    while(inputFile >> x1 >> y1 >> x2 >> y2)
    {
        if (y1 < y2 || (fabs(y1-y2) < EPS && x2 < x1)) {
            auto tmp_y = y1;
            auto tmp_x = x1;
            y1 = y2;
            x1 = x2;
            y2 = tmp_y;
            x2 = tmp_x;
        }

        duzi.emplace_back(x1, y1, x2, y2);
    }

    return duzi;
}
