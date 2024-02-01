#include "ga03_konveksniomotac.h"
#include "pomocnefunkcije.h"

#include <algorithm>
#include <QPainterPath>

KonveksniOmotac::KonveksniOmotac(QWidget *pCrtanje,
                                 int pauzaKoraka,
                                 const bool &naivni,
                                 std::string imeDatoteke,
                                 int brojTacaka)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{
    if (imeDatoteke == "")
        _tacke = generisiNasumicneTacke(brojTacaka);
    else
        _tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    _k = _tacke.size();
}

void KonveksniOmotac::pokreniAlgoritam() {
    /* Slozenost ovakvog (Gremovog) algoritma: O(nlogn).
     * Dominira sortiranje, dok su ostali koraci linearni. */
    _maxTacka = _tacke[0];

    for (auto i = 1ul; i < _tacke.size(); i++) {
        if (_tacke[i].x() > _maxTacka.x() || (_tacke[i].x() == _maxTacka.x() && _tacke[i].y() < _maxTacka.y()))
            _maxTacka = _tacke[i];
    }
    AlgoritamBaza_updateCanvasAndBlock()

    std::sort(_tacke.begin(), _tacke.end(), [&](const auto& lhs, const auto& rhs) {
        int P = pomocneFunkcije::povrsinaTrougla(_maxTacka, lhs, rhs);
        return  (P < 0) ||  (P == 0 && pomocneFunkcije::distanceKvadrat(_maxTacka, lhs)
                             < pomocneFunkcije::distanceKvadrat(_maxTacka, rhs));
    });

    _konveksniOmotac.push_back(_maxTacka);
    _konveksniOmotac.push_back(_tacke[1]);
    unsigned pom = 2;
    unsigned j = 2;

    while(j < _tacke.size()) {
        if(pomocneFunkcije::povrsinaTrougla(_konveksniOmotac[pom-2],
                                            _konveksniOmotac[pom-1],
                                            _tacke[j]) < 0)
        {
            _konveksniOmotac.push_back(_tacke[j]);
            ++pom;
            ++j;
        }
        else {
            _konveksniOmotac.pop_back();
            --pom;
            // Ne smemo da povecamo j u ovom slucaju, jer nismo zavrsili sa ovom tackom
        }
        AlgoritamBaza_updateCanvasAndBlock()
    }

    _konveksniOmotac.push_back(_maxTacka);
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

void KonveksniOmotac::crtajAlgoritam(QPainter *painter) const {
    if (!painter) return;

    auto pen = painter->pen();
    pen.setColor(Qt::red);
    painter->setPen(pen);
    for(auto &tacka: _tacke) {
        painter->drawPoint(tacka);
    }

    pen.setColor(Qt::blue);
    painter->setPen(pen);
    for(auto i = 1ul; i < _konveksniOmotac.size(); i++) {
       painter->drawLine(_konveksniOmotac[i-1], _konveksniOmotac[i]);
    }
}

void KonveksniOmotac::pokreniNaivniAlgoritam() {
    /* Slozenost naivnog algoritma: O(n^3).
     * Prolazi se kroz svaki par tacaka. */
    for (_i = 0; _i < _tacke.size(); _i++) {
        for (_j = 0; _j < _tacke.size(); _j++) {
            if (_i == _j) continue;

            /* Proverava se da li su sve povrsine sa
             * trecom tackom negativne, sto znaci da
             * je trojka negativne orijentacije */
            bool svePovrsineNegativne = true;
            for (_k = 0; _k < _tacke.size(); _k++) {
                if (_k == _i || _k == _j) continue;

                _povrsina = pomocneFunkcije::povrsinaTrougla(_tacke[_i],
                                                             _tacke[_j],
                                                             _tacke[_k]);
                AlgoritamBaza_updateCanvasAndBlock()
                if (_povrsina > 0) {
                    svePovrsineNegativne = false;
                    break;
                }
            }

            /* Ako jesu, dodaje se usmerena duz */
            _k = _tacke.size();
            if (svePovrsineNegativne) {
                _naivniOmotac.emplace_back(_tacke[_i], _tacke[_j]);
                AlgoritamBaza_updateCanvasAndBlock()
            }
        }
    }

    /* Obavestavanje pozivaoca o gotovoj animaciji */
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

void KonveksniOmotac::naglasiTrenutno(QPainter *painter, unsigned long i, const char *s) const
{
    /* Transformacija cetkice */
    painter->save();
    painter->scale(1, -1);
    painter->translate(0, -2*_tacke[i].y()-10);

    /* Oznacavanje prvog temena */
    painter->drawText(_tacke[i], s);

    /* Ponistavanje transformacija */
    painter->restore();
}

void KonveksniOmotac::crtajNaivniAlgoritam(QPainter *painter) const
{
    /* Odustajanje u slucaju greske */
    if (!painter) return;

    /* Dohvatanje olovke */
    auto pen = painter->pen();

    /* Ako je algoritam u toku */
    if (_k < _tacke.size()) {
        /* Podesavanje stila olovke */
        if (_povrsina < 0) {
            pen.setColor(Qt::green);
        } else {
            pen.setColor(Qt::yellow);
        }
        painter->setPen(pen);

        /* Crtanje tekuceg trougla */
        QPainterPath put(_tacke[_i]);
        put.lineTo(_tacke[_j]);
        put.lineTo(_tacke[_k]);
        put.lineTo(_tacke[_i]);
        painter->fillPath(put, pen.color());

        /* Crtanje tekuce stranice */
        pen.setColor(Qt::red);
        painter->setPen(pen);
        painter->drawLine(_tacke[_i], _tacke[_j]);

        /* Podesavanje stila fonta */
        auto font = painter->font();
        font.setWeight(font.Bold);
        font.setPointSizeF(1.3*font.pointSizeF());
        painter->setFont(font);

        /* Podesavanje stila olovke */
        auto olovka = painter->pen();
        olovka.setColor(Qt::darkGreen);
        painter->setPen(olovka);

        /* Naglasavanje tekucih temena */
        naglasiTrenutno(painter, _i, "i");
        naglasiTrenutno(painter, _j, "j");
        naglasiTrenutno(painter, _k, "k");
    }

    /* Podesavanje stila olovke */
    pen.setColor(Qt::red);
    painter->setPen(pen);

    /* Crtanje svih tacaka */
    for(auto &tacka: _tacke) {
        painter->drawPoint(tacka);
    }

    /* Podesavanje stila olovke */
    pen.setColor(Qt::blue);
    painter->setPen(pen);

    /* Crtanje konveksnog omotaca */
    for(auto &duz : _naivniOmotac) {
       painter->drawLine(duz);
    }
}

const std::vector<QLine> &KonveksniOmotac::getNaivniOmotac() const
{
    return _naivniOmotac;
}

const std::vector<QPoint> &KonveksniOmotac::getKonveksniOmotac() const
{
    return _konveksniOmotac;
}
