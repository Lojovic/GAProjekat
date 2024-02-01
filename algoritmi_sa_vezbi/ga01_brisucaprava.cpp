#include "ga01_brisucaprava.h"

BrisucaPrava::BrisucaPrava(QWidget *pCrtanje,
                           int pauzaKoraka,
                           const bool &naivni,
                           std::string imeDatoteke,
                           int brojTacaka)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni), _yPoz(0)
{
    if (imeDatoteke == "")
      _tacke = generisiNasumicneTacke(brojTacaka);
    else
      _tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
}

void BrisucaPrava::pokreniAlgoritam()
{
    while (_pCrtanje && _yPoz < _pCrtanje->height())
    {
        _yPoz += 5;
        AlgoritamBaza_updateCanvasAndBlock()
    }

    emit animacijaZavrsila();
}

void BrisucaPrava::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen p = painter->pen();
    p.setColor(Qt::magenta);
    p.setWidth(2);
    p.setCapStyle(Qt::RoundCap);

    painter->setPen(p);
    painter->drawLine(0, _yPoz, _pCrtanje->width(), _yPoz);

    p.setWidth(5);
    for(const QPoint &pt : _tacke)
    {
        if (pt.y() < _yPoz)
            p.setColor(Qt::red);
        else
            p.setColor(Qt::black);

        painter->setPen(p);
        painter->drawPoint(pt);
    }
}

void BrisucaPrava::pokreniNaivniAlgoritam()
{
    emit animacijaZavrsila();
}

void BrisucaPrava::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;
}
