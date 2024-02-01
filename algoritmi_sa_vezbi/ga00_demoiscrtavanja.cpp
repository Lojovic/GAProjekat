#include "ga00_demoiscrtavanja.h"

DemoIscrtavanja::DemoIscrtavanja(QWidget *pCrtanje,
                                 int pauzaKoraka,
                                 const bool &naivni,
                                 std::string imeDatoteke,
                                 int brojTacaka)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni), _n(0)
{
    if (imeDatoteke == "")
        _tacke = generisiNasumicneTacke(brojTacaka);
    else
        _tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
}

void DemoIscrtavanja::pokreniAlgoritam()
{
    while (_n != _tacke.size())
    {
        _n++;
        AlgoritamBaza_updateCanvasAndBlock()
    }

    emit animacijaZavrsila();
}

void DemoIscrtavanja::crtajAlgoritam(QPainter *painter) const
{
    if (!painter) return;

    QPen p = painter->pen();
    p.setColor(Qt::red);
    p.setWidth(2);

    painter->setPen(p);

    painter->drawPolygon(_tacke.data(), static_cast<int>(_n));
}

void DemoIscrtavanja::pokreniNaivniAlgoritam()
{
    emit animacijaZavrsila();
}

void DemoIscrtavanja::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;
}
