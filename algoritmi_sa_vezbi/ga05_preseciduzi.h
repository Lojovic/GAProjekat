#ifndef GA05_PRESECIDUZI_H
#define GA05_PRESECIDUZI_H

#include <set>

#include "algoritambaza.h"
#include "pomocnefunkcije.h"

#include "ga05_datastructures.h"

class PreseciDuzi : public AlgoritamBaza
{
public:
    PreseciDuzi(QWidget *pCrtanje,
                int pauzaKoraka,
                const bool &naivni = false,
                std::string imeDatoteke = "",
                int brojDuzi = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

private:
    void naglasiTrenutnu(QPainter *painter, unsigned long i) const;

    std::vector<QLineF> generisiNasumicneDuzi(int brojDuzi = BROJ_SLUCAJNIH_OBJEKATA) const;
    std::vector<QLineF> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;

    std::vector<QLineF> _duzi;
    double _brisucaPravaY;
    std::vector<QPointF> _preseci;

    std::set<tackaDogadjaja, poredjenjeDogadjaja> _redDogadjaja;
    std::set<QLineF*, poredjenjeDuzi> _redDuzi; /* status */

    /* Polja za naivni algoritam */
    unsigned long _i, _j;
    std::vector<QPointF> _naivniPreseci;
};

#endif // GA05_PRESECIDUZI_H
