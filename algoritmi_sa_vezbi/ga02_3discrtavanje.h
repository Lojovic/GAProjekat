#ifndef DISCRTAVANJE_H
#define DISCRTAVANJE_H

#include "algoritambaza.h"

class Discrtavanje : public AlgoritamBaza
{
public:
    Discrtavanje(QWidget *pCrtanje,
                 int pauzaKoraka,
                 const bool &naivni = false,
                 std::string imeDatoteke = "",
                 int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

private:
    std::vector<QVector3D> generisiNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;
    std::vector<QVector3D> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;

    std::vector<QVector3D> _tacke;
    unsigned int _n;
};

#endif // DISCRTAVANJE_H
