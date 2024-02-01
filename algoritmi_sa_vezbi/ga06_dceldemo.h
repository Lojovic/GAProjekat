#ifndef DCELDEMO_H
#define DCELDEMO_H

#include "algoritambaza.h"
#include "pomocnefunkcije.h"
#include "ga06_dcel.h"

class DCELDemo : public AlgoritamBaza
{
public:
    DCELDemo(QWidget *pCrtanje,
             int pauzaKoraka,
             const bool &naivni = false,
             std::string imeDatoteke = "",
             int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

private:
    DCEL _polygon;
};

#endif // DCELDEMO_H
