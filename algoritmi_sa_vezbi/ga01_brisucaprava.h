#ifndef BRISUCAPRAVA_H
#define BRISUCAPRAVA_H

#include "algoritambaza.h"

class BrisucaPrava : public AlgoritamBaza
{
public:
    BrisucaPrava(QWidget *pCrtanje,
                 int pauzaKoraka,
                 const bool &naivni = false,
                 std::string imeDatoteke = "",
                 int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

private:
    std::vector<QPoint> _tacke;
    int _yPoz;
};

#endif // BRISUCAPRAVA_H
