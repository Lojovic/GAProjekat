#ifndef KONVEKSNIOMOTAC3D_H
#define KONVEKSNIOMOTAC3D_H

#include <unordered_set>
#include <experimental/unordered_set>

#include "algoritambaza.h"
#include "ga04_konveksni3dDatastructures.h"

class KonveksniOmotac3D : public AlgoritamBaza
{
public:
    KonveksniOmotac3D(QWidget *pCrtanje,
                      int pauzaKoraka,
                      const bool &naivni = false,
                      std::string imeDatoteke = "",
                      int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA);
    virtual ~KonveksniOmotac3D() override;

public:
    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

    const std::unordered_set<Ivica *, HashIvica, EqIvica> &getIvice() const;
    const std::unordered_set<Ivica *, HashIvica, EqIvica> &getNaivneIvice() const;


    const std::vector<Stranica*> &getNaivniOmotac3d() const;
    const std::unordered_set<Ivica*, HashIvica, EqIvica> &getKonveksniOmotac3d() const;
private:
    void crtajTeme(Teme* t) const;
    void crtajStranicu(Stranica* s) const;

    /* Pomocne funkcije za generisanje konveksnog omotaca. */
    /* Glavne funkcije. */
    bool Tetraedar();
    void DodajTeme(Teme* t);
    void ObrisiVisak();

    /* Pomocne funkcije. */
    float zapremina6(Stranica* s, Teme* t) const;
    bool kolinearne(Teme* a, Teme* b, Teme* c) const;
    Stranica* napraviStranicu(Teme *i1t1, Teme *i1t2,
                              Teme *i2t1, Teme *i2t2,
                              Teme *st1, Teme *st2, Teme *st3);
    Stranica* napraviDruguStranicu(Ivica* iv, Teme* t);
    Stranica* napraviPrvuStranicu(Ivica* iv, Teme* t);
    /* Ucitavanje podataka. */
    std::vector<Teme*> generisiNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;
    std::vector<Teme*> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;

    std::vector<Teme*> _tacke;

    /* Neuredjeni skup za cuvanje ivica; nije nam sustinski vazan redosled
     * ivica u skupu, pa nema potrebe za nizovima kod kojih su pretraga i
     * brisanje reda O(n) umesto jednostavnog O(1) u proseku */
    std::unordered_set<Ivica*, HashIvica, EqIvica> _ivice;
    std::unordered_set<Ivica*, HashIvica, EqIvica> _naivneIvice;

    Stranica* _tekucaStranica = nullptr;
    //vektor stranica
    std::vector<Stranica*> _naivniOmotac;
    // ovaj vektor je potreban kako bi se memorija uredno oslobodila
    std::vector<Stranica*> _stranice;
};

#endif // KONVEKSNIOMOTAC3D_H
