#ifndef GA06_PRESEKPRAVOUGAONIKA_H
#define GA06_PRESEKPRAVOUGAONIKA_H

#include "algoritambaza.h"

#include <csetjmp>
#include <experimental/optional>
#include <ext/pb_ds/assoc_container.hpp>

/* Makro za proveru performansi algoritma */
/*#define GA06_BENCHMARK*/

/* Makro za osvezavanje crteza uz eventualni daleki
 * skok za potrebe iskakanja iz duboke rekurzije */
#define PresekPravougaonika_updateCanvasAndBlock() \
    if (updateCanvasAndBlock()) \
    { \
        longjmp(_buf, 1); \
    }

/* Enumeracija pripadnosti skupu kandidata */
enum class KandidatS {Sxx, S11, S12, S21, S22};

/* Struktura koja predstavlja pravougaonik */
struct Pravougaonik : QRect {
    /* Konstruktori i destruktor strukture */
    Pravougaonik(QPoint &, QPoint &);
    Pravougaonik(int, int, int, int);
    virtual ~Pravougaonik() = default;

    /* Cuvanje intervala za brisucu pravu */
    const int xLevo, xDesno,
              yGore, yDole;

    /* Cuvanje pripadnosti S skupovima
     * kandidata za podeli pa vladaj */
    KandidatS kS;
};

/* Struktura za poredjenje pravougaonika */
struct PravougaonikComp {
    bool operator()(const Pravougaonik *, const Pravougaonik *) const;
};

/* Definisanje preseka */
using Presek = std::pair<const Pravougaonik *, const Pravougaonik *>;

/* Struktura za poredjenje preseka */
struct PresekComp {
    bool operator()(const Presek &, const Presek &) const;
};

/* Definisanje skupa preseka */
using IntersecSet = std::set<Presek, PresekComp>;

/* Definisanje niza preseka */
using IntersecVec = std::vector<Presek>;

/* Politika azuriranja cvorova stabla; sabloni
 * su iteratori, kao i dva koja zanemarujemo:
 * funkcija za poredjenje i politika alokacije */
template <class CIterator,
          class Iterator,
          class, class>
struct IntervalUpdatePolicy {
    /* Definisanje metapodatka kao celobrojne
     * vrednosti, posto intervalno drvo cuva
     * maksimalnu desnu tacku poddrveta */
    using metadata_type = int;

    /* Neophodan virtuelni destruktor */
    virtual ~IntervalUpdatePolicy() = default;

    /* Virtuelne metode iteracije */
    virtual CIterator node_begin() const = 0;
    virtual CIterator node_end() const = 0;

    /* Provera presecanja dva intervala */
    bool imaPreseka(const Pravougaonik *i1,
                    CIterator i2) const
    {
        return i1->xLevo < (**i2)->xDesno &&
               (**i2)->xLevo < i1->xDesno;
    }

    /* Registrovanje preseka po potrebi */
    void registruj(const Pravougaonik *p,
                   const Pravougaonik *i,
                   IntersecVec &preseci) const
    {
        /* Odredjivanje adekvatnog poretka */
        PravougaonikComp pc;
        if (!pc(p, i)) {
            std::swap(p, i);
        }

        /* Dodavanje ako nema dodirivanja */
        if (p->yGore != i->yDole &&
            i->yGore != p->yDole) {
            preseci.emplace_back(p, i);
        }
    }

    /* Trazenje intervala u podstablu */
    void pretrazi(CIterator it,
                  const Pravougaonik *i,
                  IntersecVec &preseci) const
    {
        /* Baza indukcije je prazno podstablo */
        if (it == node_end()) return;

        /* Hvatanje preseka sa korenom */
        if (imaPreseka(i, it)) {
            registruj(**it, i, preseci);
        }

        /* Provera levog sina ako ima nade */
        const auto levi = it.get_l_child();
        if (levi != node_end() &&
            levi.get_metadata() > i->xLevo){
            pretrazi(levi, i, preseci);
        }

        /* Provera desnog sina ako ima nade */
        if (i->xDesno > (**it)->xLevo) {
            const auto desni = it.get_r_child();
            pretrazi(desni, i, preseci);
        }
    }

    /* Pronalazak svih preklapajucih intervala,
     * pocevsi od korena intervalnog stabla */
    IntersecVec nadjiPreseke(const Pravougaonik *i) const
    {
        IntersecVec preseci;
        pretrazi(node_begin(), i, preseci);
        return preseci;
    }

    /* Popravka stabla nakon izmene, kako
     * bi se odrzala invarijanta drveta */
    void operator()(Iterator it,
                    CIterator end_it)
    {
        /* Podrazumevano je to desna
         * granica tekuceg intervala */
        auto metadata = (**it)->xDesno;

        /* Azuriranje nalevo */
        const auto levo = it.get_l_child();
        if (levo != end_it) {
            metadata = std::max(
                metadata,
                levo.get_metadata());
        }

        /* Azuriranje nadesno */
        const auto desno = it.get_r_child();
        if (desno != end_it) {
            metadata = std::max(
                metadata,
                desno.get_metadata());
        }

        /* Azuriranje metapodatka, uz neophodno
         * uveravanje kompilatora const castom */
        const_cast<int&>(it.get_metadata()) = metadata;
    }
};

/* Definisanje intervalnog drveta */
using IntervalTree = __gnu_pbds::tree<const Pravougaonik *, /* interval */
                     __gnu_pbds::null_type, /* mapiranje (ne postoji) */
                     PravougaonikComp, /* uredjenje nad pravougaonicima */
                     __gnu_pbds::rb_tree_tag, /* crveno-crno stablo */
                     IntervalUpdatePolicy>; /* politika azuriranja */

/* Enumeracija tipa dogadjaja */
enum class TipDogadjaja {DONJA, GORNJA};

/* Struktura koja predstavlja dogadjaj */
struct Dogadjaj {
    /* Konstruktor strukture */
    Dogadjaj(const Pravougaonik *, TipDogadjaja);

    /* Cuvanje pravougaonika i tipa */
    const Pravougaonik *pravougaonik;
    const TipDogadjaja tipDogadjaja;

    /* Cuvanje visine dogadjaja */
    const int y;

    /* Dohvataci koordinata */
    int getXLevo() const;
    int getXDesno() const;
};

/* Struktura za poredjenje dogadjaja */
struct DogadjajComp {
    bool operator()(const Dogadjaj &, const Dogadjaj &) const;
};

/* Definisanje reda dogadjaja */
using EventQueue = std::set<Dogadjaj, DogadjajComp>;

/* Enumeracija tipa vertikalne ivice */
enum class TipIvice {LEVA, DESNA};

/* Definisanje vertikalne stranice */
using VertIvica = std::pair<TipIvice, Pravougaonik *>;

/* Klasa koja predstavlja algoritam */
class PresekPravougaonika : public AlgoritamBaza {
public:
    /* Konstruktor i destruktor klase */
    PresekPravougaonika(QWidget *, int,
                        const bool & = false,
                        std::string = "",
                        int = BROJ_SLUCAJNIH_OBJEKATA);
    virtual ~PresekPravougaonika() override;

    /* Dohvataci rezultata algoritama */
    IntersecSet getGlavni() const;
    IntersecSet getNaivni() const;
    IntersecSet getGruba() const;

    /* Virtuelni metodi iz natklase */
    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *) const final;

    /* Dodatni metod za grubu silu */
    void pokreniAlgoritamGrubeSile();

private:
    /* Rad sa podacima, inicijalizacija */
    bool sekuSe(const Pravougaonik *, const Pravougaonik *) const;
    void ubaciPresek(const Pravougaonik *, const Pravougaonik *, IntersecVec &);
    void generisiNasumicnePravougaonike(int);
    void ucitajPodatkeIzDatoteke(std::string);

    /* Niz pravougaonika, pokazivaca */
    Pravougaonik **_pravougaonici = nullptr;
    unsigned int _n = 0;

    /* Pomocni metod za lepsu animaciju */
    QRect uzmiPresek(IntersecVec, unsigned long) const;

    /* Pomocni metodi za strategiju podeli pa vladaj */
    bool proveriIndeks(unsigned int, unsigned int) const;
    void azurirajIndeks(unsigned int &,
                        unsigned int,
                        KandidatS) const;
    int uzmiIvicu(const VertIvica &) const;

    /* Glavni metodi za strategiju podeli pa vladaj */
    void stab(unsigned int, unsigned int,
              KandidatS, KandidatS);
    void detect(unsigned int, unsigned int);
    void report();

    /* Niz za algoritam grube sile */
    IntersecVec _preseciGruba;

    /* Strukture za metod brisuce prave */
    EventQueue _dogadjaji;
    IntervalTree _status;
    IntersecVec _preseciNaivni;

    /* Nizovi za strategiju podeli pa vladaj */
    VertIvica *_V = nullptr;
    Pravougaonik **_H = nullptr;
    Pravougaonik **_Hh = nullptr;
    IntersecVec _preseciGlavni;

    /* Pomocna polja za lepsu animaciju */
    int _brisucaPravaY;
    std::vector<double> _podele;
    std::experimental::optional<unsigned int> _pocetakNovih;

#ifndef GA06_BENCHMARK
    /* Staticki bafer za daleki skok */
    static jmp_buf _buf;
#endif
};

#endif // GA06_PRESEKPRAVOUGAONIKA_H
