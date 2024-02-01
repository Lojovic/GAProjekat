#include "ga06_presekPravougaonika.h"

#include <fstream>

/* Konstrukcija pravougaonika od dve tacke */
Pravougaonik::Pravougaonik(QPoint &topLeft, QPoint &bottomRight)
    : QRect(topLeft, bottomRight), xLevo(topLeft.x()), xDesno(bottomRight.x()),
                                   yGore(topLeft.y()), yDole(bottomRight.y())
{ }


/* Konstrukcija pravougaonika od parametara */
Pravougaonik::Pravougaonik(int x, int y, int w, int h)
    : QRect(x, y, w, h), xLevo(x), xDesno(x+w),
                         yGore(y+h), yDole(y)
{ }

/* Struktura za poredjenje pravougaonika */
bool PravougaonikComp::operator()(const Pravougaonik *l,
                                  const Pravougaonik *d) const
{
    /* Odustajanje u slucaju greske */
    if (!l && !d) {
        return false;
    /* Manji je levlji pravougaonik */
    } else if (l->xLevo != d->xLevo) {
        return l->xLevo < d->xLevo;
    /* Manji je levlji pravougaonik */
    } else if (l->xDesno != d->xDesno) {
        return l->xDesno < d->xDesno;
    /* Inace se porede pokazivaci */
    } else {
        return l < d;
    }
}

/* Struktura za poredjenje preseka */
bool PresekComp::operator()(const Presek &l, const Presek &d) const
{
    /* Komparator pravougaonika */
    PravougaonikComp pc;

    /* Jednakost nezavisno od redosleda */
    if ((l.first == d.first && l.second == d.second) ||
        (l.first == d.second && l.second == d.first)) {
        return false;
    /* Inace poredjenje pravougaonika */
    } else if (l.first != d.first) {
        return pc(l.first, d.first);
    } else {
        return pc(l.second, d.second);
    }
}

/* Konstrukcija dogadjaja */
Dogadjaj::Dogadjaj(const Pravougaonik *p, TipDogadjaja t)
    : pravougaonik(p), tipDogadjaja(t),
      y(t == TipDogadjaja::GORNJA ? p->yGore : p->yDole)
{ }

/* Dohvatanje levog intervala pravougaonika */
int Dogadjaj::getXLevo() const
{
    return pravougaonik->xLevo;
}

/* Dohvatanje desnog intervala pravougaonika */
int Dogadjaj::getXDesno() const
{
    return pravougaonik->xDesno;
}

/* Struktura za poredjenje dogadjaja */
bool DogadjajComp::operator()(const Dogadjaj &l,
                              const Dogadjaj &d) const
{
    /* Manji je dogadjaj na visoj poziciji */
    if (l.y != d.y) {
        return l.y > d.y;
    /* Odustajanje u slucaju greske */
    } else if (!l.pravougaonik || !d.pravougaonik) {
        return false;
    /* Manji je dogadjaj dodavanja od izbacivanja */
    } else if (l.tipDogadjaja != d.tipDogadjaja) {
        return l.tipDogadjaja < d.tipDogadjaja;
    /* Manji je dogadjaj na levljoj poziciji */
    } else if (l.getXLevo() != d.getXLevo()) {
        return l.getXLevo() < d.getXLevo();
    /* Manji je dogadjaj na levljoj poziciji */
    } else if (l.getXDesno() != d.getXDesno()) {
        return l.getXDesno() < d.getXDesno();
    /* Inace se porede pokazivaci */
    } else {
        return l.pravougaonik < d.pravougaonik;
    }
}

/* Konstrukcija algoritma */
PresekPravougaonika::PresekPravougaonika(QWidget *pCrtanje,
                                         int pauzaKoraka,
                                         const bool &naivni,
                                         std::string imeDatoteke,
                                         int brojPravougaonika)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni),
      _brisucaPravaY(_pCrtanje ? _pCrtanje->height()-3 : 0)
{
    /* Inicijalizacija niza pravougaonika */
    if (imeDatoteke == "")
        generisiNasumicnePravougaonike(brojPravougaonika);
    else
        ucitajPodatkeIzDatoteke(imeDatoteke);
}

/* Deinicijalizacija algoritma */
PresekPravougaonika::~PresekPravougaonika()
{
    for (auto i = 0ul; i < _n; i++) {
        delete _pravougaonici[i];
    }
    delete [] _pravougaonici;
    delete [] _V;
    delete [] _H;
    delete [] _Hh;
}

/* Dohvatanje skupa preseka algoritma podeli pa vladaj */
IntersecSet PresekPravougaonika::getGlavni() const
{
    return IntersecSet(std::cbegin(_preseciGlavni),
                       std::cend(_preseciGlavni));
}

/* Dohvatanje skupa preseka algoritma brisuce prave */
IntersecSet PresekPravougaonika::getNaivni() const
{
    return IntersecSet(std::cbegin(_preseciNaivni),
                       std::cend(_preseciNaivni));
}

/* Dohvatanje skupa preseka algoritma grube sile */
IntersecSet PresekPravougaonika::getGruba() const
{
    return IntersecSet(std::cbegin(_preseciGruba),
                       std::cend(_preseciGruba));
}

#ifndef GA06_BENCHMARK
/* Staticki bafer za daleki skok */
jmp_buf PresekPravougaonika::_buf;
#endif

/* Algoritam zasnovan na strategiji podeli pa vladaj, pogledati
 * https://www.widelands.org/~sirver/wl/141229_devail_rects.pdf */
void PresekPravougaonika::pokreniAlgoritam()
{
    /* Ciscenje skupa preseka */
    _preseciGlavni.clear();

#ifndef GA06_BENCHMARK
    /* Postavljanje dalekog skoka, za lako iskakanje iz
     * duboke rekurzije po potrebi (prekid animacije)*/
    if (setjmp(_buf)) return;
#endif

    /* Prijavljivanje svih preseka u paru */
    report();

#ifdef GA06_BENCHMARK
    /* Popunjavanje crteza konacnim rezultatom */
    if (_pCrtanje) {
        _pCrtanje->update();
    }
#endif

    /* Obavestavanje pozivaoca o gotovoj animaciji */
    emit animacijaZavrsila();
}

/* Iscrtavanje tekuceg stanja algoritma */
void PresekPravougaonika::crtajAlgoritam(QPainter *painter) const
{
    /* Odustajanje u slucaju greske */
    if (!painter) return;

    /* Iscrtavanje svih preseka */
    for (auto i = 0ul; i < _preseciGlavni.size(); i++) {
        painter->fillRect(uzmiPresek(_preseciGlavni, i), Qt::yellow);
    }

    /* Iscrtavanje novih preseka */
    if (_pocetakNovih) {
        for (auto i = _pocetakNovih.value();
             i < _preseciGlavni.size(); i++) {
            painter->fillRect(uzmiPresek(_preseciGlavni, i), Qt::red);
        }
    }

    /* Iscrtavanje svakog pravougaonika */
    for (auto i = 0ul; i < _n; i++) {
        painter->drawRect(*_pravougaonici[i]);
    }

    /* Podesavanje stila crtanja */
    auto olovka = painter->pen();
    olovka.setColor(Qt::darkGreen);
    olovka.setWidthF(olovka.widthF()/1.5);
    painter->setPen(olovka);

    /* Iscrtavanje tekucih podela */
    for (const auto s : _podele) {
        painter->drawLine(QPointF(s, 0),
                          QPointF(s, _pCrtanje->height()));
    }
}

/* Algoritam zasnovan na brisucoj (pokretnoj) pravoj */
void PresekPravougaonika::pokreniNaivniAlgoritam()
{
    /* Ciscenje skupa preseka */
    _preseciNaivni.clear();

    /* Popunjavanje niza dogadjaja */
    for (auto i = 0ul; i < _n; i++) {
        _dogadjaji.emplace(_pravougaonici[i], TipDogadjaja::GORNJA);
        _dogadjaji.emplace(_pravougaonici[i], TipDogadjaja::DONJA);
    }

    /* Prolazak redom kroz dogadjaje */
    for (const auto &dogadjaj : _dogadjaji) {
        const auto pravougaonik = dogadjaj.pravougaonik;

        /* Obrada nailaska na gornju stranicu */
        if (dogadjaj.tipDogadjaja == TipDogadjaja::GORNJA) {
#ifndef GA06_BENCHMARK
            /* Azuriranje polozaja brisuce prave po potrebi */
            if (_brisucaPravaY != pravougaonik->yGore) {
                _brisucaPravaY = pravougaonik->yGore;
                PresekPravougaonika_updateCanvasAndBlock()
            }
#endif

            /* Nalazenje svih trenutnih preseka */
            const auto preseci = _status.nadjiPreseke(pravougaonik);

            /* Ubacivanje novog pravougaonika u status */
            _status.insert(pravougaonik);

#ifndef GA06_BENCHMARK
            /* Osvezavanje crteza novim pravougaonikom */
            PresekPravougaonika_updateCanvasAndBlock()
#endif

            _preseciNaivni.insert(std::cend(_preseciNaivni),
                                  std::cbegin(preseci), std::cend(preseci));

#ifndef GA06_BENCHMARK
            /* Osvezavanje crteza ako ima novih preseka */
            if (!preseci.empty()) {
                PresekPravougaonika_updateCanvasAndBlock()
            }
#endif
        /* Obrada nailaska na donju stranicu */
        } else /*if (dogadjaj.tipDogadjaja == TipDogadjaja::DONJA)*/ {
#ifndef GA06_BENCHMARK
            /* Azuriranje polozaja brisuce prave po potrebi */
            if (_brisucaPravaY != pravougaonik->yDole) {
                _brisucaPravaY = pravougaonik->yDole;
                PresekPravougaonika_updateCanvasAndBlock()
            }
#endif

            /* Izbacivanje starog pravougaonika iz statusa */
            _status.erase(pravougaonik);

#ifndef GA06_BENCHMARK
            /* Osvezavanje crteza jednim pravougaonikom manje */
            PresekPravougaonika_updateCanvasAndBlock()
#endif
        }
    }

    /* Ciscenje reda dogadjaja */
    _dogadjaji.clear();

    /* Zavrsni frejm sa skroz spustenom pravom */
#ifndef GA06_BENCHMARK
    _brisucaPravaY = 3;
    PresekPravougaonika_updateCanvasAndBlock()
#endif

    /* Obavestavanje pozivaoca o gotovoj animaciji */
    emit animacijaZavrsila();
}

void PresekPravougaonika::crtajNaivniAlgoritam(QPainter *painter) const
{
    /* Odustajanje u slucaju greske */
    if (!painter) return;

    /* Iscrtavanje svih preseka */
    for (auto i = 0ul; i < _preseciNaivni.size(); i++) {
        painter->fillRect(uzmiPresek(_preseciNaivni, i), Qt::yellow);
    }

    /* Iscrtavanje svakog pravougaonika */
    for (auto i = 0ul; i < _n; i++) {
        painter->drawRect(*_pravougaonici[i]);
    }

    /* Podesavanje stila crtanja */
    auto olovka = painter->pen();
    olovka.setColor(Qt::red);
    painter->setPen(olovka);

    /* Crtanje pravougaonika iz statusa */
    for (const auto p : _status)
    {
        painter->drawRect(*p);
    }

    /* Podesavanje stila crtanja */
    olovka.setColor(Qt::darkGreen);
    olovka.setWidth(2*olovka.width());
    painter->setPen(olovka);

    /* Crtanje brisuce prave */
    painter->drawLine(0, _brisucaPravaY,
                      _pCrtanje->width(),
                      _brisucaPravaY);
}

/* Provera da li postoji presek para pravougaonika */
bool PresekPravougaonika::sekuSe(const Pravougaonik *p1,
                                 const Pravougaonik *p2) const
{
    return p1 && p2 && /* nisu null */
           p1->xLevo < p2->xDesno &&
           p1->xDesno > p2->xLevo &&
           p1->yDole < p2->yGore &&
           p1->yGore > p2->yDole;
}

/* Ubacivanje preseka u skup ako treba */
void PresekPravougaonika::ubaciPresek(const Pravougaonik *p1,
                                      const Pravougaonik *p2,
                                      IntersecVec &preseci)
{
    /* Odredjivanje adekvatnog poretka */
    PravougaonikComp pc;
    if (!pc(p1, p2)) {
        std::swap(p1, p2);
    }

    /* Dodavanje ako nema dodirivanja */
    if (p1->xDesno != p2->xLevo) {
        preseci.emplace_back(p1, p2);
    }
}

/* Kvadratni algoritam grube sile */
void PresekPravougaonika::pokreniAlgoritamGrubeSile()
{
    /* Ciscenje skupa preseka */
    _preseciGruba.clear();

    /* Uporedjivanje pravougaonika u parovima */
    for (auto i = 0ul; i < _n; i++) {
        const auto *p1 = _pravougaonici[i];
        for (auto j = i+1; j < _n; j++) {
            const auto *p2 = _pravougaonici[j];

            /* Odredjivanje preseka i dodavanje
             * u skup ako postoji (nije prazan) */
            if (sekuSe(p1, p2)) {
                ubaciPresek(p1, p2, _preseciGruba);
            }
        }
    }
}

/* Pravljenje slucajnih pravougaonika */
void PresekPravougaonika::generisiNasumicnePravougaonike(int brojPravougaonika)
{
    /* Alokacija potrebnog prostora */
    _n = 0;
    _pravougaonici = new Pravougaonik *[static_cast<unsigned>(brojPravougaonika)];

    /* Generisanje duplo veceg broja tacaka */
    const auto tacke = generisiNasumicneTacke(2*brojPravougaonika);

    /* Za svaki par pokusaj pravljenja pravougaonika */
    for (auto i = 0ul; i < tacke.size(); i+=2) {
        auto tacka1 = tacke[i];
        auto tacka2 = tacke[i+1];

        /* Preskakanje losih parova */
        if (tacka1.x() == tacka2.x() ||
            tacka1.y() == tacka2.y()) {
            continue;
        }

        /* Prva treba da bude levo */
        if (tacka1.x() > tacka2.x()) {
            std::swap(tacka1.rx(), tacka2.rx());
        }

        /* Prva treba da bude gore */
        if (tacka1.y() < tacka2.y()) {
            std::swap(tacka1.ry(), tacka2.ry());
        }

        /* Pravljenje pomocu sredjenih tacaka */
        _pravougaonici[_n++] = new Pravougaonik(tacka1, tacka2);
    }
}

/* Citanje pravougaonika iz datoteke */
void PresekPravougaonika::ucitajPodatkeIzDatoteke(std::string imeDatoteke)
{
    /* Otvaranje zadatog fajla */
    std::ifstream datoteka(imeDatoteke);

    /* Alokacija potrebnog prostora */
    datoteka >> _n;
    _pravougaonici = new Pravougaonik *[_n];

    /* Ucitavanje svakog pravougaonika pojedinacno */
    for (auto i = 0ul; i < _n; i++) {
        int x, y, w, h;
        datoteka >> x >> y >> w >> h;
        _pravougaonici[i] = new Pravougaonik(x, y, w, h);
    }
}

/* Racunanje tacnog preseka dva pravougaonika */
QRect PresekPravougaonika::uzmiPresek(IntersecVec preseci, unsigned long i) const
{
    return *preseci[i].first & *preseci[i].second;
}

/* Pomocna funkcija za proveru indeksa */
bool PresekPravougaonika::proveriIndeks(unsigned int i,
                                        unsigned int d) const
{
    return i < d && _H[i];
}

/* Pomocna funkcija za azuriranje indeksa */
void PresekPravougaonika::azurirajIndeks(unsigned int &i,
                                         unsigned int d,
                                         KandidatS k) const
{
    do i++; while (proveriIndeks(i, d) && _H[i]->kS != k);
}

/* Odredjivanje preseka najuzih kandidata; skupovi
 * A i B iz rada su apstrahovani enumeracijom */
void PresekPravougaonika::stab(unsigned int l, unsigned int d,
                               KandidatS A, KandidatS B)
{
    /* Odredjivanje pocetnih indeksa */
    auto i = l-1; azurirajIndeks(i, d, A);
    auto j = l-1; azurirajIndeks(j, d, B);

    /* Prolazak kroz sve kandidate */
    while (proveriIndeks(i, d) && proveriIndeks(j, d)) {
        /* Ako je prvi skup ispod drugog */
        if (_H[i]->yDole < _H[j]->yDole) {
            auto k = j;
            /* Registrovanje svih susednih */
            while (proveriIndeks(k, d) &&
                   _H[k]->yDole < _H[i]->yGore) {
                ubaciPresek(_H[i], _H[k], _preseciGlavni);
                azurirajIndeks(k, d, B);
            }
            azurirajIndeks(i, d, A);
        /* Ako je prvi skup iznad drugog */
        } else {
            auto k = i;
            /* Registrovanje svih susednih */
            while (proveriIndeks(k, d) &&
                   _H[k]->yDole < _H[j]->yGore) {
                ubaciPresek(_H[j], _H[k], _preseciGlavni);
                azurirajIndeks(k, d, A);
            }
            azurirajIndeks(j, d, B);
        }
    }
}

/* Uzimanje odgovarajuce ivice pravougaonika */
int PresekPravougaonika::uzmiIvicu(const VertIvica &ivica) const
{
    return ivica.first == TipIvice::LEVA ?
                ivica.second->xLevo : ivica.second->xDesno;
}

/* Pronalazenje preseka u potprostoru [l, d) */
void PresekPravougaonika::detect(unsigned int l, unsigned int d)
{
    /* Nema preseka ako je samo jedan pravougaonik */
    if (d-l < 2) return;

    /* Podrazumevano nijedan pravougaonik nije zanimljiv */
    for (auto i = l; i < d; i++) {
        _V[i].second->kS = KandidatS::Sxx;
    }

    /* Odredjivanje sredine intervala za podelu; ovo
     * imitira skupove V1 i V2 iz rada bez kopiranja */
    const auto s = l+(d-l)/2;

#ifndef GA06_BENCHMARK
    /* Ubacivanje podele u niz svih */
    _podele.push_back((uzmiIvicu(_V[s-1])+
                       uzmiIvicu(_V[s]))/2.);
    PresekPravougaonika_updateCanvasAndBlock()
#endif

    /* Odredjivanje skupova S11 i S12 iz rada */
    for (auto i = l; i < s; i++) {
        /* Pravougaonik je skroz levo od podele */
        if (_V[i].first == TipIvice::DESNA) {
            _V[i].second->kS = KandidatS::S11;
        /* Pravougaonik obuhvata desnu polovinu */
        } else if (_V[i].second->xDesno >= uzmiIvicu(_V[d-1])) {
            _V[i].second->kS = KandidatS::S12;
        }
    }

    /* Odredjivanje skupova S22 i S21 iz rada */
    for (auto i = s; i < d; i++) {
        /* Pravougaonik je skroz desno od podele */
        if (_V[i].first == TipIvice::LEVA) {
            _V[i].second->kS = KandidatS::S22;
        /* Pravougaonik obuhvata levu polovinu */
        } else if (_V[i].second->xLevo <= uzmiIvicu(_V[l])) {
            _V[i].second->kS = KandidatS::S21;
        }
    }

#ifndef GA06_BENCHMARK
    /* Pamcenje prvog indeksa novih preseka */
    _pocetakNovih = _preseciGlavni.size();
#endif

    /* Odredjivanje preseka najuzih kandidata */
    stab(l, d, KandidatS::S12, KandidatS::S22);
    stab(l, d, KandidatS::S21, KandidatS::S11);
    stab(l, d, KandidatS::S12, KandidatS::S21);

#ifndef GA06_BENCHMARK
    /* Osvezavanje crteza novim presecima */
    if (_pocetakNovih.value() < _preseciGlavni.size()) {
        PresekPravougaonika_updateCanvasAndBlock()
        _pocetakNovih = _preseciGlavni.size();
        PresekPravougaonika_updateCanvasAndBlock()
    }
#endif

    /* Particionisanje niza H na H1 i H2 iz rada */
    auto H1 = l, H2 = s;
    for (auto i = l; proveriIndeks(i, d); i++) {
        switch (_H[i]->kS) {
        /* Neodredjeni idu u obe polovine */
        case KandidatS::Sxx:
            _Hh[H1++] = _H[i];
            _Hh[H2++] = _H[i];
            break;
        /* Levo idu oni koju su u toj polovini */
        case KandidatS::S11:
        case KandidatS::S12:
            _Hh[H1++] = _H[i];
            break;
        /* Desno idu oni koji su u toj polovini */
        case KandidatS::S21:
        case KandidatS::S22:
            _Hh[H2++] = _H[i];
            break;
        }
    }

    /* Popunjavanje preostalih rupa */
    std::fill(_Hh+H1, _Hh+s, nullptr);
    std::fill(_Hh+H2, _Hh+d, nullptr);

    /* Prepisivanje particionisanog niza nazad, i
     * to vrlo jednostavno, zamenom pokazivaca */
    std::swap(_H, _Hh);

    /* Pronalazenje preseka u potprostorima */
    detect(s, d); detect(l, s);

    /* Vracanje pokazivaca kako bi sve bilo okej */
    std::swap(_H, _Hh);

#ifndef GA06_BENCHMARK
    /* Izbacivanje podele iz niza */
    _podele.pop_back();
    PresekPravougaonika_updateCanvasAndBlock()
#endif
}

/* Prijavljivanje svih preseka u paru */
void PresekPravougaonika::report()
{
    /* Pravljenje i sortiranje niza vertikalnih ivica */
    _V = new VertIvica[2*_n];
    for (auto i = 0ul; i < _n; i++) {
        _V[2*i] = std::make_pair(TipIvice::LEVA, _pravougaonici[i]);
        _V[2*i+1] = std::make_pair(TipIvice::DESNA, _pravougaonici[i]);
    }
    std::sort(_V, _V+2*_n, [this](const VertIvica &a, const VertIvica &b)
                                 /* Poredjenje odgovarajucih ivica */
                                 { return uzmiIvicu(a) < uzmiIvicu(b); });

    /* Pravljenje i sortiranje niza pravougaonika
     * posmatranih kao vertikalnih intervala */
    _H = new Pravougaonik *[2*_n];
    std::copy(_pravougaonici, _pravougaonici+_n, _H);
    std::sort(_H, _H+_n, [](const Pravougaonik *a, const Pravougaonik *b)
                           { return a->yDole < b->yDole; });

    /* Popunjavanje preostalih praznina */
    std::fill_n(_H+_n, _n, nullptr);

    /* Pravljenje pomocnog sortiranog niza */
    _Hh = new Pravougaonik *[2*_n];

    /* Pronalazenje preseka u (celom) potprostoru */
    detect(0, 2*_n);
}
