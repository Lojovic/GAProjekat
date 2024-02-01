#include "ga04_konveksniomotac3d.h"
#include "pomocnefunkcije.h"

#include <fstream>
#include <cfloat>
#include <cmath>

KonveksniOmotac3D::KonveksniOmotac3D(QWidget *pCrtanje,
                                     int pauzaKoraka,
                                     const bool &naivni,
                                     std::string imeDatoteke,
                                     int brojTacaka)
        : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{
    if (imeDatoteke != "")
        _tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        _tacke = generisiNasumicneTacke(brojTacaka);
}

KonveksniOmotac3D::~KonveksniOmotac3D()
{
    for (auto t: _tacke)
        delete t;

    for (auto i: _ivice)
        delete i;

    for (auto ni: _naivneIvice)
        delete ni;

    for (auto s: _stranice)
        delete s;

    for(auto s: _naivniOmotac){
        delete s;
    }
}

/*--------------------------------------------------------------------------------------------------*/
/*---------------------------------Ucitavanje podataka----------------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/
std::vector<Teme *> KonveksniOmotac3D::generisiNasumicneTacke(int brojTacaka) const
{
    srand(static_cast<unsigned>(time(nullptr)));

    std::vector<Teme*> randomPoints;

    for(int i=0; i < brojTacaka; i++)
        randomPoints.emplace_back(new Teme(
           1.f*rand()/RAND_MAX,
           1.f*rand()/RAND_MAX,
           1.f*rand()/RAND_MAX));

    return randomPoints;
}

std::vector<Teme *> KonveksniOmotac3D::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<Teme*> points;
    float x, y, z;
    while(inputFile >> x >> y >> z)
    {
        points.emplace_back(new Teme(x, y, z));
    }
    return points;
}

const std::unordered_set<Ivica *, HashIvica, EqIvica> &KonveksniOmotac3D::getIvice() const
{
    return _ivice;
}

const std::unordered_set<Ivica *, HashIvica, EqIvica> &KonveksniOmotac3D::getNaivneIvice() const
{
    return _naivneIvice;
}


/*--------------------------------------------------------------------------------------------------*/
/*---------------------------------Implementacija---------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/

void KonveksniOmotac3D::pokreniAlgoritam()
{
    /* Slozenost tetraedra: O(n), samo
     * jedan prolaz kroz temena. */
    if(!Tetraedar()){
        emit animacijaZavrsila();
        return;
    }
    AlgoritamBaza_updateCanvasAndBlock()

    /* Svako teme se dodaje u linearnom vremenu po
     * tekucem broju ivica. Prema lemi iz knjige, broj
     * ivica je najvise O(n). Brisanje viska je takodje
     * linearno, ovog puta po broju izbacenih ivica, ali
     * i to je svakako ograniceno sa O(n). Sveukupno, ova
     * petlja i samim tim ceo algoritam su reda O(n^2).
     * Pre zamene niza skupom, dodavanje temena je bilo
     * kvadratno, pa je ukupno vreme bilo reda O(n^3). */
    for(auto i=0ul; i<_tacke.size() ;i++)
    {
        if(!_tacke[i]->getObradjeno())
        {
            DodajTeme(_tacke[i]);
            ObrisiVisak();
            _tacke[i]->setObradjeno(true);
            AlgoritamBaza_updateCanvasAndBlock()
        }
    }

    emit animacijaZavrsila();
}

bool KonveksniOmotac3D::Tetraedar()
{
    auto it = _tacke.begin();
    Teme *t1 = *it;
    t1->setObradjeno(true);
    it++;
    Teme *t2 = *it;
    t2->setObradjeno(true);
    it++;
    for(;it != _tacke.end(); it++)
        if(!kolinearne(t1, t2, *it))
            break;
    if(it == _tacke.end())
        return false;

    Teme *t3 = *it;
    t3->setObradjeno(true);

    for(;it != _tacke.end(); it++)
        if(static_cast<bool>(pomocneFunkcije::zapremina(t1->koordinate(),
                                                        t2->koordinate(),
                                                        t3->koordinate(),
                                                        (*it)->koordinate())))
            break;
    if(it == _tacke.end())
        return false;

    Teme *t4 = *it;
    t4->setObradjeno(true);

    Ivica *i1 = new Ivica(t1, t2);
    Ivica *i2 = new Ivica(t2, t3);
    Ivica *i3 = new Ivica(t3, t1);
    Stranica *s1 = new Stranica(t1, t2, t3);
     _stranice.push_back(s1);
    _ivice.insert(i1);
    _ivice.insert(i2);
    _ivice.insert(i3);

    i1->postavi_stranicu(s1);
    i2->postavi_stranicu(s1);
    i3->postavi_stranicu(s1);

    if(zapremina6(s1, t4) < 0){
        s1->izmeniRedosledTemena();
        i1->izmeniRedosledTemena();
        i2->izmeniRedosledTemena();
        i3->izmeniRedosledTemena();
    }

    i1->postavi_stranicu(napraviDruguStranicu(i1, t4));
    i2->postavi_stranicu(napraviDruguStranicu(i2, t4));
    i3->postavi_stranicu(napraviDruguStranicu(i3, t4));

    return true;
}

void KonveksniOmotac3D::DodajTeme(Teme* t)
{
    bool vidljiva = false;

    /* Kroz stranice iteriramo tako sto zapravo
     * uzimamo dve po dve prolaskom kroz ivice */
    for(auto ivica : _ivice){
        Stranica *s = ivica->s1();
        float zapremina = zapremina6(s, t);

        if(zapremina <= 0){
            vidljiva = true;
            s->setVidljiva(true);
        }

        s = ivica->s2();
        zapremina = zapremina6(s, t);

        if(zapremina <= 0){
            vidljiva = true;
            s->setVidljiva(true);
        }
    }

    if(!vidljiva){
        return;
    }

    /* Skup ivica nije konstantan, pa je ovde neophodno kopirati ih u vektor, kako
     * bi se proslo samo kroz stare ivice, ne i one koje su novododate u prolazu.
     * Ovo mozda zvuci skupo, ali je zapravo vrlo jeftino i ne utice na slozenost.
     * Problem u originalnom kodu: prolazilo se i kroz novododate ivice, pri cemu
     * su one po pravilu menjane (samim sobom). U kombinaciji sa preterano velikom
     * tolerancijom na numericku gresku, ovo je po pravilu vodilo ka null stranicama,
     * te SIGSEGV (segmentation fault) pri svim iole vecim ulazima (npr. preko sto) */
    for(auto ivica: std::vector<Ivica*>(_ivice.cbegin(), _ivice.cend())){
        const auto s1Losa = ivica->s1()->getVidljiva();
        const auto s2Losa = ivica->s2()->getVidljiva();
        if(s1Losa && s2Losa)
            ivica->setObrisati(true);
        else if(s1Losa)
            ivica->zameniVidljivuStranicu(napraviPrvuStranicu(ivica, t), 0);
        else if(s2Losa)
            ivica->zameniVidljivuStranicu(napraviDruguStranicu(ivica, t), 1);
    }
}

void KonveksniOmotac3D::ObrisiVisak()
{
    std::experimental::erase_if(_ivice, [](Ivica* ivica){
        bool rez = ivica->obrisati();
        if(rez)
            delete ivica;
        return rez;
    });
}

/*--------------------------------------------------------------------------------------------------*/
/*---------------------------------Pomocni metodi---------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/
float KonveksniOmotac3D::zapremina6(Stranica *s, Teme *t) const
{
    return pomocneFunkcije::zapremina(s->t1()->koordinate(),
                                      s->t2()->koordinate(),
                                      s->t3()->koordinate(),
                                      t->koordinate());
}

bool KonveksniOmotac3D::kolinearne(Teme *a, Teme *b, Teme *c) const
{
    return pomocneFunkcije::kolinearne3D(a->koordinate(), b->koordinate(), c->koordinate());
}

Stranica* KonveksniOmotac3D::napraviStranicu(Teme *i1t1, Teme *i1t2,
                                             Teme *i2t1, Teme *i2t2,
                                             Teme *st1, Teme *st2, Teme *st3)
{
    /* Pretraga u neuredjenom skupu je u konstantnom vremenu. Pravi
     * se zeljena ivica i zatim pronalazi fjom find. Ako je uspesno
     * nadjena, uzima se taj pokazivac. Inace se prosiruje skup. */
    Ivica* i1 = new Ivica(i1t1, i1t2);
    auto iv1 = _ivice.find(i1);
    if(iv1 != _ivice.end()){
        delete i1;
        i1 = *iv1;
    } else {
        _ivice.insert(i1);
    }

    Ivica* i2 = new Ivica(i2t1, i2t2);
    auto iv2 = _ivice.find(i2);
    if(iv2 != _ivice.end()){
        delete i2;
        i2 = *iv2;
    } else {
        _ivice.insert(i2);
    }

    Stranica *s = new Stranica(st1, st2, st3);
     _stranice.push_back(s);
    i1->postavi_stranicu(s);
    i2->postavi_stranicu(s);
    return s;
}

Stranica* KonveksniOmotac3D::napraviPrvuStranicu(Ivica *iv, Teme *t)
{
    return napraviStranicu(t, iv->t1(), iv->t2(), t, iv->t1(), iv->t2(), t);
}

Stranica* KonveksniOmotac3D::napraviDruguStranicu(Ivica *iv, Teme *t)
{
    return napraviStranicu(iv->t1(), t, t, iv->t2(), iv->t2(), iv->t1(), t);
}

/*--------------------------------------------------------------------------------------------------*/
/*---------------------------------Crtanje i Naivni algoritam---------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/
void KonveksniOmotac3D::crtajTeme(Teme *t) const
{
    /* Crtanje tacaka preko koordinata */
    glColor3d(1, 0, 0.4);
    glVertex3f(t->x(), t->y(), t->z());
}

void KonveksniOmotac3D::crtajStranicu(Stranica* s) const
{
    /* Postavljanje boje */
    glColor4dv(s->boje);

    /* Crtanje stranice kao trougla */
    glBegin(GL_POLYGON);
        glVertex3f(s->t1()->x(), s->t1()->y(), s->t1()->z());
        glVertex3f(s->t2()->x(), s->t2()->y(), s->t2()->z());
        glVertex3f(s->t3()->x(), s->t3()->y(), s->t3()->z());
    glEnd();
}

void KonveksniOmotac3D::crtajAlgoritam(QPainter*) const
{
    /* Crtanje svih temena */
    glBegin(GL_POINTS);
        for(auto teme : _tacke)
            crtajTeme(teme);
    glEnd();

    /* Crtanje svih stranica prolaskom
     * dve po dve kroz skup ivica */
    for(auto ivica: _ivice)
    {
        crtajStranicu(ivica->s1());
        crtajStranicu(ivica->s2());
    }
}

void KonveksniOmotac3D::pokreniNaivniAlgoritam()
{
    /* Za svake tri nekolinearne tacke se proverava da li stranica koju obrazuju
     * pripada konveksnom omotacu.
     * Ako se sve ostale tacke nalaze sa iste strane stranice, stranica pripada
     * konveksnom omotacu, u suprotnom ne pripada.
     * Slozenost naivnog algoritma je O(n^4). */

    for (auto i = 0ul; i < _tacke.size(); i++) {
        for (auto j = i+1; j < _tacke.size(); j++) {
            for (auto k = j+1; k < _tacke.size(); k++) {
                /* if(i==j || i==k || j==k)
                    continue;*/
                if (kolinearne(_tacke[i], _tacke[j], _tacke[k]))
                    continue;

                //cuva se tekuca stranica, radi lakseg pracenja algoritma
                _tekucaStranica = new Stranica(_tacke[i], _tacke[j], _tacke[k]);

                AlgoritamBaza_updateCanvasAndBlock()

                // Pronalazenje jedne nenula zapremine
                auto zapremina = 0.f;
                for (auto tacka : _tacke) {
                    zapremina = zapremina6(_tekucaStranica, tacka);
                    if (fabsf(zapremina) > EPSf)
                        break;
                }
                // Ako ne postoje 4 nekomplanarne tacke zavrsava se algoritam
                if (fabsf(zapremina) <= EPSf)
                    return;


                int znakZapremine = zapremina > 0 ? 1 : -1;

                auto it = _tacke.begin();
                for (; it != _tacke.end(); it++) {
                    zapremina = zapremina6(_tekucaStranica, *it);
                    if (zapremina * znakZapremine < -EPSf)
                        break;
                }

               //delete stranica;

                if (it == _tacke.end()) {
                    _naivniOmotac.push_back(_tekucaStranica);
                    AlgoritamBaza_updateCanvasAndBlock()
                    /* Sve zapremine su istog znaka (dozvoljavamo da su neke i nula),
                     * pa stranica pripada konveksnom omotacu */

                     /*Postoji funkcija za iscrtavanje stranice*/
                     auto *ivica1 = new Ivica(_tacke[i], _tacke[j]);
                     auto *ivica2 = new Ivica(_tacke[i], _tacke[k]);
                     auto *ivica3 = new Ivica(_tacke[j], _tacke[k]);
                     _naivneIvice.insert(ivica1);
                     _naivneIvice.insert(ivica2);
                     _naivneIvice.insert(ivica3);
                }

            }
        }
    }
    //da ne bi iscrtalo poslednju tekucu stranicu
    _tekucaStranica = nullptr;
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}

void KonveksniOmotac3D::crtajNaivniAlgoritam(QPainter *) const
{
    // crtanje svih temena
    glBegin(GL_POINTS);
        for(auto teme : _tacke)
            crtajTeme(teme);
    glEnd();

    if(_tekucaStranica!=nullptr){
        crtajStranicu(_tekucaStranica);
    }

    for(auto stranica: _naivniOmotac){
         crtajStranicu(stranica);
    }
}

const std::vector<Stranica*> &KonveksniOmotac3D::getNaivniOmotac3d() const{
    return _naivniOmotac;
}

const std::unordered_set<Ivica*, HashIvica, EqIvica> &KonveksniOmotac3D::getKonveksniOmotac3d() const{
   return _ivice;
}
