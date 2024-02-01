#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <set>
#include <unordered_map>

#include "algoritambaza.h"
#include "pomocnefunkcije.h"
#include "ga07_datastructures.h"

class Triangulation : public AlgoritamBaza
{
public:
    Triangulation(QWidget *pCrtanje,
                int pauzaKoraka,
                const bool &naivni = false,
                std::string imeDatoteke = "",
                int brojDuzi = BROJ_SLUCAJNIH_OBJEKATA);

    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *painter) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *painter) const final;

    /* Pomocne funkcije */
private:
    std::vector<QPointF> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;
    /* Funkcija generise nasumicne tacke */
    std::vector<QPointF> generisiNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;
    /* Ucitane nasumicne tacke se pripremaju za pravljenje prostog poligona
     * (zadatog u smeru suprotno od kazaljke na satu) */
    std::vector<QPointF> ucitajNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;

    /* Monotono particionisanje */
private:
    void initialiseEventQueue();

    ///
    /// \brief MonotonePartition
    /// Od datog poligona kreira vise monotonih poligona
    void monotonePartition();

    /* Obrada temena za pravljenje monotonog poligona
     * Ove metode su opsane pseudo kodom u Computational Geometry
     * na stranama 53 i 54
     * Brojevi koji stoje u komentarima implementacije ovih metoda
     * odgovaraju brojevima u pseudo kodu u knjizi
     */
    HalfEdge *levaStranica(Vertex *v);
    void handleStartVertex(Vertex *v);
    void handleEndVertex(Vertex *v);
    void handleSplitVertex(Vertex *v);
    void handleMergeVertex(Vertex *v);
    void handleRegularVertex(Vertex *v);

    /* Pomocne funkcije */
private:
    void connectDiagonalsDCEL();
    bool sameDirectionVectors(HalfEdge* e1, HalfEdge* e2);

    /* Triangulacija */
private:
    /*
     * Ova metoda je opsana pseudo kodom u Computational Geometry
     * na strani 57
     * Brojevi koji stoje u komentarima implementacije ove metode
     * odgovara brojevima u pseudo kodu u knjizi
     */
    void triangulacija(Field *f);
    bool istiLanac(HalfEdge *e1, HalfEdge *e2);
    bool leviLanac(HalfEdge *e1, HalfEdge *e2);
    bool desniLanac(HalfEdge *e1, HalfEdge *e2);

private:
    double _brisucaPravaY;
    DCEL _polygon;

    /* MOTONE PARTITION */
    std::set<Vertex*, EventQueueCompTriangulation> _eventQueue;
    std::set<HalfEdge*, StatusQueueCompTriangulation> _statusQueue;
    /* neuredjena jer koristimo samo ubacivanje i pretragu,
     * pa bolje da budu oba u vremenu O(1) umesto O(logn) */
    std::unordered_map<HalfEdge*, Vertex*> _helpers;
    std::vector<std::pair<Vertex*, Vertex*>> _allDiagonals;
    // potrebno samo za crtanje
    bool _monotone;

    /* TRIANGULATION  */
    std::vector<HalfEdge*> _stekTriangulacije;
    std::set<HalfEdge*, EventQueueCompTriangulation2> _eventQueueTriangulation;
};


#endif // TRIANGULATION_H
