///
/// Svaki algoritam koji cemo implementirati treba da:
/// 1. Nasledi AlgoritamBaza (koji predstavlja apstrakciju koja vodi racuna o animaciji i iscrtavanju).
/// 2. Implementira metod "pokreniAlgoritam()" u kom ce biti implementacija algoritma.
/// 3. Implementira metod "crtajAlgoritam()" u kom ce biti implementirano iscrtavanje.
/// 4. Implementira metod "pokreniNaivniAlgoritam()" za potrebe poredjenja performansi.
/// 5. Implementira metod "crtajNaivniAlgoritam()" za ilustrovanje naivnog algoritma.
///
/// Svaki put kada se promeni stanje algoritma (kada je potrebno promeniti crtez),
/// potrebno je pozvati metod AlgoritamBaza_updateCanvasAndBlock();
/// Ovo je zapravo makro, koji pozove metod updateCanvasAndBlock() i u slucaju da
/// detektuje da algoritam treba da se zaustavi (jedna bool promenljiva koja se
/// postavlja na true kada kosirnik klikne "Zaustavi"), samo se zamenjuje komandom
/// return;
///
/// Razlog za ovaj mehanizam je to sto se algoritam izvrsava u okviru svoje niti i nije bezbedno
/// (sa aspekta resursa) ubiti tu nit dok se ne zavrsi. Na ovaj nacin, svaki put
/// kada se vrsi azuriranje crteza, proverava se i stanje pomenute logicke promenljive
/// i ukoliko ona signalizira da treba da se zaustavimo, na mestu azuriranja crteza ce
/// se naci komanda return; i na taj nacin smo postigli efekat koji smo zeleli
/// (izaslo se iz funkcije koja se izvrsava u okviru niti za animaciju)

#ifndef ALGORITAMBAZA_H
#define ALGORITAMBAZA_H

#include <QObject>
#include <QPainter>
#include <QSemaphore>
#include <QWidget>
#include <QOpenGLWidget>
#include <QCheckBox>

#include "animacijanit.h"
#include "config.h"

#define AlgoritamBaza_updateCanvasAndBlock() \
    if (updateCanvasAndBlock()) \
    { \
        return; \
    }

/* Nasledjuje se QObject koji omogucava rad sa signalima i slotovima.
 */

class AlgoritamBaza : public QObject
{
    Q_OBJECT
private:
    static int constexpr INVALID_TIMER_ID = -1;
    static int constexpr DRAWING_BORDER = 10;

    /* Parametri za animaciju */
    int _pauzaKoraka;
    int _timerId;
    QSemaphore _semafor;
    bool _unistiAnimaciju;

    /* Nit koja izvrsava algoritam */
    AnimacijaNit *_pNit;

    ///
    /// \brief timerEvent - funkcija koja se poziva na svakih _delayMs ms.
    ///     U njoj samo oslobadjamo semafor i na taj nacin omogucavamo da se predje na sledeci
    ///     korak algoritma
    /// \param event
    /// The QTimer class provides a high-level programming interface for timers. To use it, create a QTimer,
    /// connect its timeout() signal to the appropriate slots, and call start(). From then on, it will emit
    /// the timeout() signal at constant intervals.
    ///
    void timerEvent(QTimerEvent *event) override;

public:
    AlgoritamBaza(QWidget *pCrtanje, int pauzaKoraka, const bool &naivni);

    /* Algoritme ne treba da bude moguce kopirati */
    virtual ~AlgoritamBaza() override = default;
    AlgoritamBaza(const AlgoritamBaza &) = delete;
    AlgoritamBaza(AlgoritamBaza &&) = delete;
    AlgoritamBaza& operator=(const AlgoritamBaza &) = delete;
    AlgoritamBaza& operator=(AlgoritamBaza &&) = delete;

    ///
    /// \brief pokreniAlgoritam - funkcija koju ce svaka podklasa implementirati,
    ///         a koja predstavlja izvrsavanje konkretnog algoritma
    ///
    virtual void pokreniAlgoritam() = 0;
    virtual void crtajAlgoritam(QPainter *painter) const = 0;
    virtual void pokreniNaivniAlgoritam() = 0;
    virtual void crtajNaivniAlgoritam(QPainter *painter) const = 0;

    bool naivni() const;
    void crtaj(QPainter *painter = nullptr) const;

    ///
    /// \brief pokreniAnimaciju - funkcija za pokretanje animacije
    ///     1. Pravi i pokrece nit za konkretan algoritam
    ///     2. Pokrece tajmer
    ///
    void pokreniAnimaciju();

    ///
    /// \brief pauzirajIliNastaviAnimaciju - pauziranja, odnosno ponovno pokretanje animacije
    ///         Ako je tajmer pokrenut, zaustavlja ga
    ///         Ako tajmer nije pokrenut, pokrece ga
    ///
    void pauzirajIliNastaviAnimaciju();

    ///
    /// \brief sledeciKorakAnimacije - izvrsavanje algoritma korak po korak
    ///     1. Zaustavljanje tajmera (da bi se zaustavilo kontinuirano izvrsavanje algoritma)
    ///     2. Oslobadjanje semafora (da bi se izvrsio jedan korak algoritma)
    ///
    void sledeciKorakAnimacije();

    ///
    /// \brief zaustaviAnimaciju - zaustavljanje animacije
    ///     1. Postavljanje logicke promenljive na true
    ///     2. Pustanje semafora (da bi se animacija "nastavila", tj. da bi se doslo do updateCanvasAndBlock())
    ///         2.1. U updateCanvasAndBlock ce se detektovati da je logicka promenljiva true i algoritam ce se zaustaviti
    ///     3. Oslobadjanje odgovarajucih resursa
    ///
    void zaustaviAnimaciju();

    ///
    /// \brief promeniDuzinuPauze - funkcija za promenu brzine animacije
    /// \param duzinaPauze - pauza izmedju dva koraka, izrazena u ms
    ///
    void promeniDuzinuPauze(int duzinaPauze);

signals:
    void animacijaZavrsila();

protected:
    ///
    /// \brief _pCrtanje - oblast crtanja
    ///
    QWidget *const _pCrtanje;

    ///
    /// \brief updateCanvasAndBlock - azuriranje crteza i blokiranje dok se semafor ne oslobodi
    ///     1. Azuriranje crteza iz renderarea klase
    ///     2. Blokiranje semafora (semafor se oslobadja u tajmeru, tako da se na taj nacin
    ///         simulira animacija)
    /// \return true ako animacija treba da se zaustavi, false u suprotnom
    ///
    bool updateCanvasAndBlock();

    std::vector<QPoint> generisiNasumicneTacke(int brojTacaka = BROJ_SLUCAJNIH_OBJEKATA) const;
    std::vector<QPoint> ucitajPodatkeIzDatoteke(std::string imeDatoteke) const;

    /* Parametar naivnosti */
    const bool &_naivni;
};

#endif // ALGORITAMBAZA_H
