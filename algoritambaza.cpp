#include "algoritambaza.h"
#include <fstream>

void AlgoritamBaza::timerEvent(QTimerEvent */* unused */)
{
     _semafor.release();
}

AlgoritamBaza::AlgoritamBaza(QWidget *pCrtanje, int pauzaKoraka, const bool &naivni)
    : QObject{}, _pauzaKoraka{pauzaKoraka},
      _timerId{INVALID_TIMER_ID}, _semafor(0),
      _unistiAnimaciju(false), _pNit(nullptr),
      _pCrtanje{pCrtanje}, _naivni(naivni)
{}

bool AlgoritamBaza::naivni() const
{
    return _naivni;
}

void AlgoritamBaza::crtaj(QPainter *painter) const
{
    if (_naivni)
        crtajNaivniAlgoritam(painter);
    else
        crtajAlgoritam(painter);
}

void AlgoritamBaza::pokreniAnimaciju()
{
    delete _pNit;
    _pNit = new AnimacijaNit(this);
    _pNit->start();
    _timerId = startTimer(_pauzaKoraka);
    _pCrtanje->update();
}

void AlgoritamBaza::pauzirajIliNastaviAnimaciju()
{
    if (INVALID_TIMER_ID != _timerId)
    {
        killTimer(_timerId);
        _timerId = INVALID_TIMER_ID;
    }
    else
        _timerId = startTimer(_pauzaKoraka);
}

void AlgoritamBaza::sledeciKorakAnimacije()
{
    if (INVALID_TIMER_ID != _timerId)
    {
        killTimer(_timerId);
        _timerId = INVALID_TIMER_ID;
    }

    _semafor.release();
}

void AlgoritamBaza::zaustaviAnimaciju()
{
    _unistiAnimaciju = true;
    _semafor.release();
    _pNit->quit();
    _pNit->wait();

    delete _pNit;
    _pNit = nullptr;
    _unistiAnimaciju = false;
}

void AlgoritamBaza::promeniDuzinuPauze(int duzinaPauze)
{
    if (INVALID_TIMER_ID != _timerId)
        killTimer(_timerId);
    _pauzaKoraka = duzinaPauze;
    _timerId = startTimer(_pauzaKoraka);
}

bool AlgoritamBaza::updateCanvasAndBlock()
{
    if (_pCrtanje) {
        _pCrtanje->update();
        _semafor.acquire();
    }
    return _unistiAnimaciju;
}

std::vector<QPoint> AlgoritamBaza::generisiNasumicneTacke(int brojTacaka) const
{
    srand(static_cast<unsigned>(time(nullptr)));
    int xMax;
    int yMax;

    if (_pCrtanje)
    {
        xMax = _pCrtanje->width() - DRAWING_BORDER;
        yMax = _pCrtanje->height() - DRAWING_BORDER;
    }
    else
    {
        xMax = CANVAS_WIDTH;
        yMax = CANVAS_HEIGHT;
    }

    int xMin = DRAWING_BORDER;
    int yMin = DRAWING_BORDER;

    std::vector<QPoint> randomPoints;

    int xDiff = xMax-xMin;
    int yDiff = yMax-yMin;
    for(int i=0; i < brojTacaka; i++)
        randomPoints.emplace_back(xMin + rand()%xDiff, yMin + rand()%yDiff);

    return randomPoints;
}

std::vector<QPoint> AlgoritamBaza::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QPoint> points;
    int x, y;
    while(inputFile >> x >> y)
        points.emplace_back(x, y);
    return points;
}
