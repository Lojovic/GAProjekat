#include "ga02_3discrtavanje.h"

#include <QVector3D>
#include <fstream>

Discrtavanje::Discrtavanje(QWidget *pCrtanje,
                           int pauzaKoraka,
                           const bool &naivni,
                           std::string imeDatoteke,
                           int brojTacaka)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni), _n(0)
{
    if (imeDatoteke != "")
        _tacke = ucitajPodatkeIzDatoteke(imeDatoteke);
    else
        _tacke = generisiNasumicneTacke(brojTacaka);
}

void Discrtavanje::pokreniAlgoritam()
{
    while (_n != _tacke.size())
    {
        _n++;
        AlgoritamBaza_updateCanvasAndBlock()
    }

    emit animacijaZavrsila();
}

void Discrtavanje::crtajAlgoritam(QPainter *) const
{
    glBegin(GL_TRIANGLE_STRIP);
        for(auto i = 0u; i < _n; i++)
        {
            double drand = 3.*rand()/RAND_MAX;
            if (drand < 1)
                glColor3d(1, 0, 0);
            else if (drand < 2)
                glColor3d(0, 1, 0);
            else
                glColor3d(0, 0, 1);
            glVertex3f(_tacke[i].x(), _tacke[i].y(), _tacke[i].z());
        }
    glEnd();
}

void Discrtavanje::pokreniNaivniAlgoritam()
{
    emit animacijaZavrsila();
}

void Discrtavanje::crtajNaivniAlgoritam(QPainter *) const
{

}

std::vector<QVector3D> Discrtavanje::generisiNasumicneTacke(int brojTacaka) const
{
    srand(static_cast<unsigned>(time(nullptr)));

    std::vector<QVector3D> randomPoints;

    for(int i=0; i < brojTacaka; i++)
        randomPoints.emplace_back(
           1.f*rand()/RAND_MAX,
           1.f*rand()/RAND_MAX,
           1.f*rand()/RAND_MAX);

    return randomPoints;
}

std::vector<QVector3D> Discrtavanje::ucitajPodatkeIzDatoteke(std::string imeDatoteke) const
{
    std::ifstream inputFile(imeDatoteke);
    std::vector<QVector3D> points;
    float x, y, z;
    while(inputFile >> x >> y >> z)
    {
        points.emplace_back(x, y, z);
    }
    return points;
}
