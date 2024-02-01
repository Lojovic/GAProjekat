#include "timemeasurementthread.h"

#include "config.h"
#include "mainwindow.h"
#include "algoritambaza.h"

/* Ovde ukluciti zaglavlja novih algoritama. */
#include "./algoritmi_sa_vezbi/ga00_demoiscrtavanja.h"
#include "./algoritmi_sa_vezbi/ga01_brisucaprava.h"
#include "./algoritmi_sa_vezbi/ga02_3discrtavanje.h"
#include "./algoritmi_sa_vezbi/ga03_konveksniomotac.h"
#include "./algoritmi_sa_vezbi/ga04_konveksniomotac3d.h"
#include "./algoritmi_sa_vezbi/ga05_preseciduzi.h"

//#include "ga06_presekPravougaonika.h"

TimeMeasurementThread::TimeMeasurementThread(TipAlgoritma tipAlgoritma, int minValue, int step, int maxValue)
    : QThread(), _algorithmType(tipAlgoritma), _minValue(minValue), _step(step), _maxValue(maxValue)
{
}

void TimeMeasurementThread::run()
{
    clock_t begin, end;
    double optimalTime, naiveTime;

    AlgoritamBaza *pAlgorithm = nullptr;

    /* _i_ je broj nasumicnih tacaka.
     * Kada se radi poredjenje, onda se instancira algoritam sa _i_ nasumicnih tacaka
     * i poredi se vreme izvrsavanja efikasnog i naivnog algoritma
     * nad tih _i_ nasumicnih tacaka.
     */
    for(int i= _minValue; i <= _maxValue; i += _step)
    {

        /* Ovde kreirati instancu klase algoritma. */
        switch (_algorithmType) {
        case TipAlgoritma::DEMO_ISCRTAVANJA:
            pAlgorithm = new DemoIscrtavanja(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::BRISUCA_PRAVA:
            pAlgorithm = new BrisucaPrava(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::_3D_ISCRTAVANJE:
            pAlgorithm = new BrisucaPrava(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::KONVEKSNI_OMOTAC:
            pAlgorithm = new KonveksniOmotac(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::KONVEKSNI_OMOTAC_3D:
            pAlgorithm = new KonveksniOmotac3D(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::PRESECI_DUZI:
            pAlgorithm = new PreseciDuzi(nullptr, 0, false, "", i);
            break;
        case TipAlgoritma::TRIANGULACIJA:
            pAlgorithm = new Triangulation(nullptr, 0, false, "", i);
            break;
       /* case TipAlgoritma::PRESEK_PRAVOUGAONIKA:
            pAlgorithm = new PresekPravougaonika(nullptr, 0, false, "", i);
            break;*/
        default:
            break;
        }

        if(pAlgorithm)
        {
#ifndef SKIP_OPTIMAL
            begin = clock();
            pAlgorithm->pokreniAlgoritam();
            end = clock();
            optimalTime = double(end - begin) / CLOCKS_PER_SEC;
#else
            optimalTime = 0;
#endif

#ifndef SKIP_NAIVE
            begin = clock();
            pAlgorithm->pokreniNaivniAlgoritam();
            end = clock();
            naiveTime = double(end - begin) / CLOCKS_PER_SEC;
#else
            naiveTime = 0;
#endif
            // std::cout << "Zavrsio oba poziva, stigao do crtanja" << std::endl;
            emit updateChart(i, optimalTime, naiveTime);
            delete pAlgorithm;
            pAlgorithm = nullptr;
        }
    }

    emit finishChart();
}
