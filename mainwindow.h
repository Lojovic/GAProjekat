#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "algoritambaza.h"
#include "oblastcrtanja.h"
#include "oblastcrtanjaopengl.h"
#include "config.h"
#include "timemeasurementthread.h"

/* Ovde ukljuciti zaglavlja novih algoritma. */
#include "./algoritmi_sa_vezbi/ga00_demoiscrtavanja.h"
#include "./algoritmi_sa_vezbi/ga01_brisucaprava.h"
#include "./algoritmi_sa_vezbi/ga02_3discrtavanje.h"
#include "./algoritmi_sa_vezbi/ga03_konveksniomotac.h"
#include "./algoritmi_sa_vezbi/ga04_konveksniomotac3d.h"
#include "./algoritmi_sa_vezbi/ga05_preseciduzi.h"
#include "./algoritmi_sa_vezbi/ga06_dceldemo.h"
#include "./algoritmi_sa_vezbi/ga07_triangulation.h"

//#include "ga06_presekPravougaonika.h"
#include "./algoritmi_studentski_projekti/shortcutHull.h"

/* Enumeracija algoritama */
enum class TipAlgoritma {
    ALGORITMI_SA_VEZBI,
    DEMO_ISCRTAVANJA,
    BRISUCA_PRAVA,
    _3D_ISCRTAVANJE,
    KONVEKSNI_OMOTAC,
    KONVEKSNI_OMOTAC_3D,
    PRESECI_DUZI,
    DCEL_DEMO,
    TRIANGULACIJA,
    SEPARATOR,
    STUDENTSKI_PROJEKTI,
    SHORTCUT_HULL
};

/* Enumeracija tabova */
enum TabIndex {
    ALGORITAM_2D,
    ALGORITAM_3D,
    POREDJENJE
};

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    virtual ~MainWindow() override;

private slots:
    void on_datoteka_dugme_clicked();

    void on_Nasumicni_dugme_clicked();

    void on_Ponisti_dugme_clicked();

    void on_Zapocni_dugme_clicked();

    void on_Zaustavi_dugme_clicked();

    void on_Sledeci_dugme_clicked();

    void on_Ispocetka_dugme_clicked();

    void on_polozajKursora(int x, int y);

    /* za Chart */
    void on_merenjeButton_clicked();
    void on_lineSeriesChange(double dim, double optimal, double naive);
    void on_chartFinished();

    void on_tipAlgoritma_currentIndexChanged(int index);

    void on_naivniCheck_stateChanged(int);

    void na_krajuAnimacije();


private:
    void animacijaButtonAktivni(bool param_aktivnosti);

    void animacijaParametriButtonAktivni(bool param_aktivnosti);

    void napraviNoviAlgoritam();

private:
    Ui::MainWindow *ui;

    AlgoritamBaza *_pAlgoritamBaza;
    OblastCrtanja *_pOblastCrtanja;
    OblastCrtanjaOpenGL *_pOblastCrtanjaOpenGL;
    std::string _imeDatoteke;
    bool _naivni;
    int _duzinaPauze;
    int _brojSlucajnihObjekata;

    /* Chart deo */
    QLineSeries *const _naiveSeries = new QLineSeries();
    QLineSeries *const _optimalSeries = new QLineSeries();

    TimeMeasurementThread *_mThread;
};

#endif // MAINWINDOW_H
