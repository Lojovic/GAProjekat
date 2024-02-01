#ifndef TIMEMEASUREMENTTHREAD_H
#define TIMEMEASUREMENTTHREAD_H

#include <QThread>

/* QChart */
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>

QT_CHARTS_USE_NAMESPACE

#include "algoritambaza.h"
enum class TipAlgoritma;

class TimeMeasurementThread : public QThread
{
    Q_OBJECT

public:
    TimeMeasurementThread(TipAlgoritma algorithmType, int minValue, int step, int maxValue);

    // QThread interface
    void run() override;

signals:
    void updateChart(double dim, double optimal, double naive);
    void finishChart();

private:
    const TipAlgoritma _algorithmType;
    const int _minValue;
    const int _step;
    const int _maxValue;
};

#endif // TIMEMEASUREMENTTHREAD_H
