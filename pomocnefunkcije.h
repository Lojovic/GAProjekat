#ifndef POMOCNEFUNKCIJE_H
#define POMOCNEFUNKCIJE_H

#include <QPoint>
#include <QLineF>
#include <QVector3D>
#include <cmath>

#include <config.h>

namespace pomocneFunkcije {
int povrsinaTrougla(const QPoint& A, const QPoint& B, const QPoint& C);
int distanceKvadrat(const QPoint& A, const QPoint& B);
bool presekDuzi(const QLineF& l1, const QLineF& l2, QPointF& presek);
double distanceKvadratF(const QPointF& A, const QPointF& B);
bool ispod(const QPointF& A, const QPointF& B);
bool konveksan(const QPointF& A, const QPointF& B, const QPointF& C);

/* Pomocni metodi za 3D. */
bool kolinearne3D(const QVector3D& a, const QVector3D& b, const QVector3D& c);
float zapremina(const QVector3D& a, const QVector3D& b, const QVector3D& c, const QVector3D& d);
}

#endif // POMOCNEFUNKCIJE_H
