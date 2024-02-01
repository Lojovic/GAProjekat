#include "pomocnefunkcije.h"

#include <QtGlobal>

int pomocneFunkcije::povrsinaTrougla(const QPoint& A, const QPoint& B, const QPoint& C)
{
    /* (Dvostruka) Povrsina trougla.
     * 2P(Trougla) = |ax ay 1|
     *               |bx by 1|
     *               |cx cy 1|
     */
    return (B.x() - A.x())*(C.y() - A.y()) - (C.x() - A.x())*(B.y() - A.y());
}

int pomocneFunkcije::distanceKvadrat(const QPoint& A, const QPoint& B)
{
    return (A.x() - B.x())*(A.x() - B.x()) + (A.y() - B.y())*(A.y() - B.y());
}

bool pomocneFunkcije::kolinearne3D(const QVector3D& a, const QVector3D& b, const QVector3D& c)
{
    /* Proverava se povrsina trougla. Zaparvo da li su sva tri unakrsna proizvoda
     * vektora jednaki nuli.
     * Kako je u pitanju float vrednost, gledamo da li je apsultna vrednost manja od
     * nekog malog broja.
     * Tj. proveravamo sledece
     * |  i      j      k  |
     * |ax-bx  ay-by  az-bz| = (0, 0, 0)
     * |ax-cx  ay-cy  az-cz|
     */
    return (fabsf((c.z() - a.z()) * (b.y() - a.y()) -
              (b.z() - a.z()) * (c.y() - a.y())) < EPSf) &&
           (fabsf((b.z() - a.z()) * (c.x() - a.x()) -
              (b.x() - a.x()) * (c.z() - a.z())) < EPSf) &&
           (fabsf((b.x() - a.x()) * (c.y() - a.y()) -
              (b.y() - a.y()) * (c.x() - a.x())) < EPSf);
}


float pomocneFunkcije::zapremina(const QVector3D& a, const QVector3D& b, const QVector3D& c, const QVector3D& d)
{
    /* Zapremina se racuna kao determinanta
     * 6V(Tetraedar) = |ax ay az 1|
     *                 |bx by bz 1|
     *                 |cx cy cz 1|
     *                 |dx dy dz 1|
     */
    float bxdx = b.x() - d.x();
    float bydy = b.y() - d.y();
    float bzdz = b.z() - d.z();
    float cxdx = c.x() - d.x();
    float cydy = c.y() - d.y();
    float czdz = c.z() - d.z();
    float vol =    (a.z() - d.z()) * (bxdx*cydy - bydy*cxdx)
                +  (a.y() - d.y()) * (bzdz*cxdx - bxdx*czdz)
                +  (a.x() - d.x()) * (bydy*czdz - bzdz*cydy);

    if (fabsf(vol) < EPSf)
        return 0;
    else
        return vol;
}

bool pomocneFunkcije::presekDuzi(const QLineF& l1, const QLineF& l2, QPointF& presek)
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    return l1.intersects(l2, &presek) == QLineF::BoundedIntersection;
#else
    return l1.intersect(l2, &presek) == QLineF::BoundedIntersection;
#endif
}

double pomocneFunkcije::distanceKvadratF(const QPointF& A, const QPointF& B)
{
    return (A.x() - B.x())*(A.x() - B.x()) + (A.y() - B.y())*(A.y() - B.y());
}

// da li je A ispod B
bool pomocneFunkcije::ispod(const QPointF &A, const QPointF &B)
{
    if(A.y() < B.y())
        return true;
    else if(fabs(A.y() - B.y()) < EPS)
    {
       if(A.x() > B.x()) return true;
    }
    return false;
}

bool pomocneFunkcije::konveksan(const QPointF &A, const QPointF &B, const QPointF &C)
{
    double P = (B.x() - A.x())*(C.y() - A.y()) - (C.x() - A.x())*(B.y() - A.y());

     return (P > 0) ||
            (fabs(P) < EPS && pomocneFunkcije::distanceKvadratF(A, B)
                            < pomocneFunkcije::distanceKvadratF(A, C));
}
