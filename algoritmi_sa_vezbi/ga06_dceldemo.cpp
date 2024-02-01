#include "ga06_dceldemo.h"

DCELDemo::DCELDemo(QWidget *pCrtanje,
                   int pauzaKoraka,
                   const bool &naivni,
                   std::string imeDatoteke,
                   int) /* brojTacaka, odnosno nasumicno ne koristimo u ovom primeru */
    :AlgoritamBaza(pCrtanje, pauzaKoraka, naivni),
     _polygon(imeDatoteke, pCrtanje->height(), pCrtanje->width())
{
}

void DCELDemo::pokreniAlgoritam() {
    emit animacijaZavrsila();
}

void DCELDemo::crtajAlgoritam(QPainter *painter) const {
    if (!painter) return;

    QPen pen = painter->pen();
    for(auto i=0ul; i<_polygon.edges().size(); i++)
    {
        /* Crta se poligon */
        pen.setColor(Qt::yellow);
        painter->setPen(pen);
        painter->drawLine(_polygon.edge(i)->origin()->coordinates(),
                          _polygon.edge(i)->next()->origin()->coordinates());

        /* Crta se poligon "u suprotnom smeru", koriscenjem twin. Da bi se video efekat
         * na crtezu, koordinate su za malo pomerene u odnosu na originalnu tacku. */

        pen.setColor(Qt::red);
        painter->setPen(pen);
        painter->drawLine(QPointF(_polygon.edge(i)->twin()->origin()->x() + 5,
                                  _polygon.edge(i)->twin()->origin()->y() + 5),
                          QPointF(_polygon.edge(i)->twin()->next()->origin()->x() - 5,
                                  _polygon.edge(i)->twin()->next()->origin()->y() + 5));
   }

    int curr_num = 0;
    painter->setBrush(Qt::red);
    painter->setPen(Qt::white);
    /* Crtaju se temena, ali kao elipsa, radi lepote. */
    for(Vertex* v: _polygon.vertices())
    {
       painter->drawEllipse(v->coordinates(), 10, 10);

       /* Okretanje cetkice kako brojevi ne bi bili obrnuti */
       painter->save();
       painter->scale(1, -1);
       painter->translate(0, -2*v->y());

       painter->drawText(QPointF(v->x() - 4, v->y() + 4),
                         QString::number(curr_num));
       curr_num++;

       /* Ponistavanje transformacija */
       painter->restore();
    }

}

void DCELDemo::pokreniNaivniAlgoritam() {
    emit animacijaZavrsila();
}

void DCELDemo::crtajNaivniAlgoritam(QPainter *painter) const
{
    if (!painter) return;
}
