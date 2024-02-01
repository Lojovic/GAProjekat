#ifndef OBLASTCRTANJA_H
#define OBLASTCRTANJA_H

#include <QWidget>

#include "algoritambaza.h"

class OblastCrtanja : public QWidget
{
    Q_OBJECT
public:
    explicit OblastCrtanja(QWidget *parent = nullptr);
    void postaviAlgoritamKojiSeIzvrsava(AlgoritamBaza *pAlgoritamBaza);

    void set_obrisiSve(bool param);

protected:
    ///
    /// \brief paintEvent -
    ///         Funkcija (preciznije event handler) u kojoj je implementirano
    ///         iscrtavanje po nasoj oblasti za crtanje.
    ///
    void paintEvent(QPaintEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;

private:
    AlgoritamBaza *_pAlgoritamBaza;
    bool _obrisiSve;

signals:
    void polozajKursora(int x, int y);
};

#endif // OBLASTCRTANJA_H
