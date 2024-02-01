#ifndef ANIMACIJANIT_H
#define ANIMACIJANIT_H

#include <QThread>

class AlgoritamBaza;

class AnimacijaNit : public QThread
{
public:
    ///
    /// \brief AnimationThread
    /// \param pAlgoritam - algoritam koji se izvrsava unutar thread-a
    ///
    AnimacijaNit(AlgoritamBaza *pAlgoritam);

    ///
    /// \brief run - metod koji se izvrsava kada se thread pokrene
    ///
    virtual void run() override;

private:
    AlgoritamBaza *const _pAlgoritam;
};

#endif // ANIMACIJANIT_H
