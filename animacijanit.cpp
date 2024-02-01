#include "animacijanit.h"
#include "algoritambaza.h"

AnimacijaNit::AnimacijaNit(AlgoritamBaza *pAlgoritam)
    :QThread{}, _pAlgoritam(pAlgoritam)
{}

void AnimacijaNit::run()
{
    if (!_pAlgoritam) return;

    if (_pAlgoritam->naivni())
        _pAlgoritam->pokreniNaivniAlgoritam();
    else
        _pAlgoritam->pokreniAlgoritam();
}
