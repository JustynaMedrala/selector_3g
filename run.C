#include "selector_3g.h"
#include <TSystem.h>

void run()
{
    gROOT->ProcessLine(".L selector_3g.C++g");    

    // Tworzymy obiekt Selector3g
    Selector3g s("/data/4/users/jsowa/data/Run10/MC/out_1.root", "output.root", true);

    // Uruchamiamy analizę
    s.runAnalysis();
}
