#include "selector_3g.h"
#include <TSystem.h>

void run()
{
    Selector3g s("/data/4/users/jsowa/data/Run10/MC/2025_10_16-23_08_22.ntu.root", "output.root", true);
    
    s.runAnalysis();
}
