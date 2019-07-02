#include "in.h"
#ifndef ERN_H
#define ERN_H

class ERn
{
public:
    void ERnetwork(double avgK,int nn,iMatrix& aa);
    void ERnetworkEff(double avgK,int nn,iMatrix& aa);
    void PCNN(int Ne,int Ni,dMatrix & aa);
    void SER2EFFER(int N,iMatrix &aa,iMatrix &bb);
};

#endif // ERN_H
