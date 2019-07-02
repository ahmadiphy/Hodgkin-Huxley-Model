#include "regular_net.h"
// looped connections
void Regular_net::regT1(int N,int M,iMatrix &a)
{
    for(int i=0;i<N;++i)
    {
        for(int j=1;j<=M;++j)
        {
            int k=i+j;
            if(i+j>=N)
            {
                k=k-N;
            }
            a[i].push_back(k);
            a[k].push_back(i);
        }
    }
}
