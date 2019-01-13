#include "ern.h"
//
using namespace::std;
//
void ERn::ERnetwork(double avgK, int nn, iMatrix &aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,nn-1);
    double linkss;
    cout<<"ok 1"<<endl;
    linkss=double((avgK*nn)/2);
    cout<<"ok 2"<<endl;
    const int links=int(1*linkss);
    int i=0;
    int r1=0,r2=0;
    cout<<i<<" "<<links<<endl;
    do{
        cout<<"****"<<i<<" ";
        r1=0;
        r2=0;
        r1=distall(gen);
        r2=distall(gen);
        if(aa[r1][r2]==0 && r1 != r2)
        {
            aa[r1][r2]=1;
            aa[r2][r1]=1;
            ++i;
        }
    }while (i<links);
    cout<<"ok 3 "<<links<<endl;
}
//
void ERn::ERnetworkEff(double avgK, int nn, iMatrix &aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_int_distribution<> distall(0,nn-1);
    double linkss;
    linkss=double((avgK*nn)/2);
    int links;
    links=int(linkss);
    int i=0;
    while(i<links){
        int r1=0,r2=0;
        r1=distall(gen);
        r2=distall(gen);
        if(aa[r1][r2]==0 && r1 != r2)
        {
            aa[r1][r2]=1;
            aa[r2][r1]=1;
            ++i;
        }
    }
}
//
void ERn::PCNN(int Ne, int Ni, dMatrix & aa)
{
    int N=Ni+Ne;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist1_1(0,1);
    uniform_real_distribution<> dist1_2(0,1);
    for(int i=0;i<N;++i)
    {
        for(int j=0;j<Ne;++j)
        {
            double randN=dist1_1(gen);
            aa[i][j]=0.5*randN;
        }
        for(int j=Ne;j<N;++j)
        {
            double randN=dist1_2(gen);
            aa[i][j]=-1*randN;
        }
    }
}
//
void ERn::SER2EFFER(int N, iMatrix &aa, iMatrix &bb)
{
    for(int i=0;i<N;++i)
    {
        for(int j=0;j<N;++j)
        {
            if(aa[i][j]==1)
                bb[i].push_back(j);
        }
    }
}
