#include "in.h"
#include "global_variables.h"
#ifndef HHN_H
#define HHN_H


class HHN
{
private:
    double Ti=200;
    int state;//chek parameter for initialize
    int Spike,myNum;
    double v,n,m,h,I,cashI;//cashI is sum of
    std::vector<int> inC; //for non symetric connections
    std::vector<double> inW; //for non symetric connections
    std::vector<int> outC;
    std::vector<double> outW;
public:
    HHN();
    double fv(double v,double n,double m,double h, double I);
    double alpha_n(double v);
    double alpha_m(double v);
    double alpha_h(double v);
    double beta_n(double v);
    double beta_m(double v);
    double beta_h(double v);
    double fn(double v,double n);
    double fm(double v,double m);
    double fh(double v,double h);

    double k1v(double v, double n, double m, double h, double I);
    double k1n(double v, double n);
    double k1m(double v, double m);
    double k1h(double v, double h);

    double k2v(double v, double n, double m, double h, double I);
    double k2n(double v, double n, double m, double h, double I);
    double k2m(double v, double n, double m, double h, double I);
    double k2h(double v, double n, double m, double h, double I);

    double k3v(double v, double n, double m, double h, double I);
    double k3n(double v, double n, double m, double h, double I);
    double k3m(double v, double n, double m, double h, double I);
    double k3h(double v, double n, double m, double h, double I);

    double k4v(double v, double n, double m, double h, double I);
    double k4n(double v, double n, double m, double h, double I);
    double k4m(double v, double n, double m, double h, double I);
    double k4h(double v, double n, double m, double h, double I);

    double RKv(double v, double n, double m, double h, double I);
    double RKn(double v, double n, double m, double h, double I);
    double RKm(double v, double n, double m, double h, double I);
    double RKh(double v, double n, double m, double h, double I);

    void Initialize();// inh is time step
    void Connect(std::vector<int>  & cv);
    void UpdateI();
    void Run();
    double get(char var_str);
};

#endif // HHN_H
