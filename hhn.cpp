#include "hhn.h"
//#include "pinknumber.h"

using namespace std;

HHN::HHN()
{
    G_Na=G_Na/C_m;
    G_L=G_L/C_m;
    G_K=G_K/C_m;
    state=0;
}

//----------------------------------------------------------------------------------------
//--------------------------   HH Model v' Function   ----------------------------
//----------------------------------------------------------------------------------------
double HHN::fv(double v,double n,double m, double h,double I)
{
    double result=I-(G_Na*pow(m,3)*h*(v-v_Na))-(G_K*pow(n,4)*(v-v_K))-(G_L*(v-v_L));
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model alpha_n Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::alpha_n(double v)
{
    double result=(0.01*(10.0-v))/(exp(1.0 -(0.1*v))-1.0);
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model alpha_m Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::alpha_m(double v)
{
    double result=(0.1*(25.0-v))/(exp(2.5-(0.1*v)) -1.0);
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model alpha_h Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::alpha_h(double v)
{
    double result=0.07*exp((-1*v)/20.0);
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model beta_n  Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::beta_n(double v)
{
    double result=0.125*exp((-1*v)/80.0);
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model beta_m  Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::beta_m(double v)
{
    double result=4.0*exp((-1*v)/18.0);
    return result;
}
//----------------------------------------------------------------------------------------
//----------------------------   HH Model beta_h  Function   -----------------------------
//----------------------------------------------------------------------------------------
double HHN::beta_h(double v)
{
    double result=1.0/(exp(3.0-(0.1*v))+1.0);
    return result;
}
//----------------------------------------------------------------------------------------
//--------------------------   HH Model n' Function   ----------------------------
//----------------------------------------------------------------------------------------
double HHN::fn(double v,double n)
{
    double result=alpha_n(v)*(1-n)-beta_n(v)*n;
    return result;
}
//----------------------------------------------------------------------------------------
//--------------------------   HH Model m' Function   ----------------------------
//----------------------------------------------------------------------------------------
double HHN::fm(double v,double m)
{
    double result=alpha_m(v)*(1-m)-beta_m(v)*m;
    return result;
}
//----------------------------------------------------------------------------------------
//--------------------------   HH Model h' Function   ----------------------------
//----------------------------------------------------------------------------------------
double HHN::fh(double v,double h)
{
    double result=alpha_h(v)*(1-h)-beta_h(v)*h;
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   step 1 in runge kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double HHN::k1v(double v,double n,double m,double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fv(v,n,m,h,I);
    return result;
}
//
double HHN::k1n(double v, double n)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fn(v,n);
    return result;
}
//
double HHN::k1m(double v,double m)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fm(v,m);
    return result;
}
//
double HHN::k1h(double v, double h)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fh(v,h);
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 2 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double HHN::k2v(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fv(v+0.5*k1v(v,n,m,h,I),n+0.5*k1n(v,n),m+0.5*k1m(v,m),h+0.5*k1h(v,h),I);
    return result;
}
//
double HHN::k2n(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fn(v+0.5*k1v(v,n,m,h,I),n+0.5*k1n(v,n));
    return result;
}
//
double HHN::k2m(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fm(v+0.5*k1v(v,n,m,h,I),m+0.5*k1m(v,m));
    return result;
}
//
double HHN::k2h(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fh(v+0.5*k1v(v,n,m,h,I),h+0.5*k1h(v,h));
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 3 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double HHN::k3v(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fv(v+0.5*k2v(v,n,m,h,I),n+0.5*k2n(v,n,m,h,I),m+0.5*k2m(v,n,m,h,I),h+0.5*k2h(v,n,m,h,I),I);
    return result;
}
//
double HHN::k3n(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fn(v+0.5*k2v(v,n,m,h,I),n+0.5*k2n(v,n,m,h,I));
    return result;
}
//
double HHN::k3m(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fm(v+0.5*k2v(v,n,m,h,I),m+0.5*k2m(v,n,m,h,I));
    return result;
}
//
double HHN::k3h(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fh(v+0.5*k2v(v,n,m,h,I),h+0.5*k2h(v,n,m,h,I));
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 4 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double HHN::k4v(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fv(v+k3v(v,n,m,h,I),n+k3n(v,n,m,h,I),m+k3m(v,n,m,h,I),h+k3h(v,n,m,h,I),I);
    return result;
}
//
double HHN::k4n(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fn(v+k3v(v,n,m,h,I),n+k3n(v,n,m,h,I));
    return result;
}
//
double HHN::k4m(double v,double n,double m, double h,double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fm(v+k3v(v,n,m,h,I),m+k3m(v,n,m,h,I));
    return result;
}
//
double HHN::k4h(double v, double n, double m, double h, double I)//fr is result of f(x,...) function and h is time step length
{
    double result=dt*fh(v+k3v(v,n,m,h,I),h+k3h(v,n,m,h,I));
    return result;
}

//----------------------------------------------------------------------------------------
//-----------------------------   Runge Kutta Result   -----------------------------------
//----------------------------------------------------------------------------------------
double HHN::RKv(double v,double n,double m,double h,double I)
{
    double result=v+0.166666667*(k1v(v,n,m,h,I)+2*k2v(v,n,m,h,I)+2*k3v(v,n,m,h,I)+k4v(v,n,m,h,I));
    return result;
}
//
double HHN::RKn(double v,double n,double m,double h,double I)
{
    double result=n+0.166666667*(k1n(v,n)+2*k2n(v,n,m,h,I)+2*k3n(v,n,m,h,I)+k4n(v,n,m,h,I));
    return result;
}
//
double HHN::RKm(double v,double n,double m,double h,double I)
{
    double result=m+0.166666667*(k1m(v,m)+2*k2m(v,n,m,h,I)+2*k3m(v,n,m,h,I)+k4m(v,n,m,h,I));
    return result;
}
//
double HHN::RKh(double v,double n,double m,double h,double I)
{
    double result=h+0.166666667*(k1h(v,h)+2*k2h(v,n,m,h,I)+2*k3h(v,n,m,h,I)+k4h(v,n,m,h,I));
    return result;
}
//----------------------------------------------------------------------------------------
//-----------------------------   Initialize Neurons   -----------------------------------
//----------------------------------------------------------------------------------------
void HHN::Initialize()
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    poisson_distribution<> pdist(10);
    normal_distribution<double> distribution(10.0,6.0);
    uniform_real_distribution<double> udistribution(0.0,10.0);
    gamma_distribution<double> gdistribution(1,6);
    I=0;//pdist(gen);
    v=v_0;
    state=1;
}
//----------------------------------------------------------------------------------------
//-------------------------------   Connect Neurons   ------------------------------------
//----------------------------------------------------------------------------------------
void HHN::Connect(std::vector<int> &cv)
{
    for(int i=0;i<cv.size();++i)
    {
        outC.push_back(cv[i]);
    }
}

//----------------------------------------------------------------------------------------
//-----------------------------   Run Neural Dynamics   ----------------------------------
//----------------------------------------------------------------------------------------
void HHN::Run()
{
    if(state==1)
    {
        double tempV,tempN,tempM,tempH;
        if(v>=30)
        {
            Spike=1;
        }else
        {
            Spike=0;
        }
        I=I/C_m;
        tempV=RKv(v,n,m,h,I);//solve v equation RK4
        tempN=RKn(v,n,m,h,I);
        tempM=RKm(v,n,m,h,I);
        tempH=RKh(v,n,m,h,I);
        v=tempV;
        n=tempN;
        m=tempM;
        h=tempH;
    }else
    {
        cout<<"   ERROR! (initial value problem)  "<<endl;
        exit (EXIT_FAILURE);
    }
}
//----------------------------------------------------------------------------------------
//------------------------------   Get Random Current   ----------------------------------
//----------------------------------------------------------------------------------------
void HHN::UpdateI()
{
    //PinkNumber pink1;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    poisson_distribution<> pdist(10);
    normal_distribution<double> distribution(10.0,10.0);
    uniform_real_distribution<double> udistribution(0.0,1.0);
    gamma_distribution<double> gdistribution(1,6);
    //double D=outC.size();
    double rI=udistribution(gen);
    double n=-1.4500001;
    double x0=0.01;
    double x1=1000.0001;
    double plr = pow(((pow(x1,(n+1))-pow(x0,(n+1)))*rI+pow(x0,(n+1))),(1/(n+1)));//pink
    //cashI++;
    //double g=2.35;
//    if(cashI==Ti)
   // {
    //    I=I+(-80)*(gi)*;
    //    cashI=0;
   // }
    //cashI+0.00055;
    //double xI=10*sin(cashI);
    I=plr;//(xI+rI)/2;//udistribution(gen);//=rI+cashI/D;
    //cout<<plr<< endl;
    //cout<<I<<endl;
}
double HHN::get(string var_str)
{
    double result=0
    if(var_str=='v')
        result=v;
    if(var_str=='n')
        result=n;
    if(var_str=='m')
        result=m;
    if(var_str=='h')
        result=h;
    if(var_str=='I')
        result=I;
    return result;
}
