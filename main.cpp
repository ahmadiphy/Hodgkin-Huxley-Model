#include "in.h"
#include "global_variables.h"
#include "hhn.h"
#include "matrixf.h"
#include "regular_net.h"
#include "ern.h"
#include <iomanip>
#include <chrono>
using namespace std;
double a,b,c,d,h;
double G_Na,G_K,G_L,C_m,v_K,v_Na,v_L,dt,v_0;

int main()
{
    //
    //omp_set_num_threads(3);
    //----------------------------------------------------------------------------------------
    //-----------------------------  core dump error checking --------------------------------
    //----------------------------------------------------------------------------------------
    struct rlimit core_limit;
    core_limit.rlim_cur = RLIM_INFINITY;
    core_limit.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_CORE, &core_limit) < 0) {
        /* ERROR */
    }
    //----------------------------------------------------------------------------------------
    //---------------------------- geting system time and date -------------------------------
    //----------------------------------------------------------------------------------------
    time_t now = time(0)+16200;
    char* dt = ctime(&now);
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);
    //----------------------------------------------------------------------------------------
    //----------------------------------  make directory -------------------------------------
    //----------------------------------------------------------------------------------------
    const int dir_err = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        cout<<"Error Creating Directory"<<endl;
        cout<<"-------------- No Problem ---------------"<<endl;
    }
    //----------------------------------------------------------------------------------------
    //------------------------------------ Make Network --------------------------------------
    //----------------------------------------------------------------------------------------
    Regular_net rg1;
    int l=9;
    int m0=2,ll=pow(2,l);
    int N=m0*ll;
    Matrixf mf1;
    iMatrix mat1(N,iRow());
    rg1.regT1(N,150,mat1);
    mf1.mat_shrink_to_fit(N,mat1);//behineh saziye ram
    //----------------------------------------------------------------------------------------
    //-------------------------------- info file write ---------------------------------------
    ostringstream info;
    info<<"./data/info.txt";
    ofstream infoPrint(info.str().c_str(),ios_base::binary);
    //infoPrint<<"."<< endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    infoPrint<<"------------- Simulation of HH Neurons on HMN --------------"<<endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    //----------------------------------------------------------------------------------------
    cout<<"The Runing Start date and time is: "<<dt <<endl;
    infoPrint<<"The Runing Start date and time is: "<<dt <<endl;
    //----------------------------------------------------------------------------------------
    //---------------------------- Initialize Global Variables -------------------------------
    //----------------------------------------------------------------------------------------
    G_Na=120.0;
    G_K=36.0;
    G_L=0.3;
    C_m=1.0;
    v_K=-12;
    v_Na=150;
    v_L=10.613;
    dt=0.001;
    v_0=0;
    //----------------------------------------------------------------------------------------
    //------------------------------------  One HH Neuron  -----------------------------------
    //----------------------------------------------------------------------------------------
    HHN hhn1;
    hhn1.Initialize();
    int loopTry=200000;
    ostringstream vData;
    vData<<"./data/data_v.dat";
    ofstream vPrint(vData.str().c_str(),ios_base::binary);
    for(int i=0;i<loopTry;++i)
    {
        hhn1.Run();
        hhn1.UpdateI();
        vPrint<<hhn1.get('v')<<" "<<hhn1.get('I')<<" "<<hhn1.get('n')<<" "<<hhn1.get('m')<<" "<<hhn1.get('h')<<endl;
    }
    //----------------------------------------------------------------------------------------
    //-------------------------------------  file writing  -----------------------------------
    //----------------------------------------------------------------------------------------
    ostringstream outData,vuData;
    outData<<"./data/data_activity.dat";
    vuData<<"./data/data_vu.dat";
    ofstream outPrint(outData.str().c_str(),ios_base::binary);
    ofstream vuPrint(vuData.str().c_str(),ios_base::binary);
    //----------------------------------------------------------------------------------------
    //------------------------------------  Time printing  -----------------------------------
    //----------------------------------------------------------------------------------------
    now = time(0)+16200;
    dt = ctime(&now);
    gmtm = gmtime(&now);
    dt = asctime(gmtm);
    cout<<"The Runing End   date and time is: "<<dt<< endl;
    infoPrint<<"The Runing End   date and time is: "<<dt<< endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    //----------------------------------------------------------------------------------------
}
