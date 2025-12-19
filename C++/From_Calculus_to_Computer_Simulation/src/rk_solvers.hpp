#pragma once
#include <vector>
#include <functional>
#include <cmath>

class RK45 {
public:
    double h, tol;

    RK45(double h_init=0.01, double tol_init=1e-6)
        : h(h_init), tol(tol_init) {}

    std::vector<double> step(
        std::function<std::vector<double>(double,std::vector<double>)> f,
        double t, std::vector<double> y)
    {
        int n=y.size();
        std::vector<double> k1,k2,k3,k4,k5,k6,temp(n),y4(n),y5(n);

        k1 = f(t,y);

        for(int i=0;i<n;i++) temp[i]=y[i] + h*(1.0/4*k1[i]);
        k2 = f(t+1.0/4*h,temp);

        for(int i=0;i<n;i++) temp[i]=y[i] + h*(3.0/32*k1[i] + 9.0/32*k2[i]);
        k3 = f(t+3.0/8*h,temp);

        for(int i=0;i<n;i++) temp[i]=y[i] + h*(1932.0/2197*k1[i]-7200.0/2197*k2[i]+7296.0/2197*k3[i]);
        k4 = f(t+12.0/13*h,temp);

        for(int i=0;i<n;i++) temp[i]=y[i] + h*(439.0/216*k1[i]-8*k2[i]+3680.0/513*k3[i]-845.0/4104*k4[i]);
        k5 = f(t+h,temp);

        for(int i=0;i<n;i++) temp[i]=y[i] + h*(-8.0/27*k1[i]+2*k2[i]-3544.0/2565*k3[i]+1859.0/4104*k4[i]-11.0/40*k5[i]);
        k6 = f(t+1.0/2*h,temp);

        for(int i=0;i<n;i++)
            y4[i] = y[i] + h*(25.0/216*k1[i]+1408.0/2565*k3[i]+2197.0/4104*k4[i]-1.0/5*k5[i]);

        for(int i=0;i<n;i++)
            y5[i] = y[i] + h*(16.0/135*k1[i]+6656.0/12825*k3[i]+28561.0/56430*k4[i]-9.0/50*k5[i]+2.0/55*k6[i]);

        double err=0;
        for(int i=0;i<n;i++) err+=pow(y5[i]-y4[i],2);
        err=sqrt(err);

        double s=pow(tol/(err+1e-12),0.25);
        h *= std::min(4.0,std::max(0.1,0.9*s));

        return (err<tol? y5:y);
    }

    std::vector<double> solve(
        std::function<std::vector<double>(double,std::vector<double>)> f,
        double t0, std::vector<double> y0, double tf)
    {
        double t=t0;
        std::vector<double> y=y0;

        while(t<tf){
            auto ynew = step(f,t,y);
            if(ynew!=y){t+=h; y=ynew;}
        }
        return y;
    }
};


//----------------------------------- RK78 ----------------------------------//

class RK78 {
    public:
        double h, tol;
    
        RK78(double h_init=0.01, double tol_init=1e-9)
            : h(h_init), tol(tol_init) {}
    
        // ---- RK78 step returns {} if step rejected ----
        std::vector<double> step(
            std::function<std::vector<double>(double,std::vector<double>)> f,
            double t,
            const std::vector<double>& y,
            double &err_out)
        {
            int n=y.size();
            std::vector<std::vector<double>> k(13,std::vector<double>(n));
            std::vector<double> yt(n), y7(n), y8(n);
    
            // Coef: Fehlberg RK78 (corrigido)
            static double a[] = {
                0,2.0/27,1.0/9,1.0/6,5.0/12,1.0/2,5.0/6,
                1.0/6,2.0/3,1.0/3,1.0,0,1.0
            };
    
            static double b[13][12] = {
                {0},
                {2.0/27},
                {1.0/36,1.0/12},
                {1.0/24,0,1.0/8},
                {5.0/12,0,-25.0/16,25.0/16},
                {1.0/20,0,0,1.0/4,1.0/5},
                {-25.0/108,0,0,125.0/108,-65.0/27,125.0/54},
                {31.0/300,0,0,0,61.0/225,-2.0/9,13.0/900},
                {2.0,0,0,-53.0/6,704.0/45,-107.0/9,67.0/90,3.0},
                {-91.0/108,0,0,23.0/108,-976.0/135,311.0/54,-19.0/60,17.0/6,-1.0/12},
                {2383.0/4100,0,0,-341.0/164,4496.0/1025,-301.0/820,
                 2133.0/4100,45.0/82,45.0/164,18.0/41},
                {3.0/205,0,0,0,0,-6.0/41,-3.0/205,-3.0/41,3.0/41,6.0/41,0},
                {-1777.0/4100,0,0,-341.0/164,4496.0/1025,-289.0/820,
                 2193.0/4100,51.0/82,33.0/164,12.0/41,0,1}
            };
    
            static double c7[]={41.0/840,0,0,0,0,34.0/105,9.0/35,9.0/35,
                                9.0/280,9.0/280,0,0,41.0/840};
    
            static double c8[]={0,0,0,0,0,34.0/105,9.0/35,9.0/35,
                                9.0/280,9.0/280,41.0/840,41.0/840,0};
    
            // ------- compute k-stages --------
            for(int s=0;s<13;s++){
                for(int i=0;i<n;i++){
                    yt[i]=y[i];
                    for(int j=0;j<s;j++) yt[i]+=h*b[s][j]*k[j][i];
                }
                k[s]=f(t+a[s]*h,yt);
            }
    
            // ------- build y7, y8 --------
            double err=0;
            for(int i=0;i<n;i++){
                y7[i]=y[i]; y8[i]=y[i];
                for(int s=0;s<13;s++){
                    y7[i]+=h*c7[s]*k[s][i];
                    y8[i]+=h*c8[s]*k[s][i];
                }
                err+=(y8[i]-y7[i])*(y8[i]-y7[i]);
            }
            err_out=sqrt(err);
    
            // adapt step-size
            double s=pow(tol/(err+1e-14),1.0/8.0);
            if(s<0.1) s=0.1; if(s>4.0) s=4.0;
            h*=s;
    
            return (err<tol)?y8:std::vector<double>(); // return {} = rejected
        }
    
        // ===== SOLVER NOW INCLUDED =====
        std::vector<double> solve(
            std::function<std::vector<double>(double,std::vector<double>)> f,
            double t0,std::vector<double> y0,double tf)
        {
            double t=t0;
            auto y=y0;
    
            while(t<tf){
                double err;
                auto ytrial=step(f,t,y,err);
    
                if(!ytrial.empty()){ // accepted
                    y=ytrial;
                    t+=h;
                }
            }
            return y;
        }
    };