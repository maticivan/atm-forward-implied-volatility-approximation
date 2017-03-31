
//*************************************************************************************************
//*************************************************************************************************
//* The MIT License (MIT)                                                                         *
//* Copyright (C) 2017 Ivan Matic, Rados Radoicic, and Dan Stefanica                              *
//*                                                                                               *
//* Permission is hereby granted, free of charge, to any person obtaining a copy of this          *
//* software and associated documentation files (the "Software"), to deal in the Software         *
//* without restriction, including without limitation the rights to use, copy, modify, merge,     *
//* publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons    *
//* to whom the Software is furnished to do so, subject to the following conditions:              *
//*                                                                                               *
//* The above copyright notice and this permission notice shall be included in all copies or      *
//* substantial portions of the Software.                                                         *
//*                                                                                               *
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,           *
//* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR      *
//* PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE     *
//* FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR          *
//* OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        *
//* DEALINGS IN THE SOFTWARE.                                                                     *
//*************************************************************************************************




#include <iostream>
#include <cmath>
#include <iomanip>

typedef long int myint;
typedef double mydouble;


mydouble cubeRoot(const mydouble &x){
    mydouble absX=x;
    mydouble signX=1.0;
    if(absX<0.0){
        absX *= (-1.0);
        signX = -1.0;
    }
    return signX * (exp(log(absX)/3.0) );
}
mydouble sigmaImpApprox(const mydouble &Cm,
                        const mydouble &q,
                        const mydouble &r,
                        const mydouble &S0,
                        const mydouble &T){
    
    //*****************************************************
    //* Input:                                            *
    //* Cm = market price of the call option              *
    //* q  = constant rate of dividends                   *
    //* r  = constant interest rate                       *
    //* S0 = price of the underlying security             *
    //* T  = option maturity                              *
    //*                                                   *
    //* Output:                                           *
    //* ATM-forward implied volatility approximation      *
    //*****************************************************

    mydouble sigma;
    mydouble CmTilde= Cm/(S0*exp(-q*T));
    mydouble pi=3.141592653589793238;
    mydouble piSquare=pi*pi;
    mydouble piCube=piSquare *pi;
    mydouble gamma2= -1.0/3.0+1.0/pi;
    mydouble gamma4= 7.0/90.0 -2.0/(3.0*pi)+ 4.0/(3.0*piSquare);
    mydouble gamma6= -1.0/70.0 + 4.0/(15.0*pi) - 4.0/(3.0*piSquare) + 2.0/(piCube);
    mydouble piHalfLog= 0.5*pi* log(1.0- CmTilde*CmTilde);
    
    if(CmTilde<=0.65949631){
        mydouble gamma4OverGamma6= gamma4/gamma6;
        mydouble gamma4OverGamma6Square=gamma4OverGamma6 * gamma4OverGamma6;
        mydouble gamma4OverGamma6Cube= gamma4OverGamma6Square *gamma4OverGamma6;
        mydouble u= gamma2/gamma6- 0.375 * gamma4OverGamma6Square;
        mydouble v= 0.125* gamma4OverGamma6Cube - 0.5* (gamma2/gamma6) * gamma4OverGamma6+1.0/gamma6;
        mydouble w= piHalfLog/gamma6 - 0.01171875*gamma4OverGamma6*gamma4OverGamma6Cube - 0.25 * gamma4OverGamma6/gamma6 + 0.0625*(gamma2/gamma6)*gamma4OverGamma6Square;
        mydouble P= (-12.0*w-u*u)/12.0;
        mydouble Q= (u*(-u*u+36.0*w))/(108.0)-0.125*v*v;
        mydouble D= sqrt(0.25*Q*Q+ (P*P*P/27.0));
        
        
        mydouble m=cubeRoot(-0.5*Q+D)  + cubeRoot(-0.5*Q-D) -u/3.0;
        sigma= (2.0/ (sqrt(T))) * sqrt(0.5*(- sqrt(2*m)  + sqrt(2*v/(sqrt(2*m)) - 2*(u+m)) ) - 0.25*gamma4OverGamma6);
        
    }
    if((CmTilde>0.65949631)&&(CmTilde<0.99999971)){
        mydouble gamma4Square=gamma4*gamma4;
        mydouble P= (3.0*gamma4 - gamma2*gamma2)/ (3.0*gamma4Square);
        mydouble Q= (2.0*gamma2 *gamma2*gamma2-9.0*gamma2*gamma4)/(27.0*gamma4Square * gamma4)+ piHalfLog/gamma4;
        mydouble D= sqrt(0.25*Q*Q+ (P*P*P/27.0));
        
        sigma= (2.0/ (sqrt(T))) * sqrt(cubeRoot(-0.5*Q+D) + cubeRoot(-0.5*Q-D) - gamma2/(3*gamma4));
    }
    if(CmTilde>=0.99999971){
        sigma= (2.0/ (sqrt(T))) * sqrt(-piHalfLog);
    }
    return sigma;
}
int main(){
    mydouble Cm= 19.9;
    mydouble q= 0.0;
    mydouble r= 0.1;
    mydouble S0=60.0;
    mydouble T= 1.0;
    std::cout<<sigmaImpApprox(Cm,q,r,S0,T)<<std::endl;
    return 0;
}
