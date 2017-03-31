#################################################################################################
#################################################################################################
# The MIT License (MIT)                                                                         #
# Copyright (C) 2017 Ivan Matic, Rados Radoicic, and Dan Stefanica                              #
#                                                                                               #
# Permission is hereby granted, free of charge, to any person obtaining a copy of this          #
# software and associated documentation files (the "Software"), to deal in the Software         #
# without restriction, including without limitation the rights to use, copy, modify, merge,     #
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons    #
# to whom the Software is furnished to do so, subject to the following conditions:              #
#                                                                                               #
# The above copyright notice and this permission notice shall be included in all copies or      #
# substantial portions of the Software.                                                         #
#                                                                                               #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,           #
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR      #
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE     #
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR          #
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                                     #
#################################################################################################


cubeRoot<- function(x){
    ### This one calculates x^(1/3), even when x is negative.
    ### Original x^(1/3) does not work because R seems to implement cube roots using logarithms
    absX<- x
    signX<- 1.0
    if(absX<0){
        absX<- -1.0*x
        signX<- -1.0
    }
    return(signX*  absX^(1/3))
}


sigmaImpApprox<- function(Cm, q, r, S0, T){
    #####################################################
    # Input:                                            #
    # Cm = market price of the call option              #
    # q  = constant rate of dividends                   #
    # r  = constant interest rate                       #
    # S0 = price of the underlying security             #
    # T  = option maturity                              #
    #                                                   #
    # Output:                                           #
    # ATM-forward implied volatility approximation      #
    #####################################################

    CmTilde<- Cm/(S0*exp(-q*T))
    gamma2<- -1/3+1/pi
    gamma4<- 7/90 -2/(3*pi)+ 4/(3*pi^2)
    gamma6<- -1/70 + 4/(15*pi) - 4/(3*pi^2) + 2/(pi^3)
    piHalfLog<- 0.5*pi* log(1- CmTilde^2)
    
    if(CmTilde<=0.65949631){
        gamma4OverGamma6<- gamma4/gamma6
        gamma4OverGamma6Cube<- gamma4OverGamma6^3
        u<- gamma2/gamma6- 0.375 * gamma4OverGamma6^2
        v<- 0.125* gamma4OverGamma6Cube - 0.5* (gamma2/gamma6) * gamma4OverGamma6+1/gamma6
        w<- piHalfLog/gamma6 - (3/256)*gamma4OverGamma6*gamma4OverGamma6Cube - 0.25 * gamma4OverGamma6/gamma6 + 0.0625*(gamma2/gamma6)*gamma4OverGamma6^2
        P<- (-12*w-u*u)/12
        Q<- (-u^3+36*u*w)/(108)-v^2/8
        D<- (0.25*Q^2+ (P^3/27))^(0.5)
        
        
        m<-cubeRoot(-0.5*Q+D)  + cubeRoot(-0.5*Q-D) -u/3
        sigma<- (2/ (T^(0.5))) * (0.5*(- (2*m)^(0.5) + (2*v/((2*m)^(0.5)) - 2*(u+m))^(0.5)) - 0.25*gamma4OverGamma6)^(0.5)
        
    }
    if((CmTilde>0.65949631)&&(CmTilde<0.99999971)){
        P<- (3*gamma4 - gamma2^2)/(3*gamma4^2)
        Q<- (2*gamma2^3-9*gamma2*gamma4)/(27*gamma4^3)+ piHalfLog/gamma4
        D<-(0.25*Q^2+ (P^3/27))^(0.5)
        sigma<-(2/ (T^(0.5))) * (cubeRoot(-0.5*Q+D) + cubeRoot(-0.5*Q-D) - gamma2/(3*gamma4))^(0.5)
    }
    if(CmTilde>=0.99999971){
        sigma<-(2/ (T^(0.5))) * (-piHalfLog)^(0.5)
    }
    return(sigma)
}

