### Fits extended water retention curve models 
### T.P. Leao February 2025 

# dependencies
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12
#font = {'family' : 'normal', 'weight' : 'bold',  'size'   : 22}

import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from scipy.stats import t

# Ignore all warnings
import warnings
warnings.filterwarnings("ignore")


num = [3,5,8,10,50,100,500]
for zz in num:
    print(zz)

    nex = 0
    while (nex < 100):
        nex = nex + 1
        print(nex)

        # nint used are 3, 5, 8, 10, 20, 50, 100, 500 
        nint = zz

        # simulated range for predicted plot
        sim1 = np.arange(1, 100, (100 - 1)/nint)
        sim2 = np.arange(100, 1000, (1000 - 100)/nint)
        sim3 = np.arange(1000, 10000, (10000 - 1000)/nint)
        sim4 = np.arange(10000, 100000, (100000 - 10000)/nint)
        #sim5 = np.arange(100000, 1000000, (1000000 - 100000)/nint)

        simt = np.concatenate((sim1, sim2, sim3, sim4), axis=None)

        #simt = np.arange(1, 100000, (100000 - 1)/100)


        ##simt= np.log(simt)

        ##print(len(simt))

        ##notice you wrote 0.001, 0.005 and 0.01 in the paper and here it is 0.001, 0.010 and 0.020
        sig1 = 0.01
        sig2 = 0.025 
        sig3 = 0.05

        sigm = sig3

        errt = np.random.normal(0, sigm, len(simt))

        qs3p = 0.655
        qr3p = 0.016
        #w1 = .2
        alpha3p = 0.051
        n3p = 1.792
        alpha30p = 1.57E-05
        n30p = 2.472
        w2p = 0.349

        ##print(len(simt))

        # "true" model for the monte-carlo simulation - VGe
        def funct(x, sig, err):
            return  qr3p + (qs3p - qr3p)   * ( (1-w2p) * ( 1 + (alpha3p*x)**n3p ) ** (-( 1 - 1 / n3p)) + w2p * ( 1 + (alpha30p*x)**n30p ) ** (-( 1 - 1 / n30p)) )  + err

        yt = funct(simt, sigm, errt)

        # plot - define paramters according to user need

        #plt.plot(simt, yt, "o")
        #plt.plot(Sim, func(Sim, *popt), linestyle = 'solid', linewidth = 1.2, )


        #color = 'gray'
        #plt.legend(['EG', 'VG', 'MSN', 'FX'] )
        #plt.legend(['Degree 2', 'Degree 3', 'Degree 4', 'Degree 5', 'Degree 6', 'Degree 7', 'Degree 8', 'Degree 9', 'Degree 10'], frameon = False )

        #plt.ylabel("Volumetric water content (m$^3$ m$^{-3}$)")
        #plt.ylim([0.0, 0.8])
        #plt.xlim([0.0, 1000])
        #plt.xlabel("Water potential (|hPa|)")
        #plt.xscale('log')
        ##plt.savefig('/home/t01/1_2024/ArtigosTPL/ArtigoAnaJulia/ManuscriptAJ/OX_H_e.pdf', format='pdf', dpi=300)

        #plt.show()

        #plt.hist(yt, 60, density = True)
        #plt.show()


        #plt.hist(errt, 60, density = True)
        #plt.show()



        # macro for significance
    ##    def sigfun(x):
    ##        if x > 0.05:
    ##            return "ns"
    ##        elif (x <= 0.05 and x > 0.01):
    ##            return "*"
    ##        elif (x <= 0.01 and x > 0.001):
    ##            return "**"
    ##        else:
    ##            return "***" 


        ### Define saturated water content for all stupid equations - Soil used in MCM is the Regosol
        ##thetasat = 0.554
        ##thetares = 0.005
        #print("The saturation water content for the selected dataset is: ", qs)


        ##print("#########################################################################")
        ##print("#########################################################################")
        ##print("#########################################################################")


        ################# Groenvevelt-Grant

        ##print("Groenevelt-Grant ################################################")

        # vector of initial parameter values

        k0 = 7.89E-04
        k1 = -887.401
        epsilon = -0.301

        InitPar = [k0, k1, epsilon]

        # function to be fit - change for a custom function
        def func(x, k0, k1, epsilon):
            return  k1 * (np.exp(-k0/6.9**epsilon) - np.exp(-k0/((np.log10(x))**epsilon)))  

        # fits data

        try:
            popt, pcov = curve_fit(func, simt, yt, p0 = InitPar, method = 'lm', maxfev=5000) 
        except RuntimeError:
        ##    print("Error - curve_fit failed")
            popt = np.array([0.0,  0.0,  0.0])
            pcov = np.zeros((3, 3), dtype='d')

        # degrees of freedom n - number of parameters
        n = len(simt)
        d = len(popt)
        df = n-d


        # standard error
        se0 = np.sqrt(pcov[0,0])
        se1 = np.sqrt(pcov[1,1])
        se2 = np.sqrt(pcov[2,2])

        # t value
        t0 = popt[0]/np.sqrt(pcov[0,0])
        t1 = popt[1]/np.sqrt(pcov[1,1])
        t2 = popt[2]/np.sqrt(pcov[2,2])

        # significance
        # uses the sf function, equivalent (but more precise) than 1 - cdf
        # see scipy.stats.t documentation
        p0 = 2*(t.sf(np.abs(t0), df, loc=0, scale =1))
        p1 = 2*(t.sf(np.abs(t1), df, loc=0, scale =1))
        p2 = 2*(t.sf(np.abs(t2), df, loc=0, scale =1))

        ##print(round(popt[0],4),  round(se0,4), round(t0,4), np.format_float_scientific(p0, exp_digits=2, precision=4), sigfun(p0)) 
        ##print(round(popt[1],4),  round(se1,4), round(t1,4), np.format_float_scientific(p1, exp_digits=2, precision=4), sigfun(p1)) 
        ##print(round(popt[2],4),  round(se2,4), round(t2,4), np.format_float_scientific(p2, exp_digits=2, precision=4), sigfun(p2)) 


        # calculates the pseudo r2 from the correlation coefficient
        r = np.corrcoef(yt, func(simt, *popt))
        r2 = r[0,1]**2

        # prints R2 to terminal
        ##print("R2: ", round(r2,4))


        # calculates and prints R2adj to terminal
        r2adj = 1- (1 - r2) * (n-1)/(df)
        ##print("R2adj: ", round(r2adj ,4))



        # calculates the root mean square error - rmse
        RMSE = np.sqrt(np.sum((yt- func(simt, *popt))**2)/df)
        ##print("RMSE: ", round(RMSE, 4))


        # calculates sum of squares of deviations SSR and prints to terminal
        SSR = np.sum((yt- func(simt, *popt))**2)
        ##print("SSR: ", round(SSR, 4))

        # calculates the Akaike Information Criterion - AIC 
        AIC = n * np.log(SSR/n) + 2 * d 
        ##print("AIC: ", round(AIC, 4))

        # calculates the corrected Akaike Information Criterion - AICc 
        AICc = n * np.log(SSR/n) + 2 * d + (2*d*(d+1))/(n - d -1) 
        ##print("AICc: ", round(AICc, 4))

        # calculates the Bayesian Information Criterion - BIC 
        BIC = n * np.log(SSR/n) + d * np.log(n)
        ##print("BIC: ", round(BIC, 4))

        # calculates the Kayshap Information Criterion - KIC

        sigma = np.linalg.det((pcov))

        KIC = n * np.log(SSR/n) - d * np.log((2*np.pi)) - np.log(sigma) 
        ##print("KIC: ", round(KIC, 4))


        ################# Extended Groenvelt-Grant

        ##print("Extended Groenevelt-Grant ################################################")

        # vector of initial parameter values

        A1_0 = 0.387
        A2_0 = 0.052
        phi1_0 = 50.770
        phi2_0 = 2076.007
        k0_0 = 7.78E+09
        k1_0 = 0.237
        n_0 = 14.413

        InitPar0 = [A1_0, A2_0, phi1_0, phi2_0, k0_0, k1_0, n_0]

        # function to be fit - change for a custom function
        def func0(x, A1, A2, phi1, phi2, k0, k1, n):
            return A1 * np.exp(-x/phi1) + A2 * np.exp(-x/phi2) + k1 * (np.exp(-k0/6.9**n) - np.exp(-k0/(np.log10(x)**n)))  

        # fits data

        try:
            popt0, pcov0 = curve_fit(func0, simt, yt, p0 = InitPar0, method = 'lm', maxfev=5000) 
        except RuntimeError:
        ##    print("Error - curve_fit failed")
            popt0 = np.array([-0.0,  -0.0,  -0.0,  -0.0,  -0.0,  -0.0,  -0.0])
            pcov0 = np.zeros((7, 7), dtype='d')

        # degrees of freedom n - number of parameters
        n0 = len(simt)
        d0 = len(popt0)
        df0 = n0-d0

        # standard error
        se00 = np.sqrt(pcov0[0,0])
        se10 = np.sqrt(pcov0[1,1])
        se20 = np.sqrt(pcov0[2,2])
        se30 = np.sqrt(pcov0[3,3])
        se40 = np.sqrt(pcov0[4,4])
        se50 = np.sqrt(pcov0[5,5])
        se60 = np.sqrt(pcov0[6,6])

        # t value
        t00 = popt0[0]/np.sqrt(pcov0[0,0])
        t10 = popt0[1]/np.sqrt(pcov0[1,1])
        t20 = popt0[2]/np.sqrt(pcov0[2,2])
        t30 = popt0[3]/np.sqrt(pcov0[3,3])
        t40 = popt0[4]/np.sqrt(pcov0[4,4])
        t50 = popt0[5]/np.sqrt(pcov0[5,5])
        t60 = popt0[6]/np.sqrt(pcov0[6,6])

        # significance
        # uses the sf function, equivalent (but more precise) than 1 - cdf
        # see scipy.stats.t documentation
        p00 = 2*(t.sf(np.abs(t00), df0, loc=0, scale =1))
        p10 = 2*(t.sf(np.abs(t10), df0, loc=0, scale =1))
        p20 = 2*(t.sf(np.abs(t20), df0, loc=0, scale =1))
        p30 = 2*(t.sf(np.abs(t30), df0, loc=0, scale =1))
        p40 = 2*(t.sf(np.abs(t40), df0, loc=0, scale =1))
        p50 = 2*(t.sf(np.abs(t50), df0, loc=0, scale =1))
        p60 = 2*(t.sf(np.abs(t60), df0, loc=0, scale =1))

        ##print(round(popt0[0],4),  round(se00,4), round(t00,4), np.format_float_scientific(p00, exp_digits=2, precision=4), sigfun(p00)) 
        ##print(round(popt0[1],4),  round(se10,4), round(t10,4), np.format_float_scientific(p10, exp_digits=2, precision=4), sigfun(p10)) 
        ##print(round(popt0[2],4),  round(se20,4), round(t20,4), np.format_float_scientific(p20, exp_digits=2, precision=4), sigfun(p20)) 
        ##print(round(popt0[3],4),  round(se30,4), round(t30,4), np.format_float_scientific(p30, exp_digits=2, precision=4), sigfun(p30)) 
        ##print(round(popt0[4],4),  round(se40,4), round(t40,4), np.format_float_scientific(p40, exp_digits=2, precision=4), sigfun(p40)) 
        ##print(round(popt0[5],4),  round(se50,4), round(t50,4), np.format_float_scientific(p50, exp_digits=2, precision=4), sigfun(p50)) 
        ##print(round(popt0[6],4),  round(se60,4), round(t60,4), np.format_float_scientific(p60, exp_digits=2, precision=4), sigfun(p60)) 


        # calculates the pseudo r2 from the correlation coefficient
        r0 = np.corrcoef(yt, func0(simt, *popt0))
        r20 = r0[0,1]**2

        ### prints R2 to terminal
        ##print("R2: ", round(r20,4))
        ##
        ##
        ### calculates and prints R2adj to terminal
        r2adj0 = 1- (1 - r20) * (n0-1)/(df0)
        ##print("R2adj: ", round(r2adj0 ,4))

        # calculates the root mean square error - rmse
        RMSE0 = np.sqrt(np.sum((yt- func0(simt, *popt0))**2)/df0)
        ##print("RMSE: ", round(RMSE0, 4))


        # calculates sum of squares of deviations SSR and prints to terminal
        SSR0 = np.sum((yt- func0(simt, *popt0))**2)
        ##print("SSR: ", round(SSR0, 4))

        # calculates the Akaike Information Criterion - AIC 
        AIC0 = n0 * np.log(SSR0/n0) + 2 * d0 
        ##print("AIC: ", round(AIC0, 4))

        # calculates the corrected Akaike Information Criterion - AICc 
        AICc0 = n0 * np.log(SSR0/n0) + 2 * d0 + (2*d0*(d0+1))/(n0 - d0 -1) 
        ##print("AICc: ", round(AICc0, 4))

        # calculates the Bayesian Information Criterion - BIC 
        BIC0 = n0 * np.log(SSR0/n0) + d0 * np.log(n0)
        ##print("BIC: ", round(BIC0, 4))

        # calculates the Kayshap Information Criterion - KIC

        sigma0 = np.linalg.det((pcov0))

        KIC0 = n0 * np.log(SSR0/n0) - d0 * np.log((2*np.pi)) - np.log(sigma0) 
        ##print("KIC: ", round(KIC0, 4))




        ################# MSN


        ##print("Mehta-Shiozawa-Nakano  ################################################")
        # vector of initial parameter values

        #A11 = .15 
        #B11 = .3
        #C11 = 8000
        #alpha11 = 70
        #n11 = 3

        #A11 = .0015 
        #B11 = .3
        #C11 = 800000
        #alpha11 = 70000
        #n11 = 3

        #A11 = 4.4e-01 
        #B11 = 3.44e-01 
        #C11 = 9.9e-01 
        #alpha11 = 2
        #n11 = -4

        A11 = 0.198
        B11 = 0.488
        C11 = 15.657
        alpha11 = 0.027
        n11 = 3.379


        InitPar1 = [A11, B11, C11, alpha11, n11]

        # function to be fit - change for a custom function
        ### I ***** UP HERE IT IS (1/N - 1) AND NOT (1/N + 1)
        def func1(x, A11, B11, C11, alpha11, n11):
            return A11 * (1 + (alpha11 * x) **n11) **(-(1-1/n)) + B11 * (1 - np.log( x + 1) / C11)    

        # fits data

        try:
            popt1, pcov1 = curve_fit(func1, simt, yt, p0 = InitPar1, method = 'lm', maxfev=5000) 
        except RuntimeError:
        ##    print("Error - curve_fit failed")
            popt0 = np.array([0.0,  0.0,  0.0, 0.0, 0.0])
            pcov0 = np.zeros((5, 5), dtype='d')

        # degrees of freedom n - number of parameters
        n1 = len(simt)
        d1 = len(popt1)
        df1 = n1-d1

        # standard error
        se01 = np.sqrt(pcov1[0,0])
        se11 = np.sqrt(pcov1[1,1])
        se21 = np.sqrt(pcov1[2,2])
        se31 = np.sqrt(pcov1[3,3])
        se41 = np.sqrt(pcov1[4,4])

        # t value
        t01 = popt1[0]/np.sqrt(pcov1[0,0])
        t11 = popt1[1]/np.sqrt(pcov1[1,1])
        t21 = popt1[2]/np.sqrt(pcov1[2,2])
        t31 = popt1[3]/np.sqrt(pcov1[3,3])
        t41 = popt1[4]/np.sqrt(pcov1[4,4])

        # significance
        # uses the sf function, equivalent (but more precise) than 1 - cdf
        # see scipy.stats.t documentation
        p01 = 2*(t.sf(np.abs(t01), df1, loc=0, scale =1))
        p11 = 2*(t.sf(np.abs(t11), df1, loc=0, scale =1))
        p21 = 2*(t.sf(np.abs(t21), df1, loc=0, scale =1))
        p31 = 2*(t.sf(np.abs(t31), df1, loc=0, scale =1))
        p41 = 2*(t.sf(np.abs(t41), df1, loc=0, scale =1))

        ##print(round(popt1[0],4),  round(se01,4), round(t01,4), np.format_float_scientific(p01, exp_digits=2, precision=4), sigfun(p01)) 
        ##print(round(popt1[1],4),  round(se11,4), round(t11,4), np.format_float_scientific(p11, exp_digits=2, precision=4), sigfun(p11)) 
        ##print(round(popt1[2],4),  round(se21,4), round(t21,4), np.format_float_scientific(p21, exp_digits=2, precision=4), sigfun(p21)) 
        ##print(round(popt1[3],4),  round(se31,4), round(t31,4), np.format_float_scientific(p31, exp_digits=2, precision=4), sigfun(p31)) 
        ##print(round(popt1[4],4),  round(se41,4), round(t41,4), np.format_float_scientific(p41, exp_digits=2, precision=4), sigfun(p41)) 


        # calculates the pseudo r2 from the correlation coefficient
        r1 = np.corrcoef(yt, func1(simt, *popt1))
        r21 = r1[0,1]**2

        # prints R2 to terminal
        ##print("R2: ", round(r21,4))

        # calculates and prints R2adj to terminal
        r2adj1 = 1- (1 - r21) * (n1-1)/(df1)
        ##print("R2adj: ", round(r2adj1, 4))

        # calculates the root mean square error - rmse
        RMSE1 = np.sqrt(np.sum((yt- func1(simt, *popt1))**2)/df1)
        ##print("RMSE: ", round(RMSE1, 4))


        # calculates sum of squares of deviations SSR and prints to terminal
        SSR1 = np.sum((yt- func1(simt, *popt1))**2)
        ##print("SSR: ", round(SSR1, 4))

        # calculates the Akaike Information Criterion - AIC 
        AIC1 = n1 * np.log(SSR1/n1) + 2 * d1 
        ##print("AIC: ", round(AIC1, 4))

        # calculates the corrected Akaike Information Criterion - AICc 
        AICc1 = n1 * np.log(SSR1/n1) + 2 * d1 + (2*d1*(d1+1))/(n1 - d1 -1) 
        ##print("AICc: ", round(AICc1, 4))

        # calculates the Bayesian Information Criterion - BIC 
        BIC1 = n1 * np.log(SSR1/n1) + d1 * np.log(n1)
        ##print("BIC: ", round(BIC1, 4))

        # calculates the Kayshap Information Criterion - KIC

        sigma1 = np.linalg.det((pcov1))

        KIC1 = n1 * np.log(SSR1/n1) - d1 * np.log((2*np.pi)) - np.log(sigma1) 
        ##print("KIC: ", round(KIC1, 4))



        ################# FX

        ##print("Fredlund-Xing  ################################################")

        # vector of initial parameter values


        qs2 = 0.609
        phir2 = 7.30E+08
        phi02 = 215412.290
        alpha2 = 29.064
        n2 = 10.236
        m2 = 0.225

        InitPar2 = [qs2, phir2, phi02, alpha2, n2, m2]

        # function to be fit - change for a custom function
        def func2(x, qs2, phir2, phi02, alpha2, n2, m2):
            return qs2   * ( 1 - np.log(1 + x/phir2)/np.log(1+phi02/phir2))   * (np.log(np.exp(1) + (x/alpha2)**n2))**(-m2) 

        # fits data

        try:
            popt2, pcov2 = curve_fit(func2, simt, yt, p0 = InitPar2, method = 'lm', maxfev=5000) 
        except RuntimeError:
        ##    print("Error - curve_fit failed")
            popt0 = np.array([0.0,  0.0,  0.0, 0.0,  0.0,  0.0])
            pcov0 = np.zeros((6, 6), dtype='d')

        # degrees of freedom n - number of parameters
        n2 = len(simt)
        d2 = len(popt2)
        df2 = n2-d2

        # standard error
        se02 = np.sqrt(pcov2[0,0])
        se12 = np.sqrt(pcov2[1,1])
        se22 = np.sqrt(pcov2[2,2])
        se32 = np.sqrt(pcov2[3,3])
        se42 = np.sqrt(pcov2[4,4])
        se52 = np.sqrt(pcov2[5,5])

        # t value
        t02 = popt2[0]/np.sqrt(pcov2[0,0])
        t12 = popt2[1]/np.sqrt(pcov2[1,1])
        t22 = popt2[2]/np.sqrt(pcov2[2,2])
        t32 = popt2[3]/np.sqrt(pcov2[3,3])
        t42 = popt2[4]/np.sqrt(pcov2[4,4])
        t52 = popt2[5]/np.sqrt(pcov2[5,5])

        # significance
        # uses the sf function, equivalent (but more precise) than 1 - cdf
        # see scipy.stats.t documentation
        p02 = 2*(t.sf(np.abs(t02), df2, loc=0, scale =1))
        p12 = 2*(t.sf(np.abs(t12), df2, loc=0, scale =1))
        p22 = 2*(t.sf(np.abs(t22), df2, loc=0, scale =1))
        p32 = 2*(t.sf(np.abs(t32), df2, loc=0, scale =1))
        p42 = 2*(t.sf(np.abs(t42), df2, loc=0, scale =1))
        p52 = 2*(t.sf(np.abs(t52), df2, loc=0, scale =1))

        ##print(round(popt2[0],4),  round(se02,4), round(t02,4), np.format_float_scientific(p02, exp_digits=2, precision=4), sigfun(p02)) 
        ##print(round(popt2[1],4),  round(se12,4), round(t12,4), np.format_float_scientific(p12, exp_digits=2, precision=4), sigfun(p12)) 
        ##print(round(popt2[2],4),  round(se22,4), round(t22,4), np.format_float_scientific(p22, exp_digits=2, precision=4), sigfun(p22)) 
        ##print(round(popt2[3],4),  round(se32,4), round(t32,4), np.format_float_scientific(p32, exp_digits=2, precision=4), sigfun(p32)) 
        ##print(round(popt2[4],4),  round(se42,4), round(t42,4), np.format_float_scientific(p42, exp_digits=2, precision=4), sigfun(p42)) 
        ##print(round(popt2[5],4),  round(se52,4), round(t52,4), np.format_float_scientific(p52, exp_digits=2, precision=4), sigfun(p52)) 


        # calculates the pseudo r2 from the correlation coefficient
        r2a = np.corrcoef(yt, func2(simt, *popt2))
        r22 = r2a[0,1]**2

        # prints R2 to terminal
        ##print("R2: ", round(r22,4))


        # calculates and prints R2adj to terminal
        r2adj2 = 1- (1 - r22) * (n2-1)/(df2)
        ##print("R2adj: ", round(r2adj2, 4))


        # calculates the root mean square error - rmse
        RMSE2 = np.sqrt(np.sum((yt- func2(simt, *popt2))**2)/df2)
        ##print("RMSE: ", round(RMSE2, 4))

        # calculates sum of squares of deviations SSR and prints to terminal
        SSR2 = np.sum((yt- func2(simt, *popt2))**2)
        ##print("SSR: ", round(SSR2, 4))

        # calculates the Akaike Information Criterion - AIC 
        AIC2 = n2 * np.log(SSR2/n2) + 2 * d2 
        ##print("AIC: ", round(AIC2, 4))

        # calculates the corrected Akaike Information Criterion - AICc 
        AICc2 = n2 * np.log(SSR2/n2) + 2 * d2 + (2*d2*(d2+1))/(n2 - d2 -1) 
        ##print("AICc: ", round(AICc2, 4))

        # calculates the Bayesian Information Criterion - BIC 
        BIC2 = n2 * np.log(SSR2/n2) + d2 * np.log(n2)
        ##print("BIC: ", round(BIC2, 4))

        # calculates the Kayshap Information Criterion - KIC

        sigma2 = np.linalg.det((pcov2))

        KIC2 = n2 * np.log(SSR2/n2) - d2 * np.log((2*np.pi)) - np.log(sigma2) 
        ##print("KIC: ", round(KIC2, 4))


        ################# Extended van Genuchten

        ##print("Extended van Genuchten  ################################################")

        # vector of initial parameter values

        qs3 = 0.655
        qr3 = 0.016
        #w1 = .2
        alpha3 = 0.051
        n3 = 1.792
        alpha30 = 1.57E-05
        n30 = 2.472
        w2 = 0.349

        InitPar3 = [qs3, qr3,  alpha3, n3, alpha30, n30, w2]
        #InitPar3 = [qs3, alpha3, n3, alpha30, n30, w2]
        #InitPar3 = [qr3, alpha3, n3, alpha30, n30, w2]

        # function to be fit - change for a custom function
        def func3(x,  qs3, qr3, alpha3, n3, alpha30, n30, w2):
            return qr3 + (qs3 - qr3)   * ( (1-w2) * ( 1 + (alpha3*x)**n3 ) ** (-( 1 - 1 / n3)) + w2 * ( 1 + (alpha30*x)**n30 ) ** (-( 1 - 1 / n30)) )    

        #np.clip(,0,0.5)

        # fits data

        try:
            popt3, pcov3 = curve_fit(func3, simt, yt, p0 = InitPar3, method = 'lm', maxfev=5000) 
        except RuntimeError:
        ##    print("Error - curve_fit failed")
            popt3 = np.array([0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0])
            pcov3 = np.zeros((7, 7), dtype='d')


        # degrees of freedom n - number of parameters
        n3 = len(simt)
        d3 = len(popt3)
        df3 = n3-d3


        # standard error
        se03 = np.sqrt(pcov3[0,0])
        se13 = np.sqrt(pcov3[1,1])
        se23 = np.sqrt(pcov3[2,2])
        se33 = np.sqrt(pcov3[3,3])
        se43 = np.sqrt(pcov3[4,4])
        se53 = np.sqrt(pcov3[5,5])
        se63 = np.sqrt(pcov3[6,6])

        # t value
        t03 = popt3[0]/np.sqrt(pcov3[0,0])
        t13 = popt3[1]/np.sqrt(pcov3[1,1])
        t23 = popt3[2]/np.sqrt(pcov3[2,2])
        t33 = popt3[3]/np.sqrt(pcov3[3,3])
        t43 = popt3[4]/np.sqrt(pcov3[4,4])
        t53 = popt3[5]/np.sqrt(pcov3[5,5])
        t63 = popt3[6]/np.sqrt(pcov3[6,6])

        # significance
        # uses the sf function, equivalent (but more precise) than 1 - cdf
        # see scipy.stats.t documentation
        p03 = 2*(t.sf(np.abs(t03), df3, loc=0, scale =1))
        p13 = 2*(t.sf(np.abs(t13), df3, loc=0, scale =1))
        p23 = 2*(t.sf(np.abs(t23), df3, loc=0, scale =1))
        p33 = 2*(t.sf(np.abs(t33), df3, loc=0, scale =1))
        p43 = 2*(t.sf(np.abs(t43), df3, loc=0, scale =1))
        p53 = 2*(t.sf(np.abs(t53), df3, loc=0, scale =1))
        p63 = 2*(t.sf(np.abs(t63), df3, loc=0, scale =1))

        ##print(round(popt3[0],4),  round(se03,4), round(t03,4), np.format_float_scientific(p03, exp_digits=2, precision=4), sigfun(p03)) 
        ##print(round(popt3[1],4),  round(se13,4), round(t13,4), np.format_float_scientific(p13, exp_digits=2, precision=4), sigfun(p13)) 
        ##print(round(popt3[2],4),  round(se23,4), round(t23,4), np.format_float_scientific(p23, exp_digits=2, precision=4), sigfun(p23)) 
        ##print(round(popt3[3],4),  round(se33,4), round(t33,4), np.format_float_scientific(p33, exp_digits=2, precision=4), sigfun(p33)) 
        ##print(round(popt3[4],4),  round(se43,4), round(t43,4), np.format_float_scientific(p43, exp_digits=2, precision=4), sigfun(p43)) 
        ##print(round(popt3[5],4),  round(se53,4), round(t53,4), np.format_float_scientific(p53, exp_digits=2, precision=4), sigfun(p53)) 
        ##print(round(popt3[6],4),  round(se63,4), round(t63,4), np.format_float_scientific(p63, exp_digits=2, precision=4), sigfun(p63)) 



        # calculates the pseudo r2 from the correlation coefficient
        r3 = np.corrcoef(yt, func3(simt, *popt3))
        r23 = r3[0,1]**2

        # prints R2 to terminal
        ##print("R2: ", round(r23,4))

        # calculates and prints R2adj to terminal
        r2adj3 = 1- (1 - r23) * (n3-1)/(df3)
        ##print("R2adj: ", round(r2adj3, 4))

        # calculates the root mean square error - rmse
        RMSE3 = np.sqrt(np.sum((yt- func3(simt, *popt3))**2)/df3)
        ##print("RMSE: ", round(RMSE3, 4))

        # calculates sum of squares of deviations SSR and prints to terminal
        SSR3 = np.sum((yt- func3(simt, *popt3))**2)
        ##print("SSR: ", round(SSR3, 4))

        # calculates the Akaike Information Criterion - AIC 
        AIC3 = n3 * np.log(SSR3/n3) + 2 * d3 
        ##print("AIC: ", round(AIC3, 4))

        # calculates the corrected Akaike Information Criterion - AICc 
        AICc3 = n3 * np.log(SSR3/n3) + 2 * d3 + (2*d3*(d3+1))/(n3 - d3 -1) 
        ##print("AICc: ", round(AICc3, 4))

        # calculates the Bayesian Information Criterion - BIC 
        BIC3 = n3 * np.log(SSR3/n3) + d3 * np.log(n3)
        ##print("BIC: ", round(BIC3, 4))

        # calculates the Kayshap Information Criterion - KIC

        sigma3 = np.linalg.det((pcov3))

        KIC3 = n3 * np.log(SSR3/n3) - d3 * np.log((2*np.pi)) - np.log(sigma3) 
        ##print("KIC: ", round(KIC3, 4))

        #######################################################################################################


        # simulated range for predicted plot
    ##    Sim = np.arange(np.min(simt), np.max(simt), np.abs(np.max(simt) - np.min(simt))/100000)

        # plot - define paramters according to user need

    ##    plt.plot(simt, yt, "o")
        #plt.errorbar(simt, yt, yerr = ft['sd'], fmt = "ok", ecolor = "r", fillstyle = "none", capsize = 2.0)
    ##    plt.plot(Sim, func(Sim, *popt), linestyle = 'solid', linewidth = 1.2, )
    ##    plt.plot(Sim, func0(Sim, *popt0), linestyle = 'dotted', linewidth = 1.2)
    ##    plt.plot(Sim, func1(Sim, *popt1), linestyle = 'dashed', linewidth = 1.2)
    ##    plt.plot(Sim, func2(Sim, *popt2), linestyle = 'dashdot', linewidth = 1.2, color = "red")
    ##    plt.plot(Sim, func3(Sim, *popt3), linestyle = (0, (3,1,1,1,1)), linewidth = 1.2, color = "gray")


        #color = 'gray'
        #plt.legend(['EG', 'VG', 'MSN', 'FX'] )
    ##    plt.legend(['Sim', 'GG', 'GGe', 'MSN', 'FX', 'VGe'], frameon = False )
        
    ##    plt.ylabel("Volumetric water content (m$^3$ m$^{-3}$)")
    ##    plt.ylim([0.0, 0.8])
    ##    plt.xlabel("Water potential (|hPa|)")
    ##    plt.xscale('log')
        ##plt.savefig('/home/t01/1_2024/ArtigosTPL/ArtigoAnaJulia/ManuscriptAJ/GX_Cg_e.pdf', format='pdf', dpi=300)

    ##    plt.show()

 ###### Simplified code to check fit
        
##        Sim = np.arange(np.min(simt), np.max(simt), np.abs(np.max(simt) - np.min(simt))/100000)
##        plt.plot(simt, yt, "o")
##        plt.plot(Sim, func(Sim, *popt), linestyle = 'solid', linewidth = 1.2, )
##        plt.plot(Sim, func0(Sim, *popt0), linestyle = 'dotted', linewidth = 1.2)
##        plt.plot(Sim, func1(Sim, *popt1), linestyle = 'dashed', linewidth = 1.2)
##        plt.plot(Sim, func2(Sim, *popt2), linestyle = 'dashdot', linewidth = 1.2, color = "red")
##        plt.plot(Sim, func3(Sim, *popt3), linestyle = (0, (3,1,1,1,1)), linewidth = 1.2, color = "gray")
##
##
##        plt.legend(['Sim', 'GG', 'GGe', 'MSN', 'FX', 'VGe'], frameon = False )
##        
##        plt.ylabel("Volumetric water content (m$^3$ m$^{-3}$)")
##        plt.ylim([0.0, 0.8])
##        plt.xlabel("Water potential (|hPa|)")
##        plt.xscale('log')
##  
##        plt.show()

        ##print("Fitting parameters")
        ##print("GG", popt)
        ##print("GGe", popt0)
        ##print("MSN", popt1)
        ##print("FX", popt2)
        ##print("VGe", popt3)


        ##print("Statistics")
        ##print("GG", round(r2,4) , round(RMSE,4), round(AIC,1), round(AICc,1), round(BIC,1), round(KIC,1))
        ##print("GGe", round(r20,4), round(RMSE0,4), round(AIC0,1), round(AICc0,1), round(BIC0,1), round(KIC0,1))
        ##print("MSN", round(r21,4), round(RMSE1,4), round(AIC1,1), round(AICc1,1), round(BIC1,1), round(KIC1,1))
        ##print("FX", round(r22,4), round(RMSE2,4), round(AIC2,1), round(AICc2,1), round(BIC2,1), round(KIC2,1))
        ##print("VGe", round(r23,4), round(RMSE3,4), round(AIC3,1), round(AICc3,1), round(BIC3,1), round(KIC3,1))

        ##print("Nobs", len(simt))

        ##print("Statistics")
        print("Model", "R2" , "R2adj" ,"RMSE", "AIC", "AICc", "BIC", "KIC")
        print("FX", round(r22,4), round(r2adj2,4), round(RMSE2,4), round(AIC2,1), round(AICc2,1), round(BIC2,1), round(KIC2,1), np.log(sigma2))
        print("GG", round(r2,4), round(r2adj,4) , round(RMSE,4), round(AIC,1), round(AICc,1), round(BIC,1), round(KIC,1), np.log(sigma))
        print("GGe", round(r20,4), round(r2adj0,4), round(RMSE0,4), round(AIC0,1), round(AICc0,1), round(BIC0,1), round(KIC0,1), np.log(sigma0))
        print("MSN", round(r21,4), round(r2adj1,4), round(RMSE1,4), round(AIC1,1), round(AICc1,1), round(BIC1,1), round(KIC1,1), np.log(sigma1))
        print("VGe", round(r23,4), round(r2adj3,4), round(RMSE3,4), round(AIC3,1), round(AICc3,1), round(BIC3,1), round(KIC3,1), np.log(sigma3))



        #print(pcov1)







