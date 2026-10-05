### Fits water retention curve models to wet range 
### T.P. Leao Jan, 10, 2024 

# dependencies
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12
#font = {'family' : 'normal', 'weight' : 'bold',  'size'   : 22}


import numpy as np
from scipy.optimize import curve_fit, least_squares
from scipy import special
import pandas as pd
from scipy.stats import t
from array import array


soil = 'LVA_A'


# import external file 

if soil == "GX_A":
    ft =  pd.read_csv(r"~/code/python/artigoaj/gx_a_w.txt", delim_whitespace = True)

elif soil == "GX_Cg":
    ft =  pd.read_csv(r"~/code/python/artigoaj/gx_cg_w.txt", delim_whitespace = True)

elif soil == "LV_A":
    ft =  pd.read_csv(r"~/code/python/artigoaj/lv_a_w.txt", delim_whitespace = True)


elif soil == "LV_Bw":
    ft =  pd.read_csv(r"~/code/python/artigoaj/lv_bw_w.txt", delim_whitespace = True)


elif soil == "LVA_A":
    ft =  pd.read_csv(r"~/code/python/artigoaj/lva_a_w.txt", delim_whitespace = True)


elif soil == "LVA_Bw":
    ft =  pd.read_csv(r"~/code/python/artigoaj/lva_bw_w.txt", delim_whitespace = True)


elif soil == "OX_H":
    ft =  pd.read_csv(r"~/code/python/artigoaj/ox_h_w.txt", delim_whitespace = True)


elif soil == "RR_A":
    ft =  pd.read_csv(r"~/code/python/artigoaj/rr_a_w.txt", delim_whitespace = True)

else:
    print("Error")

ft['psi'] = ft['psi']*10


### For some reason, the segmented nonlinear regression procedure does not work with the original data, it has to be converted to array
psi = np.array(ft['psi'])
theta = np.array(ft['theta'])


#inspect data

#plt.errorbar(ft['psi'], ft['theta'], yerr = ft['sd'], fmt = "ok", ecolor = "r", fillstyle = "none", capsize = 2.0)
#plt.ylabel("Volumetric water content (v/v)")
#plt.xlabel("Water potential (kPa)")
#plt.errorbar(x, y, e, linestyle='None', marker='^')
#plt.ylim([0.0, 1.0])
#plt.xscale('log')
#plt.show()


# macro for significance
def sigfun(x):
    if x > 0.05:
        return "ns"
    elif (x <= 0.05 and x > 0.01):
        return "*"
    elif (x <= 0.01 and x > 0.001):
        return "**"
    else:
        return "***" 


### Define saturated water content for all stupid equations
if soil == 'GX_A':
    qsat = 0.675
    qair = 0.023 
elif soil == 'GX_Cg':
    qsat = 0.689
    qair = 0.014
elif soil == 'LV_A':
    qsat = 0.650
    qair = 0.015 
elif soil == 'LV_Bw':
    qsat = 0.668
    qair = 0.016 
elif soil == 'LVA_A':
    qsat = 0.718
    qair = 0.014 
elif soil == 'LVA_Bw':
    qsat = 0.690
    qair = 0.013 
elif soil == 'OX_H':
    qsat = 0.791
    qair = 0.052 
elif soil == 'RR_A':
    qsat = 0.554
    qair = 0.005 
else:
    qsat = 'aa'
    qair = 'bb' 

print(qsat, qair)




print("#########################################################################")
print("#########################################################################")
print("#########################################################################")

################# Brooks-Corey

print("Brooks-Corey ################################################")

# vector of initial parameter values

thetas = 0.5
thetar = 0.1
psib = 20
lamb = 1.8

InitPar = [thetas, thetar, psib, lamb]

# function to be fit - change for a custom function
def func(x, qs, qr, p, l):
    return np.piecewise(x, [(x <= p), (x > p)], [lambda x: qs, lambda x: qr + (qs - qr)*(p/x)**l])


# fits data

popt, pcov = curve_fit(func, psi, theta, p0 = InitPar, method = 'lm') 

# This is necessary because the piecewise function is not passing the argument to further calculations

thetavg = np.where((psi >= popt[2]), popt[1] + (popt[0] - popt[1])*(popt[2]/psi)**popt[3], popt[0])


# degrees of freedom n - number of parameters
n = len(psi)
d = len(popt)
df = n-d

#print(thetavg)

# standard error
se0 = np.sqrt(pcov[0,0])
se1 = np.sqrt(pcov[1,1])
se2 = np.sqrt(pcov[2,2])
se3 = np.sqrt(pcov[3,3])

# t value
t0 = popt[0]/np.sqrt(pcov[0,0])
t1 = popt[1]/np.sqrt(pcov[1,1])
t2 = popt[2]/np.sqrt(pcov[2,2])
t3 = popt[3]/np.sqrt(pcov[3,3])

# significance
# uses the sf function, equivalent (but more precise) than 1 - cdf
# see scipy.stats.t documentation
p0 = 2*(t.sf(np.abs(t0), df, loc=0, scale =1))
p1 = 2*(t.sf(np.abs(t1), df, loc=0, scale =1))
p2 = 2*(t.sf(np.abs(t2), df, loc=0, scale =1))
p3 = 2*(t.sf(np.abs(t3), df, loc=0, scale =1))

print(round(popt[0],4),  round(se0,4), round(t0,4), np.format_float_scientific(p0, exp_digits=2, precision=4), sigfun(p0)) 
print(round(popt[1],4),  round(se1,4), round(t1,4), np.format_float_scientific(p1, exp_digits=2, precision=4), sigfun(p1)) 
print(round(popt[2],4),  round(se2,4), round(t2,4), np.format_float_scientific(p2, exp_digits=2, precision=4), sigfun(p2)) 
print(round(popt[3],4),  round(se3,4), round(t3,4), np.format_float_scientific(p3, exp_digits=2, precision=4), sigfun(p3)) 


# calculates the pseudo r2 from the correlation coefficient
r = np.corrcoef(theta, thetavg)
r2 = r[0,1]**2


# prints R2 to terminal
print("R2: ", round(r2,4))

# calculates and prints R2adj to terminal
r2adj = 1- (1 - r2) * (n-1)/(df)
print("R2adj: ", round(r2adj ,4))



# calculates the root mean square error - rmse
RMSE = np.sqrt(np.sum((theta - thetavg)**2)/df)
print("RMSE: ", round(RMSE, 4))

# calculates sum of squares of deviations SSR and prints to terminal
SSR = np.sum((theta - thetavg)**2)
print("SSR: ", round(SSR, 4))

# calculates the Akaike Information Criterion - AIC 
AIC = n * np.log(SSR/n) + 2 * d 
print("AIC: ", round(AIC, 4))

# calculates the corrected Akaike Information Criterion - AICc 
AICc = n * np.log(SSR/n) + 2 * d + (2*d*(d+1))/(n - d -1) 
print("AICc: ", round(AICc, 4))

# calculates the Bayesian Information Criterion - BIC 
BIC = n * np.log(SSR/n) + d * np.log(n)
print("BIC: ", round(BIC, 4))

# calculates the Kayshap Information Criterion - BIC

sigma = np.linalg.det(pcov)

KIC = n * np.log(SSR/n) - d * np.log((2*np.pi)) - np.log(sigma) 
print("KIC: ", round(KIC, 4))



################# van Genuchten

print("van Genuchten  ################################################")

# vector of initial parameter values

qs = qsat
qr = .5 
alpha = .2
n = 1.2

InitPar0 = [qr, alpha, n]

# function to be fit - change for a custom function
def func0(x, qr,  alpha, n):
    return qr + (qs - qr) / ( 1 + (alpha*x)**n ) ** ( 1 - 1 / n)   

# fits data

popt0, pcov0 = curve_fit(func0, ft['psi'], ft['theta'], p0 = InitPar0, method = 'lm', maxfev=5000) 

# degrees of freedom n - number of parameters
n0 = len(ft['psi'])
d0 = len(popt0)
df0 = n0-d0

# standard error
se00 = np.sqrt(pcov0[0,0])
se10 = np.sqrt(pcov0[1,1])
se20 = np.sqrt(pcov0[2,2])
#se30 = np.sqrt(pcov0[3,3])

# t value
t00 = popt0[0]/np.sqrt(pcov0[0,0])
t10 = popt0[1]/np.sqrt(pcov0[1,1])
t20 = popt0[2]/np.sqrt(pcov0[2,2])
#t30 = popt0[3]/np.sqrt(pcov0[3,3])

# significance
# uses the sf function, equivalent (but more precise) than 1 - cdf
# see scipy.stats.t documentation
p00 = 2*(t.sf(np.abs(t00), df0, loc=0, scale =1))
p10 = 2*(t.sf(np.abs(t10), df0, loc=0, scale =1))
p20 = 2*(t.sf(np.abs(t20), df0, loc=0, scale =1))
#p30 = 2*(t.sf(np.abs(t30), df0, loc=0, scale =1))

print(round(popt0[0],4),  round(se00,4), round(t00,4), np.format_float_scientific(p00, exp_digits=2, precision=4), sigfun(p00)) 
print(round(popt0[1],4),  round(se10,4), round(t10,4), np.format_float_scientific(p10, exp_digits=2, precision=4), sigfun(p10)) 
print(round(popt0[2],4),  round(se20,4), round(t20,4), np.format_float_scientific(p20, exp_digits=2, precision=4), sigfun(p20)) 
#print(round(popt0[3],4),  round(se30,4), round(t30,4), np.format_float_scientific(p30, exp_digits=2, precision=4), sigfun(p30)) 

# calculates the pseudo r2 from the correlation coefficient
r0 = np.corrcoef(ft['theta'], func0(ft['psi'], *popt0))
r20 = r0[0,1]**2

# prints R2 to terminal
print("R2: ", round(r20,4))

# calculates and prints R2adj to terminal
r2adj0 = 1- (1 - r20) * (n0-1)/(df0)
print("R2adj: ", round(r2adj0 ,4))

# calculates the root mean square error - rmse
RMSE0 = np.sqrt(np.sum((ft['theta']- func0(ft['psi'], *popt0))**2)/df0)
print("RMSE: ", round(RMSE0, 4))

# calculates sum of squares of deviations SSR and prints to terminal
SSR0 = np.sum((ft['theta']- func0(ft['psi'], *popt0))**2)
print("SSR: ", round(SSR0, 4))

# calculates the Akaike Information Criterion - AIC 
AIC0 = n0 * np.log(SSR0/n0) + 2 * d0 
print("AIC: ", round(AIC0, 4))

# calculates the corrected Akaike Information Criterion - AICc 
AICc0 = n0 * np.log(SSR0/n0) + 2 * d0 + (2*d0*(d0+1))/(n0 - d0 -1) 
print("AICc: ", round(AICc0, 4))

# calculates the Bayesian Information Criterion - BIC 
BIC0 = n0 * np.log(SSR0/n0) + d0 * np.log(n0)
print("BIC: ", round(BIC0, 4))

# calculates the Kayshap Information Criterion - BIC

sigma0 = np.linalg.det(pcov0)

KIC0 = n0 * np.log(SSR0/n0) - d0 * np.log((2*np.pi)) - np.log(sigma0) 
print("KIC: ", round(KIC0, 4))



################## Kosugi


print("Kosugi  ################################################")

# vector of initial parameter values

qs1 = qsat
qr1 = 0.22 
hm = 2
sigmaK = 1.5

InitPar1 = [qr1, hm, sigmaK]

# function to be fit - change for a custom function
def func1(x, qr1, hm, sigmaK):
    return qr1 + (qs1 - qr1) *  1.0 / 2.0 * special.erfc( np.log( x / hm )/( np.sqrt(2) * sigmaK) )   

# fits data

popt1, pcov1 = curve_fit(func1, ft['psi'], ft['theta'], p0 = InitPar1, method = 'lm', maxfev=5000)  

# degrees of freedom n - number of parameters
n1 = len(ft['psi'])
d1 = len(popt1)
df1 = n1-d1

# standard error
se01 = np.sqrt(pcov1[0,0])
se11 = np.sqrt(pcov1[1,1])
se21 = np.sqrt(pcov1[2,2])
#se31 = np.sqrt(pcov1[3,3])

# t value
t01 = popt1[0]/np.sqrt(pcov1[0,0])
t11 = popt1[1]/np.sqrt(pcov1[1,1])
t21 = popt1[2]/np.sqrt(pcov1[2,2])
#t31 = popt1[3]/np.sqrt(pcov1[3,3])

# significance
# uses the sf function, equivalent (but more precise) than 1 - cdf
# see scipy.stats.t documentation
p01 = 2*(t.sf(np.abs(t01), df1, loc=0, scale =1))
p11 = 2*(t.sf(np.abs(t11), df1, loc=0, scale =1))
p21 = 2*(t.sf(np.abs(t21), df1, loc=0, scale =1))
#p31 = 2*(t.sf(np.abs(t31), df1, loc=0, scale =1))

print(round(popt1[0],4),  round(se01,4), round(t01,4), np.format_float_scientific(p01, exp_digits=2, precision=4), sigfun(p01)) 
print(round(popt1[1],4),  round(se11,4), round(t11,4), np.format_float_scientific(p11, exp_digits=2, precision=4), sigfun(p11)) 
print(round(popt1[2],4),  round(se21,4), round(t21,4), np.format_float_scientific(p21, exp_digits=2, precision=4), sigfun(p21)) 
#print(round(popt1[3],4),  round(se31,4), round(t31,4), np.format_float_scientific(p31, exp_digits=2, precision=4), sigfun(p31)) 


# calculates the pseudo r2 from the correlation coefficient
r1 = np.corrcoef(ft['theta'], func1(ft['psi'], *popt1))
r21 = r1[0,1]**2

# prints R2 to terminal
print("R2: ", round(r21,4))

# prints R2adj to terminal
r2adj1 = 1- (1 - r21) * (n1-1)/(df1)
print("R2adj: ", round(r2adj1 ,4))

# calculates the root mean square error - rmse
RMSE1 = np.sqrt(np.sum((ft['theta']- func1(ft['psi'], *popt1))**2)/df1)
print("RMSE: ", round(RMSE1, 4))

# calculates sum of squares of deviations SSR and prints to terminal
SSR1 = np.sum((ft['theta']- func1(ft['psi'], *popt1))**2)
print("SSR: ", round(SSR1, 4))

# calculates the Akaike Information Criterion - AIC 
AIC1 = n1 * np.log(SSR1/n1) + 2 * d1 
print("AIC: ", round(AIC1, 4))

# calculates the corrected Akaike Information Criterion - AICc 
AICc1 = n1 * np.log(SSR1/n1) + 2 * d1 + (2*d1*(d1+1))/(n1 - d1 -1) 
print("AICc: ", round(AICc1, 4))

# calculates the Bayesian Information Criterion - BIC 
BIC1 = n1 * np.log(SSR1/n1) + d1 * np.log(n1)
print("BIC: ", round(BIC1, 4))

# calculates the Kayshap Information Criterion - BIC

sigma1 = np.linalg.det(pcov1)

KIC1 = n1 * np.log(SSR1/n1) - d1 * np.log((2*np.pi)) - np.log(sigma1) 
print("KIC: ", round(KIC1, 4))


print("Rieu-Sposito  ################################################")

# vector of initial parameter values

qmax = 0.5
D = 2.8
hmin = 20
A = 1
InitPar2 = [qmax, hmin, D]

# function to be fit - change for a custom function
def func2(x, qmax, hmin, D):
    return qmax - ( 1 -  (hmin/x)**(3-D ) )   

# fits data

popt2, pcov2 = curve_fit(func2, ft['psi'], ft['theta'], p0 = InitPar2, method = 'lm') 

# degrees of freedom n - number of parameters
n2 = len(ft['psi'])
d2 = len(popt2)
df2 = n2-d2


# standard error
se02 = np.sqrt(pcov2[0,0])
se12 = np.sqrt(pcov2[1,1])
se22 = np.sqrt(pcov2[2,2])
#se32 = np.sqrt(pcov2[3,3])

# t value
t02 = popt2[0]/np.sqrt(pcov2[0,0])
t12 = popt2[1]/np.sqrt(pcov2[1,1])
t22 = popt2[2]/np.sqrt(pcov2[2,2])
#t32 = popt2[3]/np.sqrt(pcov2[3,3])

# significance
# uses the sf function, equivalent (but more precise) than 1 - cdf
# see scipy.stats.t documentation
p02 = 2*(t.sf(np.abs(t02), df2, loc=0, scale =1))
p12 = 2*(t.sf(np.abs(t12), df2, loc=0, scale =1))
p22 = 2*(t.sf(np.abs(t22), df2, loc=0, scale =1))
#p32 = 2*(t.sf(np.abs(t32), df2, loc=0, scale =1))

print(round(popt2[0],4),  round(se02,4), round(t02,4), np.format_float_scientific(p02, exp_digits=2, precision=4), sigfun(p02)) 
print(round(popt2[1],4),  round(se12,4), round(t12,4), np.format_float_scientific(p12, exp_digits=2, precision=4), sigfun(p12)) 
print(round(popt2[2],4),  round(se22,4), round(t22,4), np.format_float_scientific(p22, exp_digits=2, precision=4), sigfun(p22)) 
#print(round(popt2[3],4),  round(se32,4), round(t32,4), np.format_float_scientific(p32, exp_digits=2, precision=4), sigfun(p32)) 


# calculates the pseudo r2 from the correlation coefficient
r_2 = np.corrcoef(ft['theta'], func2(ft['psi'], *popt2))
r22 = r_2[0,1]**2

# prints R2 to terminal
print("R2: ", round(r22,4))

# prints R2adj to terminal
r2adj2 = 1- (1 - r22) * (n2-1)/(df2)
print("R2adj: ", round(r2adj2 ,4))


# calculates the root mean square error - rmse
RMSE2 = np.sqrt(np.sum((ft['theta']- func2(ft['psi'], *popt2))**2)/df2)
print("RMSE: ", round(RMSE2, 4))

# calculates sum of squares of deviations SSR and prints to terminal
SSR2 = np.sum((ft['theta']- func2(ft['psi'], *popt2))**2)
print("SSR: ", round(SSR2, 4))


# calculates the Akaike Information Criterion - AIC 
AIC2 = n2 * np.log(SSR2/n2) + 2 * d2 
print("AIC: ", round(AIC2, 4))

# calculates the corrected Akaike Information Criterion - AICc 
AICc2 = n2 * np.log(SSR2/n2) + 2 * d2 + (2*d2*(d2+1))/(n2 - d2 - 1) 
print("AICc: ", round(AICc2, 4))

# calculates the Bayesian Information Criterion - BIC 
BIC2 = n2 * np.log(SSR2/n2) + d2 * np.log(n2)
print("BIC: ", round(BIC2, 4))

# calculates the Kayshap Information Criterion - BIC

sigma2 = np.linalg.det(pcov2)

KIC2 = n2 * np.log(SSR2/n2) - d2 * np.log((2*np.pi)) - np.log(sigma2) 
print("KIC: ", round(KIC2, 4))


# simulated range for predicted plot
Sim = np.arange(np.min(ft['psi']), np.max(ft['psi']*10), np.abs(np.max(ft['psi']) - np.min(ft['psi']))/10000)

# plot - define parameters according to user need

#plt.plot(ft['psi'], ft['theta'], "o")
plt.errorbar(ft['psi'], ft['theta'], yerr = ft['sd'], fmt = "ok", ecolor = "r", fillstyle = "none", capsize = 2.0)
plt.plot(Sim, func(Sim, *popt), linestyle = 'solid', linewidth = 1.2, )
plt.plot(Sim, func0(Sim, *popt0), linestyle = 'dotted', linewidth = 1.2)
plt.plot(Sim, func1(Sim, *popt1), linestyle = 'dashed', linewidth = 1.2)
plt.plot(Sim, func2(Sim, *popt2), linestyle = (0, (3,1,1,1,1)), linewidth = 1.2)
#plt.plot(Sim, func2(Sim, *popt2), linestyle = 'dashdot', linewidth = 1.2, color = "red")


#color = 'gray'
#plt.legend(['EG', 'VG', 'MSN', 'FX'] )
plt.legend(['BC', 'VG', 'Kos', 'RS'], frameon = False)
#plt.legend(['BC', 'VG', 'RS'], frameon = False)

plt.ylabel("Volumetric water content (m$^3$ m$^{-3}$)")
plt.ylim([0.0, 0.8])
plt.xlabel("Water potential (|hPa|)")
plt.xscale('log')
###plt.savefig('/home/t01/1_2024/ArtigosTPL/ArtigoAnaJulia/ManuscriptAJ/RR_A_w.pdf', format='pdf', dpi=300)

plt.show()

print("Fitting parameters")
print("BC", popt)
print("Kos", popt1)
print("R-S", popt2)
print("VG", popt0)
#print("MSN", popt1)
#print("FX", popt2)



print("Statistics", soil)
print("Model", "R2    " , "R2adj  " ,"RMSE  ", "AIC  ", "AICc   ", "BIC     ", "KIC     ")
print("BC   ", round(r2,4) , round(r2adj,4) , round(RMSE,4),  round(AIC,1),  round(AICc,1),  round(BIC,1),  round(KIC,1))
print("Kos  ", round(r21,4), round(r2adj1,4), round(RMSE1,4), round(AIC1,1), round(AICc1,1), round(BIC1,1), round(KIC1,1))
print("R-S  ", round(r22,4), round(r2adj2,4), round(RMSE2,4), round(AIC2,1), round(AICc2,1), round(BIC2,1), round(KIC2,1))
print("VG   ", round(r20,4), round(r2adj0,4), round(RMSE0,4), round(AIC0,1), round(AICc0,1), round(BIC0,1), round(KIC0,1))















