#=====================================================
import numpy as np
import pandas as pd
import mplcursors
from matplotlib.backend_bases import key_press_handler
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#=====================================================
#read data from excel file, combine according to date
guiinput1 = np.ones(32);

dt_fulldata = pd.read_csv('sird_data.csv')
dt_fulldata['date'] = pd.to_datetime(dt_fulldata['date'])
dt_fulldata = dt_fulldata.set_index(['date'])
dt_fulldata = dt_fulldata.tail(dt_fulldata.shape[0] - 8)
Ptotal = int(3.153*10**7);
#=====================================================
#assign initial value
T = len(dt_fulldata);
t0 = 0;
N = T;
M = 500;

from datetime import timedelta as td, datetime as datetime
# initialdatadate = date(2021,1,1)
initialdatadate = dt_fulldata.index.tolist()[0];
sdate1 = initialdatadate.date();
edate1 = sdate1 + td(days = T + 1)
pastplot_t = pd.date_range(sdate1,edate1-td(days=1),freq='d')
sdate2 = edate1 + td(days = 1)
edate2 = sdate2 + td(days = T + 1)
futplot_t = pd.date_range(sdate2,edate2-td(days=1),freq='d')
fullplot_t = pd.date_range(sdate1,edate2-td(days=1),freq='d')

def datetoint(date):
    try:
        datedate = datetime.strptime(date,'%Y-%m-%d');
        strdate = datedict[datedate];
        intdate = int(strdate);
        return intdate
    except:
        return 0;
    
def NPIoptiontoint(guientry):
    if guientry == "Strict NPI":
        return  2;
    elif guientry == "Loose NPI":
        return  1;
    else:
        return  0;

datedict = {}
diff = edate2 - sdate1
for i in range(diff.days + 1):
    day = sdate1 + td(days=i)
    datedict[datetime.strptime(str(day),"%Y-%m-%d")] = i

Izero = dt_fulldata['cases_active'].iloc[0];
Rzero = dt_fulldata['cases_recovered'].iloc[0]; 
Dzero = dt_fulldata['deaths_total'].iloc[0]; 

Szero = (Ptotal - Izero - Rzero - Dzero);
dt = T/N;
t = np.linspace(0,N,N+1);
data_t = np.linspace(0,T,T+1);
futureplot = np.linspace(T,2*T,T+1)
dt_mcmcresult = pd.read_csv('mcmc_result_2.csv')
dt_mcmcresult = dt_mcmcresult.dropna();
mcmcsize = np.size(dt_mcmcresult, axis=0);

#=====================================================
#parameter calculations for future approximations
segment_NPI_durations = dt_mcmcresult['t'].diff().dropna();

# mean_beta1 = sum(dt_mcmcresult['beta1'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_beta2 = sum(dt_mcmcresult['beta2'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_beta3 = sum(dt_mcmcresult['beta3'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

# mean_r1 = sum(dt_mcmcresult['r1'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_r2 = sum(dt_mcmcresult['r2'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_r3 = sum(dt_mcmcresult['r3'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

# mean_d1 = sum(dt_mcmcresult['d1'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_d2 = sum(dt_mcmcresult['d2'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_d3 = sum(dt_mcmcresult['d3'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

# mean_gamma = sum(dt_mcmcresult['gamma'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_sigma = sum(dt_mcmcresult['sigma'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

# mean_p1 = sum(dt_mcmcresult['p1'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_p2 = sum(dt_mcmcresult['p2'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];
# mean_p3 = sum(dt_mcmcresult['p3'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

#Frequency grouping algorithm (Infectious rate)
#=====================================================
betaperday = np.zeros((T));
for i in np.linspace(1,T-1,T-1):
    for j in np.linspace(0,mcmcsize-1,mcmcsize-1):
        if i > dt_mcmcresult['t'].iloc[j.astype(int)]:
            betaperday[int(i)] = dt_mcmcresult['beta'].iloc[j.astype(int)]
# finding beta for each day
hist_beta = np.histogram(betaperday, bins = 25);
# finding frequency of certain value of beta
binsize = np.size(hist_beta[1])-1;
midpoint_beta = np.zeros((binsize,1));
for i in np.linspace(1,binsize,binsize):
    midpoint_beta[i.astype(int)-1] = (hist_beta[1][i.astype(int)]+ hist_beta[1][i.astype(int)-1])/2;
beta_freq = hist_beta[0].astype(int);
# plt.subplot(311)
# plt.plot(midpoint_beta,beta_freq)
beta_peaks = (np.diff(np.sign(np.diff(beta_freq))) < 0).nonzero()[0] + 1
# plt.plot(midpoint_beta[beta_peaks], beta_freq[beta_peaks], "o", color='b')
# plt.ylim([0, 1])
#=====================================================
#Trend approximation (Recovery and Fatality rate)
#=====================================================
def exponential(x, a, b):
    return a*np.exp(b*x)
rperday = np.zeros((T));
for i in np.linspace(1,T-1,T-1):
    for j in np.linspace(0,mcmcsize-1,mcmcsize-1):
        if i > dt_mcmcresult['t'].iloc[j.astype(int)]:
            rperday[int(i)] = dt_mcmcresult['r'].iloc[j.astype(int)]
# plt.subplot(312)
# plt.plot(np.linspace(1,T,T),rperday)
pars, cov = curve_fit(f=exponential, xdata=np.linspace(1,T,T), ydata=rperday, p0=[0, 0], bounds=(-np.inf, np.inf))
# plt.plot(np.linspace(1,T,T), exponential(np.linspace(1,T,T), *pars), linestyle='--', linewidth=1, color='black')
#=====================================================
dperday = np.zeros((T));
for i in np.linspace(1,T-1,T-1):
    for j in np.linspace(0,mcmcsize-1,mcmcsize-1):
        if i > dt_mcmcresult['t'].iloc[j.astype(int)]:
            dperday[int(i)] = dt_mcmcresult['d'].iloc[j.astype(int)]
# plt.subplot(313)
# plt.plot(np.linspace(1,T,T),dperday)
pars2, cov2 = curve_fit(f=exponential, xdata=np.linspace(1,T,T), ydata=dperday, p0=[0, 0], bounds=(-np.inf, np.inf))
# plt.plot(np.linspace(1,T,T), exponential(np.linspace(1,T,T), *pars2), linestyle='--', linewidth=1, color='black')
#=====================================================
mean_beta = sum(dt_mcmcresult['beta'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1];

mean_r = dt_mcmcresult['r'].iloc[-1];

mean_d = dt_mcmcresult['d'].iloc[-1];

mean_gamma = dt_mcmcresult['gamma'].iloc[-1];
mean_sigma = sum(dt_mcmcresult['sigma'][1:]*segment_NPI_durations)/dt_mcmcresult['t'].iloc[-1]+0.1;

#=====================================================
#read parameter from mcmc, assign according to step in SRK4
def parameter_func(j,guiinput):
    Beta = 0;
    r = 0;
    d = 0;
    gamma,sigma = 0, 0;
    
    [display_entry1,display_entry2,display_entry3,display_entry4,display_entry5,display_entry6,display_entry7,settings_entry1,userinput_entry11,userinput_entry12,userinput_entry13,userinput_entry21,userinput_entry22,userinput_entry23,userinput_entry31,userinput_entry32,userinput_entry33,userinput_entry41,userinput_entry42,userinput_entry43,userinput_entry51,userinput_entry52,userinput_entry53,userinput_entry61,userinput_entry62,userinput_entry63,userinput_entry71,userinput_entry72,userinput_entry73,userinput_entry81,userinput_entry82,userinput_entry83] = guiinput
    
    if j < T:
        Beta = dt_mcmcresult['beta'].iloc[j.astype(int)-1];
        r = dt_mcmcresult['r'].iloc[j.astype(int)-1];
        d = dt_mcmcresult['d'].iloc[j.astype(int)-1];
        gamma = dt_mcmcresult['gamma'].iloc[j.astype(int)-1];
        sigma = dt_mcmcresult['sigma'].iloc[j.astype(int)-1];
        if settings_entry1 == 1:
            sigma = sigma*0;
    else:
        npidatestart1 = datetoint(userinput_entry11);
        npidateend1   = datetoint(userinput_entry12);
        npidatestart2 = datetoint(userinput_entry21);
        npidateend2   = datetoint(userinput_entry22);
        npidatestart3 = datetoint(userinput_entry31);
        npidateend3   = datetoint(userinput_entry32);
        npidatestart4 = datetoint(userinput_entry41);
        npidateend4   = datetoint(userinput_entry42);
        npidatestart5 = datetoint(userinput_entry51);
        npidateend5   = datetoint(userinput_entry52);
        npidatestart6 = datetoint(userinput_entry61);
        npidateend6   = datetoint(userinput_entry62);
        npidatestart7 = datetoint(userinput_entry71);
        npidateend7   = datetoint(userinput_entry72);
        npidatestart8 = datetoint(userinput_entry81);
        npidateend8   = datetoint(userinput_entry82);
        
        futbetavalues = [[midpoint_beta[beta_peaks][2],0.09495,0.1052],
                         [midpoint_beta[beta_peaks][2],0.032846,0.0363214],
                         [midpoint_beta[beta_peaks][1],0.027115,0.029983]];

        if settings_entry1 == 1:
            sigma = mean_sigma*0;
        else:
            sigma = mean_sigma;
               
        if j > npidatestart1 and j <= npidateend1:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry13)]
        if j > npidatestart2 and j <= npidateend2:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry23)]
        if j > npidatestart3 and j <= npidateend3:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry33)]
        if j > npidatestart4 and j <= npidateend4:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry43)]
        if j > npidatestart5 and j <= npidateend5:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry53)]
        if j > npidatestart6 and j <= npidateend6:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry63)]
        if j > npidatestart7 and j <= npidateend7:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry73)]
        if j > npidatestart8 and j <= npidateend8:
            Beta, Beta2, Beta3 = futbetavalues[NPIoptiontoint(userinput_entry83)]
            
        r, d, gamma, sigma = mean_r, mean_d, mean_gamma, sigma
    return Beta, r, d, gamma, sigma
    
#=====================================================
#SRK-4 Tableus form
A = np.array([[0  ,0  ,0,0],
              [1/2,0  ,0,0],
              [0  ,1/2,0,0],
              [0  ,0  ,1,0]]);
B_1 = np.array([[0            ,0            ,0          ,0],
                [-0.7242916356,0            ,0          ,0],
                [0.4237353406 ,-0.1994437050,0          ,0],
                [-1.578475506 ,0.840100343  ,1.738375163,0]]); 
B_2 = np.array([[0           ,0,0,0],
                [2.702000410 ,0,0,0],
                [1.757261649 ,0,0,0],
                [-2.918524118,0,0,0]]); 
G_1 = np.array([[-.7800788474],
                [0.07363768240],
                [1.486520013],
                [0.2199211524]]); 
G_2 = np.array([[1.693950844],
                [1.636107882],
                [-3.024009558],
                [-0.3060491602]]); 
# alpha = np.array([1/6, 1/3, 1/3, 1/6]);

def SIRD_func(t0,dt,N,M,Szero,Izero,Rzero,Dzero,Ptotal,guiinput):
    #=====================================================
    fS  = lambda S,I,R,D,Beta,r,d,gamma,sigma:0-S/Ptotal*Beta*I+gamma*R   
    gS  = lambda S,I,R,D,Beta,r,d,gamma,sigma:-sigma*I*S/Ptotal
    #=====================================================
    fI  = lambda S,I,R,D,Beta,r,d,gamma,sigma:S/Ptotal*Beta*I-(r+d)*I
    gI  = lambda S,I,R,D,Beta,r,d,gamma,sigma:sigma*I*S/Ptotal
    #=====================================================
    fR  = lambda S,I,R,D,Beta,r,d,gamma,sigma:r*I-gamma*R
    gR  = lambda S,I,R,D,Beta,r,d,gamma,sigma:0
    #=====================================================
    fD  = lambda S,I,R,D,Beta,r,d,gamma,sigma:d*I
    gD  = lambda S,I,R,D,Beta,r,d,gamma,sigma:0
    #=====================================================
        
    SdW1 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    SdW2 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    IdW1 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    IdW2 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    RdW1 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    RdW2 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    DdW1 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
    DdW2 = np.random.standard_normal((N+1,M+1))*np.sqrt(dt);
        
    Stemp = Szero*np.ones((1,M+1));
    Itemp = Izero*np.ones((1,M+1));
    Rtemp = Rzero*np.ones((1,M+1));
    Dtemp = Dzero*np.ones((1,M+1));
    
    Ssrk4_2 = np.ones((N+1,M+1));
    Isrk4_2 = np.ones((N+1,M+1));
    Rsrk4_2 = np.ones((N+1,M+1));
    Dsrk4_2 = np.ones((N+1,M+1));
    #=====================================================
        
    for j in np.linspace(t0,t0+N,N+1):
        Beta,r,d,gamma,sigma = parameter_func(j,guiinput);
        SWinc1 = SdW1[int(j) - t0];
        IWinc1 = IdW1[int(j) - t0];
        RWinc1 = RdW1[int(j) - t0];
        DWinc1 = DdW1[int(j) - t0];
        SWinc2 = SdW2[int(j) - t0];
        IWinc2 = IdW2[int(j) - t0];
        RWinc2 = RdW2[int(j) - t0];
        DWinc2 = DdW2[int(j) - t0];
        
        SJ10 = (0.5*(dt**(1/2))*(SWinc1+(1/np.sqrt(3))*SWinc2));
        IJ10 = (0.5*(dt**(1/2))*(IWinc1+(1/np.sqrt(3))*IWinc2));
        RJ10 = (0.5*(dt**(1/2))*(RWinc1+(1/np.sqrt(3))*RWinc2));
        DJ10 = (0.5*(dt**(1/2))*(DWinc1+(1/np.sqrt(3))*DWinc2));
        
        S1 = Stemp;
        I1 = Itemp;
        R1 = Rtemp;
        D1 = Dtemp;
        
        f1_S = fS(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        g1_S = gS(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        f1_I = fI(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        g1_I = gI(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        f1_R = fR(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        g1_R = gR(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        f1_D = fD(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        g1_D = gD(S1,I1,R1,D1,Beta,r,d,gamma,sigma);
        
        S2 = Stemp + dt*(A[0,1]*f1_S) + (B_1[0,1]*SWinc1 + B_2[0,1]*SJ10)*g1_S;
        I2 = Itemp + dt*(A[0,1]*f1_I) + (B_1[0,1]*IWinc1 + B_2[0,1]*IJ10)*g1_I;
        R2 = Rtemp + dt*(A[0,1]*f1_R) + (B_1[0,1]*RWinc1 + B_2[0,1]*RJ10)*g1_R;
        D2 = Dtemp + dt*(A[0,1]*f1_D) + (B_1[0,1]*DWinc1 + B_2[0,1]*DJ10)*g1_D;
        
        f2_S = fS(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        g2_S = gS(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        f2_I = fI(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        g2_I = gI(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        f2_R = fR(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        g2_R = gR(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        f2_D = fD(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        g2_D = gD(S2,I2,R2,D2,Beta,r,d,gamma,sigma);
        
        S3 = Stemp + dt*(A[1,2]*f2_S) + (B_1[0,2]*SWinc1 + B_2[0,2]*SJ10)*g2_S + (B_1[1,2]*SJ10*g2_S);
        I3 = Itemp + dt*(A[1,2]*f2_I) + (B_1[0,2]*IWinc1 + B_2[0,2]*IJ10)*g2_I + (B_1[1,2]*IJ10*g2_I);
        R3 = Rtemp + dt*(A[1,2]*f2_R) + (B_1[0,2]*RWinc1 + B_2[0,2]*RJ10)*g2_R + (B_1[1,2]*RJ10*g2_R);
        D3 = Dtemp + dt*(A[1,2]*f2_D) + (B_1[0,2]*DWinc1 + B_2[0,2]*DJ10)*g2_D + (B_1[1,2]*DJ10*g2_D);
        # V3 = Vtemp + dt*(A[1,2]*f2_V) + (B_1[0,2]*VWinc1 + B_2[0,2]*VJ10)*g2_V + (B_1[1,2]*VJ10*g2_V);
        
        f3_S = fS(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        g3_S = gS(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        f3_I = fI(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        g3_I = gI(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        f3_R = fR(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        g3_R = gR(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        f3_D = fD(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        g3_D = gD(S3,I3,R3,D3,Beta,r,d,gamma,sigma);
        
        S4 = Stemp + dt*f3_S + (B_1[0,3]*SWinc1 + B_2[0,3]*SJ10)*g3_S + (B_1[1,3]*SWinc1*g3_S + B_1[2,3]*SWinc1*g3_S);
        I4 = Itemp + dt*f3_I + (B_1[0,3]*IWinc1 + B_2[0,3]*IJ10)*g3_I + (B_1[1,3]*IWinc1*g3_I + B_1[2,3]*IWinc1*g3_I);
        R4 = Rtemp + dt*f3_R + (B_1[0,3]*RWinc1 + B_2[0,3]*RJ10)*g3_R + (B_1[1,3]*RWinc1*g3_R + B_1[2,3]*RWinc1*g3_R);
        D4 = Dtemp + dt*f3_D + (B_1[0,3]*DWinc1 + B_2[0,3]*DJ10)*g3_D + (B_1[1,3]*DWinc1*g3_D + B_1[2,3]*DWinc1*g3_D);
        # V4 = Vtemp + dt*f3_V + (B_1[0,3]*VWinc1 + B_2[0,3]*VJ10)*g3_V + (B_1[1,3]*VWinc1*g3_V + B_1[2,3]*VWinc1*g3_V);
        
        f4_S = fS(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        g4_S = gS(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        f4_I = fI(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        g4_I = gI(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        f4_R = fR(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        g4_R = gR(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        f4_D = fD(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        g4_D = gD(S4,I4,R4,D4,Beta,r,d,gamma,sigma);
        
        Stemp = Stemp + dt*((1/6)*f1_S + (1/3)*f2_S + (1/3)*f3_S + (1/6)*f4_S) + (G_1[0]*SWinc1 + G_2[0]*SJ10)*g1_S + (G_1[1]*SWinc1 + G_2[1]*SJ10)*g2_S + (G_1[2]*SWinc1 + G_2[2]*SJ10)*g3_S + (G_1[3]*SWinc1 + G_2[3]*SJ10)*g4_S;
        Itemp = Itemp + dt*((1/6)*f1_I + (1/3)*f2_I + (1/3)*f3_I + (1/6)*f4_I) + (G_1[0]*IWinc1 + G_2[0]*IJ10)*g1_I + (G_1[1]*IWinc1 + G_2[1]*IJ10)*g2_I + (G_1[2]*IWinc1 + G_2[2]*IJ10)*g3_I + (G_1[3]*IWinc1 + G_2[3]*IJ10)*g4_I;
        Rtemp = Rtemp + dt*((1/6)*f1_R + (1/3)*f2_R + (1/3)*f3_R + (1/6)*f4_R) + (G_1[0]*RWinc1 + G_2[0]*RJ10)*g1_R + (G_1[1]*RWinc1 + G_2[1]*RJ10)*g2_R + (G_1[2]*RWinc1 + G_2[2]*RJ10)*g3_R + (G_1[3]*RWinc1 + G_2[3]*RJ10)*g4_R;
        Dtemp = Dtemp + dt*((1/6)*f1_D + (1/3)*f2_D + (1/3)*f3_D + (1/6)*f4_D) + (G_1[0]*DWinc1 + G_2[0]*DJ10)*g1_D + (G_1[1]*DWinc1 + G_2[1]*DJ10)*g2_D + (G_1[2]*DWinc1 + G_2[2]*DJ10)*g3_D + (G_1[3]*DWinc1 + G_2[3]*DJ10)*g4_D;
        
        # for array in np.linspace(1,mcmcsize-1,mcmcsize):
        #     if j == dt_mcmcresult['t'].iloc[array.astype(int)]:
        #         Stemp  = Stemp  - (Stemp.mean(axis=1) -dt_fulldata['susceptible'].iloc[j.astype(int)])
        #         Itemp  = Itemp  - (Itemp.mean(axis=1) -dt_fulldata['cases_active'].iloc[j.astype(int)])
        #         Rtemp  = Rtemp  - (Rtemp.mean(axis=1) -dt_fulldata['cases_recovered'].iloc[j.astype(int)])
        #         Dtemp  = Dtemp  - (Dtemp.mean(axis=1) -dt_fulldata['deaths_total'].iloc[j.astype(int)])
        
        Ssrk4_2[int(j) - t0] = Stemp;
        Isrk4_2[int(j) - t0] = Itemp;
        Rsrk4_2[int(j) - t0] = Rtemp;
        Dsrk4_2[int(j) - t0] = Dtemp;
        
    return Ssrk4_2, Isrk4_2, Rsrk4_2, Dsrk4_2

Ssrk4_2, Isrk4_2, Rsrk4_2, Dsrk4_2 = SIRD_func(t0,dt,N,M,Szero,Izero,Rzero,Dzero,Ptotal,guiinput1)

#=====================================================

Save_srk4 = Ssrk4_2.mean(axis=1);
Iave_srk4 = Isrk4_2.mean(axis=1);
Rave_srk4 = Rsrk4_2.mean(axis=1);
Dave_srk4 = Dsrk4_2.mean(axis=1);

Idiff = np.zeros(T+1);

def rmse(y,data):
    relativeerr = y - data;
    return ((sum(relativeerr ** 2))/len(y))** .5

rmseI = rmse((Iave_srk4[:-1]),(dt_fulldata['cases_active']))
rmseR = rmse((Rave_srk4[:-1]),(dt_fulldata['cases_recovered']))
rmseD = rmse((Dave_srk4[:-1]),(dt_fulldata['deaths_total']))

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk
 
try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

def init(top, gui, *args, **kwargs):
    global w, top_level, root
    w = gui
    top_level = top
    root = top
 
def destroy_window():
    # Function which closes the window.
    global top_level
    top_level.destroy()
    top_level = None
 
class App:
    def __init__(self, master):
         
        global top_level
        top_level = root
         
        root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
        root.state('zoomed')
        root.title("Pandemic Preparedness Simulator")
        ##########################################################
        ## Main Graphs ###########################################
        self.Frame1 = tk.Frame(root)
        self.Frame1.place(relx=0.15, rely=0.058, relheight=0.69, relwidth=0.71)
        self.Frame1.configure(relief='groove')
        self.Frame1.configure(borderwidth="2")
        self.Frame1.configure(relief='groove')
        self.Frame1.configure(highlightcolor="black")
        self.Frame1.configure(width=805)
        
        self.TLabel1 = ttk.Label(root)
        self.TLabel1.place(relx=0.16, rely=0.045, height=24, width=83)
        self.TLabel1.configure(relief='flat')
        self.TLabel1.configure(text='''Epidemic Plot''')
        
        self.TLabelSignature = ttk.Label(root)
        self.TLabelSignature.place(relx=0.012, rely=0.85, height=140, width=210)
        self.TLabelSignature.configure(relief='flat',font=("Courier", 7),foreground="gray")
        self.TLabelSignature.configure(text='''Pandemic Preparedness Simulator v3.2
made by: 
System Developer, Model & Simulation:
    Muhammad Fahmi
Model Development & System Improvise:
    Norhayati Rosli
System Validation:
    Noryanti Muhammad
Data source: MoH_Malaysia_covid19_public
Parameter Estimate: MCMC_Toolbox_Matlab
GUI package: Tkinter''')
        ##########################################################
        ## Reproductive Number plot ##############################
        self.Frame2 = tk.Frame(root)
        self.Frame2.place(relx=0.15, rely=0.75, relheight=0.195, relwidth=0.71)
        self.Frame2.configure(relief='groove')
        self.Frame2.configure(borderwidth="2")
        self.Frame2.configure(relief='groove')
        self.Frame2.configure(highlightcolor="black")
        self.Frame2.configure(width=805)
        
        self.TLabel2 = ttk.Label(root)
        self.TLabel2.place(relx=0.16, rely=0.755, height=24, width=210)
        self.TLabel2.configure(relief='flat')
        self.TLabel2.configure(text='''Real-time Reproductive Number Plot''')
        
        ##########################################################
        ## Display ###############################################
        self.Frame3 = tk.Frame(root)
        self.Frame3.place(relx=0.012, rely=0.058, relheight=0.26, relwidth=0.133)
        self.Frame3.configure(relief='groove')
        self.Frame3.configure(borderwidth="2")
        self.Frame3.configure(relief='groove')
        self.Frame3.configure(highlightcolor="black")
        self.Frame3.configure(width=805)
        
        self.TLabel31 = ttk.Label(root)
        self.TLabel31.place(relx=0.06, rely=0.06, height=24, width=110)
        self.TLabel31.configure(relief='flat')
        self.TLabel31.configure(text='''DISPLAY''')
        
        self.TLabel32 = ttk.Label(root)
        self.TLabel32.place(relx=0.02, rely=0.09, height=22, width=110)
        self.TLabel32.configure(relief='flat')
        self.TLabel32.configure(text='''Susceptible''')
        
        self.TLabel33 = ttk.Label(root)
        self.TLabel33.place(relx=0.02, rely=0.115, height=22, width=110)
        self.TLabel33.configure(relief='flat')
        self.TLabel33.configure(text='''Infected''')
        
        # self.TLabel34 = ttk.Label(root)
        # self.TLabel34.place(relx=0.04, rely=0.14, height=22, width=110)
        # self.TLabel34.configure(relief='flat')
        # self.TLabel34.configure(text='''Total Infected''')
        
        self.TLabel35 = ttk.Label(root)
        self.TLabel35.place(relx=0.04, rely=0.165, height=22, width=110)
        self.TLabel35.configure(relief='flat')
        self.TLabel35.configure(text='''Quarantine''')
        
        # self.TLabel36 = ttk.Label(root)
        # self.TLabel36.place(relx=0.04, rely=0.190, height=22, width=110)
        # self.TLabel36.configure(relief='flat')
        # self.TLabel36.configure(text='''ICU''')
        
        # self.TLabel37 = ttk.Label(root)
        # self.TLabel37.place(relx=0.04, rely=0.215, height=22, width=110)
        # self.TLabel37.configure(relief='flat')
        # self.TLabel37.configure(text='''Ventilators''')
        
        self.TLabel38 = ttk.Label(root)
        self.TLabel38.place(relx=0.02, rely=0.240, height=22, width=110)
        self.TLabel38.configure(relief='flat')
        self.TLabel38.configure(text='''Recovered''')
        
        self.TLabel39 = ttk.Label(root)
        self.TLabel39.place(relx=0.02, rely=0.265, height=22, width=110)
        self.TLabel39.configure(relief='flat')
        self.TLabel39.configure(text='''Death''')
        
        self.CButtonvar3_1 = tk.IntVar()
        self.CButton3_1 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_1, command = self.updatePlot)
        self.CButton3_1.place(relx=0.9, rely=0.14,height=10, relwidth=0.1)
        self.CButton3_1.configure(highlightcolor="black")
        # self.CButton3_1.select();
        self.CButton3_1.bind("<Return>", self.updatePlot)
        
        self.CButtonvar3_2 = tk.IntVar()
        self.CButton3_2 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_2, command = self.updatePlot)
        self.CButton3_2.place(relx=0.9, rely=0.33,height=10, relwidth=0.1)
        self.CButton3_2.configure(highlightcolor="black")
        # self.CButton3_2.select();
                
        self.CButtonvar3_3 = tk.IntVar()
        self.CButton3_3 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_3, command = self.updatePlot)
        self.CButton3_3.place(relx=0.9, rely=0.43,height=10, relwidth=0.1)
        self.CButton3_3.configure(highlightcolor="black")
        self.CButton3_3.select();
        self.CButton3_3.bind("<Return>", self.updatePlot)
        
        self.CButtonvar3_4 = tk.IntVar()
        self.CButton3_4 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_4, command = self.updatePlot)
        self.CButton3_4.place(relx=0.9, rely=0.53,height=10, relwidth=0.1)
        self.CButton3_4.configure(highlightcolor="black")
        self.CButton3_4.select();
        self.CButton3_4.bind("<Return>", self.updatePlot)
        
        self.CButtonvar3_5 = tk.IntVar()
        self.CButton3_5 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_5, command = self.updatePlot)
        self.CButton3_5.place(relx=0.9, rely=0.63,height=10, relwidth=0.1)
        self.CButton3_5.configure(highlightcolor="black")
        self.CButton3_5.select();
        self.CButton3_5.bind("<Return>", self.updatePlot)
        
        self.CButtonvar3_6 = tk.IntVar()
        self.CButton3_6 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_6, command = self.updatePlot)
        self.CButton3_6.place(relx=0.9, rely=0.72,height=10, relwidth=0.1)
        self.CButton3_6.configure(highlightcolor="black")
        self.CButton3_6.select();
        self.CButton3_6.bind("<Return>", self.updatePlot)
        
        self.CButtonvar3_7 = tk.IntVar()
        self.CButton3_7 = tk.Checkbutton(self.Frame3, variable = self.CButtonvar3_7, command = self.updatePlot)
        self.CButton3_7.place(relx=0.9, rely=0.82,height=10, relwidth=0.1)
        self.CButton3_7.configure(highlightcolor="black")
        self.CButton3_7.select();
        self.CButton3_7.bind("<Return>", self.updatePlot)
        
        ##########################################################
        ## Settings ##############################################
        
        self.Frame4 = tk.Frame(root)
        self.Frame4.place(relx=0.012, rely=0.326, relheight=0.13, relwidth=0.133)
        self.Frame4.configure(relief='groove')
        self.Frame4.configure(borderwidth="2")
        self.Frame4.configure(relief='groove')
        self.Frame4.configure(highlightcolor="black")
        self.Frame4.configure(width=805)
        
        self.TLabel41 = ttk.Label(root)
        self.TLabel41.place(relx=0.06, rely=0.33, height=24, width=110)
        self.TLabel41.configure(relief='flat')
        self.TLabel41.configure(text='''SETTINGS''')
        
        self.TLabel42 = ttk.Label(root)
        self.TLabel42.place(relx=0.02, rely=0.355, height=24, width=110)
        self.TLabel42.configure(relief='flat')
        self.TLabel42.configure(text='''Model''')
        
        self.TLabel43 = ttk.Label(root)
        self.TLabel43.place(relx=0.04, rely=0.38, height=24, width=110)
        self.TLabel43.configure(relief='flat')
        self.TLabel43.configure(text='''Deterministic SIRD''')
        
        self.TLabel44 = ttk.Label(root)
        self.TLabel44.place(relx=0.04, rely=0.405, height=24, width=110)
        self.TLabel44.configure(relief='flat')
        self.TLabel44.configure(text='''Stochastic SIRD''')

        self.varModelOption = tk.StringVar();
        self.RButton4_1 = tk.Radiobutton(self.Frame4, variable = self.varModelOption, value = '1', command = self.updatePlot)
        self.RButton4_1.place(relx=0.9, rely=0.50,height=10, relwidth=0.1)
        self.RButton4_1.configure(highlightcolor="black")
        self.RButton4_1.select();
        self.RButton4_1.bind("<Return>", self.updatePlot)
        
        self.RButton4_2 = tk.Radiobutton(self.Frame4, variable = self.varModelOption, value = '2', command = self.updatePlot)
        self.RButton4_2.place(relx=0.9, rely=0.70,height=10, relwidth=0.1)
        self.RButton4_2.configure(highlightcolor="black")
        self.RButton4_2.select();
        self.RButton4_2.bind("<Return>", self.updatePlot)
        
        ##########################################################
        ## Rt ####################################################
        
        self.Frame5 = tk.Frame(root)
        self.Frame5.place(relx=0.012, rely=0.465, relheight=0.385, relwidth=0.133)
        self.Frame5.configure(relief='groove')
        self.Frame5.configure(borderwidth="2")
        self.Frame5.configure(relief='groove')
        self.Frame5.configure(highlightcolor="black")
        self.Frame5.configure(width=805)
        
        self.TLabel51 = ttk.Label(root)
        self.TLabel51.place(relx=0.025, rely=0.47, height=24, width=160)
        self.TLabel51.configure(relief='flat')
        self.TLabel51.configure(text='''Reproduction Number, R(t)''')
        
        # self.TLabel52 = ttk.Label(root)
        # self.TLabel52.place(relx=0.124, rely=0.605, height=24, width=10)
        # self.TLabel52.configure(relief='flat')
        # self.TLabel52.configure(text='''t''')
        
        self.TLabel53 = ttk.Label(root)
        self.TLabel53.place(relx=0.02, rely=0.50, height=24, width=160)
        self.TLabel53.configure(relief='flat')
        self.TLabel53.configure(text='''Without NPI''')
        
        self.TLabel53 = ttk.Label(root)
        self.TLabel53.place(relx=0.04, rely=0.525, height=24, width=100)
        self.TLabel53.configure(relief='flat')
        self.TLabel53.configure(text='''Average R(t)''')
        
        withoutNPIvalue = midpoint_beta[beta_peaks][2]/(mean_r+mean_d);
        self.TLabel54 = ttk.Label(root)
        self.TLabel54.place(relx=0.04, rely=0.55, height=24, width=100)
        self.TLabel54.configure(relief='flat')
        self.TLabel54.configure(text=withoutNPIvalue)
        
        self.TLabel53 = ttk.Label(root)
        self.TLabel53.place(relx=0.02, rely=0.575, height=24, width=160)
        self.TLabel53.configure(relief='flat')
        self.TLabel53.configure(text='''Loose NPI''')
        
        self.TLabel53 = ttk.Label(root)
        self.TLabel53.place(relx=0.04, rely=0.60, height=24, width=100)
        self.TLabel53.configure(relief='flat')
        self.TLabel53.configure(text='''Average R(t)''')
        
        averagelooseNPIvalue = midpoint_beta[beta_peaks][2]/(mean_r+mean_d);
        self.TLabel54 = ttk.Label(root)
        self.TLabel54.place(relx=0.04, rely=0.625, height=24, width=100)
        self.TLabel54.configure(relief='flat')
        self.TLabel54.configure(text=averagelooseNPIvalue)
        
        # self.TLabel55 = ttk.Label(root)
        # self.TLabel55.place(relx=0.04, rely=0.65, height=24, width=100)
        # self.TLabel55.configure(relief='flat')
        # self.TLabel55.configure(text='''Average duration''')
        
        # averagelooseNPIduration = '''~ 80 days''';
        # self.TLabel56 = ttk.Label(root)
        # self.TLabel56.place(relx=0.04, rely=0.675, height=24, width=100)
        # self.TLabel56.configure(relief='flat')
        # self.TLabel56.configure(text=averagelooseNPIduration)
        
        self.TLabel57 = ttk.Label(root)
        self.TLabel57.place(relx=0.02, rely=0.70, height=24, width=160)
        self.TLabel57.configure(relief='flat')
        self.TLabel57.configure(text='''Strict NPI''')
        
        self.TLabel58 = ttk.Label(root)
        self.TLabel58.place(relx=0.04, rely=0.725, height=24, width=100)
        self.TLabel58.configure(relief='flat')
        self.TLabel58.configure(text='''Average R(t)''')
        
        averagestrictNPIvalue = midpoint_beta[beta_peaks][1]/(mean_r+mean_d);
        self.TLabel59 = ttk.Label(root)
        self.TLabel59.place(relx=0.04, rely=0.75, height=24, width=100)
        self.TLabel59.configure(relief='flat')
        self.TLabel59.configure(text=averagestrictNPIvalue)
        
        # self.TLabel510 = ttk.Label(root)
        # self.TLabel510.place(relx=0.04, rely=0.775, height=24, width=100)
        # self.TLabel510.configure(relief='flat')
        # self.TLabel510.configure(text='''Average duration''')
        
        # averagestrictNPIduration = '''~ 40 days''';
        # self.TLabel511 = ttk.Label(root)
        # self.TLabel511.place(relx=0.04, rely=0.80, height=24, width=100)
        # self.TLabel511.configure(relief='flat')
        # self.TLabel511.configure(text=averagestrictNPIduration)
        
        ##########################################################
        ## User input ############################################
        
        self.Frame6 = tk.Frame(root)
        self.Frame6.place(relx=0.863, rely=0.058, relheight=0.922, relwidth=0.133)
        self.Frame6.configure(relief='groove')
        self.Frame6.configure(borderwidth="2")
        self.Frame6.configure(relief='groove')
        self.Frame6.configure(highlightcolor="black")
        self.Frame6.configure(width=805)
        
        self.TLabel61 = ttk.Label(root)
        self.TLabel61.place(relx=0.905, rely=0.06, height=24, width=100)
        self.TLabel61.configure(relief='flat')
        self.TLabel61.configure(text='''USER INPUT''')
        
        self.TLabel61 = ttk.Label(root)
        self.TLabel61.place(relx=0.88, rely=0.1, height=24, width=150)
        self.TLabel61.configure(relief='flat')
        self.TLabel61.configure(text='''Start date''')
        
        self.TLabel62 = ttk.Label(root)
        self.TLabel62.place(relx=0.927, rely=0.1, height=24, width=50)
        self.TLabel62.configure(relief='flat')
        self.TLabel62.configure(text='''-''')
        
        self.TLabel62 = ttk.Label(root)
        self.TLabel62.place(relx=0.945, rely=0.1, height=24, width=70)
        self.TLabel62.configure(relief='flat')
        self.TLabel62.configure(text='''End date''')
        
        self.Entry6_11 = tk.Entry(root)
        self.Entry6_11.place(relx=0.875, rely=0.15,height=30, relwidth=0.045)
        self.Entry6_11.configure(highlightcolor="black")
        self.Entry6_11.configure(selectforeground="black")
        self.Entry6_11.insert(0,str(sdate2))
        self.Entry6_11.bind("<Return>", self.updatePlot)
        
        self.Entry6_12 = tk.Entry(root)
        self.Entry6_12.place(relx=0.94, rely=0.15,height=30, relwidth=0.045)
        self.Entry6_12.configure(highlightcolor="black")
        self.Entry6_12.configure(selectforeground="black")
        self.Entry6_12.insert(0,str(sdate2 + td(days=30)))
        self.Entry6_12.bind("<Return>", self.updatePlot)
        
        self.Entry6_13 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_13.place(relx=0.875, rely=0.2,height=30, relwidth=0.109);
        self.Entry6_13.current(1)
        self.Entry6_13.bind("<Return>", self.updatePlot)
        
        self.Entry6_21 = tk.Entry(root)
        self.Entry6_21.place(relx=0.875, rely=0.25,height=30, relwidth=0.045)
        self.Entry6_21.configure(highlightcolor="black")
        self.Entry6_21.configure(selectforeground="black")
        self.Entry6_21.insert(0,str('YYYY-MM-DD'))
        self.Entry6_21.bind("<Return>", self.updatePlot)
        
        self.Entry6_22 = tk.Entry(root)
        self.Entry6_22.place(relx=0.94, rely=0.25,height=30, relwidth=0.045)
        self.Entry6_22.configure(highlightcolor="black")
        self.Entry6_22.configure(selectforeground="black")
        self.Entry6_22.insert(0,str('YYYY-MM-DD'))
        self.Entry6_22.bind("<Return>", self.updatePlot)
        
        self.Entry6_23 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_23.place(relx=0.875, rely=0.3,height=30, relwidth=0.109);
        self.Entry6_23.current(0)
        self.Entry6_23.bind("<Return>", self.updatePlot)
        
        self.Entry6_31 = tk.Entry(root)
        self.Entry6_31.place(relx=0.875, rely=0.35,height=30, relwidth=0.045)
        self.Entry6_31.configure(highlightcolor="black")
        self.Entry6_31.configure(selectforeground="black")
        self.Entry6_31.insert(0,str('YYYY-MM-DD'))
        self.Entry6_31.bind("<Return>", self.updatePlot)
        
        self.Entry6_32 = tk.Entry(root)
        self.Entry6_32.place(relx=0.94, rely=0.35,height=30, relwidth=0.045)
        self.Entry6_32.configure(highlightcolor="black")
        self.Entry6_32.configure(selectforeground="black")
        self.Entry6_32.insert(0,str('YYYY-MM-DD'))
        self.Entry6_32.bind("<Return>", self.updatePlot)
        
        self.Entry6_33 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_33.place(relx=0.875, rely=0.4,height=30, relwidth=0.109);
        self.Entry6_33.current(0)
        self.Entry6_33.bind("<Return>", self.updatePlot)
        
        self.Entry6_41 = tk.Entry(root)
        self.Entry6_41.place(relx=0.875, rely=0.45,height=30, relwidth=0.045)
        self.Entry6_41.configure(highlightcolor="black")
        self.Entry6_41.configure(selectforeground="black")
        self.Entry6_41.insert(0,str('YYYY-MM-DD'))
        self.Entry6_41.bind("<Return>", self.updatePlot)
        
        self.Entry6_42 = tk.Entry(root)
        self.Entry6_42.place(relx=0.94, rely=0.45,height=30, relwidth=0.045)
        self.Entry6_42.configure(highlightcolor="black")
        self.Entry6_42.configure(selectforeground="black")
        self.Entry6_42.insert(0,str('YYYY-MM-DD'))
        self.Entry6_42.bind("<Return>", self.updatePlot)
        
        self.Entry6_43 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_43.place(relx=0.875, rely=0.5,height=30, relwidth=0.109);
        self.Entry6_43.current(0)
        self.Entry6_43.bind("<Return>", self.updatePlot)
        
        self.Entry6_51 = tk.Entry(root)
        self.Entry6_51.place(relx=0.875, rely=0.55,height=30, relwidth=0.045)
        self.Entry6_51.configure(highlightcolor="black")
        self.Entry6_51.configure(selectforeground="black")
        self.Entry6_51.insert(0,str('YYYY-MM-DD'))
        self.Entry6_51.bind("<Return>", self.updatePlot)
        
        self.Entry6_52 = tk.Entry(root)
        self.Entry6_52.place(relx=0.94, rely=0.55,height=30, relwidth=0.045)
        self.Entry6_52.configure(highlightcolor="black")
        self.Entry6_52.configure(selectforeground="black")
        self.Entry6_52.insert(0,str('YYYY-MM-DD'))
        self.Entry6_52.bind("<Return>", self.updatePlot)
        
        self.Entry6_53 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_53.place(relx=0.875, rely=0.6,height=30, relwidth=0.109);
        self.Entry6_53.current(0)
        self.Entry6_53.bind("<Return>", self.updatePlot)
        
        self.Entry6_61 = tk.Entry(root)
        self.Entry6_61.place(relx=0.875, rely=0.65,height=30, relwidth=0.045)
        self.Entry6_61.configure(highlightcolor="black")
        self.Entry6_61.configure(selectforeground="black")
        self.Entry6_61.insert(0,str('YYYY-MM-DD'))
        self.Entry6_61.bind("<Return>", self.updatePlot)
        
        self.Entry6_62 = tk.Entry(root)
        self.Entry6_62.place(relx=0.94, rely=0.65,height=30, relwidth=0.045)
        self.Entry6_62.configure(highlightcolor="black")
        self.Entry6_62.configure(selectforeground="black")
        self.Entry6_62.insert(0,str('YYYY-MM-DD'))
        self.Entry6_62.bind("<Return>", self.updatePlot)
        
        self.Entry6_63 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_63.place(relx=0.875, rely=0.7,height=30, relwidth=0.109);
        self.Entry6_63.current(0)
        self.Entry6_63.bind("<Return>", self.updatePlot)
        
        self.Entry6_71 = tk.Entry(root)
        self.Entry6_71.place(relx=0.875, rely=0.75,height=30, relwidth=0.045)
        self.Entry6_71.configure(highlightcolor="black")
        self.Entry6_71.configure(selectforeground="black")
        self.Entry6_71.insert(0,str('YYYY-MM-DD'))
        self.Entry6_71.bind("<Return>", self.updatePlot)
        
        self.Entry6_72 = tk.Entry(root)
        self.Entry6_72.place(relx=0.94, rely=0.75,height=30, relwidth=0.045)
        self.Entry6_72.configure(highlightcolor="black")
        self.Entry6_72.configure(selectforeground="black")
        self.Entry6_72.insert(0,str('YYYY-MM-DD'))
        self.Entry6_72.bind("<Return>", self.updatePlot)
        
        self.Entry6_73 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_73.place(relx=0.875, rely=0.8,height=30, relwidth=0.109);
        self.Entry6_73.current(0)
        self.Entry6_73.bind("<Return>", self.updatePlot)
        
        self.Entry6_81 = tk.Entry(root)
        self.Entry6_81.place(relx=0.875, rely=0.85,height=30, relwidth=0.045)
        self.Entry6_81.configure(highlightcolor="black")
        self.Entry6_81.configure(selectforeground="black")
        self.Entry6_81.insert(0,str('YYYY-MM-DD'))
        self.Entry6_81.bind("<Return>", self.updatePlot)
        
        self.Entry6_82 = tk.Entry(root)
        self.Entry6_82.place(relx=0.94, rely=0.85,height=30, relwidth=0.045)
        self.Entry6_82.configure(highlightcolor="black")
        self.Entry6_82.configure(selectforeground="black")
        self.Entry6_82.insert(0,str('YYYY-MM-DD'))
        self.Entry6_82.bind("<Return>", self.updatePlot)
        
        self.Entry6_83 = ttk.Combobox(master, values = ['Without NPI','Loose NPI','Strict NPI']);
        self.Entry6_83.place(relx=0.875, rely=0.9,height=30, relwidth=0.109);
        self.Entry6_83.current(0)
        self.Entry6_83.bind("<Return>", self.updatePlot)
        
        ##########################################################
        ## create plot with default values from Entry widget
        
        fig = Figure()
        self.subplot1 = fig.add_subplot(111)
        line = self.subplot1.plot(pastplot_t,Iave_srk4)
        self.canvas = FigureCanvasTkAgg(fig, master=self.Frame1)
        mplcursors.cursor(line)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        
        self.Frame7 = tk.Frame(root)
        self.Frame7.place(relx=0.15, rely=0.94, relheight=0.04, relwidth=0.71)
        self.Frame7.configure(relief='groove')
        self.Frame7.configure(borderwidth="2")
        self.Frame7.configure(relief='groove')
        self.Frame7.configure(highlightcolor="black")
        self.Frame7.configure(width=805)
        fig2 = Figure()
        self.subplot2 = fig2.add_subplot(111)
        line2 = self.subplot2.plot(pastplot_t,Idiff)
        self.canvas2 = FigureCanvasTkAgg(fig2, master=self.Frame2)
        mplcursors.cursor(line2)
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(side='top', fill='both', expand=1)
                
        toolbar = NavigationToolbar2Tk(self.canvas, self.Frame1)
        toolbar2 = NavigationToolbar2Tk(self.canvas2, self.Frame7)
        
        def on_key_press(event):
            print("you pressed {}".format(event.key))
            key_press_handler(event, self.canvas, toolbar)
            key_press_handler(event, self.canvas2, toolbar2)
            
        self.canvas.mpl_connect("key_press_event", on_key_press)
        self.canvas2.mpl_connect("key_press_event", on_key_press)
        
        def _quit():
            root.quit()
            root.destroy()  
       
        self.updatePlot()
        
    #=====================================================           
    def updatePlot(self, *args):
        
        display_entry1 = int(self.CButtonvar3_1.get());
        display_entry2 = int(self.CButtonvar3_2.get());
        display_entry3 = int(self.CButtonvar3_3.get());
        display_entry4 = int(self.CButtonvar3_4.get());
        display_entry5 = int(self.CButtonvar3_5.get());
        display_entry6 = int(self.CButtonvar3_6.get());
        display_entry7 = int(self.CButtonvar3_7.get());
        
        settings_entry1 = int(float(self.varModelOption.get()));
        
        userinput_entry11 = str(self.Entry6_11.get());
        userinput_entry12 = str(self.Entry6_12.get());
        userinput_entry13 = str(self.Entry6_13.get());
        userinput_entry21 = str(self.Entry6_21.get());
        userinput_entry22 = str(self.Entry6_22.get());
        userinput_entry23 = str(self.Entry6_23.get());
        userinput_entry31 = str(self.Entry6_31.get());
        userinput_entry32 = str(self.Entry6_32.get());
        userinput_entry33 = str(self.Entry6_33.get());
        userinput_entry41 = str(self.Entry6_41.get());
        userinput_entry42 = str(self.Entry6_42.get());
        userinput_entry43 = str(self.Entry6_43.get());
        userinput_entry51 = str(self.Entry6_51.get());
        userinput_entry52 = str(self.Entry6_52.get());
        userinput_entry53 = str(self.Entry6_53.get());
        userinput_entry61 = str(self.Entry6_61.get());
        userinput_entry62 = str(self.Entry6_62.get());
        userinput_entry63 = str(self.Entry6_63.get());
        userinput_entry71 = str(self.Entry6_71.get());
        userinput_entry72 = str(self.Entry6_72.get());
        userinput_entry73 = str(self.Entry6_73.get());
        userinput_entry81 = str(self.Entry6_81.get());
        userinput_entry82 = str(self.Entry6_82.get());
        userinput_entry83 = str(self.Entry6_83.get());
        
        guiinput2 = [display_entry1,display_entry2,display_entry3,display_entry4,display_entry5,display_entry6,display_entry7,
        settings_entry1,
        userinput_entry11,userinput_entry12,userinput_entry13,
        userinput_entry21,userinput_entry22,userinput_entry23,
        userinput_entry31,userinput_entry32,userinput_entry33,
        userinput_entry41,userinput_entry42,userinput_entry43,
        userinput_entry51,userinput_entry52,userinput_entry53,
        userinput_entry61,userinput_entry62,userinput_entry63,
        userinput_entry71,userinput_entry72,userinput_entry73,
        userinput_entry81,userinput_entry82,userinput_entry83]
        
        Ssrk4_2, Isrk4_2, Rsrk4_2, Dsrk4_2 = SIRD_func(t0,dt,N,M,Szero,Izero,Rzero,Dzero,Ptotal,guiinput2)
        
        # Save_srk4 = Ssrk4_2.mean(axis=1);
        # Iave_srk4 = Isrk4_2.mean(axis=1);
        # Rave_srk4 = Rsrk4_2.mean(axis=1);
        # Dave_srk4 = Dsrk4_2.mean(axis=1);
        
        percentileSsrk4 = np.ones((N+1,3));
        percentileIsrk4 = np.ones((N+1,3));
        percentileRsrk4 = np.ones((N+1,3));
        percentileDsrk4 = np.ones((N+1,3));
        
        lowper = 25;
        mean = 50;
        highper = 75;
            
        for j in np.linspace(t0,t0+N,N+1):
            percentileSsrk4[j.astype(int)] = np.percentile(np.transpose(Ssrk4_2[j.astype(int),:]),[lowper,mean,highper]);
            percentileIsrk4[j.astype(int)] = np.percentile(np.transpose(Isrk4_2[j.astype(int),:]),[lowper,mean,highper]);
            percentileRsrk4[j.astype(int)] = np.percentile(np.transpose(Rsrk4_2[j.astype(int),:]),[lowper,mean,highper]);
            percentileDsrk4[j.astype(int)] = np.percentile(np.transpose(Dsrk4_2[j.astype(int),:]),[lowper,mean,highper]);

        Save_srk4 = percentileSsrk4[:,1];
        Iave_srk4 = percentileIsrk4[:,1];
        Rave_srk4 = percentileRsrk4[:,1];
        Dave_srk4 = percentileDsrk4[:,1];

        Sfut_2, Ifut_2, Rfut_2, Dfut_2 = SIRD_func(T,dt,N,M,dt_fulldata['susceptible'].iloc[-1],dt_fulldata['cases_active'].iloc[-1],dt_fulldata['cases_recovered'].iloc[-1],dt_fulldata['deaths_total'].iloc[-1],Ptotal,guiinput2)

        Sfut = Sfut_2.mean(axis=1);
        Ifut = Ifut_2.mean(axis=1);
        Rfut = Rfut_2.mean(axis=1);
        Dfut = Dfut_2.mean(axis=1);
        
        percentileSfut = np.ones((N+1,2));
        percentileIfut = np.ones((N+1,2));
        percentileRfut = np.ones((N+1,2));
        percentileDfut = np.ones((N+1,2));
            
        for j in np.linspace(t0,t0+N,N+1):
            percentileSfut[j.astype(int)] = np.percentile(np.transpose(Sfut_2[j.astype(int),:]),[lowper,highper]);
            percentileIfut[j.astype(int)] = np.percentile(np.transpose(Ifut_2[j.astype(int),:]),[lowper,highper]);
            percentileRfut[j.astype(int)] = np.percentile(np.transpose(Rfut_2[j.astype(int),:]),[lowper,highper]);
            percentileDfut[j.astype(int)] = np.percentile(np.transpose(Dfut_2[j.astype(int),:]),[lowper,highper]);
        
        futpltxend = max(datetoint(userinput_entry12),datetoint(userinput_entry22),datetoint(userinput_entry32),datetoint(userinput_entry42),datetoint(userinput_entry52),datetoint(userinput_entry62),datetoint(userinput_entry72),datetoint(userinput_entry82))
        futpltxlim = futpltxend-T;
        
        dt_futdata = dt_fulldata[['cases_recovered']].copy();
        dt_futdata['S_fut'] = Sfut[1:]
        dt_futdata['I_fut'] = Ifut[1:]
        dt_futdata['R_fut'] = Rfut[1:]
        dt_futdata['D_fut'] = Dfut[1:]
        dt_futdata = dt_futdata[['S_fut','I_fut','R_fut','D_fut']]
        dt_futdata = dt_futdata.shift(T, freq='D')
        
        self.subplot1.cla()
        if display_entry1 == 1:
            self.subplot1.plot(pastplot_t,Save_srk4,'slategray')
            self.subplot1.plot(futplot_t[:futpltxlim],Sfut[:futpltxlim],'darkgray')
            self.subplot1.plot(pastplot_t[:T],dt_fulldata['susceptible'],color = 'slategray',linestyle = 'dotted')
            self.subplot1.fill_between(pastplot_t,percentileSsrk4[:,0],percentileSsrk4[:,2],facecolor='slategray' ,alpha=0.2)
            self.subplot1.fill_between(futplot_t[:futpltxlim],percentileSfut[:futpltxlim,0],percentileSfut[:futpltxlim,1],facecolor='slategray' ,alpha=0.2)
            
        if display_entry3 == 1:
            self.subplot1.plot(pastplot_t,Iave_srk4,'darkviolet')
            self.subplot1.plot(futplot_t[:futpltxlim],Ifut[:futpltxlim],'purple')
            self.subplot1.plot(pastplot_t[:T],dt_fulldata['cases_active'],color = 'darkviolet',linestyle = 'dotted')
            self.subplot1.fill_between(pastplot_t,percentileIsrk4[:,0],percentileIsrk4[:,2],facecolor='darkviolet' ,alpha=0.2)
            self.subplot1.fill_between(futplot_t[:futpltxlim],percentileIfut[:futpltxlim,0],percentileIfut[:futpltxlim,1],facecolor='purple' ,alpha=0.2)
            
        if display_entry6 == 1:
            self.subplot1.plot(pastplot_t,Rave_srk4,'lime')
            self.subplot1.plot(futplot_t[:futpltxlim],Rfut[:futpltxlim],'green')
            self.subplot1.plot(pastplot_t[:T],dt_fulldata['cases_recovered'],color = 'lime',linestyle = 'dotted')
            self.subplot1.fill_between(pastplot_t,percentileRsrk4[:,0],percentileRsrk4[:,2],facecolor='lime' ,alpha=0.2)
            self.subplot1.fill_between(futplot_t[:futpltxlim],percentileRfut[:futpltxlim,0],percentileRfut[:futpltxlim,1],facecolor='green' ,alpha=0.2)
            
        if display_entry7 == 1:
            self.subplot1.plot(pastplot_t,Dave_srk4,'salmon')
            self.subplot1.plot(futplot_t[:futpltxlim],Dfut[:futpltxlim],'crimson')
            self.subplot1.plot(pastplot_t[:T],dt_fulldata['deaths_total'],color = 'salmon',linestyle = 'dotted')
            self.subplot1.fill_between(pastplot_t,percentileDsrk4[:,0],percentileDsrk4[:,2],facecolor='salmon' ,alpha=0.2)
            self.subplot1.fill_between(futplot_t[:futpltxlim],percentileDfut[:futpltxlim,0],percentileDfut[:futpltxlim,1],facecolor='crimson' ,alpha=0.2)
                    
        lines = [Line2D([0], [0], color='slategray', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='darkgray', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='slategray', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='blue', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='royalblue', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='blue', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='darkviolet', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='purple', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='darkviolet', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='fuchsia', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='darkmagenta', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='fuchsia', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='aquamarine', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='turquoise', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='aquamarine', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='lime', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='green', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='lime', linewidth=3, linestyle='dotted'),
                  Line2D([0], [0], color='salmon', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='crimson', linewidth=3, linestyle='-'),
                  Line2D([0], [0], color='salmon', linewidth=3, linestyle='dotted')]
        
        labels = ['Susceptible case','Susceptible approximation','Susceptible data',
                  'Total Infected case','Total Infected approximation','Total Infected data',
                  'Quarantine case','Quarantine approximation','Quarantine data',
                  'ICU case','ICU approximation','ICU data',
                  'Ventilators case','Ventilators approximation','Ventilators data',
                  'Recovered case','Recovered approximation','Recovered data',
                  'Cummulative Death case','Cummulative Death approximation','Cummulative Death data']
        
        if display_entry1 == 0:
            lines[0:3] = [0,0,0];
            labels[0:3] = [0,0,0];
        if display_entry2 == 0:
            lines[3:6] = [0,0,0];
            labels[3:6] = [0,0,0];
        if display_entry3 == 0:
            lines[6:9] = [0,0,0];
            labels[6:9] = [0,0,0];
        if display_entry4 == 0:
            lines[9:12] = [0,0,0];
            labels[9:12] = [0,0,0];
        if display_entry5 == 0:
            lines[12:15] = [0,0,0];
            labels[12:15] = [0,0,0];
        if display_entry6 == 0:
            lines[15:18] = [0,0,0];
            labels[15:18] = [0,0,0];
        if display_entry7 == 0:
            lines[18:21] = [0,0,0];
            labels[18:21] = [0,0,0];
        
        lines = [line for line in lines if line != 0];
        labels = [label for label in labels if label != 0];
        
        self.subplot1.legend(lines,labels,loc="upper left")
        self.subplot1.grid(True)
        self.subplot1.set_xlabel('time (day)')
        self.subplot1.set_ylabel('Case')
        self.canvas.draw()
        
        ##########################################################
        Idatalist = dt_fulldata['cases_active'].values.tolist() + dt_futdata['I_fut'].values.tolist();
        Idiff = (np.diff(Idatalist)/Idatalist[1:])+1
        from scipy.signal import savgol_filter
        Idiff = savgol_filter(Idiff, 21, 3) # window size 51, polynomial order 3
        self.subplot2.cla()
        self.subplot2.plot(fullplot_t[:T+futpltxlim],Idiff[:T+futpltxlim],'b')
        self.subplot2.legend(['Real-time Reproductive Number, R(t)'],loc="upper left")
        self.subplot2.grid(True)
        self.subplot2.autoscale(enable=True)
        self.canvas2.draw()
        
if __name__ == '__main__':
    root = tk.Tk()
    app = App(root)
    init(root, app)
    root.mainloop()



