from numpy import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker
import sys,getopt,csv

#[time]=ms
#[length]=cm
#Data acquisition (simulation) parameters
dt=.02 #timestep=50us
N=2000 #how many axons  
SWEEP_TIME=5 #ms
#Nerve physiological quantities
SLOW_V_MEAN     = 1.5 #cm/ms    
SLOW_V_STD      = .2 #cm/ms     2
SLOW_THRESHOLD_MEAN     =1.5#norm. potential 
SLOW_THRESHOLD_STD      =0.1 #norm. potential 
FAST_V_MEAN=3.5 #cm/s
FAST_V_STD      = .9 #cm/ms     
FAST_THRESHOLD_MEAN     =0.5 #norm. potential 
FAST_THRESHOLD_STD      =0.2 #norm. potential   
PROB_SLOW=0.3
TAU_MEAN=0.1
TAU_STD=0.05
LENGTH_STD=2
LENGTH_MEAN=5


#Electrodes distances
AB=1            #Distance between recording electrodes, cm
DIST=2          #Stimulus-to-recording dinstance, cm
STIM_AMPLITUDE=1.2#V
STIM_DURATION=.05 #ms
SHOW_AB_SEP=False

#OUTPUT
SAVEFLAG=False
SAVEPATH=''

def randnormpos(mu,sigma):
        x=-1
        while x<0:
                x=random.randn()*sigma+mu
        return x
        
class axon():
        def __init__(self,):
                self.tau=randnormpos(TAU_MEAN,TAU_STD)
                if random.rand()<PROB_SLOW:
                        self.v=randnormpos(SLOW_V_MEAN,SLOW_V_STD)
                        self.th=randnormpos(SLOW_THRESHOLD_MEAN,SLOW_THRESHOLD_STD)
                else:
                        self.v=randnormpos(FAST_V_MEAN,FAST_V_STD)
                        self.th=randnormpos(FAST_THRESHOLD_MEAN,FAST_THRESHOLD_STD)
                self.length=randnormpos(LENGTH_MEAN,LENGTH_STD)
                
        def apply_stimulus(self,stim_amp,stim_dur,distance,tspace):
                X=0*tspace
                #print stim_amp*(1-exp(-stim_dur/self.tau))
                if stim_amp*(1-exp(-stim_dur/self.tau))>self.th and distance<=self.length:
                        t=arange(0,stim_dur+dt,dt)
                        self.latency= [t1 for t1 in t if stim_amp*(1-exp(-t1/self.tau))>self.th][0]
                        t0=self.latency+distance/self.v
                        if t0<tspace[-1]:
                                i0=[i for i in range(len(tspace)) if tspace[i]-t0>=0][0]
                                newt=tspace[i0:]                        
                                X[i0:]=-exp(-(newt-t0)/self.tau)
                return X
                
def response(): 
        
        t=arange(-1,SWEEP_TIME,dt)
        A=[axon() for k in range(N)]


        d=DIST
        S_A=sum(matrix([a.apply_stimulus(STIM_AMPLITUDE,STIM_DURATION,DIST,t) for a in A])/N*50,axis=0).T+(random.rand(len(t),1)/8.0-1.0/16.0)
        S_B=sum(matrix([a.apply_stimulus(STIM_AMPLITUDE,STIM_DURATION,DIST+AB,t) for a in A])/N*50,axis=0).T+(random.rand(len(t),1)/8.0-1.0/16.0)
        S=S_A-S_B
        return t,S

def measure_response(t,X):
        xs=sorted(X)
        R=dict()
        R['Amplitude_B']=mean(xs[-5:])
        R['Amplitude_A']=-mean(xs[:5])
        nsmp=len(X)
        X_FW_I =array([ i for i in range(nsmp) if X[i] <= -R['Amplitude_A']*0.5])
        R['Latency']=mean(X_FW_I[:2])*dt-1
        R['FWHM']=mean(X_FW_I[-2:])-R['Latency']
        return R


def main(argv): 
        try:
                opts,args=getopt.getopt(argv,"x",['distance_AB=','distance_record=',
                                                                 'tau_mean=','tau_std=','length_mean=','length_std=',
                                                                'speed_mean=','speed_std=',
                                                                'speed2_mean=','speed2_std=',
                                                                'two_populations','stimulus_amplitude=','stimulus_duration=',
                                                                'show_AB','output='])
        except:
                print 'One of the arguments is wrong!'
                exit (1)
        for o,a in opts:
                 if o=='--output':
                         global SAVEFLAG
                         global SAVEPATH
                         SAVEFLAG=True
                         SAVEPATH=a
                         #print a
                 if o=='--speed_mean':
                         global FAST_V_MEAN
                         FAST_V_MEAN=float(a)
                 elif o=='--speed_std':
                         global FAST_V_STD
                         FAST_V_STD=float(a)
                 elif o=='--speed2_mean':
                         global SLOW_V_MEAN
                         SLOW_V_MEAN=float(a)
                 elif o=='--speed2_std':
                         global SLOW_V_STD
                         SLOW_V_STD=float(a)
                 elif o=='--stimulus_duration':
                         global STIM_DURATION
                         STIM_DURATION=float(a)
                 elif o=='--stimulus_amplitude':
                         global STIM_AMPLITUDE
                         STIM_AMPLITUDE=float(a)        
                 elif o=='--distance_record':
                         global DIST
                         DIST=float(a)
                 elif o=='--distance_AB':
                         global AB
                         AB=float(a)
                 elif o=='--tau_mean':
                         global TAU_MEAN
                         TAU_MEAN=float(a)
                 elif o=='--tau_std':
                         global TAU_STD
                         TAU_STD=float(a)
                 elif o=='--length_std':
                         global LENGTH_STD
                         LENGTH_STD=float(a)
                 elif o=='--length_mean':
                         global LENGTH_MEAN
                         LENGTH_MEAN=float(a)
                 elif o=='--two_populations':
                         global PROB_SLOW
                         PROB_SLOW=0.3
                         print PROB_SLOW
                 elif o=='--show_AB':
                         global SHOW_AB_SEP
                         SHOW_AB_SEP=True
        t=arange(-1,SWEEP_TIME,dt)
        A=[axon() for k in range(N)]
 
 
        d=DIST
        S_A=sum(matrix([a.apply_stimulus(STIM_AMPLITUDE,STIM_DURATION,DIST,t) for a in A])/N*50,axis=0).T+(random.rand(len(t),1)/8.0-1.0/16.0)
        S_B=sum(matrix([a.apply_stimulus(STIM_AMPLITUDE,STIM_DURATION,DIST+AB,t) for a in A])/N*50,axis=0).T+(random.rand(len(t),1)/8.0-1.0/16.0)
        S=S_A-S_B
 
 
        CLR_BGR=(.1,.1,.1)
        CLR_FRG=(.9,.9,.9)
        CLR_BGR=(1,1,1)
        CLR_FRG=(0.1,0.1,0.1)
        CLR_SIGNAL=(.2,0,.2)
        fig,ax=subplots()
        fig.set_facecolor(CLR_BGR)
        ax.set_axis_bgcolor(CLR_BGR)
        if SHOW_AB_SEP:
                 plot(t,S_A,'-',color=(.5,1,.5,.6),linewidth=2)
                 plot(t,-S_B,'-',color=(1,.5,.5,.6),linewidth=2)
         
        plot(t,S,color=CLR_SIGNAL,linewidth=3)
        xlabel('t(ms)')
        ylabel('V(mV)')
 
        ax.xaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
        ax.xaxis.set_minor_locator(plticker.MultipleLocator(base=0.25))
        ax.yaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
        ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=0.25))
        ax.tick_params(axis='x',colors=CLR_FRG)
        ax.xaxis.label.set_color(CLR_FRG)
        ax.tick_params(axis='y',colors=CLR_FRG)
        ax.yaxis.label.set_color(CLR_FRG)
 
        ax.grid (True,which='major',color=(0.9,0.9,0.9),linestyle='-')
        ax.grid (True,which='minor',color=(0.5,0.5,0.5),linestyle=':')
 
        if SAVEFLAG:
                #print shape(array(t)),shape(S),shape(S_A),shape(S_B)
                M=matrix(hstack([S,S_A,S_B]))
                tt=matrix(t).T
                M=hstack([tt,M])
                print shape(M),shape(tt)
                savetxt(SAVEPATH+'.csv',M,fmt='%.3f', delimiter=', ',header='time(ms), Signal(mV), SignalA(mV), SignalB(mV)')
                PTABLE=matrix([ STIM_AMPLITUDE,STIM_DURATION,DIST])
                print PTABLE
                savetxt(SAVEPATH+'_params.csv',PTABLE,delimiter=', ',fmt='%.3f',header='STIM_AMP(mA), STIM_DUR(ms), DISTANCE(cm)')
        show()
if __name__=="__main__":
         main(sys.argv[1:])

