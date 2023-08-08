sigma=5.67051e-5
Lsun=3.828e33
Rsun=6.96e10
Msun=1.989e33
c=2.99792458e10
G=6.67e-8
yr=3600*24.*365.25

#ZAMS quantities

aaaa=0.39704170
bbbb=8.52762600
cccc=0.00025546
dddd=5.43288900
eeee=5.56357900
ffff=0.78866060
gggg=0.00586685

hhhh=1.71535900
iiii=6.59778800
jjjj=10.08855000
kkkk=1.012490166
llll=0.07490166
mmmm=0.01077422
nnnn=3.08223400
oooo=17.84778000
pppp=0.00022582

Lzams(M)=(aaaa*M**5.5+bbbb*M**11) \
	/ (cccc+M**3+dddd*M**5+eeee*M**7+ffff*M**8+gggg*M**9.5)

Rzams(M)=(hhhh*M**2.5+iiii*M**6.5+jjjj*M**11+kkkk*M**19+llll*M**19.5) \
	/ (mmmm+nnnn*M**2+oooo*M**8.5+M**18.5+pppp*M**19.5)

Tzams(M)=(Lzams(M)*Lsun/(4*pi*(Rzams(M)*Rsun)**2*sigma))**0.25

#TAMS quantities

ssss=log(0.02)
a17_1=0.1451
a17_2=0.1461+0.1237*(ssss+2)
a17_12=(a17_1<a17_2?a17_1:a17_2)

a17_3=0.097
a17_13=(a17_12>a17_3?a17_12:a17_3)

a17_4=0.097-0.1072*(ssss+3)
a17_14=(a17_13>a17_4?a17_13:a17_4)

a17=10**a17_14#max(0.097-0.1072*(ssss+3),max(0.097,min(0.1451,0.1461+0.1237*(ssss+2))))

a20=2.652091e1
a18=a20*2.187715e-1
a19=a20*1.46644
a21=1.472103
a22=3.071408
a23=2.61789e0
a24=1.075567e-2
a25=1.476246e0
a26=5.502535e0
c1=-8.672073e-2

Rtams(M)=(M<a17?(a18+a19*M**a21)/(a20+M**a22):(c1*M**3+a23*M**a26+a24*M**(a26+1.5))/(a25+M**5))

a13=7.859573e2
a14=3.858911e3
a15=2.888720e2
a16=7.196580
a11=1.031538*a14 
a12=1.043715*a14 

Ltams(M)=(a11*M**3+a12*M**4+a13*M**(a16+1.8))/(a14+a15*M**5+M**a16)

Ttams(M)=(Ltams(M)*Lsun/(4*pi*(Rtams(M)*Rsun)**2*sigma))**0.25

# MS time
a_1=1.593890e3 
a_2=2.706708e3
a_3=1.466143e2
a_4=4.141960e-2
a_5=3.426349e-1
a_6=1.949814e1
a_7=4.903830e0
a_8=5.212154e-2
a_9=1.312179e0
a10=8.073972e-1

t_bgb(M)=(a_1+a_2*M**4+a_3*M**5.5+M**7)/(a_4*M**2+a_5*M**7)

xxxx=0.95 # Only for Z=0.02
insidemumu1(M)=a_6/M**a_7
insidemumu2(M)=a_8+a_9/M**a10
insidemumu(M)=(insidemumu1(M)>insidemumu2(M)?insidemumu1(M):insidemumu2(M))
mumu(M)=(-0.5>-0.01*insidemumu(M)?0.5:1-0.01*insidemumu(M))
t_hook(M)=mumu(M)*t_bgb(M)

t_MS(M)=(mumu(M)>xxxx?t_hook(M):xxxx*t_bgb(M)) # in Myr


# Luminosity and Radius time evolution
tau_1(tt,M)=(1.<tt/t_hook(M)?1.:tt/t_hook(M))
epseps=0.01
tau_2_1(tt,M)=(tt-(1.-epseps)*t_hook(M))/(epseps*t_hook(M))
tau_2_2(tt,M)=(1.<tau_2_1(tt,M)?1:tau_2_1(tt,M))
tau_2(tt,M)=(0.>tau_2_2(tt,M)?0:tau_2_2(tt,M))

etaeta(M)=(M<1.?10:(M>1.1?20:(100*M-90)))

M_hook=1.0185 # only for Z=0.02

a33=1.4 # only for Z=0.02
a34=1.910302e-1
a35=3.931056e-1
a36=3.267776e-1
a37=5.990212e-1
a38=7.330122e-1
a39=1.172768e0
a40=3.982622e-1
a41=3.571038e0
a42=1.25 # only for Z=0.02
a43=6.300e-2
a44=1.200e0

delta_L_31(M)=a34/M**a35
delta_L_32(M)=a36/M**a37
delta_L_3(M)=(delta_L_31(M)<delta_L_32(M)?delta_L_31(M):delta_L_32(M))
delta_L_2(M)=delta_L_3(a33)*((M-M_hook)/(a33-M_hook))**0.4

delta_L(M)=(M<=M_hook?0.:(M<a33?delta_L_2(M):delta_L_3(M)))

delta_R_2(M)=a43*sqrt((M-M_hook)/(a42-M_hook))
delta_R_4(M)=(a38+a39*M**3.5)/(a40*M**3+M**a41)-1.
delta_R_3(M)=a43+(delta_R_4(2.)-a43)*((M-a42)/(2.-a42))**a44

delta_R(M)=(M<=M_hook?0.:(M<=a42?delta_R_2(M):(M<2.?delta_R_3(M):delta_R_4(M))))

a45=2.321400e-1 
a46=1.163659e-2
a47=1.048020e-2
a48=1.555590e0
a49=9.7700e-2
a50=2.4000e-1
a51=3.3000e-1
a52=1.1064e0
a53=1.1900e0
a54=3.855707e-1
a55=3.579064e-1
a56=9.587587e-1
a57=1.4 # only for Z=0.02

alpha_L_1(M)=(a45+a46*M**a48)/(M**0.4+a47*M**1.9)
alpha_L_2(M)=a49
alpha_L_3(M)=a49+5.*(0.3-a49)*(M-0.5)
alpha_L_4(M)=0.3+(a50-0.3)*(M-0.7)/(a52-0.7)
alpha_L_5(M)=a50+(a51-a50)*(M-a52)/(a53-a52)
alpha_L_6(M)=a51+(alpha_L_1(2.)-a51)*(M-a53)/(2.-a53)

alpha_L(M)=(M>=2.?alpha_L_1(M):\
           (M<0.5?alpha_L_2(M):\
           (M<0.7?alpha_L_3(M):\
	   (M<a52?alpha_L_4(M):\
	   (M<a53?alpha_L_5(M):alpha_L_6(M))))))

beta_L_1(M)=(0.>a54-a55*M**a56?0.:a54-a55*M**a56)
beta_L_2(M)=beta_L_1(a57)-10.*(M-a57)*beta_L_1(a57)
beta_L(M)=(M<=a57?beta_L_1(M):(beta_L_1(M)<=0.?beta_L_1(M):(beta_L_2(M)>0.?beta_L_2(M):0.)))

a58=4.907546e-1
a59=4.537070e0
a60=1.796220e0
a61=2.256216e0
a62=8.4300e-2 # only for Z=0.02 
a63=7.3600e-2 # only for Z=0.02 
a64=1.2100e-1 # only for Z=0.02 
a65=1.564231e-3
a66=0.8 # only for Z=0.02
a67=5.210157e0
a68=0.8 # only for Z=0.02
a69=1.071489e0
a70=7.108492e-1
a71=3.478514e0
a72=0.95e0 # only for Z=0.02
a73=3.969331e-3
a74=1.600e0
a75=1. # only for Z=0.02
a76=1.192334e-2
a77=-0.3868776e0
a78=7.615495e-1
a79=2. # only for Z=0.02
a80=0.0585542 # only for Z=0.02
a81=1.5 # only for Z=0.02

alpha_R_1(M)=a58*M**a60/(a59+M**a61)
a64=alpha_R_1(a66) # if a68>a66
alpha_R_2(M)=a62
alpha_R_3(M)=a62+(a63-a62)*(M-0.5)/0.15
alpha_R_4(M)=a63+(a64-a63)*(M-0.65)/(a68-0.65)
alpha_R_5(M)=a64+(alpha_R_1(a66)-a64)*(M-a68)/(a66-a68)
alpha_R_6(M)=alpha_R_1(a67)+a65*(M-a67)

alpha_R(M)=(M<0.5?alpha_R_2(M):(M<0.65?alpha_R_3(M):(M<a68?alpha_R_4(M):(M<a66?alpha_R_5(M):(M<=a67?alpha_R_1(M):alpha_R_6(M))))))

betap_R_1(M)=a69*M**3.5/(a70+M**a71)
betap_R_2(M)=1.06
betap_R_3(M)=1.06+(a72-1.06)*(M-1.)/(a74-1.06)
betap_R_4(M)=a72+(betap_R_1(2.)-a72)*(M-a74)/(2.-a74)
betap_R_5(M)=betap_R_1(16.)+a73*(M-16.)

betap_R(M)=(M<=1.?betap_R_2(M):(M<a74?betap_R_3(M):(M<2.?betap_R_4(M):(M<=16.?betap_R_1(M):betap_R_5(M)))))
beta_R(M)=betap_R(M)-1.

gamagama_1(M)=a76+a77*(M-a78)**a79
gamagama_2(M)=gamagama_1(1.)+(a80-gamagama_1(1.))*((M-1.)/(a75-1.))**a81
gamagama_3(M)=a80-10.*(M-a75)*a80

gamagama(M)=(M<=1.?gamagama_1(M):(M<=a75?gamagama_2(M):(M<a75+0.1?gamagama_3(M):0.)))


logLms_on_Lzams(tt,M)=alpha_L(M)*(tt/t_MS(M))+beta_L(M)*(tt/t_MS(M))**etaeta(M)+(log10(Ltams(M)/Lzams(M))-alpha_L(M)-beta_L(M))*(tt/t_MS(M))**2 - delta_L(M)*(tau_1(tt,M)**2-tau_2(tt,M)**2)

Lms(tt,M)=(tt<t_MS(M)?Lzams(M)*10**logLms_on_Lzams(tt,M):NaN)

logRms_on_Rzams(tt,M)=alpha_R(M)*(tt/t_MS(M))+beta_R(M)*(tt/t_MS(M))**10+gamagama(M)*(tt/t_MS(M))**40+(log10(Rtams(M)/Rzams(M))-alpha_R(M)-beta_R(M)-gamagama(M))*(tt/t_MS(M))**3 - delta_R(M)*(tau_1(tt,M)**3-tau_2(tt,M)**3)

Rms(tt,M)=(tt<t_MS(M)?Rzams(M)*10**logRms_on_Rzams(tt,M):NaN)

Tms(tt,M)=(Lms(tt,M)*Lsun/(4*pi*(Rms(tt,M)*Rsun)**2*sigma))**0.25

RL_on_a(q)=0.49*q**(2./3.)/(0.6*q**(2./3.)+log(1+q**(1./3.)))

#set xr [1:1000]
#set log x
#set samples 10000
#pl for [i=1:20] Tms(x,i)

#set xr [5:3.5]
#set yr [3:6]
#pl for [i=1:50:3] [t=0:500:0.01] '+' u (log10(Tms($1,i))):(log10(Lms($1,i))):(i) w l lw 2 lc pal not
