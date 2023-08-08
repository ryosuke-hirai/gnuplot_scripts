arad=7.5646e-15
G  = 6.67428e-8
msun = 1.989e33
rsun = 6.963e10
lsun = 3.9e33
kbol = 1.38064852e-16
amu = 1.6605402e-24
m_p = 1.6726231e-24
m_n = 1.6749286e-24
m_e = 9.1093897e-28
N_A = 6.0221367e23
sigma = 5.67051e-5
year = 3600*24*365.25
h = 6.6260755e-27
au = 1.5e13
clight = 2.99792458e10

eion(e,rho,mu,T,v1,v2)=e/rho-0.5*(v1**2+v2**2)-1.5*kbol*T/(mu*amu)-arad*T**4/rho
entropy(rho,T,mu)=1/mu*log(T**1.5/rho)+4*arad/3*T**3/rho/kbol*amu
