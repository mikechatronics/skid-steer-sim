"""Headless Python reference for skidSteerPro.m (same equations, no GUI).
Run:  python skid_steer_reference.py     (needs numpy)
Used to cross-check the MATLAB model; also handy for batch parameter sweeps."""
import numpy as np
g=9.81; Jw=0.01; Jm=6e-5; Kc=0.2; Ki=1.5; tau_s=0.05; rho_air=1.2
TER={ 'hard':dict(mx=0.80,my=0.70,K=0.004,crr=0.02,bek=None),
 'gravel':dict(mx=0.55,my=0.50,K=0.010,crr=0.04,bek=None),
 'sand':dict(mx=0.55,my=0.50,K=0.010,crr=0.01,bek=(1.1,0.99,1528.43)),
 'loam':dict(mx=0.50,my=0.45,K=0.025,crr=0.01,bek=(0.7,5.27,1515.04)),
 'snow':dict(mx=0.25,my=0.20,K=0.025,crr=0.01,bek=(1.6,4.37,196.72)) }
def P0(**kw):
    p=dict(v=0.8,w=0.8,L=0.6,b=0.5,r=0.1,bw=0.06,xcg=0,h=0.12,m=20,kz=1.0,
      mx=0.8,my=0.7,K=0.004,crr=0.02,bek=None,slope=0.0,psi=0.0,CdA=0.1,
      Vbat=24,Rbat=0.05,Rm=0.6,kt=0.05,Ng=10,eta=0.9,Imax=20,rho=0.5)
    p.update(kw); 
    if 'ter' in kw: p.update(TER[kw['ter']])
    return p
def sim(p,T=8,dt=5e-4,log=False):
    L,b,r=p['L'],p['b'],p['r']
    xcg=np.clip(p['xcg'],-L/2,L/2)
    xg=np.array([L/2,L/2,-L/2,-L/2]);yg=np.array([b/2,-b/2,b/2,-b/2])
    xr=xg-xcg;yr=yg
    m=p['m'];W=m*g*np.cos(p['slope'])
    Iz=p['kz']*m*(L**2+b**2)/12+m*xcg**2
    wl=(p['v']-p['w']*b/2)/r; wr=(p['v']+p['w']*b/2)/r
    omc=np.array([wl,wr,wl,wr])
    Jeff=Jw+Jm*p['Ng']**2
    a0=(L/2+xcg)/L
    N=W*np.array([a0,a0,1-a0,1-a0])/2
    X=Y=th=vx=vy=w=0.0
    om=np.zeros(4);zi=np.zeros(4);zx=np.zeros(4);zy=np.zeros(4)
    I=np.zeros(4);Ib=0.0;lastFx=np.zeros(4);lastFy=np.zeros(4)
    mu_x,mu_y,K=p['mx'],p['my'],p['K']
    vs=0.0
    for k in range(int(T/dt)):
        u=vx-w*yr; v=vy+w*xr
        sx=u-r*om; sy=v
        smag=np.sqrt(sx**2+sy**2)
        zx=zx+dt*(sx-smag*zx/K); zy=zy+dt*(sy-smag*zy/K)
        kk=mu_x*N/K; mi=np.maximum(N/g,0.2)
        c=2*0.7*np.sqrt(np.maximum(kk,0)*mi)
        fxn=-(zx/K)-c*sx/np.maximum(mu_x*N,1e-6)   # normalised by mu N
        fyn=-(zy/K)-c*sy/np.maximum(mu_y*N,1e-6)
        mag=np.hypot(fxn,fyn); sc=np.where(mag>1,1/np.maximum(mag,1e-9),1)
        Fx=mu_x*N*fxn*sc; Fy=mu_y*N*fyn*sc
        # rolling resistance
        Rr=p['crr']*N
        if p['bek'] is not None:
            n,kc,kp=p['bek']; kpp=kc*1000+p['bw']*kp*1000
            zs=(3*N/((3-n)*kpp*np.sqrt(2*r)))**(2/(2*n+1))
            Rr=Rr+kpp*zs**(n+1)/(n+1)
        Trr=r*Rr*np.tanh(om/0.3)
        # motor
        e=omc-om; zi=np.clip(zi+dt*e, -1/Ki, 1/Ki)
        duty=np.clip(Kc*e+Ki*zi,-1,1)
        Vt=np.clip(p['Vbat']-p['Rbat']*Ib,0,1.3*p['Vbat'])
        omm=p['Ng']*om
        I=np.clip((duty*Vt-p['kt']*omm)/p['Rm'],-p['Imax'],p['Imax'])
        taum=p['kt']*I
        tauw=np.where(taum*omm>=0,p['Ng']*p['eta']*taum,p['Ng']*taum/p['eta'])
        Ib=np.sum(duty*I)
        # body
        gx=-g*np.sin(p['slope'])*np.cos(th-p['psi']); gy=g*np.sin(p['slope'])*np.sin(th-p['psi'])
        V=np.hypot(vx,vy); Dx=-0.5*rho_air*p['CdA']*V*vx; Dy=-0.5*rho_air*p['CdA']*V*vy
        Msc=-(2/3)*0.5*(mu_x+mu_y)*N*(p['bw']/2)*np.tanh(w/0.1)
        ax=(Fx.sum()+Dx)/m+gx+w*vy
        ay=(Fy.sum()+Dy)/m+gy-w*vx
        aw=((xr*Fy-yr*Fx).sum()+Msc.sum())/Iz
        # load transfer target
        SFx=Fx.sum(); SFy=Fy.sum(); h=p['h']
        Nf=(-h*SFx+W*(L/2+xcg))/L; Nr=W-Nf
        T_=-h*SFy/b
        Nt=np.array([Nf/2+p['rho']*T_, Nf/2-p['rho']*T_, Nr/2+(1-p['rho'])*T_, Nr/2-(1-p['rho'])*T_])
        lift=Nt.min()<=0
        N=N+dt*(np.maximum(Nt,0)-N)/tau_s
        om=om+dt*(tauw-r*Fx-Trr)/Jeff
        c_,s_=np.cos(th),np.sin(th)
        X+=dt*(vx*c_-vy*s_);Y+=dt*(vx*s_+vy*c_);th+=dt*w
        vx+=dt*ax;vy+=dt*ay;w+=dt*aw
        lastFx,lastFy=Fx,Fy
    return dict(vx=vx,vy=vy,w=w,om=om,N=N,I=I,Ib=Ib,Vt=Vt,lift=lift,th=th,X=X,Y=Y,zx=zx)
def rep(name,p,**kw):
    s=sim(p,**kw);w=s['w']
    chi=np.nan
    if abs(w)>0.02:
        vL=p['r']*s['om'][[0,2]].mean();vR=p['r']*s['om'][[1,3]].mean()
        chi=((s['vx']-vL)/w-(s['vx']-vR)/w)/p['b']
    print(f"{name:28s} vx={s['vx']:+.3f} vy={s['vy']:+.3f} w={s['w']:+.3f} chi={chi:.2f} N={np.round(s['N'],1)} I={np.round(s['I'],1)} Ib={s['Ib']:.1f} Vt={s['Vt']:.1f} lift={s['lift']}")
if __name__=='__main__':
    rep('hard straight',P0(w=0))
    rep('hard turn',P0())
    rep('hard spin',P0(v=0,w=2))
    rep('sand straight',P0(w=0,ter='sand'))
    rep('sand turn',P0(ter='sand'))
    rep('snow turn',P0(ter='snow'))
    rep('standstill',P0(v=0,w=0))
    rep('slope 15deg hold',P0(v=0,w=0,slope=np.radians(15)),T=4)
    rep('slope 15deg climb',P0(w=0,slope=np.radians(15)))
    rep('high CG fast turn',P0(v=2,w=2,h=0.35,L=0.4,b=0.3),T=5)
    rep('weak battery',P0(Imax=4,v=0,w=2))
