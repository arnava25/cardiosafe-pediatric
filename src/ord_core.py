"""
Canonical O'Hara-Rudy 2011 endocardial ventricular myocyte.
Faithful INaK (4-state), two-compartment INaCa, SR Ca cycling with buffering.
Drug hooks: GKr_mult, GCaL_mult, GNa_mult, GKs_mult (1.0 = no effect).
Units: mV, ms, uA/uF, mM.
"""
import numpy as np

R=8314.0; T=310.0; F=96485.0
nao=140.0; cao=1.8; ko=5.4
L=0.01; rad=0.0011
vcell=1000*np.pi*rad**2*L
Ageo=2*np.pi*rad**2+2*np.pi*rad*L
Acap=2*Ageo
vmyo=0.68*vcell; vnsr=0.0552*vcell; vjsr=0.0048*vcell; vss=0.02*vcell

# published endo steady-state-ish initial conditions
Y0=np.array([
 -87.0, 7.0, 7.0, 145.0, 145.0, 1.0e-4, 1.0e-4, 1.2, 1.2, 0.0,   #0-9  V nai nass ki kss cai cass cansr cajsr m
 0.6, 0.6, 0.6, 0.6, 0.6,                                        #10-14 hf hs j hsp jp
 0.0, 0.5, 0.5,                                                  #15-17 mL hL hLp
 0.0, 1.0, 1.0, 0.0, 1.0, 1.0,                                   #18-23 a iF iS ap iFp iSp
 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,                    #24-32 d ff fs fcaf fcas jca nca ffp fcafp
 0.0, 0.0, 0.0, 0.0, 1.0,                                        #33-37 xrf xrs xs1 xs2 xk1
 0.0, 0.0, 0.0,                                                  #38-40 Jrelnp Jrelp CaMKt
], dtype=float)

def rhs(y,t,GKr_mult=1.0,GCaL_mult=1.0,GNa_mult=1.0,GKs_mult=1.0,
        stim_amp=-80.0,stim_dur=0.5,CL=1000.0):
    (V,nai,nass,ki,kss,cai,cass,cansr,cajsr,m,hf,hs,j,hsp,jp,
     mL,hL,hLp,a,iF,iS,ap,iFp,iSp,d,ff,fs,fcaf,fcas,jca,nca,ffp,fcafp,
     xrf,xrs,xs1,xs2,xk1,Jrelnp,Jrelp,CaMKt)=y
    cai=max(cai,1e-9); cass=max(cass,1e-9); cajsr=max(cajsr,1e-9); cansr=max(cansr,1e-9)
    nai=max(nai,1e-6); nass=max(nass,1e-6); ki=max(ki,1e-6); kss=max(kss,1e-6)
    vfrt=V*F/(R*T); vffrt=V*F*F/(R*T)

    ENa=(R*T/F)*np.log(nao/nai)
    EK=(R*T/F)*np.log(ko/ki)
    PKNa=0.01833
    EKs=(R*T/F)*np.log((ko+PKNa*nao)/(ki+PKNa*nai))

    # CaMK
    KmCaMK=0.15; aCaMK=0.05; bCaMK=0.00068; CaMKo=0.05; KmCaM=0.0015
    CaMKb=CaMKo*(1-CaMKt)/(1+KmCaM/cass); CaMKa=CaMKb+CaMKt
    dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt

    # INa
    mss=1/(1+np.exp(-(V+39.57)/9.871))
    tm=1/(6.765*np.exp((V+11.64)/34.77)+8.552*np.exp(-(V+77.42)/5.955))
    dm=(mss-m)/tm
    hss=1/(1+np.exp((V+82.90)/6.086))
    thf=1/(1.432e-5*np.exp(-(V+1.196)/6.285)+6.149*np.exp((V+0.5096)/20.27))
    ths=1/(0.009794*np.exp(-(V+17.95)/28.05)+0.3343*np.exp((V+5.730)/56.66))
    Ahf=0.99; Ahs=0.01; dhf=(hss-hf)/thf; dhs=(hss-hs)/ths
    h=Ahf*hf+Ahs*hs
    jss=hss; tj=2.038+1/(0.02136*np.exp(-(V+100.6)/8.281)+0.3052*np.exp((V+0.9941)/38.45))
    dj=(jss-j)/tj
    hssp=1/(1+np.exp((V+89.1)/6.086)); thsp=3*ths; dhsp=(hssp-hsp)/thsp; hp=Ahf*hf+Ahs*hsp
    tjp=1.46*tj; djp=(jss-jp)/tjp
    GNa=75.0*GNa_mult; fINap=1/(1+KmCaMK/CaMKa)
    INa=GNa*(V-ENa)*m**3*((1-fINap)*h*j+fINap*hp*jp)

    # INaL
    mLss=1/(1+np.exp(-(V+42.85)/5.264)); tmL=tm; dmL=(mLss-mL)/tmL
    hLss=1/(1+np.exp((V+87.61)/7.488)); thL=200.0; dhL=(hLss-hL)/thL
    hLssp=1/(1+np.exp((V+93.81)/7.488)); thLp=3*thL; dhLp=(hLssp-hLp)/thLp
    GNaL=0.0075*GNa_mult; fINaLp=1/(1+KmCaMK/CaMKa)
    INaL=GNaL*(V-ENa)*mL*((1-fINaLp)*hL+fINaLp*hLp)

    # Ito (endo)
    ass=1/(1+np.exp(-(V-14.34)/14.82))
    ta=1.0515/(1/(1.2089*(1+np.exp(-(V-18.4099)/29.3814)))+3.5/(1+np.exp((V+100)/29.3814)))
    da=(ass-a)/ta
    iss=1/(1+np.exp((V+43.94)/5.711)); delta_epi=1.0
    tiF=(4.562+1/(0.3933*np.exp(-(V+100)/100)+0.08004*np.exp((V+50)/16.59)))*delta_epi
    tiS=(23.62+1/(0.001416*np.exp(-(V+96.52)/59.05)+1.780e-8*np.exp((V+114.1)/8.079)))*delta_epi
    AiF=1/(1+np.exp((V-213.6)/151.2)); AiS=1-AiF
    diF=(iss-iF)/tiF; diS=(iss-iS)/tiS; i_=AiF*iF+AiS*iS
    assp=1/(1+np.exp(-(V-24.34)/14.82)); dap=(assp-ap)/ta
    dti_dev=1.354+1e-4/(np.exp((V-167.4)/15.89)+np.exp(-(V-12.23)/0.2154))
    dti_rec=1-0.5/(1+np.exp((V+70)/20))
    tiFp=dti_dev*dti_rec*tiF; tiSp=dti_dev*dti_rec*tiS
    diFp=(iss-iFp)/tiFp; diSp=(iss-iSp)/tiSp; ip=AiF*iFp+AiS*iSp
    Gto=0.02; fItop=1/(1+KmCaMK/CaMKa)
    Ito=Gto*(V-EK)*((1-fItop)*a*i_+fItop*ap*ip)

    # ICaL
    dss=1/(1+np.exp(-(V+3.940)/4.230))
    td=0.6+1/(np.exp(-0.05*(V+6))+np.exp(0.09*(V+14))); dd=(dss-d)/td
    fss=1/(1+np.exp((V+19.58)/3.696))
    tff=7+1/(0.0045*np.exp(-(V+20)/10)+0.0045*np.exp((V+20)/10))
    tfs=1000+1/(3.5e-5*np.exp(-(V+5)/4)+3.5e-5*np.exp((V+5)/6))
    Aff=0.6; Afs=0.4; dff=(fss-ff)/tff; dfs=(fss-fs)/tfs; f=Aff*ff+Afs*fs
    fcass=fss
    tfcaf=7+1/(0.04*np.exp(-(V-4)/7)+0.04*np.exp((V-4)/7))
    tfcas=100+1/(0.00012*np.exp(-V/3)+0.00012*np.exp(V/7))
    Afcaf=0.3+0.6/(1+np.exp((V-10)/10)); Afcas=1-Afcaf
    dfcaf=(fcass-fcaf)/tfcaf; dfcas=(fcass-fcas)/tfcas; fca=Afcaf*fcaf+Afcas*fcas
    tjca=75.0; djca=(fcass-jca)/tjca
    ktaup=2.5; tffp=ktaup*tff; dffp=(fss-ffp)/tffp; fp=Aff*ffp+Afs*fs
    tfcafp=ktaup*tfcaf; dfcafp=(fcass-fcafp)/tfcafp; fcap=Afcaf*fcafp+Afcas*fcas
    Kmn=0.002; k2n=1000; km2n=jca*1.0
    anca=1/(k2n/km2n+(1+Kmn/cass)**4); dnca=anca*k2n-nca*km2n
    PhiCaL=4*vffrt*(cass*np.exp(2*vfrt)-0.341*cao)/(np.exp(2*vfrt)-1) if abs(vfrt)>1e-7 else 4*F*(cass-0.341*cao)
    PhiCaNa=1.0*vffrt*(0.75*nass*np.exp(vfrt)-0.75*nao)/(np.exp(vfrt)-1) if abs(vfrt)>1e-7 else F*(0.75*nass-0.75*nao)
    PhiCaK=1.0*vffrt*(0.75*kss*np.exp(vfrt)-0.75*ko)/(np.exp(vfrt)-1) if abs(vfrt)>1e-7 else F*(0.75*kss-0.75*ko)
    PCa=0.0001*GCaL_mult; PCap=1.1*PCa
    PCaNa=0.00125*PCa; PCaK=3.574e-4*PCa; PCaNap=0.00125*PCap; PCaKp=3.574e-4*PCap
    fICaLp=1/(1+KmCaMK/CaMKa)
    ICaL=(1-fICaLp)*PCa*PhiCaL*d*(f*(1-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1-nca)+jca*fcap*nca)
    ICaNa=(1-fICaLp)*PCaNa*PhiCaNa*d*(f*(1-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1-nca)+jca*fcap*nca)
    ICaK=(1-fICaLp)*PCaK*PhiCaK*d*(f*(1-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1-nca)+jca*fcap*nca)

    # IKr
    xrss=1/(1+np.exp(-(V+8.337)/6.789))
    txrf=12.98+1/(0.3652*np.exp((V-31.66)/3.869)+4.123e-5*np.exp(-(V-47.78)/20.38))
    txrs=1.865+1/(0.06629*np.exp((V-34.70)/7.355)+1.128e-5*np.exp(-(V-29.74)/25.94))
    Axrf=1/(1+np.exp((V+54.81)/38.21))
    dxrf=(xrss-xrf)/txrf; dxrs=(xrss-xrs)/txrs; Xr=Axrf*xrf+(1-Axrf)*xrs
    rkr=1/(1+np.exp((V+55)/75))*1/(1+np.exp((V-10)/30))
    GKr=0.046*GKr_mult; IKr=GKr*np.sqrt(ko/5.4)*Xr*rkr*(V-EK)

    # IKs
    xs1ss=1/(1+np.exp(-(V+11.60)/8.932))
    txs1=817.3+1/(2.326e-4*np.exp((V+48.28)/17.80)+0.001292*np.exp(-(V+210)/230))
    dxs1=(xs1ss-xs1)/txs1; xs2ss=xs1ss
    txs2=1/(0.01*np.exp((V-50)/20)+0.0193*np.exp(-(V+66.54)/31)); dxs2=(xs2ss-xs2)/txs2
    KsCa=1+0.6/(1+(3.8e-5/cai)**1.4)
    GKs=0.0034*GKs_mult; IKs=GKs*KsCa*xs1*xs2*(V-EKs)

    # IK1
    xk1ss=1/(1+np.exp(-(V+2.5538*ko+144.59)/(1.5692*ko+3.8115)))
    txk1=122.2/(np.exp(-(V+127.2)/20.36)+np.exp((V+236.8)/69.33)); dxk1=(xk1ss-xk1)/txk1
    rk1=1/(1+np.exp((V+105.8-2.6*ko)/9.493))
    GK1=0.1908; IK1=GK1*np.sqrt(ko)*rk1*xk1*(V-EK)

    # INaCa (two compartment)
    kna1=15.0; kna2=5.0; kna3=88.12; kasymm=12.5
    wna=6e4; wca=6e4; wnaca=5e3; kcaon=1.5e6; kcaoff=5e3; qna=0.5224; qca=0.1670
    hca=np.exp(qca*vfrt); hna=np.exp(qna*vfrt)
    def ncx(na_c,ca_c,frac):
        h1=1+na_c/kna3*(1+hna); h2=(na_c*hna)/(kna3*h1); h3=1/h1
        h4=1+na_c/kna1*(1+na_c/kna2); h5=na_c*na_c/(h4*kna1*kna2); h6=1/h4
        h7=1+nao/kna3*(1+1/hna); h8=nao/(kna3*hna*h7); h9=1/h7
        h10=kasymm+1+nao/kna1*(1+nao/kna2); h11=nao*nao/(h10*kna1*kna2); h12=1/h10
        k1=h12*cao*kcaon; k2=kcaoff; k3p=h9*wca; k3pp=h8*wnaca; k3=k3p+k3pp
        k4p=h3*wca/hca; k4pp=h2*wnaca; k4=k4p+k4pp
        k5=kcaoff; k6=h6*ca_c*kcaon; k7=h5*h2*wna; k8=h8*h11*wna
        x1=k2*k4*(k7+k6)+k5*k7*(k2+k3); x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
        x3=k1*k3*(k7+k6)+k8*k6*(k2+k3); x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
        s=x1+x2+x3+x4; E1=x1/s; E2=x2/s; E3=x3/s; E4=x4/s
        KmCaAct=150e-6; allo=1/(1+(KmCaAct/ca_c)**2)
        JncxNa=3*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp; JncxCa=E2*k2-E1*k1
        Gncx=0.0008
        return frac*Gncx*allo*(1*JncxNa+2*JncxCa)
    INaCa_i=ncx(nai,cai,0.8); INaCa_ss=ncx(nass,cass,0.2)

    # INaK (4-state)
    k1p=949.5;k1m=182.4;k2p=687.2;k2m=39.4;k3p=1899.0;k3m=79300.0;k4p=639.0;k4m=40.0
    Knai0=9.073;Knao0=27.78;delta=-0.1550
    Knai=Knai0*np.exp(delta*vfrt/3); Knao=Knao0*np.exp((1-delta)*vfrt/3)
    Kki=0.5;Kko=0.3582;MgADP=0.05;MgATP=9.8;Kmgatp=1.698e-7;Hh=1e-7;eP=4.2;Khp=1.698e-7;Knap=224.0;Kxkur=292.0
    P=eP/(1+Hh/Khp+nai/Knap+ki/Kxkur)
    denom_i=(1+nai/Knai)**3+(1+ki/Kki)**2-1
    denom_o=(1+nao/Knao)**3+(1+ko/Kko)**2-1
    a1=(k1p*(nai/Knai)**3)/denom_i; b1=k1m*MgADP
    a2=k2p; b2=(k2m*(nao/Knao)**3)/denom_o
    a3=(k3p*(ko/Kko)**2)/denom_o; b3=(k3m*P*Hh)/(1+MgATP/Kmgatp)
    a4=(k4p*MgATP/Kmgatp)/(1+MgATP/Kmgatp); b4=(k4m*(ki/Kki)**2)/denom_i
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1
    s=x1+x2+x3+x4; E1=x1/s;E2=x2/s;E3=x3/s;E4=x4/s
    JnakNa=3*(E1*a3-E2*b3); JnakK=2*(E4*b1-E3*a1)
    Pnak=30.0; INaK=Pnak*(1*JnakNa+1*JnakK)

    # background
    xkb=1/(1+np.exp(-(V-14.48)/18.34)); IKb=0.003*xkb*(V-EK)
    PNab=3.75e-10
    INab=PNab*vffrt*(nai*np.exp(vfrt)-nao)/(np.exp(vfrt)-1) if abs(vfrt)>1e-7 else PNab*F*(nai-nao)
    PCab=2.5e-8
    ICab=PCab*4*vffrt*(cai*np.exp(2*vfrt)-0.341*cao)/(np.exp(2*vfrt)-1) if abs(vfrt)>1e-7 else PCab*4*F*(cai-0.341*cao)
    IpCa=0.0005*cai/(0.0005+cai)

    # SR / Ca fluxes
    JdiffNa=(nass-nai)/2; JdiffK=(kss-ki)/2; Jdiff=(cass-cai)/0.2
    bt=4.75; a_rel=0.5*bt
    Jrel_inf=a_rel*(-ICaL)/(1+(1.5/cajsr)**8)
    tau_rel=max(bt/(1+0.0123/cajsr),0.001)
    dJrelnp=(Jrel_inf-Jrelnp)/tau_rel
    btp=1.25*bt; a_relp=0.5*btp
    Jrel_infp=a_relp*(-ICaL)/(1+(1.5/cajsr)**8)
    tau_relp=max(btp/(1+0.0123/cajsr),0.001)
    dJrelp=(Jrel_infp-Jrelp)/tau_relp
    fJrelp=1/(1+KmCaMK/CaMKa); Jrel=(1-fJrelp)*Jrelnp+fJrelp*Jrelp
    Jupnp=0.004375*cai/(cai+0.00092)
    Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017)
    fJupp=1/(1+KmCaMK/CaMKa); Jleak=0.0039375*cansr/15
    Jup=(1-fJupp)*Jupnp+fJupp*Jupp-Jleak
    Jtr=(cansr-cajsr)/100

    # stimulus
    Istim=stim_amp if (t%CL)<stim_dur else 0.0

    # buffers
    cmdnmax=0.05; kmcmdn=0.00238; trpnmax=0.07; kmtrpn=0.0005
    BSRmax=0.047; KmBSR=0.00087; BSLmax=1.124; KmBSL=0.0087; csqnmax=10.0; kmcsqn=0.8
    Bcai=1/(1+cmdnmax*kmcmdn/(kmcmdn+cai)**2+trpnmax*kmtrpn/(kmtrpn+cai)**2)
    Bcass=1/(1+BSRmax*KmBSR/(KmBSR+cass)**2+BSLmax*KmBSL/(KmBSL+cass)**2)
    Bcajsr=1/(1+csqnmax*kmcsqn/(kmcsqn+cajsr)**2)

    dnai=-(INa+INaL+3*INaCa_i+3*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo
    dnass=-(ICaNa+3*INaCa_ss)*Acap/(F*vss)-JdiffNa
    dki=-(Ito+IKr+IKs+IK1+IKb+Istim-2*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo
    dkss=-(ICaK)*Acap/(F*vss)-JdiffK
    dcai=Bcai*(-(IpCa+ICab-2*INaCa_i)*Acap/(2*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo)
    dcass=Bcass*(-(ICaL-2*INaCa_ss)*Acap/(2*F*vss)+Jrel*vjsr/vss-Jdiff)
    dcansr=Jup-Jtr*vjsr/vnsr
    dcajsr=Bcajsr*(Jtr-Jrel)

    Itot=INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim
    dV=-Itot

    return [dV,dnai,dnass,dki,dkss,dcai,dcass,dcansr,dcajsr,dm,dhf,dhs,dj,dhsp,djp,
            dmL,dhL,dhLp,da,diF,diS,dap,diFp,diSp,dd,dff,dfs,dfcaf,dfcas,djca,dnca,dffp,dfcafp,
            dxrf,dxrs,dxs1,dxs2,dxk1,dJrelnp,dJrelp,dCaMKt]
