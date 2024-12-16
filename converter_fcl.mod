#define INSTABILITY

// -----------------------------
//  CONVERTER CIRCUITS:
//
//  CONV_PCC
//  AVG_MODEL
//  AVG_CAP_LEG
//  COMPL_MODEL
//  AVG_CONV
//  COMPL_CONV
// ----------------------------


// -----------------------------
//  CONV_PCC
// ----------------------------

subckt CONV_PCC in_a in_b in_c a1 b1 c1 a3 b3 c3 ia ib ic on
parameters V1=1 V2=1 L=1

;TRANFORMER
T1  in_a gnd  a3 pqs transformer t1=V1 t2=V2
T2  in_b gnd  b3 pqs transformer t1=V1 t2=V2
T3  in_c gnd  c3 pqs transformer t1=V1 t2=V2

Rpqs pqs gnd resistor r=10k

Lta a2 a3 inductor l=L ic=0
Ltb b2 b3 inductor l=L ic=0
Ltc c2 c3 inductor l=L ic=0

Rta a1 a2 on  gnd  ONSW
Rtb b1 b2 on  gnd  ONSW
Rtc c1 c2 on  gnd  ONSW

; Current sensing
Ia0 ia gnd ccvs sensedev="Lta" gain1=1
Ib0 ib gnd ccvs sensedev="Ltb" gain1=1
Ic0 ic gnd ccvs sensedev="Ltc" gain1=1

ends

// -----------------------------
//  AVG_CAP_LEG
// ----------------------------

subckt AVG_CAP_LEG pos neg phase i_up i_dw drv_up drv_dw blocked
parameters L=1 C=1 R=1 LEVELS=1 VDC_NOM=1

GATE gtc gnd blocked gnd vcvs func=1-v(blocked,gnd)

D1up pos  pos2 DIO
Sup  pos  pos1  gtc gnd SW
D2up pos2 pos1 DIO
Vup  pos2 pos3 cap_up drv_up blocked vcvs func=v(cap_up)*((v(drv_up)*LEVELS)*(1-v(blocked))+LEVELS*v(blocked))
D3up pos3 pos  DIO_BLK
Rup  pos3 pos4 resistor r=R
Lup  pos4 phase inductor l=L

Rdw  phase neg4 resistor r=R
Ldw  neg4  neg3 inductor l=L
D1dw neg3  neg1 DIO
Sdw  neg3  neg2 gtc gnd SW
D2dw neg1  neg2 DIO
Vdw  neg1  neg cap_dw drv_dw blocked vcvs func=v(cap_dw)*((v(drv_dw)*LEVELS)*(1-v(blocked))+LEVELS*v(blocked))
D3dw neg   neg3 DIO_BLK

I_up  i_up gnd  ccvs sensedev="Lup" gain1=1
I_dw  i_dw gnd  ccvs sensedev="Ldw" gain1=1
I_vup i_vup gnd ccvs sensedev="Vup" gain1=1
I_vdw i_vdw gnd ccvs sensedev="Vdw" gain1=1

I_cup cap_up gnd i_vup drv_up blocked vccs func=-v(i_vup)*((v(drv_up)*LEVELS)*(1-v(blocked))+LEVELS*v(blocked))
Cup cap_up gnd capacitor c=C ic=VDC_NOM/LEVELS
I_cdw cap_dw gnd i_vdw drv_dw blocked vccs func=-v(i_vdw)*((v(drv_dw)*LEVELS)*(1-v(blocked))+LEVELS*v(blocked))
Cdw cap_dw gnd capacitor c=C ic=VDC_NOM/LEVELS

ends


// -----------------------------
//  AVG_CAP
// ----------------------------
subckt AVG_CAP_MODEL a1 b1 c1 dc_pos dc_neg a3 b3 c3 za zb zc vdp vqp idp iqp idp_ref iqp_ref vdc_mea idc_pos idc_neg omega blocked acon
parameters L=1 C=1 R=1 VDIG=1 LEVELS=1 VDC_NOM=1 TCROSS=1 CNTR_L=1 R_GND=1 W_NOM=1
parameter CONV_NAME="Unknown"

Leg1 dc_pos dc_neg a1 ia_up ia_lo drv_ua drv_la blocked AVG_CAP_LEG L=L C=C R=R LEVELS=LEVELS VDC_NOM=VDC_NOM
Leg2 dc_pos dc_neg b1 ib_up ib_lo drv_ub drv_lb blocked AVG_CAP_LEG L=L C=C R=R LEVELS=LEVELS VDC_NOM=VDC_NOM
Leg3 dc_pos dc_neg c1 ic_up ic_lo drv_uc drv_lc blocked AVG_CAP_LEG L=L C=C R=R LEVELS=LEVELS VDC_NOM=VDC_NOM

Idc_pos   idc_pos gnd ia_up ib_up ic_up vcvs func=v(ia_up)+v(ib_up)+v(ic_up)
Idc_neg   idc_neg gnd ia_lo ib_lo ic_lo vcvs func=v(ia_lo)+v(ib_lo)+v(ic_lo)

Rgnda a1  gnd   resistor r=R_GND
Rgndb b1  gnd   resistor r=R_GND
Rgndc c1  gnd   resistor r=R_GND

Conv a1 b1 c1 dc_pos dc_neg a3 b3 c3 za zb zc \
     vdp vqp idp iqp idp_ref iqp_ref \
     ia_up ib_up ic_up ia_lo ib_lo ic_lo \
     drv_ua drv_ub drv_uc drv_la drv_lb drv_lc \
     omega idc_pos vdc_mea blocked acon COMPL_CONV L=CNTR_L
ends


// -----------------------------
//  COMPL_CONV
// ----------------------------

// Generic converter controller and driver with circulating current controller
subckt COMPL_CONV a1 b1 c1 dc_pos dc_neg a3 b3 c3 za zb zc \
            vdp vqp idp iqp   \
            idp_ref iqp_ref   \
	    ia_up ib_up ic_up \
	    ia_lo ib_lo ic_lo \
	    drv_ua drv_ub drv_uc \
	    drv_la drv_lb drv_lc \
	    omega idc_pos vdc_mea blocked on

parameter L=1 
parameter L_PU=L/(VAC_CONV^2/(P_NOM*W_NOM))

; AC CURRENT MEASUREMENTS
Ia1 ia gnd za gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Ib1 ib gnd zb gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Ic1 ic gnd zc gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes

; AC VOLTAGE MEASUREMENTS
Va va gnd a3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0
Vb vb gnd b3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0
Vc vc gnd c3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0

; POWER MEASUREMENTS
Pwr_mod pow_mod gnd ia ib ic va vb vc vcvs func=v(ia)*v(va)+v(ib)*v(vb)+v(ic)*v(vc)
Pwr pow gnd pow_mod gnd svcvs numer=[WN^2] denom=[WN^2,2*XI*WN,1] ic=0 skipck=yes


; DC VOLTAGE MEASUREMENTS
VdcMea vdc_mea  gnd dc_pos  dc_neg  svcvs dcgain=1/VDC_NOM numer=[WN_1^2] \
                                    denom=[WN_1^2,2*XI_1*WN_1,1] ic=1 skipck=yes


; PARK TRANSFORM OF PCC VOLTAGE AND CURRENT (ABC TO DQ0)
Vdq0 vdp vqp va vb vc omega ABCDQ0 GAIN=(2/3)/(VAC_CONV) TETA_GAIN=1
Idq0 idp iqp ia ib ic omega ABCDQ0 GAIN=(2/3)/(I_NOM*sqrt(2)) TETA_GAIN=1

; PLL (PCC VOLTAGE IS USED)
Ipll pllin     gnd  vdp  vqp  vcvs  func=v(vqp)/max(0.1,(v(vdp)^2+v(vqp)^2))
Dw      dw     gnd  pllin  gnd  svcvs numer=[122.5k,490] denom=[0,1] \
                                      maxdcgain=MAXDCGAIN  ic=0 skipck=yes
Wnew  omega    gnd  dw       vcvs   func=v(dw) + W_NOM

; FILTERS
Vd_filt  vdp_filt gnd  vdp gnd svcvs dcgain=1 numer=[1] denom=[1, 1/(2*pi*11)]
Vq_filt  vqp_filt gnd  vqp gnd svcvs dcgain=1 numer=[1] denom=[1, 1/(2*pi*11)]

; INTERNAL CONTROL LOOP
Delta_ed    d_ed gnd idp_ref  idp svcvs numer=[KI_I,KP_I] denom=[0,1] \
                                      maxdcgain=MAXDCGAIN ic=0
Ud   udp gnd d_ed iqp omega vdp_filt vcvs func=v(vdp_filt) + v(d_ed) + v(iqp)*v(omega)/W_NOM*L_PU

Delta_eq    d_eq gnd iqp_ref  iqp svcvs numer=[KI_I,KP_I] denom=[0,1] \
                                      maxdcgain=MAXDCGAIN ic=0
Uq   uqp gnd d_eq idp omega vqp_filt vcvs func=v(vqp_filt) + v(d_eq) - v(idp)*v(omega)/W_NOM*L_PU

; LIMITER
Il   udp_drv  uqp_drv  udp uqp  vdc_mea  blocked U_LIMITER_S

; Inverse Park transform of converter reference voltage (DQ0 to DQ0)
Uabc uf_a uf_b uf_c udp_drv uqp_drv gnd omega DQ0ABC GAIN=VAC_CONV TETA_GAIN=1

; COMPUTATION OF CIRCULATING CURRENTS
Iza iza gnd ia_up ia_lo on vcvs func=0.5*(v(ia_up)+v(ia_lo))*((v(on) > 0.5*VDIG) ? 1:0)
Izb izb gnd ib_up ib_lo on vcvs func=0.5*(v(ib_up)+v(ib_lo))*((v(on) > 0.5*VDIG) ? 1:0)
Izc izc gnd ic_up ic_lo on vcvs func=0.5*(v(ic_up)+v(ic_lo))*((v(on) > 0.5*VDIG) ? 1:0)

; Park transform of circulating currents
Izdq0 izdn izqn iza izb izc omega ABCDQ0 GAIN=(2/3) TETA_GAIN=-2

; Circulating current control (izd_ref = izq_ref = 0)
UzdReg uzd_ref_mod gnd izdn blocked PI_REG MAX=1M/10 MIN=-1M/10 KP=10*5 KI=1000*10 CF=10k
UzqReg uzq_ref_mod gnd izqn blocked PI_REG MAX=1M/10 MIN=-1M/10 KP=10*5 KI=1000*10 CF=10k
UzdRef uzd_ref gnd uzd_ref_mod izqn omega vcvs func=v(uzd_ref_mod)+2*v(omega)*L_ARM*v(izqn)
UzqRef uzq_ref gnd uzq_ref_mod izdn omega vcvs func=v(uzq_ref_mod)-2*v(omega)*L_ARM*v(izdn)

; Inverse Park transform of circulanting current control signals
Uzabc uz_a uz_b uz_c uzd_ref uzq_ref gnd omega DQ0ABC GAIN=1 TETA_GAIN=-2

; Upper and lower arm drivers
DrvUa drv_ua gnd uf_a uz_a vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_a) - v(uz_a))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)
DrvUb drv_ub gnd uf_b uz_b vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_b) - v(uz_b))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)
DrvUc drv_uc gnd uf_c uz_c vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_c) - v(uz_c))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)

DrvLa drv_la gnd uf_a uz_a vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_a) - v(uz_a))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)
DrvLb drv_lb gnd uf_b uz_b vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_b) - v(uz_b))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)
DrvLc drv_lc gnd uf_c uz_c vdc_mea vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_c) - v(uz_c))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)

ends

