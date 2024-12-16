
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
     omega idc_pos vdc_mea blocked acon COMPL_CONV L=CNTR_L W_NOM=W_NOM
ends


// -----------------------------
//  COMPL_CONV
// ----------------------------

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
parameter W_NOM=1

// AC CURRENT AND VOLTAGE MEASUREMENTS (WITH FILTER)
Ia ia gnd za gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Ib ib gnd zb gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Ic ic gnd zc gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Va va gnd a3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Vb vb gnd b3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes
Vc vc gnd c3 gnd svcvs numer=[WN_1^2] denom=[WN_1^2,2*XI_1*WN_1,1] ic=0 skipck=yes

// POWER MEASUREMENT (WITH FILTER)
Pwr_mod pow_mod gnd ia ib ic va vb vc vcvs func=v(ia)*v(va)+v(ib)*v(vb)+v(ic)*v(vc)
Pwr pow gnd pow_mod gnd svcvs numer=[WN^2] denom=[WN^2,2*XI*WN,1] ic=0 skipck=yes

// DC VOLTAGE MEASUREMENT (WITH FILTER)
VdcMea vdc_mea  gnd dc_pos  dc_neg  svcvs dcgain=1/VDC_NOM numer=[(2*pi*70)^2] \
                                    denom=[(2*pi*70)^2,2*XI_1*(2*pi*70),1] ic=1 skipck=yes

// CIRCULATING CURRENT MEASUREMENT
Iza iza gnd ia_up ia_lo on vcvs func=0.5*(v(ia_up)+v(ia_lo))*((v(on) > 0.5*VDIG) ? 1:0)
Izb izb gnd ib_up ib_lo on vcvs func=0.5*(v(ib_up)+v(ib_lo))*((v(on) > 0.5*VDIG) ? 1:0)
Izc izc gnd ic_up ic_lo on vcvs func=0.5*(v(ic_up)+v(ic_lo))*((v(on) > 0.5*VDIG) ? 1:0)

// PLL (PCC VOLTAGE IS USED)
Ipll w_err    gnd  vdp   vqp  vcvs  func=v(vqp)/max(0.1,(v(vdp)^2+v(vqp)^2))
Dw      dw    gnd  w_err gnd  \
	svcvs numer=[0.8615G, 9.324M, 1515, 16.4] denom=[197136,452,1] maxdcgain=MAXDCGAIN  ic=0 skipck=yes
Wnew  omega   gnd  dw         vcvs  func=v(dw) + W_NOM 

// PARK TRANSFORM OF PCC VOLTAGE AND CURRENT (ABC TO DQ0) 
Vdq0_pos vdp_raw vqp_raw va vb vc omega ABCDQ0 GAIN=(2/3)/(VAC_CONV)      TETA_GAIN=+1
Idq0_pos idp_raw iqp_raw ia ib ic omega ABCDQ0 GAIN=(2/3)/(I_NOM*sqrt(2)) TETA_GAIN=+1

// ADDITIONAL FILTERS FOR D-Q, P-N VOLTAGES
Vdp_mod  vdp_mod gnd  vdp_raw gnd svcvs dcgain=1 numer=[1] denom=[1, 1/(2*pi*11)]
Vqp_mod  vqp_mod gnd  vqp_raw gnd svcvs dcgain=1 numer=[1] denom=[1, 1/(2*pi*11)]

// NOTCH FILTERS (REQUIRED TO ELIMINATE 2*W COMPONENT IN D-Q VOLTAGES AND CURRENTS) (POSITIVE SEQUENCE)
Vdp_filt vdp gnd vdp_mod gnd svcvs dcgain=1 numer=[(2*W_NOM)^2,0,1] denom=[(2*W_NOM)^2,2*W_NOM/10,1]
Vqp_filt vqp gnd vqp_mod gnd svcvs dcgain=1 numer=[(2*W_NOM)^2,0,1] denom=[(2*W_NOM)^2,2*W_NOM/10,1]
Idp_filt idp gnd idp_raw gnd svcvs dcgain=1 numer=[(2*W_NOM)^2,0,1] denom=[(2*W_NOM)^2,2*W_NOM/10,1]
Iqp_filt iqp gnd iqp_raw gnd svcvs dcgain=1 numer=[(2*W_NOM)^2,0,1] denom=[(2*W_NOM)^2,2*W_NOM/10,1]

// INNER CURRENT LOOP
LOOP_DP udp idp vdp iqp idp_ref omega INNER_LOOP KP=KP_I KI=KI_I L=L_PU GAIN=-1 W=W_NOM MAXGAIN=MAXDCGAIN
LOOP_QP uqp iqp vqp idp iqp_ref omega INNER_LOOP KP=KP_I KI=KI_I L=L_PU GAIN=+1 W=W_NOM MAXGAIN=MAXDCGAIN

// INVERSE PARK TRANSFORM OF REFERENCE VOLTAGES
UABCP uf_a uf_b uf_c udp uqp gnd omega DQ0ABC GAIN=VAC_CONV TETA_GAIN=+1

//  CIRCULATING CURRENT CONTROL
Izdq0n izdn izqn iza izb izc omega ABCDQ0 GAIN=(2/3) TETA_GAIN=-2

UzdReg uzd_mod gnd izdn blocked PI_REG MAX=1M/10 MIN=-1M/10 KP=KP_CIRC KI=KI_CIRC CF=10k
UzqReg uzq_mod gnd izqn blocked PI_REG MAX=1M/10 MIN=-1M/10 KP=KP_CIRC KI=KI_CIRC CF=10k 
UzdRef uzd_ref gnd uzd_mod izqn omega vcvs func=v(uzd_mod)+2*v(omega)*L_ARM*v(izqn)
UzqRef uzq_ref gnd uzq_mod izdn omega vcvs func=v(uzq_mod)-2*v(omega)*L_ARM*v(izdn)

Uzabc uz_a uz_b uz_c uzd_ref uzq_ref gnd omega DQ0ABC GAIN=1 TETA_GAIN=-2

//  COMPUTATION OF DRIVING SIGNALS FOR EACH ARM
DrvUa drv_ua gnd uf_a uz_a vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_a) - v(uz_a))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))
DrvUb drv_ub gnd uf_b uz_b vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_b) - v(uz_b))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))
DrvUc drv_uc gnd uf_c uz_c vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) - v(uf_c) - v(uz_c))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))
DrvLa drv_la gnd uf_a uz_a vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_a) - v(uz_a))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))
DrvLb drv_lb gnd uf_b uz_b vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_b) - v(uz_b))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))
DrvLc drv_lc gnd uf_c uz_c vdc_mea blocked vcvs func=limit((VDC_NOM/2*v(vdc_mea) + v(uf_c) - v(uz_c))/(max(0.1,v(vdc_mea))*VDC_NOM),0,1)*(1-v(blocked))

ends

