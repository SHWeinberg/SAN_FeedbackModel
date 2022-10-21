clear
tic


%% constants
Rc = 8314.472;
T = 310;
F = 96485.3415;
C = 3.2e-5;
clamp_mode = 0;
t_holding = 0.5;
t_test = 0.5;

V_test = -35;
V_holding = -45;
Nao = 140;
Ki = 140;
Ko = 5.4;
Kc = Ko;
Cao = 1.8;
Nai_clamp = 0;
g_Na_L = 0;
Km_f = 45;
alpha = 0.5927;
blockade = 0;
y_shift = 0;
Km_Kp = 1.4;
Km_Nap = 14;
Qci = 0.1369;
Qn = 0.4315;
Qco = 0;
K3ni = 26.44;
Kci = 0.0207;
K1ni = 395.3;
K2ni = 2.289;
Kcni = 26.44;
K3no = 4.663;
K1no = 1628;
K2no = 561.4;
Kco = 3.663;
blockade_NaCa = 0;
delta_m = 1e-5;
k_dL = 4.3371;
V_dL = -16.4508;
shift_fL = 0;
k_fL = 0;
alpha_fCa = 0.0075;
Km_fCa = 0.000338;
offset_fT = 0;
MaxSR = 15;
MinSR = 1;
EC50_SR = 0.45;
HSR = 2.5;
koCa = 10000;
kiCa = 500;
kim = 5;
kom = 60;
tau_dif_Ca = 4e-5;
tau_tr = 0.04;
K_up = 0.0006;
slope_up = 5e-5;
TC_tot = 0.031;
TMC_tot = 0.062;
CM_tot = 0.045;
CQ_tot = 10;
kf_TC = 88800;
kf_TMM = 2277;
kf_TMC = 227700;
kf_CM = 227700;
kf_CQ = 534;
kb_TC = 446;
kb_TMC = 7.51;
kb_TMM = 751;
kb_CM = 542;
kb_CQ = 445;
Mgi = 2.5;
V_jsr_part = 0.0012;
V_i_part = 0.46;
V_nsr_part = 0.0116;
R_cell = 4;
L_cell = 70;
L_sub = 0.02;
RTONF = ( Rc.*T)./F;

%volumes
V_cell =  1e-09.* pi.*power(R_cell, 2).*L_cell;
V_sub =  1e-09.*2.0* pi.*L_sub.*(R_cell - L_sub./2).*L_cell;
V_jsr =  V_jsr_part.*V_cell;
V_i =  V_i_part.*V_cell - V_sub;
V_nsr =  V_nsr_part.*V_cell;

%% initial conditions
V = -40;
Nai = 5;
y = 0.181334538702451;
m = 0.440131579215766;
h = 1.3676940140066e-5;
dL = 0.001921;
fL = 0.846702;
fCa = 0.844449;
dT = 0.268909;
fT = 0.020484;
R = 0.912317231017262;
O = 1.7340201253e-7;
I = 7.86181717518e-8;
RI = 0.211148145512825;
Ca_sub = 1.0779e-04*0;
Ca_jsr = 0.2605*0;
Ca_nsr = 0.2605*0;
Cai = 5.8397e-05*0;
fTMM=0.501049376634;
fCMi=0.0373817991524254;
fCMs=0.054381370046;
fTC=0.0180519400676086;
fTMC=0.281244308217086;
fCQ=0.299624275428735;
q = 0.430836;
r = 0.014523;
paS = 0.322999177802891;
paF = 0.0990510403258968;
piy = 0.705410877258545;
n = 0.1162;
a = 0.00277;

%% time + ca_tgt
time = 3600; %s
dt = 1e-4; %s
t_n = time/dt; 

Ca_tgt = 2e-04;



%% conductances
%ryr and serca
ks = 250000000;
P_up_basal_init = 12;

%conductance initial conditions
start_scale = .00001;%.025;
g_Na = 0.0125*start_scale;     
g_f = 0.03*start_scale;
g_to = 0.002*start_scale;
g_Kr = 0.0021637*start_scale;
g_Ks = 0.0016576*start_scale;
i_NaK_max = 0.063*start_scale;
K_NaCa = 4*start_scale;
P_CaL = 0.2*start_scale;
P_CaT = 0.02*start_scale;
P_up_basal = P_up_basal_init * 2;

%m's
m_g_Na = 0.0125*start_scale;
m_g_f = 0.03*start_scale;
m_g_to = 0.002*start_scale;
m_g_Kr = 0.0021637*start_scale;
m_g_Ks = 0.0016576*start_scale;
m_i_NaK_max = 0.063*start_scale;
m_K_NaCa = 4*start_scale;
m_P_CaL = 0.2*start_scale;
m_P_CaT = 0.02*start_scale; 
m_P_serca = P_up_basal_init * 2;

%conductance time constants
tau_g = 1;% 1
tau_g_x_scale = .5; %.5
tau_g_Na = 1*tau_g_x_scale;
tau_g_f = 0.4167*tau_g_x_scale;
tau_g_to = 6.25*tau_g_x_scale;
tau_g_Kr = 5.7771*tau_g_x_scale;
tau_g_Ks = 7.5410*tau_g_x_scale;
tau_i_NaK_max = 0.1984*tau_g_x_scale;
tau_K_NaCa = 0.0031*tau_g_x_scale;
tau_P_CaL = 0.0625*tau_g_x_scale;
tau_P_CaT = 0.6250*tau_g_x_scale;
tau_P_serca = -0.001041*tau_g_x_scale;






%% save data
k_save = 1;
data_step = 10;
save_size = fix(t_n/data_step);
data_V = zeros(save_size,1);
data_Cai = zeros(save_size,1);
data_g_Na = zeros(save_size,1);
data_g_f = zeros(save_size,1);
data_g_to = zeros(save_size,1);
data_g_Kr = zeros(save_size,1);
data_g_Ks = zeros(save_size,1);
data_i_NaK_max = zeros(save_size,1);
data_K_NaCa = zeros(save_size,1);
data_P_CaL = zeros(save_size,1);
data_P_CaT = zeros(save_size,1);
data_P_serca = zeros(save_size,1);

%% ACh and iso
g_KACh = 0.00864;
ACh = 0.00001;
Iso_1_uM = 0;
ACh_on = 0;

g_leak = 0;
bool_cond = 0;

g_f_scale = 1;
K_NaCa_scale = 1;

d_Ca = 0;
i_stim = 0;
min_chan = 1/40;
max_chan = 40;
k_stim = 0;

%% main for
for t = 1:t_n
    if rem(t,1000000) == 0
        disp(string(t) + "/"+string(t_n)+" " + string(t/t_n))
    end

    %KOs
    if t*dt >1000 && bool_cond == 0
        K_NaCa_scale = .25; tau_K_NaCa = 10000000000; %d_Ca = 25;     
        g_f_scale = 0.25; tau_g_f = 10000000000; %
        bool_cond = 1;
    end
    
    
   %pacing
%     if rem(t*dt,0.18) == 0 && (t*dt)>1000 &&(t*dt<2000)  && V<-30
%         i_stim = -1;
%         k_stim = 0;
%     elseif k_stim<5 && (t*dt)>750 &&(t*dt<3000)
%         i_stim = -1;
%         k_stim = k_stim+1;
%     else
%         i_stim = 0;
%     end
    

    %% ion conc/current functions
    %ion potentials
    E_Na =  RTONF.*log(Nao./Nai);
    E_K =  RTONF.*log(Ko./Ki);
    E_mh =  RTONF.*log((Nao+ 0.12*Ko)./(Nai+ 0.12*Ki));
    E_Ks =  RTONF.*log((Ko+ 0.12*Nao)./(Ki+ 0.12*Nai));
    
    %i_KACh
    if ACh_on
        alpha_a = (3.59880 - 0.0256410)./(1+1.21550e-06./power( 1.0*ACh, 1.69510))+0.0256410;
        beta_a =  10.*exp( 0.0133000.*(V+40.0000));
        a_infinity = alpha_a./(alpha_a+beta_a);
        tau_a = 1.0/(alpha_a+beta_a);
        a = a + dt*((a_infinity - a)./tau_a);
        i_KACh = ACh_on.*g_KACh.*(V - E_K).*(1+exp((V+20.0000)./20.0000)).*a;
    else
        i_KACh = 0;
    end
    
    %i_f 
    %i_f_y_gate
    
    if ACh_on
        ACh_shift = -1-9.898*power(1*ACh, 0.618)/(power(1*ACh, 0.618)+0.00122423);
    else
        ACh_shift = 0;
    end
    
    if Iso_1_uM
        Iso_shift_i_f = 7.5;
    else
        Iso_shift_i_f = 0;
    end
    
    tau_y = 0.7166529/(0.0708*exp(-(V+5-ACh_shift-Iso_shift_i_f)/20.2791)+10.6*exp((V-ACh_shift-Iso_shift_i_f)/18));
    y_infinity = 1/(1+exp((V+52.5-ACh_shift-Iso_shift_i_f)/9));
    y = y + dt * (y_infinity-y)/tau_y;

    i_fNa = (y^2)*Kc/(Kc+Km_f)*g_f_scale*g_f*(V-E_Na);%*ICs_on_Icontrol;
    i_fK = (y^2)*Kc/(Kc+Km_f)*g_f_scale*g_f*(V-E_K);%*ICs_on_Icontrol;
    i_f = (i_fNa+i_fK);
    
    %i_NaK
    if Iso_1_uM
        Iso_increase_i_NaK = 1.2;
    else
        Iso_increase_i_NaK =1;
    end

    i_NaK = Iso_increase_i_NaK *i_NaK_max*power(1+power(Km_Kp/Kc, 1.2), -1)*power(1+power(Km_Nap/Nai, 1.3), -1)*power(1+exp(-(V-E_Na+110)/20), -1);
    
    %i_Na
    E0_m = V+41.0000;
    if abs(E0_m)<delta_m
        alpha_m = 2000;
    else
        alpha_m = ( 200.000.*E0_m)./(1 - exp(- 0.10*E0_m));
    end
    
    
    beta_m =  8000.00.*exp(- 0.056*(V+66.0000));
    m = m + dt*(alpha_m*(1-m)-beta_m*m);
    
    alpha_h = 20*exp(-0.125*(V+75));
    beta_h = 2000/(320*exp(-0.1*(V+75))+1);
    h = h +dt*( alpha_h*(1-h)-beta_h*h);
    
    i_Na = g_Na*(m^3)*h*(V-E_mh);
    
    %i_to
    q_infinity = 1.0/(1+exp((V+49.0000)./13.0000));
    tau_q =  0.0010*0.60*(65.1700./( 0.57*exp(- 0.08*(V+44.0000))+ 0.065*exp( 0.10*(V+45.9300)))+10.1000);
    q = q+dt*((q_infinity - q)./tau_q);
    r_infinity = 1.0/(1+exp( - (V - 19.3000)./15.0000));
    tau_r =  0.0010*0.66*1.4*(15.5900./( 1.03700.*exp( 0.090*(V+30.6100))+ 0.369000.*exp(  - 0.12*(V+23.8400)))+2.98000);
    r = r + dt*((r_infinity - r)./tau_r);
    i_to = g_to*(V-E_K)*q*r;
    
    %i_Kr
    pa_infinity = 1/(1+exp(-(V+14.8)/8.5));
    tau_paS = 0.84655354/(4.2*exp(V/17)+0.15*exp(-V/21.6));
    tau_paF = 1/(30*exp(V/10)+exp(-V/12));

    paS = paS + dt*((pa_infinity - paS)./tau_paS);
    paF = paF + dt*((pa_infinity - paF)./tau_paF);
    
    tau_pi = 1/(100*exp(-V/54.645)+656*exp(V/106.157));
    pi_infinity = 1/(1+exp((V+28.6)/17.1));
    piy = piy + dt *((pi_infinity - piy)./tau_pi);
    i_Kr =  g_Kr.*(V - E_K).*( 0.90*paF+ 0.10*paS).*piy;
    
    %i_Ks
    if Iso_1_uM
        Iso_increase_g_Ks = 1.2;
        Iso_shift_Ks = -14;
    else
        Iso_shift_Ks = 0;
        Iso_increase_g_Ks = 1;
    end
    shift = 0;
    n_infinity = 14/(1+exp(-(V-40-Iso_shift_Ks)/12))/(14/(1+exp(-(V-40-Iso_shift_Ks)/12))+1*exp(-(V-Iso_shift_Ks)/45));
    alpha_n = 28/(1+exp(-(V-40-Iso_shift_Ks)/3));
    beta_n = 1*exp(-(V-Iso_shift_Ks-shift-5)/25);
    tau_n = 1/(alpha_n+beta_n);

    n = n+ dt*( (n_infinity - n)./tau_n);    
    i_Ks =  Iso_increase_g_Ks* g_Ks.*(V - E_Ks).*power(n, 2);
     
    %% Calcium
    %Ca buffering
    delta_fTMM =  kf_TMM.*Mgi.*(1 - (fTMC+fTMM)) -  kb_TMM.*fTMM;
    delta_fCMs =  kf_CM.*Ca_sub.*(1 - fCMs) -  kb_CM.*fCMs;
    delta_fTC =  kf_TC.*Cai.*(1 - fTC) -  kb_TC.*fTC;
    delta_fTMC =  kf_TMC.*Cai.*(1 - (fTMC+fTMM)) -  kb_TMC.*fTMC;
    delta_fCQ =  kf_CQ.*Ca_jsr.*(1 - fCQ) -  kb_CQ.*fCQ;
    delta_fCMi =  kf_CM.*Cai.*(1 - fCMi) -  kb_CM.*fCMi;
    
    fTMM = fTMM + dt .*delta_fTMM;
    fCMs = fCMs + dt .*delta_fCMs;
    fTC = fTC + dt.*delta_fTC;
    fTMC = fTMC + dt.*delta_fTMC;
    fCQ = fCQ + dt.*delta_fCQ;
    fCMi = fCMi + dt.*delta_fCMi;

    
    %i_CaL
    if Iso_1_uM
        Iso_shift_CaL = -8;
        Iso_slope_CaL = 0.69;
        Iso_increase_i_CaL = 1.23;
    else
        Iso_shift_CaL = 0;
        Iso_slope_CaL = 1;
        Iso_increase_i_CaL = 1;
    end
    
    if V== - 41.8000
        adVm = - 41.8001;
    elseif V==0
        adVm = - 41.8000 == 0;
    elseif V== - 6.80000
        adVm = - 6.80001;
    else
        adVm = V;
    end

    if V== - 1.80000
        bdVm =- 1.80001;
    else
        bdVm = V;
    end

    dL_infinity = 1/(1+exp(-(V+20.3-Iso_shift_CaL)/(Iso_slope_CaL*4.2)));
    alpha_dL = -0.02839*(adVm+41.8-Iso_shift_CaL)/(exp(-(adVm+41.8-Iso_shift_CaL)/2.5)-1)-0.0849*(adVm+6.8-Iso_shift_CaL)/(exp(-(adVm+6.8-Iso_shift_CaL)/4.8)-1);
    beta_dL = 0.01143*(bdVm+1.8-Iso_shift_CaL)/(exp((bdVm+1.8-Iso_shift_CaL)/2.5)-1);
    tau_dL = 0.001/(alpha_dL+beta_dL);

    dL = dL + dt*((dL_infinity - dL)./tau_dL);
    fL_infinity = 1/(1+exp((V+37.4)/5.3));
    tau_fL = 0.001*(44.3+230*exp(-((V+36)/10)^2));
    fL = fL + dt*((fL_infinity - fL)./tau_fL);

    fCa_infinity = Km_fCa/(Km_fCa+Ca_sub);
    tau_fCa = 0.001*fCa_infinity/alpha_fCa;
    fCa = fCa+dt*((fCa_infinity - fCa)./tau_fCa);

    i_siCa = 2*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)*2/RTONF)))*(Ca_sub-Cao*exp(-2*(V-0)/RTONF))*dL*fL*fCa;
    i_siK = 0.000365*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)/RTONF)))*(Ki-Kc*exp(-1*(V-0)/RTONF))*dL*fL*fCa;
    i_siNa = 0.0000185*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)/RTONF)))*(Nai-Nao*exp(-1*(V-0)/RTONF))*dL*fL*fCa;
    ACh_block = ACh_on * 0.31*ACh/(ACh+0.00009);
    i_CaL = (i_siCa+i_siK+i_siNa)*(1-ACh_block)*Iso_increase_i_CaL;
    
    %i_CaT
    fT_infinity = 1.0/(1+exp((V+58.7000)./3.80000));
    tau_fT = 1.0/( 16.6700.*exp( - (V+75.0000)./83.3000)+ 16.6700.*exp((V+75.0000)./15.3800));
    fT = fT + dt*((fT_infinity - fT)./tau_fT);

    dT_infinity = 1.0/(1+exp( - (V+38.3000)./5.50000));
    tau_dT = 0.001/( 1.06800.*exp((V+38.3000)./30.0000)+ 1.06800.*exp( - (V+38.3000)./30.0000));
    dT = dT + dt*((dT_infinity - dT)./tau_dT);
    i_CaT =  (( 2.0*P_CaT.*V)./( RTONF.*(1 - exp(( - 1.0*V.*2)./RTONF)))).*(Ca_sub -  Cao.*exp((  - 2.0*V)./RTONF)).*dT.*fT;


    %operator split for Ca and i_NaCa
    exp_Qco = exp(Qco.*V./RTONF);
    exp_mQn = exp(-Qn.*V./(2.*RTONF));
    exp_Qn = exp(Qn.*V./(2.*RTONF));
    exp_mQci = exp(-Qci.*V./RTONF);
        
    %operator splitting and variable time steps for Ca 
    max_g_Na = max(g_Na,[],'all');
    max_g_Na_scale = max_g_Na/0.0125;
    %>3, 20, >7, 30; >11, 40, >15, 50; >19,60; >23, 70; >27,80; >31, 90; >35,100; >39, 11
    Ca_steps = 10;
    
    if max_g_Na_scale > 3
        Ca_steps = 20;
    end
    
    if max_g_Na_scale > 7
        Ca_steps = 30;
    end

    if max_g_Na_scale > 11
        Ca_steps = 40;
    end

    if max_g_Na_scale > 15
        Ca_steps = 50;
    end

    if max_g_Na_scale > 19
        Ca_steps = 60;
    end
    
    if max_g_Na_scale > 23
        Ca_steps = 70;
    end
    
    if max_g_Na_scale > 27
        Ca_steps = 80;
    end

    if max_g_Na_scale > 31
        Ca_steps = 90;
    end

    if max_g_Na_scale > 35
        Ca_steps = 100;
    end

    if max_g_Na_scale > 39
        Ca_steps = 110;
    end
   
    
    dt_Ca = dt./Ca_steps;
    for k_Ca = 1:Ca_steps   
        %Ca_SR_release
        Ca_jsr_pow = Ca_jsr.*Ca_jsr.*sqrt(Ca_jsr);
        kCaSR = MaxSR - (MaxSR - MinSR)./(1+(EC50_SR.^HSR./(Ca_jsr_pow)));
        koSRCa = koCa./kCaSR;
        kiSRCa =  kiCa.*kCaSR;
        
        delta_R = ( kim.*RI -  kiSRCa.*Ca_sub.*R) - ( koSRCa.*power(Ca_sub, 2).*R -  kom.*O);
        delta_O = ( koSRCa.*power(Ca_sub, 2).*R -  kom.*O) - ( kiSRCa.*Ca_sub.*O -  kim.*I);
        delta_I = ( kiSRCa.*Ca_sub.*O -  kim.*I) - ( kom.*I -  koSRCa.*power(Ca_sub, 2).*RI);
        delta_RI = ( kom.*I -  koSRCa.*power(Ca_sub, 2).*RI) - ( kim.*RI -  kiSRCa.*Ca_sub.*R);
        
        R = R + dt_Ca .* delta_R;
        O = O + dt_Ca .* delta_O ;
        I = I + dt_Ca .* delta_I;
        RI = RI + dt_Ca .* delta_RI;
        
        j_SRCarel =  ks.*O.*(Ca_jsr - Ca_sub);
        
        
        %i_NaCa
        k34 = Nao./(K3no+Nao);
        k41 = exp_mQn;
        do = 1+Cao./Kco.*(1+exp_Qco)+Nao./K1no.*(1+Nao./K2no.*(1+Nao./K3no));
        k23 = Nao./K1no.*Nao./K2no.*(1+Nao./K3no).*exp_mQn./do;
        k21 = Cao./Kco.*exp_Qco./do;
        k32 = exp_Qn;
        k43 = Nai./(K3ni+Nai);
        x1 = k41.*k34.*(k23+k21)+k21.*k32.*(k43+k41);
        di = 1+Ca_sub./Kci.*(1+exp_mQci+Nai./Kcni)+Nai./K1ni.*(1+Nai./K2ni.*(1+Nai./K3ni));
        k12 = Ca_sub./Kci.*exp_mQci./di;
        k14 = Nai./K1ni.*Nai./K2ni.*(1+Nai./K3ni).*exp_Qn./di;
        x2 = k32.*k43.*(k14+k12)+k41.*k12.*(k34+k32);
        x3 = k14.*k43.*(k23+k21)+k12.*k23.*(k43+k41);
        x4 = k23.*k34.*(k14+k12)+k14.*k21.*(k34+k32);
        i_NaCa = K_NaCa_scale * K_NaCa.*(x2.*k21-x1.*k12)./(x1+x2+x3+x4);
        
        %Cai fluxes
        b_up = 0;
        if Iso_1_uM
            b_up = -0.25;
        end
        if ACh_on
           b_up =  0.7.*ACh./(0.00009+ACh);
        end
        
        P_up =  P_up_basal.*(1 - b_up); 
        j_Ca_dif = (Ca_sub - Cai)./tau_dif_Ca;
        j_up = P_up./(1+K_up./Cai);
        j_tr = (Ca_nsr - Ca_jsr)./tau_tr;

        %Ca conc
        Cai = Cai + dt_Ca.*( ( 1.0.*( j_Ca_dif.*V_sub -  j_up.*V_nsr))./V_i - ( CM_tot.*delta_fCMi+ TC_tot.*delta_fTC+ TMC_tot.*delta_fTMC));
        Ca_sub = Ca_sub + dt_Ca.*(( j_SRCarel.*V_jsr)./V_sub - ((i_siCa + i_CaT - 2.0.*i_NaCa)./( 2.0.*F.*V_sub)+j_Ca_dif+ CM_tot.*delta_fCMs));
        Ca_nsr = Ca_nsr + dt_Ca.*(j_up - (j_tr.*V_jsr)./V_nsr);
        Ca_jsr = Ca_jsr + dt_Ca.*(j_tr - (j_SRCarel+ CQ_tot.*delta_fCQ));
        Nai = Nai + dt_Ca.*(( (1 - Nai_clamp).* - 1.0.*(i_Na+i_fNa+i_siNa+ 3.0.*i_NaK+ 3.0.*i_NaCa))./( 1.0.*(V_i+V_sub).*F));
    end
            
    Cai(Cai<0) = 0;


    
    if ~isreal(Cai)
        disp('Complex Cai');
        break
    end
    
    if isnan(Cai)
        disp('Complex Cai');
        break
    end
    
    if Cai<0
        Cai = 0;
    end
    
    %% current and voltages
    i_tot = i_f + i_Kr + i_Ks + i_to + i_NaK + i_NaCa + i_Na + i_CaL + i_CaT + i_KACh + i_stim;
    V =  V + dt*(- i_tot./C);
    
    
    if ~isreal(V)
        disp('Complex V');
        break
    end

    %% variable conductances 
    [g_Na,m_g_Na] = variable_cond(g_Na,m_g_Na,tau_g_Na,dt,Ca_tgt,Cai,tau_g);
    [g_f,m_g_f] = variable_cond(g_f,m_g_f,tau_g_f,dt,Ca_tgt,Cai,tau_g);
    [g_to,m_g_to] = variable_cond(g_to,m_g_to,tau_g_to,dt,Ca_tgt,Cai,tau_g);
    [g_Kr,m_g_Kr] = variable_cond(g_Kr,m_g_Kr,tau_g_Kr,dt,Ca_tgt,Cai,tau_g);
    [g_Ks,m_g_Ks] = variable_cond(g_Ks,m_g_Ks,tau_g_Ks,dt,Ca_tgt,Cai,tau_g);
    [i_NaK_max,m_i_NaK_max] = variable_cond(i_NaK_max,m_i_NaK_max,tau_i_NaK_max,dt,Ca_tgt,Cai,tau_g);
    [K_NaCa,m_K_NaCa] = variable_cond(K_NaCa,m_K_NaCa,tau_K_NaCa,dt,Ca_tgt,Cai,tau_g);
    [P_CaL,m_P_CaL] = variable_cond(P_CaL,m_P_CaL,tau_P_CaL,dt,Ca_tgt,Cai,tau_g);
    [P_CaT,m_P_CaT] = variable_cond(P_CaT,m_P_CaT,tau_P_CaT,dt,Ca_tgt,Cai,tau_g);
    [P_up_basal,m_P_serca] = variable_cond(P_up_basal,m_P_serca,tau_P_serca,dt,Ca_tgt,Cai,tau_g);
    if m_g_Na<0.0125*min_chan m_g_Na=0.0125*min_chan; end
    if m_g_f<0.03*min_chan m_g_f = 0.03*min_chan; end
    if m_g_to<0.002*min_chan m_g_to = 0.002*min_chan; end
    if m_g_Kr<0.0021637*min_chan m_g_Kr  = 0.0021637*min_chan; end
    if m_g_Ks<0.0016576*min_chan m_g_Ks = 0.0016576*min_chan; end
    if m_i_NaK_max<0.063*min_chan m_i_NaK_max = 0.063*min_chan; end
    if m_K_NaCa<4*min_chan m_K_NaCa = 4*min_chan; end
    if m_P_CaL<0.2*min_chan m_P_CaL = 0.2*min_chan; end
    if m_P_CaT<0.02*min_chan m_P_CaT = 0.02*min_chan; end

    if m_P_serca<P_up_basal_init./10 m_P_serca = P_up_basal_init./10; end
    if m_g_Na>0.0125*max_chan m_g_Na=0.0125*max_chan; end
    if m_g_f>0.03*max_chan m_g_f = 0.03*max_chan; end
    if m_g_to>0.002*max_chan m_g_to = 0.002*max_chan; end
    if m_g_Kr>0.0021637*max_chan m_g_Kr  = 0.0021637*max_chan; end
    if m_g_Ks>0.0016576*max_chan m_g_Ks = 0.0016576*max_chan; end
    if m_i_NaK_max>0.063*max_chan m_i_NaK_max = 0.063*max_chan; end
    if m_K_NaCa>4*max_chan m_K_NaCa = 4*max_chan; end
    if m_P_CaL>0.2*max_chan m_P_CaL = 0.2*max_chan; end
    if m_P_CaT>0.02*max_chan m_P_CaT = 0.02*max_chan; end
    if m_P_serca>P_up_basal_init*10 m_P_serca = P_up_basal_init*10; end

    %% save data
    if rem(t,data_step) == 0 
        data_V(k_save) = V;
        data_Cai(k_save) = Cai;
        
        data_g_Na(k_save) = g_Na;
        data_I_stim(k_save) = i_stim;
        data_g_f(k_save) =  g_f;
        data_g_to(k_save) = g_to;
        data_g_Kr(k_save) =  g_Kr;
        data_g_Ks(k_save) = g_Ks;
        data_i_NaK_max(k_save) =  i_NaK_max;
        data_K_NaCa(k_save) =  K_NaCa;
        data_P_CaL(k_save) = P_CaL;
        data_P_CaT(k_save) = P_CaT;
        data_P_serca(k_save) = P_up_basal;
        k_save = k_save+1;
    end

end


%% feedback model function
function [g_x,m_g_x] = variable_cond(g_x,m_g_x,tau_g_x,dt,Ca_tgt,Cai,tau_g) 
    m_g_x = m_g_x + (Ca_tgt - Cai)*dt/tau_g_x;
    g_x = g_x + (m_g_x - g_x)*dt/tau_g;
    if g_x<0
        g_x = 0;
    end
end

