function [dY, out] = bondarenko2004(time, Y)
% --------------------------- DESCRIPTION -----------------------------
% Bondarenko's Et al. mouse ventricular myocyte action potential model.
% Source: https://doi.org/10.1152/ajpheart.00185.2003
%
% Code Author: Diego Gazzoni
% Date: October, 2023
% Project: Master's Degree Thesis
% ---------------------------------------------------------------------
    type = "a"; % s = septum, a = apex

    % * Cell Geometry *************************** (T2)
    Acap = 1.534e-4;   % cm^2
    Vmyo = 25.84e-6;   % uL
    Vjsr = 0.12e-6;    % uL
    Vnsr = 2.098e-6;   % uL
    Vss =  1.485e-9;   % uL

    % * Extracellular Ion Concentrations ******** (T3)
    Ko =  5400;        % uM
    Nao = 140000;      % uM
    Cao = 1800;        % uM

    % * General ********************************* (T7)
    C = 1.0;           % uF/cm^2
    F = 96.5;          % C/mM;
    T = 25 + 273.15;   % K
    R = 8.314;         % J/M*K

    %% * Variables ******************************* (order from T8)
    V     =     Y(1);  % Membrane voltage
    Cai   =     Y(2);  % Intracellular calcium
    Cass  =     Y(3);  % Subspace calcium
    Cajsr =     Y(4);  % Junctional calcium
    Cansr =     Y(5);  % Network calcium
    LTRPNCa =   Y(6);  % Concentration of Ca2+ bound to low affinity troponin-binding sites
    HTRPNCa =   Y(7);  % Concentration of Ca2+ bound to high affinity troponing-binding sites
    O =         Y(8);  % L-Ca opening probability
    C1 =        Y(9);  % L-Ca closed proability                                     --> phantom variable
    C2 =        Y(10); % L-Ca closed probability
    C3 =        Y(11); % L-Ca closed probability
    C4 =        Y(12); % L-Ca closed probability
    I1 =        Y(13); % L-Ca inactivation probability 
    I2 =        Y(14); % L-Ca inactivation probability
    I3 =        Y(15); % L-Ca inactivation probability 
    PC1 =       Y(16); % Fraction of RyR channels in state PC1                      --> phantom variable
    PC2 =       Y(17); % Fraction of RyR channels in state PC2
    PO1 =       Y(18); % Fraction of RyR channels in state PO1
    PO2 =       Y(19); % Fraction of RyR channels in state PO2
    PRyR =      Y(20); % RyR modulation factor
    CNa3 =      Y(21); % Closed state of fast Na channel                            --> phantom variable
    CNa2 =      Y(22); % Closed state of fast Na channel
    CNa1 =      Y(23); % Closed state of fast Na channel
    ONa =       Y(24); % Open state of fast Na channel
    IFNa =      Y(25); % Fast inactivated state of fast Na channel
    I1Na =      Y(26); % Slow inactivated state 1 of fast Na channel
    I2Na =      Y(27); % Slow inactivated state 2 of fast Na channel
    ICNa2 =     Y(28); % Closed-Inactivated state of fast Na channel
    ICNa3 =     Y(29); % Closed-Inaactivated state of fast Na channel
    Nai =       Y(30); % Myoplasmic Na [uM]
    Ki =        Y(31); % Myoplasmic K  [uM]
    atof =      Y(32); % Gating variable for transient out K+ current
    itof =      Y(33); % Gating variable for transient out K+ current
    nKs  =      Y(34); % Gating variable for slow delayed-rectified K+ current
    atos =      Y(35); % Gating variable for transient out K+ current
    itos =      Y(36); % Gating variable for transient out K+ current
    aur  =      Y(37); % Gating variable for ultrarapidly activating delayed-rect K+ current
    iur =       Y(38); % Gating variable for ultrarapidly activating delayed-rect K+ current
    akss =      Y(39); % Gating variable for noninactivating steady-state K+ current
    ikss =      Y(40); % Gating variable for noninactivating steady-state K+ current --> phantom variable
    CK0 =       Y(41); % mERG closed                                                 --> phantom variable
    CK1 =       Y(42); % mERG closed
    CK2 =       Y(43); % mERG closed
    OK  =       Y(44); % mERG open
    IK  =       Y(45); % mERG inactivated state

    %% * L-Type calcium current (ICaL) **********************
    % This current is described thorugh a Markov kinetics.
    GCaL = 0.1729;    % mS/uF Specific conductance for L-Type channel.
    ECaL = 63.0;      % mV Resting potential for L-type channel.
    Kpcmax = 0.23324; % 1/ms Maximum time constant for Ca-induced inactivation
    Kpchalf = 20.0;   % uM Half-saturation constant for Ca2-induced inactivation
    Kpcb = 0.0005;    % 1/ms Voltage-insensitive rate constant for inactivation
    ICaLmax = 7.0;    % pA/pF Normalization constant for L-type Ca current

    alpha = (0.4*exp((V+12.0)/10.0)*(1 + 0.7*exp(-((V+40.0)^2)/10.0) - 0.75*exp(-((V+20.0)^2)/400.0))) / (1 + 0.12*exp((V+12.0)/10.0));
    beta = 0.05*exp(-(V+12.0)/13.0);
    gamma = (Kpcmax * Cass) / (Kpchalf + Cass);
    Kpcf = 13.0 * (1 - exp(-(V+14.5)^2)/100.0);
    
    dO = alpha*C4 - 4*beta*O + Kpcb*I1 - gamma*O + 0.001*(alpha*I2 - Kpcf*O);
    C1 = 1 - (O+C2+C3+C4+I1+I2+I3);
    dC2 = 4.0*alpha*C1 - beta*C2 + 2.0*beta*C3 - 3.0*alpha*C2;
    dC3 = 3.0*alpha*C2 - 2.0*beta*C3 + 3.0*beta*C4 - 2.0*alpha*C3;
    dC4 = 2.0*alpha*C3 - 3.0*beta*C4 + 4.0*beta*O - alpha*C4 + 0.01*(4*Kpcb*beta*I1 - alpha*gamma*C4) + 0.002*(4*beta*I2 - Kpcf*C4) + 4*beta*Kpcb*I3 - gamma*Kpcf*C4;
    dI1 = gamma*O - Kpcb*I1 + 0.001*(alpha*I3 - Kpcf*I1) + 0.01*(alpha*gamma*C4 - 4*beta*Kpcb*I1);
    dI2 = 0.001*(Kpcf*O - alpha*I2) + Kpcb*I3 - gamma*I2 + 0.002*(Kpcf*C4 - 4*beta*I2);
    dI3 = 0.001*(Kpcf*I1 - alpha*I3) + gamma*I2 - Kpcb*I3 + gamma*Kpcf*C4 - 4*beta*Kpcb*I3;
    
    dY(8,1) = dO; dY(9,1) = 0; dY(10,1) = dC2; dY(11,1) = dC3; dY(12,1) = dC4; dY(13, 1) = dI1; dY(14,1) = dI2; dY(15,1) = dI3;
    
    ICaL = GCaL * O * (V - ECaL);

    dPRyR = -0.04*PRyR - 0.1 * (ICaL/ICaLmax)*exp(-((V-5.0)^2) / 648.0);
    dY(20,1) = dPRyR;

    %% * Ca++ pump (IpCa) ***********************************
    IpCamax = 1.0; % pA/pF Maximum ICa pump current
    KmpCa = 0.5;   % uM half sat constant for Ca++ pump current
    IpCa = IpCamax * Cai^2 / (KmpCa^2 + Cai^2);

    %% * Na/Ca exchanger (INaCa) ****************************
    kNaCa = 292.8; % pA/pF Scaling factor for Na/Ca exchanger
    KmNa = 87500;  % uM Na+ half saturation constant
    KmCa = 1380;   % uM Ca++ half saturation constant
    ksat = 0.1;    % saturation factor at every negative potential
    ni = 0.35;     % Controls voltage dependence of Na/Ca exchanger
    
    INaCa = kNaCa * 1.0 / (KmNa^3 + Nao^3) * 1.0 / (KmCa + Cao) * 1.0 / (1.0 + ksat*exp((ni-1)*V*F/(R*T))) * (exp(ni*V*F/(R*T))* Nai^3 * Cao - exp((ni-1)*V*F/(R*T))*(Nao^3)*Cai);

    %% * Ca++ background current (ICab) **********************
    GCab = 0.000367; % mS/uF Maximum backgroud Ca current conductance
    ECaN = (R*T/(2*F))*log(Cao/Cai);

    ICab = GCab*(V - ECaN);

    %% * Na+ fast current (INa) ******************************
    GNa = 13.0; % mS/uF
    nn = 0.9*Nao + 0.1*Ko;
    qq = 0.9*Nai + 0.1*Ki;
    ENa = ((R*T)/F) * log(nn/qq); % mV

    alphaNa11 = 3.802 / (0.1027*exp(-(V+2.5)/17.0) + 0.20*exp(-(V+2.5)/150.0));
    alphaNa12 = 3.802 / (0.1027*exp(-(V+2.5)/15.0) + 0.23*exp(-(V+2.5)/150.0));
    alphaNa13 = 3.802 / (0.1027*exp(-(V+2.5)/12.0) + 0.25*exp(-(V+2.5)/150.0));
    betaNa11 = 0.1917*exp(-(V+2.5)/20.3);
    betaNa12 = 0.20*exp(-(V-2.5)/20.3);
    betaNa13 = 0.22*exp(-(V-7.5)/20.3);

    alphaNa3 = (7.0e-7) * exp(-(V + 7.0) / 7.7);
    betaNa3 =  0.0084 + 0.00002*(V + 7.0);
    alphaNa2 = 1.0 / (0.188495 * exp(-(V + 7.0) / 16.6) + 0.393956);
    betaNa2 =  alphaNa13 * alphaNa2 * alphaNa3 / (betaNa13 * betaNa3);
    alphaNa4 = alphaNa2 / 1000.0;
    betaNa4 =  alphaNa3;
    alphaNa5 = alphaNa2 / 95000.0;
    betaNa5 =  alphaNa3 / 50.0;

    CNa3 = 1.0 - (ONa + CNa1 + CNa2 + IFNa + I1Na + I2Na + ICNa2 + ICNa3);
    dCNa2 =  alphaNa11*CNa3 - betaNa11*CNa2 + betaNa12*CNa1 - alphaNa12*CNa2 + alphaNa3*ICNa2 - betaNa3*CNa2;
    dCNa1 =  alphaNa12*CNa2 - betaNa12*CNa1 + betaNa13*ONa - alphaNa13*CNa1 + alphaNa3*IFNa - betaNa3*CNa1;
    dONa =   alphaNa13*CNa1 - betaNa13*ONa + betaNa2*IFNa - alphaNa2*ONa;
    dIFNa =  alphaNa2*ONa - betaNa2*IFNa + betaNa3*CNa1 - alphaNa3*IFNa + betaNa4*I1Na - alphaNa4*IFNa + alphaNa12*ICNa2 - betaNa12*IFNa;
    dI1Na =  alphaNa4*IFNa - betaNa4*I1Na + betaNa5*I2Na - alphaNa5*I1Na;
    dI2Na =  alphaNa5*I1Na - betaNa5*I2Na;
    dICNa2 = alphaNa11*ICNa3 - betaNa11*ICNa2 + betaNa12*IFNa - alphaNa12*ICNa2 + betaNa3*CNa2 - alphaNa3*ICNa2;
    dICNa3 = betaNa11*ICNa2 - alphaNa11*ICNa3 + betaNa3*CNa3 - alphaNa3*ICNa3;
    
    dY(21,1)=0; dY(22,1)=dCNa2; dY(23,1)=dCNa1; dY(24,1)=dONa; dY(25,1)=dIFNa; dY(26,1)=dI1Na; dY(27,1)=dI2Na; dY(28,1)=dICNa2; dY(29,1)=dICNa3;

    INa = GNa * ONa * (V - ENa);

    %% * Na+ background **************************************
    GNab = 0.0026;      % mS/uF maximum sodium background conductance

    INab = GNab * (V - ENa);

    %% * K+ transient outward (IKto,f and IKto,s) ************
    EK = ((R*T)/F) * log(Ko / Ki);

    if type == "s"
        GKtof = 0.0798; % mS/uF
        GKtos = 0.0629; % mS/uF
    elseif type == "a"
        GKtof = 0.4067; % mS/uF
        GKtos = 0;      % mS/uF
    end

    % --- Fast component ---
    alphaa = 0.18064*exp(0.03577*(V+30.0));
    betaa = 0.3956*exp(-0.06237*(V+30.0));
    alphai = (0.000152*exp(-(V+13.0)/7.0)) / (1+0.067083*exp(-(V+33.5)/7.0));
    betai = (0.00095*exp((V+33.5)/7.0)) / (1+0.051335*exp((V+33.5)/7.0));

    datof = alphaa*(1-atof) - betaa*atof;
    ditof = alphai*(1-itof) - betai*itof;

    dY(32,1) = datof; dY(33,1) = ditof;

    IKtof = GKtof * (atof^3) * itof * (V - EK);

    % --- Slow component ---
    ass = 1.0/(1.0+exp(-(V+22.5)/7.7));
    iss = 1.0/(1.0+exp((V+45.2)/5.7));
    tautas = 0.493*exp(-(0.0629*V)) + 2.058;
    tautis = 270.0 + 1050.0 / (1.0 + exp((V+45.2)/5.7));

    datos = (ass - atos)/tautas;
    ditos = (iss - itos)/tautis;
    
    dY(35,1) = datos; dY(36,1) = ditos;

    IKtos = GKtos * atos * itos*  (V - EK);

    %% * Time-independent K+ current (IK1) *****************
    IK1 = 0.2938 * (Ko/(Ko+210.0)) * ((V-EK)/(1.0 + exp(0.0896*(V-EK))));

    %% * Slow delayed rectifier K+ (IKs) *******************
    GKs = 0.00575;     % mS/uF maximum normalized conductance

    alphan = (0.00000481333 * (V+26.5)) / (1.0 - exp(-0.128*(V+26.5)));
    betan = 0.0000953333*exp(-0.038*(V+26.5));

    dnKs = alphan * (1.0 - nKs) - betan * nKs;

    dY(34,1) = dnKs;

    IKs = GKs * (nKs^2) * (V - EK);

    %% * Ultrarapid activanting delayed rect K+ (IKur) *****
    if type == "s"
        GKur = 0.0975; % mS/uF maximum normalized conductance (septum)
    elseif type == "a"
        GKur = 0.160;  % mS/uF maximum normalized conductance (apex)
    end

    tauaur = 0.493*exp(-0.0629*V) + 2.058;
    tauiur = 1200.0 - 170.0/(1+exp((V+45.2)/5.7));

    daur = (ass - aur)/tauaur;
    diur = (iss - iur)/tauiur;

    dY(37,1) = daur; dY(38,1) = diur;

    IKur = (GKur) * aur * iur * (V - EK);

    %% * Noninactivating steady state K+ (IKss) ***********
    if type == "s"
        GKss = 0.0324; % mS/uF maximum normalized conductance (septum)
    elseif type == "a"
        GKss = 0.050;  % mS/uF maximum normalized conductance (apex)
    end
    
    tauKss = 39.3*exp(-0.0862*V) + 13.17;
    
    dakss = (ass - akss)/tauKss;
    dikss = 0;
    
    dY(39, 1) = dakss; dY(40,1) = dikss;

    IKss = (GKss)*akss*ikss*(V-EK);

    %% * Delayed rectifier K+ (mERG, IKr) ***********
    GKr = 0.078;   % mS/uF maximum normalized conductance
    kf = 0.02761;  % 1/ms Rate constant for rapid delayed-rectifier K+ current
    kb = 0.036778; % 1/ms Rate constant for rapid delayed-rectifier K+ current

    alphaa0 = 0.022348*exp(0.01176*V);
    betaa0 = 0.047002*exp(-0.0631*V);
    alphaa1 = 0.013733*exp(0.038198*V);
    betaa1 = 0.0000689*exp(-0.04178*V);
    alphaii = 0.090821*exp(0.023391*(V+5.0));
    betaii = 0.006497*exp(-0.03268*(V+5.0));

    CK0 = 1.0 - (CK1 + CK2 + OK + IK);
    dCK1 = alphaa0*CK0 - betaa0*CK1 + kb*CK2 - kf*CK1;
    dCK2 = kf*CK1 - kb*CK2 + betaa1*OK - alphaa1*CK2;
    dOK = alphaa1*CK2 - betaa1*OK + betaii*IK - alphaii*OK;
    dIK = alphaii*OK - betaii*IK;

    dY(41,1) = 0; dY(42,1) = dCK1; dY(43,1) = dCK2; dY(44,1) = dOK; dY(45,1) = dIK;
    
    nn = 0.98*Ko + 0.02*Nao;
    qq = 0.98*Ki + 0.02*Nai;

    IKr = OK * GKr * (V - ((R*T)/F)*log(nn/qq));

    %% * Na/K pump (INaK) ****************************
    INaKmax = 0.88; % pA/pF; Maximum exchange current
    KmNai = 21000;  % uM; Half saturation constant for Na/K exchanger
    KmKo = 1500;    % uM; Half saturation constant for Na/K exchanger
    
    sigma = 0.1429 * (exp(Nao/67300) - 1.0);
    fNaK = 1.0 / (1.0 + 0.1245*exp(-0.1*V*F/(R*T)) + 0.0365*sigma*exp(-V*F/(R*T)));

    INaK = INaKmax * (fNaK / (1.0 + (KmNai/Nai)^1.5)) * (Ko / (Ko + KmKo));

    %% * Ca-activated Cl- current (ICaCl) ************
    GClCa = 10.0; % mS/uF
    KmCl = 10.0;  % uM;
    ECl = -40.0;  % mV;
    OClCa = 0.2 / (1.0 + exp(-(V-46.7)/7.8));
    
    IClCa = GClCa * OClCa * (V - ECl) * (Cai / (Cai + KmCl));

    %% * Intracellular Ionic Dynamics ****************
    % Parameters
    LTRPNtot = 70.0;    % uM;
    HTRPNtot = 140.0;   % uM;
    CMDNtot = 50.0;     % uM
    CSQNtot = 15000.0;  % uM
    KmCMDN = 0.238;     % uM
    KmCSQN = 800.0;     % uM
    kp_htrpn = 0.00237; % 1/uM*ms
    km_htrpn = 3.2e-5;  % 1/ms
    kp_ltrpn = 0.0327;  % 1/uM*ms
    km_ltrpn = 0.0196;  % 1/ms

    v1 = 4.5;           % 1/ms
    v2 = 1.74e-5;       % 1/ms
    v3 = 0.45;          % uM/ms
    KmUp = 0.5;         % uM
    tautr = 20.0;       % ms
    tauxfer = 8.0;      % ms
    kp_a = 0.006075;    % RyR PC1 – PO1 rate constant
    km_a = 0.07125;     % RyR PO1 – PC1 rate constant
    kp_b = 0.00405;     % RyR PO1 – PO2 rate constant
    km_b = 0.965;       % RyR PO2 – PO1 rate constant
    kp_c = 0.009;       % RyR PO1 – PC2 rate constant
    km_c = 0.0008;      % RyR PC2 – PO1 rate constant
    n = 4;              % RyR Ca2+ cooperativity parameter PC1 – PO1
    m = 3;              % RyR Ca2+ cooperativity parameter PC1 – PO1

    % Calcium Buffers
    dLTRPNCa = kp_ltrpn * Cai * (LTRPNtot - LTRPNCa) - km_ltrpn * LTRPNCa;
    dHTRPNCa = kp_htrpn * Cai * (HTRPNtot - HTRPNCa) - km_htrpn * HTRPNCa;
    
    dY(6,1) = dLTRPNCa; dY(7,1) = dHTRPNCa;

    % Ryanodine Receptors
    PC1  = 1.0 - (PC2 + PO1 + PO2);
    dPO1 = kp_a*(Cass^n)*PC1 - km_a*PO1 - kp_b*(Cass^m)*PO1 + km_b*PO2 - kp_c*PO1 + km_c*PC2;
    dPO2 = kp_b*(Cass^m)*PO1 - km_b*PO2;
    dPC2 = kp_c*PO1 - km_c*PC2;

    dY(16, 1) = 0; dY(17,1) = dPC2; dY(18,1) = dPO1; dY(19,1) = dPO2; 

    % Calcium Fluxes
    Jrel  = v1 * (PO1 + PO2) * (Cajsr - Cass) * PRyR;
    Jtr   = (Cansr - Cajsr) / tautr;
    Jxfer = (Cass  - Cai) / tauxfer;
    Jleak = v2 * (Cansr - Cai);
    Jup   = v3 * ((Cai^2.0) / (KmUp^2 + Cai^2.0));
    Jtrpn = dLTRPNCa + dHTRPNCa; % kp_htrpn * Cai * (HTRPNtot - HTRPNCa) - km_htrpn * HTRPNCa + kp_ltrpn * Cai * (LTRPNtot - LTRPNCa) - km_ltrpn*LTRPNCa;

    % Calcium Concentration
    Bi =   (1.0 + CMDNtot*KmCMDN / (KmCMDN + Cai)^2.0)^-1.0;
    Bss =  (1.0 + CMDNtot*KmCMDN / (KmCMDN + Cass)^2.0)^-1.0;
    Bjsr = (1.0 + CSQNtot*KmCSQN / (KmCSQN + Cajsr)^2.0)^-1.0;
    
    % Formulation of dCass and dCai from Wang&Sobie Model
    % dCai = Bi * (Jleak + Jxfer - Jup - Jtrpn );
    % dCass = (Bss * (Jrel*(Vjsr/Vss) - Jxfer*(Vmyo/Vss) - (ICab + ICaL - 2*INaCa + IpCa)*(Acap*C/(2*Vmyo*F)))); 
    
    Cmyo = Acap*C/(2.0*Vmyo*F);
    Css = Acap*C/(2.0*Vss*F);

    dCai =   Bi * (Jleak + Jxfer - Jup - Jtrpn - (ICab - 2.0*INaCa + IpCa) * Cmyo);
    dCass =  Bss * (Jrel*Vjsr/Vss - Jxfer*Vmyo/Vss - ICaL * Css);
    dCajsr = Bjsr * (Jtr - Jrel);
    dCansr = (Jup - Jleak)*Vmyo/Vnsr - Jtr*Vjsr/Vnsr;
    
    dY(2, 1) = dCai; dY(3,1) = dCass; dY(4,1) = dCajsr; dY(5,1) = dCansr;

    % Sodium Concentration
    dNai = -(INa + INab + 3*INaCa + 3*INaK) * Acap * C/(Vmyo*F);
    dY(30,1) = dNai;

    % Potassium Concentration
    dKi = -(IKtof + IKtos + IK1 + IKs + IKss + IKur + IKr - 2*INaK) * Acap*C/(Vmyo*F);
    dY(31,1) = dKi;

    %% * Stimulus Current ***************************
    stim_amplitude = -80.0;   % picoA_per_picoF (in membrane)
    stim_duration = 0.5;      % millisecond (in membrane)
    stim_end = 100000.0;      % millisecond (in membrane)
    stim_period = 1000;       % millisecond (in membrane)
    stim_start = 20.0;        % millisecond (in membrane)

    if ((time >= stim_start) & (time <= stim_end) & (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
        Istim = stim_amplitude;
    else
        Istim = 0.0;
    end

    %% * Membrane voltage ****************************
    dY(1,1) = -(ICaL + IpCa + INaCa + ICab + INa + INab + INaK + IKtos + IKtof + IK1 + IKs + IKur + IKss + IKr + IClCa + Istim);
    
    out = [ICaL, IpCa, INaCa, ICab, INa, INab, INaK, IKtos, IKtof, IK1, IKs, IKur, IKss, IKr, IClCa, Cass, Cajsr, Cansr, Jleak, Jxfer, Jtr, Jrel, Jup, Jtrpn, Cai, Nai, Ki, Istim];

end