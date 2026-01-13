function [M_A,M_B,M_C] = chem_exchange_sim(flip_angle,Tacq,pop,nu,kex,R1,R2)

% All quantities are in Hz and seconds

% Unpacking offset, exchange and population arrays
nu_A0 = nu(1);nu_B0 = nu(2);nu_C0 = nu(3);
kex_AB = kex(1);kex_BC = kex(2);kex_AC = kex(3);
MA0 = pop(1);MB0 = pop(2);MC0 = 1-MA0-MB0;

% Exchange and recovery rates
kBA=kex_AB*MB0; % exchange from A to B
kAB=kex_AB*MA0; % exchange from B to A

kBC=kex_BC*MB0; % exchange from B to C
kCB=kex_BC*MC0; % exchange from C to B

kAC=kex_AC*MA0; % exchange from A to C
kCA=kex_AC*MC0; % exchange from C to A

% exchange matrix
K(1,1:10)=[0 0 0 0 0 0 0 0 0 0];
K(2,1:10)=[0 -kBA-kCA 0 0 kAB 0 0 kAC 0 0];
K(3,1:10)=[0 0 -kBA-kCA 0 0 kAB 0 0 kAC 0];
K(4,1:10)=[0 0 0 -kBA-kCA 0 0 kAB 0 0 kAC];
K(5,1:10)=[0 kBA 0 0 -kAB-kCB 0 0 kBC 0 0];
K(6,1:10)=[0 0 kBA 0 0 -kAB-kCB 0 0 kBC 0];
K(7,1:10)=[0 0 0 kBA 0 0 -kAB-kCB 0 0 kBC];
K(8,1:10)=[0 kCA 0 0 kCB 0 0 -kAC-kBC 0 0];
K(9,1:10)=[0 0 kCA 0 0 kCB 0 0 -kAC-kBC 0];
K(10,1:10)=[0 0 0 kCA 0 0 kCB 0 0 -kAC-kBC];

R1A=R1;R1B=R1; R1C=R1;
R2A=R2;R2B=R2;R2C=R2;


nu1=100000; % RF amplitude
nu1=-nu1;
omega1=2*pi*(nu1);T90=1/abs(nu1)/4;
Tpulse=T90*flip_angle/90; % pulse flip angle and duration

phase=0;
omega1X=omega1*cos(phase);omega1Y=omega1*sin(phase);

% RF matrix
Rrf(1,1:10)=[0 0 0 0 0 0 0 0 0 0];
Rrf(2,1:10)=[0 0 0 +omega1Y 0 0 0 0 0 0];
Rrf(3,1:10)=[0 0 0 -omega1X 0 0 0 0 0 0];
Rrf(4,1:10)=[0 -omega1Y +omega1X 0 0 0 0 0 0 0];
Rrf(5,1:10)=[0 0 0 0 0 0 +omega1Y 0 0 0];
Rrf(6,1:10)=[0 0 0 0 0 0 -omega1X 0 0 0];
Rrf(7,1:10)=[0 0 0 0 -omega1Y  omega1X 0 0 0 0];
Rrf(8,1:10)=[0 0 0 0 0 0 0 0 0 +omega1Y];
Rrf(9,1:10)=[0 0 0 0 0 0 0 0 0 -omega1X];
Rrf(10,1:10)=[0 0 0 0 0 0 0 -omega1Y omega1X 0];

% RF propagator
U_rf=expm(Rrf*Tpulse);

% pre allocating magnetization arrays
MXA_null = zeros(size(Tacq));MYA_null = zeros(size(Tacq)); MZA_null = zeros(size(Tacq));
MXB_null = zeros(size(Tacq));MYB_null = zeros(size(Tacq)); MZB_null = zeros(size(Tacq));
MXC_null = zeros(size(Tacq));MYC_null = zeros(size(Tacq)); MZC_null = zeros(size(Tacq));

I10 = eye(10);

for q=1:numel(Tacq)
    TR=Tacq(q);

    % if you work with arbitrary offsets
    nuA=-1/2/TR+nu_A0;%the -1/2/Tacq addition is to mimim the x -x phases you have in experiment
    nuB=-1/2/TR+nu_B0;
    nuC=-1/2/TR+nu_C0;
    %%%%%%%%%%% maybe later we can find a way to do the null without adding -1/2/Tacq

    nuA=-nuA;nuB=-nuB;nuC=-nuC;

    %%% here all sites have the same relaxation rates R1 and R2
    %%% however you can set them differently with this numerical approach

    omegaA=2*pi*nuA; % from frequency to angular frequency
    omegaB=2*pi*nuB;
    omegaC=2*pi*nuC;

    % relaxation and chemical shifts Liouville matrix
    R0(1,1:10)=-[0 0 0 0 0 0 0 0 0 0];
    R0(2,1:10)=-[0 R2A +omegaA 0 0 0 0 0 0 0];
    R0(3,1:10)=-[0 -omegaA R2A 0 0 0 0 0 0 0];
    R0(4,1:10)=-[-2*R1A*MA0 0 0 R1A 0 0 0 0 0 0];
    R0(5,1:10)=-[0 0 0 0 R2B +omegaB 0 0 0 0];
    R0(6,1:10)=-[0 0 0 0 -omegaB R2B 0 0 0 0];
    R0(7,1:10)=-[-2*R1B*MB0 0 0 0 0 0 R1B 0 0 0];
    R0(8,1:10)=-[0 0 0 0 0 0 0 R2C +omegaC 0];
    R0(9,1:10)=-[0 0 0 0 0 0 0 -omegaC R2C 0];
    R0(10,1:10)=-[-2*R1C*MC0 0 0 0 0 0 0 0 0 R1C];

    R0=R0+K; % Full Liouville matrix for free evolution (no pulses)
    Ufree=expm(R0*TR); % propagator for free evolution

    %%%%%%%%%%%%  finding the steady state
    U=U_rf*Ufree;
    A_sys = U-I10;
    A_sys(1,:) = [1,zeros(1,9)];
    b_sys = [0.5;zeros(9,1)];
    V = A_sys\b_sys;
    % V=null(U-eye(10,10));V=V/V(1)/2;
    
    
    MXA_null(q)=V(2);MYA_null(q)=V(3);MZA_null(q)=V(4);
    MXB_null(q)=V(5);MYB_null(q)=V(6);MZB_null(q)=V(7);
    MXC_null(q)=V(8);MYC_null(q)=V(9);MZC_null(q)=V(10);
end

M_A = abs(MXA_null+1i*MYA_null);
M_B = abs(MXB_null+1i*MYB_null);
M_C = abs(MXC_null+1i*MYC_null);


