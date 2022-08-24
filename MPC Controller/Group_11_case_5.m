% MPC Control
% Four Tank System
close all;
clear all; clc;

% Define the System Matrices
 % Given parameters
 A1=28; A2=32; A3=28; A4=32;             
 a1=0.071; a2=0.057; a3=0.071;a4=0.057;
 kc=1.0;
 g=981;
 k1=3.33; k2=3.35;
 z1=0.70; z2=0.60;
 T1=62; T2=90; T3=23; T4=30;
 
 % Matrices of the system
 Am1=[-1/T1 0 A3/(T3*A1) 0;0 -1/T2 0 A4/(A2*T4);0 0 -1/T3 0;0 0 0 -1/T4];   
 Bm1=[z1*k1/A1 0;0 z2*k2/A2 ;0 (1-z2)*k2/A3 ;(1-z1)*k1/A4 0];
 Cm1=[0 kc 0 0 ;0 0 kc 0]; %Change according to the controller (To control)
 Co1=[kc 0 0 0 ;0 0 0 kc]; %Change according to the observer (To measure)
 Dm1=[0 0;0 0];

 % Continuous State space model
 sys_cp=ss(Am1,Bm1,Cm1,Dm1); % (For MPC)
 sysK=ss(Am1,Bm1,Co1,Dm1);   % (For Kalman Filter)

 %Step response of system
 step(sys_cp);

 % Discrete State space model
 Ts=0.1;
 sys_d=c2d(sys_cp,Ts);       % (For MPC)
 sysK_d=c2d(sysK,Ts);        % (For Kalman Filter)
 
 Am=sys_d.A;        
 Bm=sys_d.B;
 Cm=sys_d.C;  %(To control)
 Co=sysK_d.C; %(To measure)
 Dm=sys_d.D;

%For MPC implementation

 %Prediction Horizon
 Np=20;
 %Control Horizon
 Nc=5;
 %Calculate the matrices
 [Phi_Phi,Phi_F,Phi_R,Phi,F,BarRs,A_e,B_e,C_e]=mpcgain_MIMO(Am,Bm,Cm,Nc,Np);
 q=size(Cm,1);
 [n,m]=size(B_e);
 nm=size(Am,1);
 %xm=zeros(nm,1);
 
 Xf=zeros(n,1); %Augmented state system

 % Define the length of the simulation
 nlen=800;

 r=ones(nlen,q);
 u=zeros(m,1); %u(k-1)=0
 y=zeros(q,1);
 uwt=[1 1]; %input weight R matrix
 ysp=[13.7 2.8];% Set point
 Rs=repmat(ysp',Np,1);
 R=zeros(m*Nc,m*Nc);
 R(1:m,1:m)=diag(uwt);
 for i=2:Nc
    R((i-1)*m+1:i*m,(i-1)*m+1:i*m)=diag(uwt);
 end

 % Hessian Matrix
 H=Phi_Phi+R;

 %Input Constraints
 DUmin=5*[-1;-1];  
 DUmax=5*[1;1];
 Umin= 0*[-1;-1];
 Umax= 20*[1;1];

%For Kalman Filter Implementation
  Am_tr = transpose(Am);
  Co_tr = transpose(Co);
  %Checking for Controllability and Observability
  %Since this is fourth order system, Controllability matrix is [B AB A^2*B  A^3*B]

  CO=[Bm(:,1) Bm(:,2) Am*Bm(:,1) Am*Bm(:,2) Am^2*Bm(:,1) Am^2*Bm(:,2) Am^3*Bm(:,1) Am^3*Bm(:,2)];
  R_co = rank(CO);
  if R_co==4
   disp('System is controllable')
  else
   disp('System is uncontrollable')
  end

  %Similarly Observability matrix is [H;HA;HA^2;HA^3]
  OB = [Co(1,:); Co(2,:); Co(1,:)*Am; Co(2,:)*Am; Co(1,:)*Am^2; Co(2,:)*Am^2; Co(1,:)*Am^3;Co(2,:)*Am^3];
  R_ob = rank(OB);
  if R_ob ==4
   disp('System is observable')
  else
   disp('System is unobservable')
  end
  
  %Define the noise covariances.
  %We can tune the Q and R for the better estimates of the states
  Q = 20.02*eye(4);
  R = 1*eye(2);
  
  %Preallocating the memory
  x_pri=zeros(4,nlen);
  x_post=zeros(4,nlen);
  E=zeros(4,nlen);
  P_pri=ones(4,4,nlen);
  P_post=ones(4,4,nlen);
  K=zeros(4,4,nlen);
  %Tolerance=zeros(1,nlen);
  trpapriori = zeros(1,nlen);
  trposteriori=zeros(1,nlen);

  %Initial condition on the state, x.
  xm=[12.4;12.7;1.8;1.4];

  %Initial guesses for state and a posteriori covariance.
  x_postm= [1;1;1;1];
  P_postm= (10^5)*eye(4);
  
  %Calculate the process and measurement noise.
  w= 0.005*rand(4,1);
  v= 0.05*rand(2,1);
  
  deltau1=zeros(2,nlen);
  u1=zeros(2,nlen);
  xm1=zeros(4,nlen);
  y1=zeros(2,nlen);
  x_postm1=zeros(4,nlen);
  x_pri1=zeros(4,nlen);

  ukprev=u;
  for kk=1:nlen;
     %pred_err=(Rs-F*Xf);
     %DeltaU=inv(Phi_Phi+R)*Phi'*pred_err;
     f=-(Phi')*(Rs-F*Xf); 
     [M,gamma]=mpc_constraint_MIMO(Umin,Umax,DUmin,DUmax,ukprev,Nc); 
     DeltaU=quadprog(H,f,M,gamma); %H is Haessian matrix
     deltau=DeltaU(1:m,1); %First move of m inputs
     deltau1(:,kk)=deltau;

     u=u+deltau;
     ukprev=u;
     u1(:,kk)=u;

     xm_old=xm;
     xm=Am*xm+Bm*u+w; %delxm(k+1)
     xm1(:,kk)=xm;

     y=Co*xm+v; %y(k) and random noise
     y_mpc=Cm*xm;
     y1(:,kk)=y_mpc;
     
     
     
     %Kalman Filter implementation
      %Predictor equations
        %Compute X_prior
       x_pri=Am*x_postm + Bm*u;
       x_pri1(:,kk)=x_pri;
       %Compute variance P prior
       P_pri=Am*P_postm*Am_tr+Q;
      %Calculating the trace of P prior
      %trpapriori(1,k) = trace(P_pri(:,:,k));

 %Corrector equations
 %Compute Kalman Gain
 K=P_pri*Co_tr/(Co*P_pri*Co_tr+R);
 
 Z_est = Co*x_pri;
 Z_true=y;
 E=Z_true-Z_est;

 %Compute X posterior
 x_postm=x_pri+K*E;
 x_postm1(:,kk)=x_postm;
 %Compute variance P posterior
 P_postm=P_pri - K*Co*P_pri;
 %Calculating the trace of P posterior
 %trposteriori(1,k) = trace(P_post(:,:,k));
        

  Xf=[xm-xm_old;y_mpc]; %full state feedback in delta
    end
    k=0:Ts:(nlen-1)*Ts;
    figure
                         %Plots all the tank measurements with legend
    for i = 1:4
        plot(xm1(i,:),'displayname',strcat('h',num2str(i)));
        hold on;
    end
    xlabel('Time')
    ylabel('All the water levels')
    legend;
  figure
    subplot(311)                          %Plots the controlled heights only.
        plot(transpose(y1(1,:)),'displayname',strcat('h2'));
        hold on;
        plot(transpose(y1(2,:)),'displayname',strcat('h3'));
    xlabel('Time')
    ylabel('Water level controlled')
    legend;
    subplot(312)
    plot(k,u1')    %Plots the pump voltages
    xlabel('Time')
    ylabel('pump volatge')
    legend('Control')
    subplot(313)
    stairs(k,deltau1')
    xlabel('Time')
    ylabel('deltau')
    legend('First incremental control')
    
  disp('set points are=');
  disp(ysp)
  X=['when Np=',num2str(Np), ' Nc=',num2str(Nc)]; 
  disp(X);
  disp('and uwt=')
  disp(uwt);
  disp('Set points acheived are=')
  disp(y_mpc);
                      
  
   
   
    
 

