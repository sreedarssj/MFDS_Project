clc
clear
close all

A1=28; A2=32; A3=28; A4=32;              %Given parametrs from lecture
a1=0.071; a2=0.057; a3=0.071;a4=0.057;
kc=0.50;
g=981;
                                         %and

h10=12.4; h20=12.7; h30=1.8; h40=1.4;   %Given parameters from the Reference paper
v10=3.00; v20=3.00;
k1=3.33; k2=3.35;
y1=0.70; y2=0.60;
T1=62; T2=90; T3=23; T4=30;


A=[-1/T1 0 A3/(T3*A1) 0;0 -1/T2 0 A4/(A2*T4);0 0 -1/T3 0;0 0 0 -1/T4];   %Matrix A B C D
B=[y1*k1/A1 0;0 y2*k2/A2 ;0 (1-y2)*k2/A3 ;(1-y1)*k1/A4 0];
H=[kc 0 0 0 ;0 kc 0 0];
D=0;

sys=ss(A,B,H,D);
sysd=c2d(sys,0.1); %convert from continous time to discrete time


Q=10*eye(4);            %Noise matrix
R=2*eye(2); 
P_postk=10^5*eye(4);
Uk=0;

hValues = xlsread('Measurements.xlsx');  %Reading from the Excel Sheet

X_postk=[1;1;1;1];
% Z_truek=hValues(1,1:2)-[h10;h20];
m= size(hValues);


X_post_store =[];    %Introducing matrix to store the values
X_pri_store=[];
P_post_store=[];
P_pri_store=[];
K_store=[];
Innov_store=[];
Resid_store=[];
Trace_P_pri_store=[];
Trace_P_post_store=[];


  for i = 1:150                 %taking the first 150 iterations
    X_prik=sysd.A*X_postk;       %Prediction equations         
    P_prik=sysd.A*P_postk*transpose(sysd.A)+Q;
    Trace_P_prik=sum(diag(P_prik));
    
    Kk=P_prik*transpose(sysd.C)*inv(sysd.C*P_prik*transpose(sysd.C)+R);             %Update equations
    Z_estk=sysd.C*X_prik;
    Z_truek = hValues(i,1:2)'-[h10;h20];
    Ek=Z_truek-Z_estk;
    X_postk=X_prik+Kk*Ek;
    P_postk=P_prik-Kk*sysd.C*P_prik;
    Trace_P_postk=sum(diag(P_postk));
    resid=Z_truek-H*X_postk;
   
    X_post_store = [X_post_store , X_postk];
    X_pri_store=[X_pri_store,X_prik];
    P_post_store=[P_post_store,P_postk];
    Trace_P_post_store=[Trace_P_post_store,Trace_P_postk];
    
    P_pri_store=[P_pri_store,P_prik];
    Trace_P_pri_store=[Trace_P_pri_store,Trace_P_prik];
    K_store=[K_store,Kk];
    Innov_store=[Innov_store,Ek];
    Resid_store=[Resid_store,resid];
    
 end

 figure                  %Plotting the 8 graphs
 plot(1:150,X_post_store)
 xlabel('x(first 150 iterations)')
 ylabel('X postk')
 grid on
 title('Variation of X postk with iterations')
 
 figure
 plot(1:150,X_pri_store)
  xlabel('x(first 150 iterations)')
 ylabel('X prik')
 grid on
 title('Variation of X prik with iterations')
 
 figure
 plot(2:150, Trace_P_pri_store(2:150),'.r')
 xlabel('x(first 150 iterations)')
 ylabel('P prik')
 legend('P prik')
 grid on
 title('Variation of P prik with iterations')
 
 figure
 plot(2:150, Trace_P_post_store(2:150),'.c')
 xlabel('x(first 150 iterations)')
 ylabel('P postk')
 legend('P postk')
 grid on
 title('Variation of P postk with iterations')
 
 figure
 hold on
 plot(2:150, Trace_P_post_store(2:150),'.c')
 plot(2:150, Trace_P_pri_store(2:150),'.r')
 xlabel('x(first 150 iterations)')
 ylabel('P postk and P prik')
 title('Variation of P postk and P prik with iterations')
 legend('P postk','P prik')
 grid on
  
 figure
plot(1:300,K_store)
 xlabel('x(first 150 iterations)')
 ylabel('Kk')
 grid on
 title('Variation of Kalman gain with iterations')
 
 figure
 plot(1:150,Innov_store)
  xlabel('x(first 150 iterations)')
 ylabel('Innovation')
 grid on
 title('Variation of Innovation with iterations')
 
 figure
 plot(1:150,Resid_store)
  xlabel('x(first 150 iterations)')
 ylabel('Residual')
 grid on
 title('Variation of Residual with iterations')
 
 
 
 
 