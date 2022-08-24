%Compute the constarint matrix
%function[M,gamma]=mpc_constraint_MIMO(phi,F,Umin,Umax,DUmin,DUmax,Ymin,Ymax,Ukprev,Nc,Np)
function [M,gamma]=mpc_constraint_MIMO(Umin,Umax,DUmin,DUmax,Ukprev,Nc)
m=size(Ukprev,1); %inputs
%q=size(Ymin,1);%outputs

%Max u
C1=repmat(Ukprev,Nc,1);
C2i=repmat(eye(m,m),Nc,1);
C2=zeros(m*Nc,m*Nc);
C2(:,1:m)=C2i;
UMIN=repmat(Umin,Nc,1);
UMAX=repmat(Umax,Nc,1);
for i=2:Nc
    C2((i-1)*m+1:Nc*m,(i-1)*m+1:i*m)=C2i((i-1)*m+1:Nc*m,1:m);
end
M1=[-C2;C2];
N1=[-UMIN+C1;UMAX-C1]; %Input amplitude

%Max rate of change of deltaU constraint
deltaMin=repmat(DUmin,Nc,1);
deltaMax=repmat(DUmax,Nc,1);

M2=[-eye(m*Nc,m*Nc);eye(m*Nc,m*Nc)];
N2=[-deltaMin;deltaMax]; %input rate of change

M=[M1;M2];
gamma=[N1;N2];


% %We need M3 also
% M3=[-Phi;Phi];
% N3=
    