function [Phi_Phi,Phi_F,Phi_R,Phi,F,BarRs,A_e,B_e,C_e]=mpcgain_MIMO(Am,Bm,Cm,Nc,Np);
[q,n1]=size(Cm);
[n1,m]=size(Bm);
A_e=eye(n1+q,n1+q);
A_e(1:n1,1:n1)=Am;
A_e(n1+1:n1+q,1:n1)=Cm*Am; %Transformed Matrix in notes
B_e=zeros(n1+q,m);
B_e(1:n1,:)=Bm
B_e(n1+1:n1+q,:)=Cm*Bm;
C_e=zeros(q,n1+q);
C_e(:,n1+1:n1+q)=eye(q,q);
n=n1+q;
h(1:q,:)=C_e;
F=zeros(q*Np,n1+q);
F(1:q,:)=C_e*A_e; %MAtrix CA;CA^2...
for kk=2:Np
    h((kk-1)*q+1:kk*q,:)=h((kk-2)*q+1:(kk-1)*q,:)*A_e;
    F((kk-1)*q+1:kk*q,:)=F((kk-2)*q+1:(kk-1)*q,:)*A_e;
end
v=h*B_e;
Phi=zeros(q*Np,m*Nc); %declare the dimension of phi
Phi(:,1:m)=v %first coloumn of phi
for i=2:Nc
    Phi(:,(i-1)*m+1:i*m)=[zeros((i-1)*(q),m);v(1:q*(Np-i+1),1:m)]; %Toeplitz Matrix
end
BarRs=ones(q*Np,1);
Phi_Phi=Phi'*Phi;
Phi_F=Phi'*F;
Phi_R=Phi'*BarRs;


