% A=[9 -2;-2 6]
% eig(A)
% [V,D] = eig(A)
% 
% H=[9 -2;-2 6]
% -inv(H)*[11;-8]
% 
% 
% ’%Optimization problem
% clear
% clc
% a=280; b=380;
% Di=[]; Dj=[];
% res=zeros(1,100);
% x1=a+(b-a)/3;
% x2=a+2*(b-a)/3;
% while true
% if f(x1)<f(x2)
% b=x2;
% else
% a=x1;
% end
% x1=a+(b-a)/3;
% x2=a+2*(b-a)/3;
% Dii=a; Djj=b;
% Table=[]
% Di=[Di,Dii];
% Dj=[Dj,Djj];
% if abs(a-b)<=0.001
% break
% end
% end
% matrix=zeros(100,1);
% for i=1:100
% res(i)=f(i+280);
% end
% plot(281:380,res);
% grid on;
% hold on;
% [~,k]=size(Di);
% temp=zeros(1,k);
% for i=1:k
% temp(i)=f(Di(i));
% end
% [minimum,argmin]=min(temp);
% plot(Di(argmin),minimum,'o','MarkerSize',2);
% xlabel('D1(mm)'); ylabel('Cost');
% Table=[Di',Dj'];
% function y=f(D1)
% L1= 300; L2=500 ; L3=400;
% Q1= 9 ; Q2= 3 ;Q3=2;
% HB= 75; HC= 89 ;HD=83;
% Ho=100;
% Opt = 4.457*10^8*L1*Q1^1.85/(D1^4.87);
% y = 1.2654*(L1*D1^1.327 + L2*(4.45*10^8*L2*Q2^1.85/(Ho-Opt-HC))^0.27+
% L3*(4.457*10^8*Q3^1.85/(Ho-Opt-HD))^0.27);
% end