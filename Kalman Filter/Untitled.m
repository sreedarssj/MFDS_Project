clc
clear all
syms x1 x2;
f(x1,x2)=100*(x2-x1^2)^2+(1-x1)^2;
P(x1,x2)=[diff(f,x1);diff(f,x2)];
initialguess=[2.2;1.9];
LearningParameter=1;
% hessian(f,[x1,x2])
z_store=[];
A_store = [];
for i=1:2
    double(initialguess);
    A=initialguess-LearningParameter*P(initialguess(1),initialguess(2));
    initialguess=A;
    z=double(f(A(1,1),A(2,1)));
    z_store=[z_store;z]; 
    A_store = [A_store double(A)];
end