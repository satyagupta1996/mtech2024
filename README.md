clc
clear all
syms d l g w d2 den e
a=pi*d^2/4;
I=pi*d^4/64;
disp('............FINITE ELEMENT METHOD.........')


%element 1
m1=(den*a*l/420)*[156 22*l 54 -12*l;22*l 4*l^2 13*l -3*l^2;54 13*l 156 -22*l;-12*l -3*l^2 -22*l 4*l^2];
m1=subs(m1,{d,l,den,e},{0.05,0.5,7800,2.1*10^11});
k1=(e*I/l^3)*[12 6*l -12 6*l;6*l 4*l^2 -6*l 2*l^2;-12 -6*l 12 -6*l;6*l 2*l^2 -6*l 4*l^2];
k1=subs(k1,{d,l,den,e},{0.05,0.5,7800,2.1*10^11});
k1=double(k1);
m1=double(m1);


%element 2
m2=(den*a*l/420)*[156 22*l 54 -12*l;22*l 4*l^2 13*l -3*l^2;54 13*l 156 -22*l;-12*l -3*l^2 -22*l 4*l^2];
m2=subs(m2,{d,l,den,e},{0.05,0.5,7800,2.1*10^11});
k2=(e*I/l^3)*[12 6*l -12 6*l;6*l 4*l^2 -6*l 2*l^2;-12 -6*l 12 -6*l;6*l 2*l^2 -6*l 4*l^2];
k2=subs(k2,{d,l,den,e},{0.05,0.5,7800,2.1*10^11});
k2=double(k2);
m2=double(m2);

%global Mass matrix
global_Mass_Matrix=zeros(6,6);
global_Mass_Matrix(1:4,1:4)=m1;
global_Mass_Matrix(3:6,3:6)=global_Mass_Matrix(3:6,3:6)+m2

%global stiffness matrix
global_Stifness_matrix=zeros(6,6);
global_Stifness_matrix(1:4,1:4)=k1;
global_Stifness_matrix(3:6,3:6)=global_Stifness_matrix(3:6,3:6)+k2

%reduced mass and stifness mastrix
global_Mass_Matrix=global_Mass_Matrix(3:4,3:4);
global_Stifness_matrix=global_Stifness_matrix(3:4,3:4);

%boundary condition
global_Mass_Matrix=(-w^2*global_Mass_Matrix);
A2=(global_Mass_Matrix+global_Stifness_matrix);

%solving for frequency
A2=det(A2);
frequency=solve(A2,w);
frequency=vpa(frequency([2,4]))

disp('------------------------------------------')



