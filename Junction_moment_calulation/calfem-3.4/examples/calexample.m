clear all
close all
clc

E = 2.10e11; % Pa
A = 0.00537673773; % mm^2 Cross section area 
I = 2.23193131e-05; % m^4 Moment of inertia 

elements = 2100;  
node = elements + 1;
fixed = [700, 1400];


  Edof = zeros(elements,7);
  Edof(:,1) = 1:elements;
  n = 1;
  for i = 1:elements
      Edof(i,2:end) = [n:n+5];
      n = n + 3;
  end

K=zeros(node*3); f=zeros(node*3,1);

f(1227)=864.881500;
f(1635)=475.685000;
f(1752)=778.393000;

ep=[E A I];

ex=[0 elements]; ey=[0 0];

Ke=beam2e(ex,ey,ep);        % element stiffness matrix

K=assem(Edof,K,Ke);

bc=[(fixed(1,1)*3 + 2) 0; (fixed(1,1)*3 + 1) 0; (fixed(1,1)*3) 0;...
    (fixed(1,2)*3 + 2) 0; (fixed(1,2)*3 + 1) 0; (fixed(1,2)*3) 0];

[a,r]=solveq(K,f,bc);

Ed=extract(Edof,a);

es1=beam2s(ex,ey,ep,Ed(700,:));
es2=beam2s(ex,ey,ep,Ed(1400,:));
es3=beam2s(ex,ey,ep,Ed(1050,:));




