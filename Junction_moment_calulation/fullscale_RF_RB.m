clc
clear 
close all
%% Material data
E = 159e9;                       % Youngs Modulus in Pa

A1 = 1.0704178e-3;               % Cross section area of RF in m^2
I1 = 4.510529e-5;                % Moment of inertia for RF in m^4

A2 = 6.850674e-4;                % Cross section area of CB in m^2
I2 = 1.8475127e-5;               % Moment of inertia for CB in m^4

%% Element data
elements = 8;
node = elements + 1;
Coord=[0 0
       0.25 0
       0.3334 0
       0.5 0
       0.5 0.275
       0.5 1.1
       0.6666 0
       0.75 0
       1 0];
   
 Dof=[1 2 3
      4 5 6
      7 8 9
      10 11 12
      13 14 15
      16 17 18
      19 20 21
      22 23 24
      25 26 27];
  
 %% Stiffness matrix K and load vector f

K=zeros(node*3);
f=zeros(node*3,1);
f(5)=199.275;
f(8)=121.78;
f(13)=-500;
f(20)=121.78;
f(23)=199.275;

%% Element stiffness matrices and Assemble Ke into K

ep1=[E A1 I1]; ep2=[E A2 I2];
set1 = [1 2 3 6 7 8 ];set2=[4 5];
Edof = [1 1 2 3 4 5 6
        2 4 5 6 7 8 9
        3 7 8 9 10 11 12
        4 10 11 12 13 14 15
        5 13 14 15 16 17 18
        6 10 11 12 19 20 21
        7 19 20 21 22 23 24
        8 22 23 24 25 26 27]; %Elements for RF-RB

[Ex,Ey]=coordxtr(Edof,Coord,Dof,2);
 
for i=1:elements
    if ismember(Edof(i,1),set1)
        Ke1=beam2e(Ex(i,:),Ey(i,:),ep1);
    else
        Ke1=beam2e(Ex(i,:),Ey(i,:),ep2);
    end
    K=assem(Edof(i,:),K,Ke1);
end

%% Boundary condition to solve the system of equations and compute support forces
bc=[16 0;
    17 0;
    18 0];            
[a,r]=solveq(K,f,bc);

%% Section forces
Ed=extract(Edof,a);

%bottom junction
es3 = beam2s(Ex(3, :), Ey(3, :), ep1, Ed(3, :)); 
es4 = beam2s(Ex(4, :), Ey(4, :), ep2, Ed(4, :));
es6 = beam2s(Ex(6, :), Ey(6, :), ep1, Ed(6, :));

M_3 = es3(2, 3)
M_4 = es4(1, 3)
M_6 = es6(1, 3)

M_equilibrum_L = M_3-M_4-M_6

%top junction
es5 = beam2s(Ex(5, :), Ey(5, :), ep1, Ed(5, :));
M_5 = es5(2, 3)

[sfac]=scalfact2(Ex,Ey,Ed,0.1); 
eldisp2(Ex,Ey,Ed,[2 1 1],sfac);


