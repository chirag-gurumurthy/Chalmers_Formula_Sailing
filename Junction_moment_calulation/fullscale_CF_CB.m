clc
clear 
close all
%% Material data
E = 159e9;                       % Youngs Modulus in Pa

A1 = 4.3e-3;                     % Cross section area of CF in m^2
I1 = 1.03514e-3;                 % Moment of inertia for CF in m^4

A2 = 1.5414e-3;                  % Cross section area of CB in m^2
I2 = 9.35303e-5;                 % Moment of inertia for CB in m^4
%% Element data
elements = 13;
node = elements + 1;
Coord=[0 0
    0.35 0
    0.467 0
    0.7 0
    0.7 0.22
    0.7 1
    0.875 0
    1.225 0
    1.4 0
    1.4 0.28
    1.4 1
    1.633 0
    1.75 0
    2.1 0];

Dof=[1 2 3
    4 5 6
    7 8 9
    10 11 12
    13 14 15
    16 17 18
    19 20 21
    22 23 24
    25 26 27
    28 29 30
    31 32 33
    34 35 36
    37 38 39
    40 41 42];
%% Stiffness matrix K and load vector f

K=zeros(node*3);
f=zeros(node*3,1);
f(5)=430;
f(8)=262.2;
f(13)=-550;
f(20)=477;
f(23)= 1062;
f(28)=-700;
f(35)= 585;
f(38)=956.2;

%% Element stiffness matrices and Assemble Ke into K

ep1=[E A1 I1]; ep2=[E A2 I2];
set1 = [1 2 3 6 7 8 11 12 13];set2=[4 5 9 10];
Edof = [1 1 2 3 4 5 6
        2 4 5 6 7 8 9
        3 7 8 9 10 11 12
        4 10 11 12 13 14 15
        5 13 14 15 16 17 18
        6 10 11 12 19 20 21
        7 19 20 21 22 23 24
        8 22 23 24 25 26 27
        9 25 26 27 28 29 30
        10 28 29 30 31 32 33
        11 25 26 27 34 35 36
        12 34 35 36 37 38 39
        13 37 38 39 40 41 42];   %Elements for CF

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
    18 0;
    31 0;
    32 0;
    33 0;
    10 0;
    11 0;
    25 0;
    26 0];            
[a,r]=solveq(K,f,bc);

%% Section forces
Ed=extract(Edof,a);

% Changed so that the coordinates and ep are correct for each element
%left junction
es3 = beam2s(Ex(3, :), Ey(3, :), ep1, Ed(3, :)); 
es4 = beam2s(Ex(4, :), Ey(4, :), ep2, Ed(4, :));
es6 = beam2s(Ex(6, :), Ey(6, :), ep1, Ed(6, :));
es5 = beam2s(Ex(5, :), Ey(5, :), ep2, Ed(5, :));

M_3 = es3(2, 3)
M_4 = es4(1, 3)
M_6 = es6(1, 3)
M_5 = es5(2, 3)

% The second moment in element 3 should equal the first moment in element 4 PLUS 6
% (Equilibrium for the node due to different direction of the sectional moments)  

M_equilibrum_L = M_3-M_4-M_6 % Should be 0

%right junction
es8 = beam2s(Ex(8, :), Ey(8, :), ep1, Ed(8, :)); 
es9 = beam2s(Ex(9, :), Ey(9, :), ep2, Ed(9, :));
es11 = beam2s(Ex(11, :), Ey(11, :), ep1, Ed(11, :));
es10 = beam2s(Ex(10, :), Ey(10, :), ep2, Ed(10, :));

M_8 = es8(2, 3)
M_9 = es9(1, 3)
M_11 = es11(1, 3)
M_10 = es10(2, 3)

M_equilibrum_R = M_8-M_9-M_11 % Should be 0

[sfac]=scalfact2(Ex,Ey,Ed,0.1); 
eldisp2(Ex,Ey,Ed,[2 1 1],sfac); 
