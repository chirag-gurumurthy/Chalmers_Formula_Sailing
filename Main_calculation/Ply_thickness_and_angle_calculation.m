clc
clear all
close all

%% Material data
% Design
theta= [90 45 0 0 0 0 45 90];         % lamina orientation (write out all)
sym= 1;                 % sym=1 if symetric laminate, sym=2 if not symetric laminate
E_L= 203e9;             % Youngs modulus (Pa) in longitudinal direction
E_T= 11.2e9;            % Youngs modulus (Pa) in transverse direction
nu_LT= 0.32;            % Major poisson's ratio
G_LT= 8.4e9;            % Shear modulus (Pa)
sigma_LU= 3500e6;       % Tensile ultimate strength (Pa) of fiber in longitudinal direction
sigma_LUp= 1540e6;      % Compressive ultimate strength (Pa) of fiber in longitudinal direction
sigma_TU= 56e6;         % Tensile ultimate strength (Pa) of fiber in transverse direction
sigma_TUp= 150e6;       % Compressive ultimate strength (Pa) of fiber in longitudinal direction
tau_LTU= 98e6;          % Shear strength (Pa) of fiber

%% Force and laminate thickness
th      = 3e-3 ;      % define hight of laminate [unit= m]
N       = [0 0 0]';   % Added mechanical force [unit= N]
M       = [0 58.1634 0]'; % Added mechanical moment [unit= N-m]

%% Stiffness matrix
[Q] = laminadata(E_L, E_T, nu_LT, G_LT);

%% create h vector 
h=zeros(1,length(theta)+1);
m=1;
for n=th/2:-th/length(theta):-th/2
    h(m)=n;
    m=m+1;
end

%% setup zero matrices 
A          =   zeros(size(Q));
B          =   A;
D          =   A;

%% Solve extensional stiffness matrix, coupling stiffness matrix, and bending stiffness matrix
for i=1:length(theta)
    % = create Q_bar =
    [ T1 , T2]  =   CMTd(theta(i));
    Q_bar       =   T1\Q*T2;
    % = create A, B, D matrices =
    A           =   A + Q_bar*(h(i)-h(i+1));
    if sym==2
        B       =   B + 1/2*Q_bar*(h(i)^2-h(i+1)^2);
    end
    D           =   D + 1/3*Q_bar*(h(i)^3-h(i+1)^3);
end

%% Stress and strain calculation 

[ep0] = A\N;     % mid plane strain
[k] = D\M;       % plate curvature

p=0;    
for i=1:length(theta)
    [ T1 , T2]  =   CMTd(theta(i));
    Q_bar       =   T1\Q*T2;
        
    for o=0:1
        p=p+1;
        z(p)=h(i+o);
        ep          =   ep0 + z(p)*k;
           
    % = compute stress =
        sigma(:,p)  =   Q_bar*ep;
    end
end

p=0;
for i=1:length(theta)
    [T1] = CMTd(theta(i));
    for j=1:2
        p=p+1;
        sigma_LT(:,p)=T1*sigma(:,p);
    end
end

odd = 1:2:length(sigma_LT(1,:));
even = 2:2:length(sigma_LT(1,:));
for jj = 1:length(sigma_LT(:,1))
    
    for kk = 1:length(even)
        if (sigma_LT(jj,odd(kk))<0 | sigma_LT(jj,even(kk))<0)
            if abs(sigma_LT(jj,odd(kk)))> abs(sigma_LT(jj,even(kk))) 
                if sigma_LT(jj,odd(kk))> 0
                    sigma_local(jj,kk) = max(abs(sigma_LT(jj,odd(kk))),abs(sigma_LT(jj,even(kk))));
                else
                    sigma_local(jj,kk) = -max(abs(sigma_LT(jj,odd(kk))),abs(sigma_LT(jj,even(kk))));
                end
            else
                if sigma_LT(jj,odd(kk))> 0
                    sigma_local(jj,kk) = max(abs(sigma_LT(jj,odd(kk))),abs(sigma_LT(jj,even(kk))));
                else
                    sigma_local(jj,kk) = -max(abs(sigma_LT(jj,odd(kk))),abs(sigma_LT(jj,even(kk))));
                end

            end
        else
            sigma_local(jj,kk) =  max(abs(sigma_LT(jj,odd(kk))),abs(sigma_LT(jj,even(kk))));
        end
    end
end

%% Failure criteria 

for q=1:length(sigma_local(1,:))
    sigL=sigma_local(1,:);
    sigT=sigma_local(2,:);
    tauLT=sigma_local(3,:);
end 

for r=1:length(sigma_local(1,:))
    
    if sigL(r)<0
    aa=(abs(sigL(r))/sigma_LUp)^2;
    else aa=(sigL(r)/sigma_LU)^2;
    end 
    
    if sigL(r)<0
    bb=(abs(sigL(r))/sigma_LUp);
    else bb=(sigL(r)/sigma_LU);
    end 
    
    if sigT(r)<0
    cc= (abs(sigT(r))/sigma_LUp);
    else  cc= (sigT(r)/sigma_LU);
    end 
        
     if sigT(r)<0
    dd= (abs(sigT(r))/sigma_TUp)^2;
    else  dd= (sigT(r)/sigma_TU)^2;
     end
     
    if (aa - (bb*cc) + dd + ((abs(tauLT(r))/tau_LTU)^2) >1)
        fprintf('Ply %d failed',r);
       fprintf('\n')
    end
end 

