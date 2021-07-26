function [Edof,Ex,Ey,B1,B2,B3,B4]=quad_mesh(p1,p2,nelx,nely,ndofs)
%        [Edof,Ex,Ey,B1,B2,B3,B4]=quad_mesh(p1,p2,nelx,nely,ndofs)
%
% Purpose:
% Generates a 2D rectangular mesh
%
% Input:
% p1    - Lower left point of rectangle [x1, y1]
% p2    - Upper right point of rectangle [x2, y2]
% nelx  - Number of elments in x direction
% nely  - Number of elments in y direction
% ndofs - Number of degrees of freedom per node
%
% Output:
% Edof  - Connectivity matrix for mesh, cf. Calfem Toolbox
% Ex    - Elementwise x-coordinates, cf. Calfem Toolbox
% Ey    - Elementwise y-coordinates, cf. Calfem Toolbox
% Bi    - Matrix containing boundary dofs for segment i (i=1,2,3,4)
%         First column -> 1st dofs, second column -> 2nd dofs and so on   
%         size = (num boundary nodes on segment) x ndofs
%         B1 = Bottom side     B2 = Right side
%         B3 = Upper side      B4 = Left side         
% Written by
% Jim Brouzoulis, February 15, 2011

% Element degrees of freedom
a = 1:(nelx+1)*(nely);
Elcon = [a; a+1; a+nelx+2; a+nelx+1];
Elcon(:,nelx+1:nelx+1:(nelx+1)*(nely))=[];
dofs = reshape(1:(nelx+1)*(nely+1)*ndofs,ndofs,(nelx+1)*(nely+1))';
Edof = [(1:nelx*nely)' reshape( dofs(Elcon(:),:)',ndofs*4,nelx*nely)'];

% Element coordinates
[y,x]=meshgrid(linspace(p1(2),p2(2),nely+1),linspace(p1(1),p2(1),nelx+1));
Ex=reshape(x(Elcon(:)),4,nelx*nely)';
Ey=reshape(y(Elcon(:)),4,nelx*nely)';

% Boundary dofs
B1 = dofs(1:nelx+1,:);
B2 = dofs((nelx+1)*(1:(nely+1)),:);
B3 = dofs((nelx+1)*(nely+1):-1:(nelx+1)*(nely)+1,:);
B4 = dofs((nelx+1)*(nely:-1:0)+1,:);

