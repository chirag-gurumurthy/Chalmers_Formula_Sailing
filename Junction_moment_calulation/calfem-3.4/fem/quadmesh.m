function [ Coord, Dof, Enode, Edof, Ex, Ey, DofRight, DofTop, ...
 DofLeft, DofBottom, cdof ] = quadmesh( Lx, Ly, nx, ny, plotMesh )
% function [ Coord, Dof, Enode, Edof, Ex, Ey, DofRight, DofTop, ...
% DofLeft, DofBottom, cdof ] = quadmesh( Lx, Ly, nx, ny )
%------------------------------------------------------------------
% Purpose: Generate a mesh a 4-node rectangular elements, together
%          with the relevant "CALFEM" discretization matrices
%------------------------------------------------------------------
% Input: Lx        Length in x-direction
%        Ly        Length in y-direction
%        nx        Number of elements in x-direction
%        ny        Number of elements in y-direction
%        plotMesh  'yes' if you want to plot the undeformed mesh,
%                  'no' if you dont want to plot it
%------------------------------------------------------------------
% Output: Coord     Nodal coordinates matrix, size = (NoNodes x 2)
%         Dof       Nodal dofs matrix, size = (NoNodes x 5)
%         Enode     Topology matrix in terms of nodes
%         Edof      Topology matrix in terms of dofs, size = (NoElem x 21)
%                   In each row in Edof, the dofs are arranged 
%                   anti-clockwise starting from the bottom left node
%                   of each element (5 dofs/node). For each node, first
%                   the three translational dofs are counted (ux,uy,uz) 
%                   and then the two rotational (theta_x,theta_y)
%         Ex        Element nodal x-coordinates, where the nodes are
%                   counted anti-clockwise starting from the bottom left
%                   node of each element.
%         Ey        Element nodal y-coordinates, where the nodes are
%                   counted anti-clockwise starting from the bottom left
%                   node of each element.
%         DofRight  Translational dofs (ux,uy,uz) at the right boundary,
%                   size = (3xNoNodesRightBoundary, 1)
%         DofTop    The same as DofRight for the top boundary
%         DofLeft   The same as DofRight for the left boundary
%         DofBottom The same as DofRight for the bottom boundary
%         cdof      Translational dofs (ux,uy,uz) at the center node
%                   of the rectangular domain (only for even number
%                   of elements in x- and y-direction)
%-----------------------------------------------------------------------
% Created by: Dimosthenis Floros, 20160214
%-----------------------------------------------------------------------

% Discretization in finite elements and nodes

NoElem = nx*ny;

NoNodes = ( nx + 1 )*( ny + 1 );

% Create a basic grid of coordinates in x- and y-direction

CoordX = linspace( 0, Lx, nx + 1 )';
CoordY = linspace( 0, Ly, ny + 1 )';

% Compute the nodal coordinates matrix

Coord = zeros( NoNodes, 2 );

for rowIndex = 1:( ny + 1 )
 
 for columnIndex = 1:( nx + 1 )
  
  Coord( ( rowIndex - 1 )*( nx +1 ) + columnIndex, 1 ) = ...
   CoordX( columnIndex );
  Coord( ( rowIndex - 1 )*( nx +1 ) + columnIndex, 2 ) = ...
   CoordY( rowIndex );
  
 end
 
end

% Compute the nodal dofs matrix

Dof = zeros( NoNodes, 5 );

for nodeIndex = 1:NoNodes
 
 Dof( nodeIndex, : ) = ( nodeIndex - 1 )*5 + 1: nodeIndex*5;
 
end

% Compute the topology matrix in terms of nodes

Enode = zeros( NoElem, 4 );

for elYindex = 1:ny
 
 for elXindex = 1:nx
 
 Enode( (elYindex-1)*nx + elXindex, 1:2 ) = ...
  (elYindex-1)*(nx+1)+elXindex:(elYindex-1)*(nx+1)+elXindex+1;
 
 Enode( (elYindex-1)*nx + elXindex, 3:4 ) = ...
  elYindex*(nx+1)+elXindex:elYindex*(nx+1)+elXindex+1;
  
 end
 
end

% Compute the topology matrix in terms of dofs

Edof = zeros( NoElem, 1 + 20 );
Edof( :, 1 ) = linspace( 1, NoElem, NoElem )';

for elemIndex = 1:NoElem
 
 Edof( elemIndex, 2:end ) = [ Dof( Enode( elemIndex, 1 ), : ) ...
                              Dof( Enode( elemIndex, 2 ), : ) ...
                              Dof( Enode( elemIndex, 4 ), : ) ...
                              Dof( Enode( elemIndex, 3 ), : )];
 
end

% Compute the element nodal coordinates in x- and y-direction

Ex = zeros( NoElem, 4 );
Ey = zeros( NoElem, 4 );

for elemIndex = 1:NoElem
 
 Ex( elemIndex, [ 1 2 3 4 ] ) = [ Coord( Enode( elemIndex, 1 ), 1 ) ...
                                  Coord( Enode( elemIndex, 2 ), 1 ) ...
                                  Coord( Enode( elemIndex, 4 ), 1 ) ...
                                  Coord( Enode( elemIndex, 3 ), 1 )];
 Ey( elemIndex, [ 1 2 3 4 ] ) = [ Coord( Enode( elemIndex, 1 ), 2 ) ...
                                  Coord( Enode( elemIndex, 2 ), 2 ) ...
                                  Coord( Enode( elemIndex, 4 ), 2 ) ...
                                  Coord( Enode( elemIndex, 3 ), 2 )];

end

% Compute the translational dofs at the boundaries

indicesRight = find( Coord( :, 1 ) == Lx );
DofRight = sort( reshape( Dof( indicesRight, [ 1 2 3 ] ), ...
                 length( indicesRight )*3, 1 ) );

indicesTop = find( Coord( :, 2 ) == Ly );
DofTop = sort( reshape( Dof( indicesTop, [ 1 2 3 ] ), ...
               length( indicesTop )*3, 1 ) );

indicesLeft = find( Coord( :, 1 ) == 0 );
DofLeft = sort( reshape( Dof( indicesLeft, [ 1 2 3 ] ), ...
                length( indicesLeft )*3, 1 ) );

indicesBottom = find( Coord( :, 2 ) == 0 );
DofBottom = sort( reshape( Dof( indicesBottom, [ 1 2 3 ] ), ...
                  length( indicesBottom )*3, 1 ) );

% Plot the undeformed mesh

if strcmp( plotMesh, 'yes' ) == 1

figure
eldraw2( Ex, Ey, [ 1 1 0 ], linspace( 1, NoElem, NoElem )' );

end

% Find the dofs at the center node

if mod(nx,2) == 0 && mod(ny,2) == 0
 
 centerElement = ny/2*nx - nx/2 - 1;
 cdof = Edof( centerElement, [ 12 13 14 ] );
 
else disp('Number of elements in x and y direction not even.');
     disp('Please choose even numbers, in order to get');
     disp('the translational dofs at the center of the shell element.');
     
     cdof = 0;
     
end

end