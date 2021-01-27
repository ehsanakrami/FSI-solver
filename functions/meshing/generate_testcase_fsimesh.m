function fsimesh = generate_testcase_fsimesh(Lx,Ly,Lz,Nx,Ny,Nz,Ne_beam_x,Ne_beam_y,Ne_beam_z,isave)

% Computes a test-case mesh, both with 8-node and 20-node blocks.
%
% Inputs: 
% - Lx,Ly,Lz: Dimensions of the box (in meters) where the fluid and part of the solid domains will fit
% - Nx,Ny,Nz: Number of nodes in each direction for that box.
% - Ne_beam_x,Ne_beam_y,Ne_beam_z: Number of elements in the solid domain in each direction
% - isave: If it is equal to 1, the mesh will be saved in a .mat file.
%
% Output:
% - fsimesh: structure containing the following fields:
%   - T_X: a (2,1)-cell containing the Nx3 node coordinates for the 8-node and 20-node meshes
%   - nodelist: a (2,2)-cell containing lists of the nodes in the solid and fluid domains (rows) for both meshes (columns)
%   - T_E_volume: a (2,2)-cell containing connectivity tables for the solid and fluid domains (rows) for both meshes (columns)
%   - T_E_surface: a (2,1)-cell containing cells for the solid and fluid domains. Each cell is (6,2) containing 2D connectivity tables for each outer face of the domain (rows) for each mesh (columns). 
%   - T_E_fsi: a (2,1)-cell containing 2D connectivity tables for the fluid-structure interface
%   - T_DOF_volume,T_DOF_surface,T_DOF_fsi: similar to the connectivity tables, but written in terms of numbers of degrees of freedom assuming 3 degrees of freedom per node

% Number of elements and nodes

Ne_x = Nx-1 ;
Ne_y = Ny-1 ;
Ne_z = Nz-1 ;

Ne = Ne_x * Ne_y * Ne_z ; % Total number of elements

Nnodes_line_edgex   = 2*Ne_x + 1 ;
Nnodes_line_middle  = 1*Ne_x + 1 ;

Nnodes_edgex = Nnodes_line_edgex * (Ne_y+1) * (Ne_z+1) ;
Nnodes_middle = Nnodes_line_middle * (Ne_y+1) * Ne_z + Nnodes_line_middle * Ne_y * (Ne_z+1) ;
Nnodes = Nnodes_edgex + Nnodes_middle ; 

Nnodes_corner = (Ne_x+1)*(Ne_y+1)*(Ne_z+1) ;

% Dimensions of elements

Le = Lx/Ne_x ; % Length of an element
We = Ly/Ne_y  ; % Width of an element
He = Lz/Ne_z ; % Height of an element

% Indexing of nodes in an element (scheme) : 

% node1  : (0,0,0)
% node2  : (2,0,0)
% node3  : (2,2,0)
% node4  : (0,2,0)
% node5  : (0,0,2)
% node6  : (2,0,2)
% node7  : (2,2,2)
% node8  : (0,2,2)
% node9  : (1,0,0)
% node10 : (2,1,0)
% node11 : (1,2,0)
% node12 : (0,1,0)
% node13 : (0,0,1)
% node14 : (2,0,1)
% node15 : (2,2,1)
% node16 : (0,2,1)
% node17 : (1,0,2)
% node18 : (2,1,2)
% node19 : (1,2,2)
% node20 : (0,1,2)

% Indexing of elements

% Elements are indexed starting from the one whose first node is at (0,0,0)
% and then going along length, then increasing width and then increasing
% height

% Indexing of nodes

% Nodes are indexed in the same way as elements are for corner nodes (1 to
% 8). Then we index nodes on the middle of edges of direction x, then nodes
% on the middle of edges of direction y , then nodes on the middle of edges
% of direction z

% Connectivity table  and node coordinates for the box

T_E_20 = zeros(Ne,20) ;

for mm = 1:Ne_x
    for nn = 1:Ne_y
        for pp = 1:Ne_z
                        
            T_E_20(mm+((nn-1)*Ne_x)+((pp-1)*Ne_x*Ne_y),1:8)  = [ mm+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))                              ,...
                                                             mm+1+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))                            ,... 
                                                             mm+(Ne_x+2)+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))                     ,...
                                                             mm+(Ne_x+1)+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))                     ,...
                                                             mm+((Ne_x+1)*(Ne_y+1))+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))          ,...
                                                             mm+1+((Ne_x+1)*(Ne_y+1))+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))        ,...
                                                             mm+(Ne_x+2)+((Ne_x+1)*(Ne_y+1))+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1)) ,...
                                                             mm+(Ne_x+1)+((Ne_x+1)*(Ne_y+1))+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1))  ] ;
                                                        
            T_E_20(mm+((nn-1)*Ne_x)+((pp-1)*Ne_x*Ne_y),9:20) = [ Nnodes_corner + mm + (nn-1)*Ne_x + (pp-1)*Ne_x*(Ne_y+1)                                                                          ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*Ne_y + 1                                         ,...
                                                             Nnodes_corner + mm + (nn-1)*Ne_x + (pp-1)*Ne_x*(Ne_y+1) + Ne_x                                                                   ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*Ne_y                                             ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + (Ne_x+1)*Ne_y*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*(Ne_y+1)                ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + (Ne_x+1)*Ne_y*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*(Ne_y+1) + 1            ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + (Ne_x+1)*Ne_y*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*(Ne_y+1) + (Ne_x+1) + 1 ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + (Ne_x+1)*Ne_y*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*(Ne_y+1) + (Ne_x+1)     ,...
                                                             Nnodes_corner + mm + (nn-1)*Ne_x + (pp-1)*Ne_x*(Ne_y+1) + Ne_x*(Ne_y+1)                                                          ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*Ne_y + (Ne_x+1)*Ne_y + 1                         ,...
                                                             Nnodes_corner + mm + (nn-1)*Ne_x + (pp-1)*Ne_x*(Ne_y+1) + Ne_x*(Ne_y+1) + Ne_x                                                   ,...
                                                             Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*Ne_y + (Ne_x+1)*Ne_y                              ] ;
                                                             
                                   
        end    
    end
end

T_E_8 = T_E_20(:,1:8) ;

T_X_20 = zeros(Nnodes,3) ;
Nnodes_8 = Nx*Ny*Nz ;
T_X_8 = zeros(Nnodes_8,3) ;

for mm = 1:(Ne_x+1)
    for nn = 1:(Ne_y+1)
        for pp = 1:(Ne_z+1)
            
            T_X_20(mm+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1)),:) = [ (mm-1)*Le , (nn-1)*We , (pp-1)*He ] ;
            T_X_8(mm+((nn-1)*(Ne_x+1))+((pp-1)*(Ne_x+1)*(Ne_y+1)),:)  = [ (mm-1)*Le , (nn-1)*We , (pp-1)*He ] ;

        end
    end
end

for mm = 1:Ne_x
    for nn = 1:(Ne_y+1)
        for pp = 1:(Ne_z+1)
            
            T_X_20( Nnodes_corner + mm + (nn-1)*Ne_x + (pp-1)*Ne_x*(Ne_y+1) , : ) = [ (mm-1)*Le + Le/2 , (nn-1)*We , (pp-1)*He ] ;
            
        end
    end
end

for mm = 1:(Ne_x+1)
    for nn = 1:Ne_y
        for pp = 1:(Ne_z+1)
            
            T_X_20( Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*Ne_y , : ) = [ (mm-1)*Le , (nn-1)*We + We/2 , (pp-1)*He ] ;
            
        end
    end
end

for mm = 1:(Ne_x+1)
    for nn = 1:(Ne_y+1)
        for pp = 1:Ne_z
            
            T_X_20( Nnodes_corner + Ne_x*(Ne_y+1)*(Ne_z+1) + (Ne_x+1)*Ne_y*(Ne_z+1) + mm + (nn-1)*(Ne_x+1) + (pp-1)*(Ne_x+1)*(Ne_y+1) , : ) = [ (mm-1)*Le , (nn-1)*We , (pp-1)*He + He/2 ] ;
            
        end
    end
end

Ne3D = size(T_E_20,1) ;

% Identifying elements in the box belonging to the solid domain

Lbeam_x = Ne_beam_x*Le ;
Lbeam_y = Ne_beam_y*We ;
Lbeam_z = Ne_beam_z*He ;

ls_beam_elts = [] ;

eps1 = min([Le He We])*1e-3 ;

for ii = 1:size(T_E_8,1)
    
    X_ii = T_X_8(T_E_8(ii,:),:) ;
    
    if    (sum(abs(X_ii(:,1) - Lx/2) <= (Lbeam_x/2 + eps1)) == 8) ...
       && (sum(abs(X_ii(:,2) - Ly/2) <= (Lbeam_y/2 + eps1)) == 8) ...
       && (sum(X_ii(:,3) <= (Lbeam_z + eps1)) == 8)
        
        ls_beam_elts = [ls_beam_elts , ii] ; %#ok<AGROW>
                
    end
    
end

% Connectivity tables for the solid and fluid domains

T_E_beam_8 = T_E_8(ls_beam_elts,:) ;
T_E_beam_20 = T_E_20(ls_beam_elts,:) ;

T_E_fluid_8 = T_E_8(setdiff(1:Ne3D,ls_beam_elts),:) ;
T_E_fluid_20 = T_E_20(setdiff(1:Ne3D,ls_beam_elts),:) ;

% 2D connectivity table of the fluid-structure interface

T_E_fsi_8 = [] ;
T_E_fsi_20 = [] ;

for ii = 1:size(T_E_beam_8,1)
    
    ls_faces_fsi = 1:6 ;
    
    nodes_face1_20 = T_E_beam_20(ii,[1 2 3 4 9 10 11 12]) ;
    nodes_face2_20 = T_E_beam_20(ii,[1 2 6 5 9 13 14 17]) ;
    nodes_face3_20 = T_E_beam_20(ii,[1 4 8 5 12 13 16 20]) ;
    nodes_face4_20 = T_E_beam_20(ii,[5 6 7 8 17 18 19 20]) ;
    nodes_face5_20 = T_E_beam_20(ii,[3 4 8 7 11 15 16 19]) ;
    nodes_face6_20 = T_E_beam_20(ii,[2 3 7 6 10 14 15 18]) ;
    
    nodes_faces_20 = [nodes_face1_20;nodes_face2_20;nodes_face3_20;nodes_face4_20;nodes_face5_20;nodes_face6_20] ;
    
    nodes_face1_8 = T_E_beam_8(ii,[1 2 3 4]) ;
    nodes_face2_8 = T_E_beam_8(ii,[1 2 6 5]) ;
    nodes_face3_8 = T_E_beam_8(ii,[1 4 8 5]) ;
    nodes_face4_8 = T_E_beam_8(ii,[5 6 7 8]) ;
    nodes_face5_8 = T_E_beam_8(ii,[3 4 8 7]) ;
    nodes_face6_8 = T_E_beam_8(ii,[2 3 7 6]) ;
    
    nodes_faces_8 = [nodes_face1_8;nodes_face2_8;nodes_face3_8;nodes_face4_8;nodes_face5_8;nodes_face6_8] ;
    
    if sum(ismember(nodes_face1_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,1) ;
    end
    if sum(ismember(nodes_face2_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,2) ;
    end
    if sum(ismember(nodes_face3_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,3) ;
    end
    if sum(ismember(nodes_face4_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,4) ;
    end
    if sum(ismember(nodes_face5_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,5) ;
    end
    if sum(ismember(nodes_face6_8,T_E_fluid_8))==4
        ls_faces_fsi = setdiff(ls_faces_fsi,6) ;
    end
    
    ls_faces_on_fsi = setdiff(1:6,ls_faces_fsi) ;
    
    T_E_fsi_20 = [T_E_fsi_20 ; nodes_faces_20(ls_faces_on_fsi,:)] ; %#ok<AGROW>
    
    T_E_fsi_8 = [T_E_fsi_8 ; nodes_faces_8(ls_faces_on_fsi,:)] ; %#ok<AGROW>
                   
end

T_E_fsi = cell(1,2) ;
T_E_fsi{1,1} = T_E_fsi_8 ;
T_E_fsi{1,2} = T_E_fsi_20 ;

% Adding elements for the solid domain in the vertical direction if the
% structure's height is greater than the fluid domain's height.

Ne_beam_out_z = Ne_beam_z - (Nz-1) ;

if Ne_beam_out_z > 0
    
    top_beam_elts = find( (T_X_8(T_E_beam_8(:,5),3) == Lz) & (T_X_8(T_E_beam_8(:,6),3) == Lz) & (T_X_8(T_E_beam_8(:,7),3) == Lz) & (T_X_8(T_E_beam_8(:,8),3) == Lz) ) ;
    
    top_beam_nodes_20 = T_E_beam_20(top_beam_elts,[5 6 7 8 13 14 15 16 17 18 19 20]) ;
    top_beam_nodes_20 = unique(top_beam_nodes_20(:)) ;
    top_beam_nodes_8 = T_E_beam_8(top_beam_elts,[5 6 7 8]) ;
    top_beam_nodes_8 = unique(top_beam_nodes_8(:)) ;
    
    sub_T_E_8 = T_E_beam_8(top_beam_elts,:) ;
    sub_T_E_20 = T_E_beam_20(top_beam_elts,:) ;
    
    Ne_sub = size(sub_T_E_8,1) ;
    Np_sub_20 = length(top_beam_nodes_20) ;
    Np_sub_8 = length(top_beam_nodes_8) ;
    
    mat_dz_20 = repmat([0 0 He],[length(top_beam_nodes_20) 1]) ;
    mat_dz_8 = repmat([0 0 He],[length(top_beam_nodes_8) 1]) ;
    
    %
    
    old_T_E_8 = sub_T_E_8 ;
    old_T_E_20 = sub_T_E_20 ;
    
    old_T_X_20 = T_X_20(top_beam_nodes_20,:) ;
    old_T_X_8 = T_X_8(top_beam_nodes_8,:) ;
    
    new_T_X_20 = old_T_X_20 + mat_dz_20 ;
    new_T_X_8 = old_T_X_8 + mat_dz_8 ;
    
    T_X_20 = [T_X_20 ; new_T_X_20] ;
    T_X_8 = [T_X_8 ; new_T_X_8] ;
    
    new_T_E_8 = zeros(Ne_sub,8) ;
    new_T_E_8(:,[1 2 3 4]) = old_T_E_8(:,[5 6 7 8]) ;
    for ii = 1:Ne_sub
        for jj = 1:4
            node_ij = old_T_E_8(ii,jj+4) ;
            new_node_ij = find(T_X_8(:,3)==T_X_8(node_ij,3)+He & T_X_8(:,1)==T_X_8(node_ij,1) & T_X_8(:,2)==T_X_8(node_ij,2)) ;
            new_T_E_8(ii,jj+4) = new_node_ij ;
        end
    end
    
    new_T_E_20 = zeros(Ne_sub,20) ;
    new_T_E_20(:,[1 2 3 4 9 10 11 12]) = old_T_E_20(:,[5 6 7 8 17 18 19 20]) ;
    for ii = 1:Ne_sub
        for jj = 1:4
            node_ij = old_T_E_20(ii,jj+4) ;
            new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
            new_T_E_20(ii,jj+4) = new_node_ij ;
        end
        for jj = 1:4
            node_ij = old_T_E_20(ii,jj+12) ;
            new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
            new_T_E_20(ii,jj+12) = new_node_ij ;
        end
        for jj = 1:4
            node_ij = old_T_E_20(ii,jj+16) ;
            new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
            new_T_E_20(ii,jj+16) = new_node_ij ;
        end
    end
    
    T_E_beam_8 = [T_E_beam_8 ; new_T_E_8] ;
    T_E_beam_20 = [T_E_beam_20 ; new_T_E_20] ;
    
    old_T_E_8 = new_T_E_8 ;
    old_T_E_20 = new_T_E_20 ;
    old_T_X_20 = new_T_X_20 ;
    old_T_X_8 = new_T_X_8 ;
    
    if Ne_beam_out_z > 1
        
        new_T_X_20 = old_T_X_20 + mat_dz_20 ;
        new_T_X_8 = old_T_X_8 + mat_dz_8 ;
        
        T_X_20 = [T_X_20 ; new_T_X_20] ;
        T_X_8 = [T_X_8 ; new_T_X_8] ;
        
        new_T_E_8 = zeros(Ne_sub,8) ;
        new_T_E_8(:,[1 2 3 4]) = old_T_E_8(:,[5 6 7 8]) ;
        for ii = 1:Ne_sub
            for jj = 1:4
                node_ij = old_T_E_8(ii,jj+4) ;
                new_node_ij = find(T_X_8(:,3)==T_X_8(node_ij,3)+He & T_X_8(:,1)==T_X_8(node_ij,1) & T_X_8(:,2)==T_X_8(node_ij,2)) ;
                new_T_E_8(ii,jj+4) = new_node_ij ;
            end
        end
        
        new_T_E_20 = zeros(Ne_sub,20) ;
        new_T_E_20(:,[1 2 3 4 9 10 11 12]) = old_T_E_20(:,[5 6 7 8 17 18 19 20]) ;
        for ii = 1:Ne_sub
            for jj = 1:4
                node_ij = old_T_E_20(ii,jj+4) ;
                new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
                new_T_E_20(ii,jj+4) = new_node_ij ;
            end
            for jj = 1:4
                node_ij = old_T_E_20(ii,jj+12) ;
                new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
                new_T_E_20(ii,jj+12) = new_node_ij ;
            end
            for jj = 1:4
                node_ij = old_T_E_20(ii,jj+16) ;
                new_node_ij = find(T_X_20(:,3)==T_X_20(node_ij,3)+He & T_X_20(:,1)==T_X_20(node_ij,1) & T_X_20(:,2)==T_X_20(node_ij,2)) ;
                new_T_E_20(ii,jj+16) = new_node_ij ;
            end
        end
        
        T_E_beam_8 = [T_E_beam_8 ; new_T_E_8] ;
        T_E_beam_20 = [T_E_beam_20 ; new_T_E_20] ;
        
        old_T_E_8 = new_T_E_8 ;
        old_T_E_20 = new_T_E_20 ;
        old_T_X_20 = new_T_X_20 ;
        old_T_X_8 = new_T_X_8 ;
        
        %
        
        for ii = 1:(Ne_beam_out_z-2)
            
            new_T_E_8 = old_T_E_8 + Np_sub_8 ;
            new_T_E_20 = old_T_E_20 + Np_sub_20 ;
            
            new_T_X_20 = old_T_X_20 + mat_dz_20 ;
            new_T_X_8 = old_T_X_8 + mat_dz_8 ;
            
            T_E_beam_8 = [T_E_beam_8 ; new_T_E_8] ; %#ok<AGROW>
            T_E_beam_20 = [T_E_beam_20 ; new_T_E_20] ; %#ok<AGROW>
            
            T_X_20 = [T_X_20 ; new_T_X_20] ; %#ok<AGROW>
            T_X_8 = [T_X_8 ; new_T_X_8] ; %#ok<AGROW>
            
            old_T_E_8 = new_T_E_8 ;
            old_T_E_20 = new_T_E_20 ;
            old_T_X_20 = new_T_X_20 ;
            old_T_X_8 = new_T_X_8 ;
            
        end
        
    else
        
        for ii = 1:(Ne_beam_out_z-1)
            
            new_T_E_8 = old_T_E_8 + Np_sub_8 ;
            new_T_E_20 = old_T_E_20 + Np_sub_20 ;
            
            new_T_X_20 = old_T_X_20 + mat_dz_20 ;
            new_T_X_8 = old_T_X_8 + mat_dz_8 ;
            
            T_E_beam_8 = [T_E_beam_8 ; new_T_E_8] ; %#ok<AGROW>
            T_E_beam_20 = [T_E_beam_20 ; new_T_E_20] ; %#ok<AGROW>
            
            T_X_20 = [T_X_20 ; new_T_X_20] ; %#ok<AGROW>
            T_X_8 = [T_X_8 ; new_T_X_8] ; %#ok<AGROW>
            
            old_T_E_8 = new_T_E_8 ;
            old_T_E_20 = new_T_E_20 ;
            old_T_X_20 = new_T_X_20 ;
            old_T_X_8 = new_T_X_8 ;
            
        end
        
    end
    
end

% Final 3D connectivity tables

T_E_volume = cell(2,2) ;
T_E_volume{1,1} = T_E_beam_8 ;
T_E_volume{1,2} = T_E_beam_20 ;
T_E_volume{2,1} = T_E_fluid_8 ;
T_E_volume{2,2} = T_E_fluid_20 ;

% Coordinates are centered, with respect to the x- and y-directions, on the
% center of the box

T_X_20(:,1) = T_X_20(:,1) - Lx/2 ;
T_X_20(:,2) = T_X_20(:,2) - Ly/2 ;
T_X_8(:,1) = T_X_8(:,1) - Lx/2 ;
T_X_8(:,2) = T_X_8(:,2) - Ly/2 ;

% 2D connectivity tables for the 6 outer faces of the fluid domain

ls_face1_8 = find(T_X_8(:,1)==0) ;
ls_face2_8 = find(T_X_8(:,1)==Lx) ;
ls_face3_8 = find(T_X_8(:,2)==0) ;
ls_face4_8 = find(T_X_8(:,2)==Ly) ;
ls_face5_8 = find(T_X_8(:,3)==0) ;
ls_face6_8 = find(T_X_8(:,3)==Lz) ;

ls_face1_20 = find(T_X_20(:,1)==0) ;
ls_face2_20 = find(T_X_20(:,1)==Lx) ;
ls_face3_20 = find(T_X_20(:,2)==0) ;
ls_face4_20 = find(T_X_20(:,2)==Ly) ;
ls_face5_20 = find(T_X_20(:,3)==0) ;
ls_face6_20 = find(T_X_20(:,3)==Lz) ;

nodes_fluid_faces_8 = cell(6,1) ;
nodes_fluid_faces_8{1} = ls_face1_8 ;
nodes_fluid_faces_8{2} = ls_face2_8 ;
nodes_fluid_faces_8{3} = ls_face3_8 ;
nodes_fluid_faces_8{4} = ls_face4_8 ;
nodes_fluid_faces_8{5} = ls_face5_8 ;
nodes_fluid_faces_8{6} = ls_face6_8 ;

nodes_fluid_faces_20 = cell(6,1) ;
nodes_fluid_faces_20{1} = ls_face1_20 ;
nodes_fluid_faces_20{2} = ls_face2_20 ;
nodes_fluid_faces_20{3} = ls_face3_20 ;
nodes_fluid_faces_20{4} = ls_face4_20 ;
nodes_fluid_faces_20{5} = ls_face5_20 ;
nodes_fluid_faces_20{6} = ls_face6_20 ;

T_E_fluid_face1_20 = [] ;
T_E_fluid_face1_8 = [] ;
T_E_fluid_face2_20 = [] ;
T_E_fluid_face2_8 = [] ;
T_E_fluid_face3_20 = [] ;
T_E_fluid_face3_8 = [] ;
T_E_fluid_face4_20 = [] ;
T_E_fluid_face4_8 = [] ;
T_E_fluid_face5_20 = [] ;
T_E_fluid_face5_8 = [] ;
T_E_fluid_face6_20 = [] ;
T_E_fluid_face6_8 = [] ;

for ii = 1:size(T_E_fluid_8,1)
    
    ls_faces_face1 = 1:6 ;
    ls_faces_face2 = 1:6 ;
    ls_faces_face3 = 1:6 ;
    ls_faces_face4 = 1:6 ;
    ls_faces_face5 = 1:6 ;
    ls_faces_face6 = 1:6 ;
    
    nodes_face1_20 = T_E_fluid_20(ii,[1 2 3 4 9 10 11 12]) ;
    nodes_face2_20 = T_E_fluid_20(ii,[1 2 6 5 9 13 14 17]) ;
    nodes_face3_20 = T_E_fluid_20(ii,[1 4 8 5 12 13 16 20]) ;
    nodes_face4_20 = T_E_fluid_20(ii,[5 6 7 8 17 18 19 20]) ;
    nodes_face5_20 = T_E_fluid_20(ii,[3 4 8 7 11 15 16 19]) ;
    nodes_face6_20 = T_E_fluid_20(ii,[2 3 7 6 10 14 15 18]) ;
    
    nodes_faces_20 = [nodes_face1_20;nodes_face2_20;nodes_face3_20;nodes_face4_20;nodes_face5_20;nodes_face6_20] ;
    
    nodes_face1_8 = T_E_fluid_8(ii,[1 2 3 4]) ;
    nodes_face2_8 = T_E_fluid_8(ii,[1 2 6 5]) ;
    nodes_face3_8 = T_E_fluid_8(ii,[1 4 8 5]) ;
    nodes_face4_8 = T_E_fluid_8(ii,[5 6 7 8]) ;
    nodes_face5_8 = T_E_fluid_8(ii,[3 4 8 7]) ;
    nodes_face6_8 = T_E_fluid_8(ii,[2 3 7 6]) ;
    
    nodes_faces_8 = [nodes_face1_8;nodes_face2_8;nodes_face3_8;nodes_face4_8;nodes_face5_8;nodes_face6_8] ;
    
    if sum(ismember(nodes_face1_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face1_8))==8
        ls_faces_face1 = setdiff(ls_faces_face1,6) ;
    end
    if sum(ismember(nodes_face1_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face2_8))==8
        ls_faces_face2 = setdiff(ls_faces_face2,6) ;
    end
    if sum(ismember(nodes_face1_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face3_8))==8
        ls_faces_face3 = setdiff(ls_faces_face3,6) ;
    end
    if sum(ismember(nodes_face1_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face4_8))==8
        ls_faces_face4 = setdiff(ls_faces_face4,6) ;
    end
    if sum(ismember(nodes_face1_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face5_8))==8
        ls_faces_face5 = setdiff(ls_faces_face5,6) ;
    end
    if sum(ismember(nodes_face1_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,2) ;
    end
    if sum(ismember(nodes_face3_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,3) ;
    end
    if sum(ismember(nodes_face4_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,4) ;
    end
    if sum(ismember(nodes_face5_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,5) ;
    end
    if sum(ismember(nodes_face6_8,ls_face6_8))==8
        ls_faces_face6 = setdiff(ls_faces_face6,6) ;
    end
    
    ls_faces_on_face1 = setdiff(1:6,ls_faces_face1) ;
    ls_faces_on_face2 = setdiff(1:6,ls_faces_face2) ;
    ls_faces_on_face3 = setdiff(1:6,ls_faces_face3) ;
    ls_faces_on_face4 = setdiff(1:6,ls_faces_face4) ;
    ls_faces_on_face5 = setdiff(1:6,ls_faces_face5) ;
    ls_faces_on_face6 = setdiff(1:6,ls_faces_face6) ;
    
    T_E_fluid_face1_20 = [T_E_fluid_face1_20 ; nodes_faces_20(ls_faces_on_face1,:)] ; %#ok<AGROW>
    T_E_fluid_face2_20 = [T_E_fluid_face2_20 ; nodes_faces_20(ls_faces_on_face2,:)] ; %#ok<AGROW>
    T_E_fluid_face3_20 = [T_E_fluid_face3_20 ; nodes_faces_20(ls_faces_on_face3,:)] ; %#ok<AGROW>
    T_E_fluid_face4_20 = [T_E_fluid_face4_20 ; nodes_faces_20(ls_faces_on_face4,:)] ; %#ok<AGROW>
    T_E_fluid_face5_20 = [T_E_fluid_face5_20 ; nodes_faces_20(ls_faces_on_face5,:)] ; %#ok<AGROW>
    T_E_fluid_face6_20 = [T_E_fluid_face6_20 ; nodes_faces_20(ls_faces_on_face6,:)] ; %#ok<AGROW>
    
    T_E_fluid_face1_8 = [T_E_fluid_face1_8 ; nodes_faces_8(ls_faces_on_face1,:)] ; %#ok<AGROW>
    T_E_fluid_face2_8 = [T_E_fluid_face2_8 ; nodes_faces_8(ls_faces_on_face2,:)] ; %#ok<AGROW>
    T_E_fluid_face3_8 = [T_E_fluid_face3_8 ; nodes_faces_8(ls_faces_on_face3,:)] ; %#ok<AGROW>
    T_E_fluid_face4_8 = [T_E_fluid_face4_8 ; nodes_faces_8(ls_faces_on_face4,:)] ; %#ok<AGROW>
    T_E_fluid_face5_8 = [T_E_fluid_face5_8 ; nodes_faces_8(ls_faces_on_face5,:)] ; %#ok<AGROW>
    T_E_fluid_face6_8 = [T_E_fluid_face6_8 ; nodes_faces_8(ls_faces_on_face6,:)] ; %#ok<AGROW>
                   
end

T_E_fluid_faces = cell(6,2) ;
T_E_fluid_faces{1,1} = T_E_fluid_face1_8 ;
T_E_fluid_faces{1,2} = T_E_fluid_face1_20 ;
T_E_fluid_faces{2,1} = T_E_fluid_face2_8 ;
T_E_fluid_faces{2,2} = T_E_fluid_face2_20 ;
T_E_fluid_faces{3,1} = T_E_fluid_face3_8 ;
T_E_fluid_faces{3,2} = T_E_fluid_face3_20 ;
T_E_fluid_faces{4,1} = T_E_fluid_face4_8 ;
T_E_fluid_faces{4,2} = T_E_fluid_face4_20 ;
T_E_fluid_faces{5,1} = T_E_fluid_face5_8 ;
T_E_fluid_faces{5,2} = T_E_fluid_face5_20 ;
T_E_fluid_faces{6,1} = T_E_fluid_face6_8 ;
T_E_fluid_faces{6,2} = T_E_fluid_face6_20 ;

% Tables of degrees of freedom

T_DOF_fluid_20 = make_T_DOF(T_E_fluid_20) ;
T_DOF_fluid_8  = make_T_DOF(T_E_fluid_8) ;

T_DOF_beam_20 = make_T_DOF(T_E_beam_20) ;
T_DOF_beam_8  = make_T_DOF(T_E_beam_8) ;

T_DOF_volume = cell(2,2) ;
T_DOF_volume{1,1} = T_DOF_beam_8 ;
T_DOF_volume{1,2} = T_DOF_beam_20 ;
T_DOF_volume{2,1} = T_DOF_fluid_8 ;
T_DOF_volume{2,2} = T_DOF_fluid_20 ;

% 2D connectivity tables for the 6 outer faces of the solid domain

ls_beam_face1_8 = find((T_X_8(:,1)+Lbeam_x/2==0) & (abs(T_X_8(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face2_8 = find((T_X_8(:,1)-Lbeam_x/2==0) & (abs(T_X_8(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face3_8 = find((abs(T_X_8(:,1)-Lx/2)<=Lbeam_x/2) & (T_X_8(:,2)+Lbeam_y/2==0)) ;
ls_beam_face4_8 = find((abs(T_X_8(:,1)-Lx/2)<=Lbeam_x/2) & (T_X_8(:,2)-Lbeam_y/2==0)) ;
ls_beam_face5_8 = find((T_X_8(:,3)==0)            & (abs(T_X_8(:,1)-Lx/2)<=Lbeam_x/2) & (abs(T_X_8(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face6_8 = find((T_X_8(:,3)==Ne_beam_z*He) & (abs(T_X_8(:,1)-Lx/2)<=Lbeam_x/2) & (abs(T_X_8(:,2)-Ly/2)<=Lbeam_y/2)) ;

ls_beam_face1_20 = find((T_X_20(:,1)+Lbeam_x/2==0) & (abs(T_X_20(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face2_20 = find((T_X_20(:,1)-Lbeam_x/2==0) & (abs(T_X_20(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face3_20 = find((abs(T_X_20(:,1)-Lx/2)<=Lbeam_x/2) & (T_X_20(:,2)+Lbeam_y/2==0)) ;
ls_beam_face4_20 = find((abs(T_X_20(:,1)-Lx/2)<=Lbeam_x/2) & (T_X_20(:,2)-Lbeam_y/2==0)) ;
ls_beam_face5_20 = find((T_X_20(:,3)==0)            & (abs(T_X_20(:,1)-Lx/2)<=Lbeam_x/2) & (abs(T_X_20(:,2)-Ly/2)<=Lbeam_y/2)) ;
ls_beam_face6_20 = find((T_X_20(:,3)==Ne_beam_z*He) & (abs(T_X_20(:,1)-Lx/2)<=Lbeam_x/2) & (abs(T_X_20(:,2)-Ly/2)<=Lbeam_y/2)) ;

nodes_beam_faces_8 = cell(6,1) ;
nodes_beam_faces_8{1} = ls_beam_face1_8 ;
nodes_beam_faces_8{2} = ls_beam_face2_8 ;
nodes_beam_faces_8{3} = ls_beam_face3_8 ;
nodes_beam_faces_8{4} = ls_beam_face4_8 ;
nodes_beam_faces_8{5} = ls_beam_face5_8 ;
nodes_beam_faces_8{6} = ls_beam_face6_8 ;

nodes_beam_faces_20 = cell(6,1) ;
nodes_beam_faces_20{1} = ls_beam_face1_20 ;
nodes_beam_faces_20{2} = ls_beam_face2_20 ;
nodes_beam_faces_20{3} = ls_beam_face3_20 ;
nodes_beam_faces_20{4} = ls_beam_face4_20 ;
nodes_beam_faces_20{5} = ls_beam_face5_20 ;
nodes_beam_faces_20{6} = ls_beam_face6_20 ;

nodes_fluid_faces_8{1} = setdiff(nodes_fluid_faces_8{1},nodes_beam_faces_8{1}) ;
nodes_fluid_faces_8{2} = setdiff(nodes_fluid_faces_8{2},nodes_beam_faces_8{2}) ;
nodes_fluid_faces_8{3} = setdiff(nodes_fluid_faces_8{3},nodes_beam_faces_8{3}) ;
nodes_fluid_faces_8{4} = setdiff(nodes_fluid_faces_8{4},nodes_beam_faces_8{4}) ;
nodes_fluid_faces_8{5} = setdiff(nodes_fluid_faces_8{5},nodes_beam_faces_8{5}) ;
nodes_fluid_faces_8{6} = setdiff(nodes_fluid_faces_8{6},nodes_beam_faces_8{6}) ;

nodes_fluid_faces_20{1} = setdiff(nodes_fluid_faces_20{1},nodes_beam_faces_20{1}) ;
nodes_fluid_faces_20{2} = setdiff(nodes_fluid_faces_20{2},nodes_beam_faces_20{2}) ;
nodes_fluid_faces_20{3} = setdiff(nodes_fluid_faces_20{3},nodes_beam_faces_20{3}) ;
nodes_fluid_faces_20{4} = setdiff(nodes_fluid_faces_20{4},nodes_beam_faces_20{4}) ;
nodes_fluid_faces_20{5} = setdiff(nodes_fluid_faces_20{5},nodes_beam_faces_20{5}) ;
nodes_fluid_faces_20{6} = setdiff(nodes_fluid_faces_20{6},nodes_beam_faces_20{6}) ;

T_E_beam_face1_20 = [] ;
T_E_beam_face1_8 = [] ;
T_E_beam_face2_20 = [] ;
T_E_beam_face2_8 = [] ;
T_E_beam_face3_20 = [] ;
T_E_beam_face3_8 = [] ;
T_E_beam_face4_20 = [] ;
T_E_beam_face4_8 = [] ;
T_E_beam_face5_20 = [] ;
T_E_beam_face5_8 = [] ;
T_E_beam_face6_20 = [] ;
T_E_beam_face6_8 = [] ;

for ii = 1:size(T_E_beam_8,1)
    
    ls_face1 = 1:6 ;
    ls_face2 = 1:6 ;
    ls_face3 = 1:6 ;
    ls_face4 = 1:6 ;
    ls_face5 = 1:6 ;
    ls_face6 = 1:6 ;
    
    nodes_face1_20 = T_E_fluid_20(ii,[1 2 3 4 9 10 11 12]) ;
    nodes_face2_20 = T_E_fluid_20(ii,[1 2 6 5 9 13 14 17]) ;
    nodes_face3_20 = T_E_fluid_20(ii,[1 4 8 5 12 13 16 20]) ;
    nodes_face4_20 = T_E_fluid_20(ii,[5 6 7 8 17 18 19 20]) ;
    nodes_face5_20 = T_E_fluid_20(ii,[3 4 8 7 11 15 16 19]) ;
    nodes_face6_20 = T_E_fluid_20(ii,[2 3 7 6 10 14 15 18]) ;
    
    nodes_faces_20 = [nodes_face1_20;nodes_face2_20;nodes_face3_20;nodes_face4_20;nodes_face5_20;nodes_face6_20] ;
    
    nodes_face1_8 = T_E_fluid_8(ii,[1 2 3 4]) ;
    nodes_face2_8 = T_E_fluid_8(ii,[1 2 6 5]) ;
    nodes_face3_8 = T_E_fluid_8(ii,[1 4 8 5]) ;
    nodes_face4_8 = T_E_fluid_8(ii,[5 6 7 8]) ;
    nodes_face5_8 = T_E_fluid_8(ii,[3 4 8 7]) ;
    nodes_face6_8 = T_E_fluid_8(ii,[2 3 7 6]) ;
    
    nodes_faces_8 = [nodes_face1_8;nodes_face2_8;nodes_face3_8;nodes_face4_8;nodes_face5_8;nodes_face6_8] ;
    
    if sum(ismember(nodes_face1_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face1_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face1_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face1_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face1_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face1_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face2_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face3_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face4_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face5_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face1_8))==8
        ls_face1 = setdiff(ls_face1,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face2_8))==8
        ls_face2 = setdiff(ls_face2,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face3_8))==8
        ls_face3 = setdiff(ls_face3,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face4_8))==8
        ls_face4 = setdiff(ls_face4,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face5_8))==8
        ls_face5 = setdiff(ls_face5,1) ;
    end
    if sum(ismember(nodes_face6_8,ls_beam_face6_8))==8
        ls_face6 = setdiff(ls_face6,1) ;
    end
    
    ls_faces_on_face1 = setdiff(1:6,ls_face1) ;
    ls_faces_on_face2 = setdiff(1:6,ls_face2) ;
    ls_faces_on_face3 = setdiff(1:6,ls_face3) ;
    ls_faces_on_face4 = setdiff(1:6,ls_face4) ;
    ls_faces_on_face5 = setdiff(1:6,ls_face5) ;
    ls_faces_on_face6 = setdiff(1:6,ls_face6) ;
    
    T_E_beam_face1_20 = [T_E_beam_face1_20 ; nodes_faces_20(ls_faces_on_face1,:)] ; %#ok<AGROW>
    T_E_beam_face2_20 = [T_E_beam_face2_20 ; nodes_faces_20(ls_faces_on_face2,:)] ; %#ok<AGROW>
    T_E_beam_face3_20 = [T_E_beam_face3_20 ; nodes_faces_20(ls_faces_on_face3,:)] ; %#ok<AGROW>
    T_E_beam_face4_20 = [T_E_beam_face4_20 ; nodes_faces_20(ls_faces_on_face4,:)] ; %#ok<AGROW>
    T_E_beam_face5_20 = [T_E_beam_face5_20 ; nodes_faces_20(ls_faces_on_face5,:)] ; %#ok<AGROW>
    T_E_beam_face6_20 = [T_E_beam_face6_20 ; nodes_faces_20(ls_faces_on_face6,:)] ; %#ok<AGROW>
    
    T_E_beam_face1_8 = [T_E_beam_face1_8 ; nodes_faces_8(ls_faces_on_face1,:)] ; %#ok<AGROW>
    T_E_beam_face2_8 = [T_E_beam_face2_8 ; nodes_faces_8(ls_faces_on_face2,:)] ; %#ok<AGROW>
    T_E_beam_face3_8 = [T_E_beam_face3_8 ; nodes_faces_8(ls_faces_on_face3,:)] ; %#ok<AGROW>
    T_E_beam_face4_8 = [T_E_beam_face4_8 ; nodes_faces_8(ls_faces_on_face4,:)] ; %#ok<AGROW>
    T_E_beam_face5_8 = [T_E_beam_face5_8 ; nodes_faces_8(ls_faces_on_face5,:)] ; %#ok<AGROW>
    T_E_beam_face6_8 = [T_E_beam_face6_8 ; nodes_faces_8(ls_faces_on_face6,:)] ; %#ok<AGROW>
                   
end

T_E_beam_faces = cell(6,2) ;
T_E_beam_faces{1,1} = T_E_beam_face1_8 ;
T_E_beam_faces{1,2} = T_E_beam_face1_20 ;
T_E_beam_faces{2,1} = T_E_beam_face2_8 ;
T_E_beam_faces{2,2} = T_E_beam_face2_20 ;
T_E_beam_faces{3,1} = T_E_beam_face3_8 ;
T_E_beam_faces{3,2} = T_E_beam_face3_20 ;
T_E_beam_faces{4,1} = T_E_beam_face4_8 ;
T_E_beam_faces{4,2} = T_E_beam_face4_20 ;
T_E_beam_faces{5,1} = T_E_beam_face5_8 ;
T_E_beam_faces{5,2} = T_E_beam_face5_20 ;
T_E_beam_faces{6,1} = T_E_beam_face6_8 ;
T_E_beam_faces{6,2} = T_E_beam_face6_20 ;

% Final 2D connectivity tables

T_E_surface = cell(2,1) ;
T_E_surface{1,1} = T_E_beam_faces ;
T_E_surface{2,1} = T_E_fluid_faces ;

% Tables of degrees of freedom for the 2D faces

T_DOF_beam_face1_20 = make_T_DOF(T_E_beam_face1_20) ;
T_DOF_beam_face1_8  = make_T_DOF(T_E_beam_face1_8) ;
T_DOF_beam_face2_20 = make_T_DOF(T_E_beam_face2_20) ;
T_DOF_beam_face2_8  = make_T_DOF(T_E_beam_face2_8) ;
T_DOF_beam_face3_20 = make_T_DOF(T_E_beam_face3_20) ;
T_DOF_beam_face3_8  = make_T_DOF(T_E_beam_face3_8) ;
T_DOF_beam_face4_20 = make_T_DOF(T_E_beam_face4_20) ;
T_DOF_beam_face4_8  = make_T_DOF(T_E_beam_face4_8) ;
T_DOF_beam_face5_20 = make_T_DOF(T_E_beam_face5_20) ;
T_DOF_beam_face5_8  = make_T_DOF(T_E_beam_face5_8) ;
T_DOF_beam_face6_20 = make_T_DOF(T_E_beam_face6_20) ;
T_DOF_beam_face6_8  = make_T_DOF(T_E_beam_face6_8) ;

T_DOF_fluid_face1_20 = make_T_DOF(T_E_fluid_face1_20) ;
T_DOF_fluid_face1_8  = make_T_DOF(T_E_fluid_face1_8) ;
T_DOF_fluid_face2_20 = make_T_DOF(T_E_fluid_face2_20) ;
T_DOF_fluid_face2_8  = make_T_DOF(T_E_fluid_face2_8) ;
T_DOF_fluid_face3_20 = make_T_DOF(T_E_fluid_face3_20) ;
T_DOF_fluid_face3_8  = make_T_DOF(T_E_fluid_face3_8) ;
T_DOF_fluid_face4_20 = make_T_DOF(T_E_fluid_face4_20) ;
T_DOF_fluid_face4_8  = make_T_DOF(T_E_fluid_face4_8) ;
T_DOF_fluid_face5_20 = make_T_DOF(T_E_fluid_face5_20) ;
T_DOF_fluid_face5_8  = make_T_DOF(T_E_fluid_face5_8) ;
T_DOF_fluid_face6_20 = make_T_DOF(T_E_fluid_face6_20) ;
T_DOF_fluid_face6_8  = make_T_DOF(T_E_fluid_face6_8) ;

T_DOF_fsi_20 = make_T_DOF(T_E_fsi_20) ;
T_DOF_fsi_8  = make_T_DOF(T_E_fsi_8) ;

T_DOF_fsi = cell(2,1) ;
T_DOF_fsi{1} = T_DOF_fsi_8 ;
T_DOF_fsi{2} = T_DOF_fsi_20 ;

T_DOF_beam_faces = cell(6,2) ;
T_DOF_beam_faces{1,1} = T_DOF_beam_face1_8 ;
T_DOF_beam_faces{1,2} = T_DOF_beam_face1_20 ;
T_DOF_beam_faces{2,1} = T_DOF_beam_face2_8 ;
T_DOF_beam_faces{2,2} = T_DOF_beam_face2_20 ;
T_DOF_beam_faces{3,1} = T_DOF_beam_face3_8 ;
T_DOF_beam_faces{3,2} = T_DOF_beam_face3_20 ;
T_DOF_beam_faces{4,1} = T_DOF_beam_face4_8 ;
T_DOF_beam_faces{4,2} = T_DOF_beam_face4_20 ;
T_DOF_beam_faces{5,1} = T_DOF_beam_face5_8 ;
T_DOF_beam_faces{5,2} = T_DOF_beam_face5_20 ;
T_DOF_beam_faces{6,1} = T_DOF_beam_face6_8 ;
T_DOF_beam_faces{6,2} = T_DOF_beam_face6_20 ;

T_DOF_fluid_faces = cell(6,2) ;
T_DOF_fluid_faces{1,1} = T_DOF_fluid_face1_8 ;
T_DOF_fluid_faces{1,2} = T_DOF_fluid_face1_20 ;
T_DOF_fluid_faces{2,1} = T_DOF_fluid_face2_8 ;
T_DOF_fluid_faces{2,2} = T_DOF_fluid_face2_20 ;
T_DOF_fluid_faces{3,1} = T_DOF_fluid_face3_8 ;
T_DOF_fluid_faces{3,2} = T_DOF_fluid_face3_20 ;
T_DOF_fluid_faces{4,1} = T_DOF_fluid_face4_8 ;
T_DOF_fluid_faces{4,2} = T_DOF_fluid_face4_20 ;
T_DOF_fluid_faces{5,1} = T_DOF_fluid_face5_8 ;
T_DOF_fluid_faces{5,2} = T_DOF_fluid_face5_20 ;
T_DOF_fluid_faces{6,1} = T_DOF_fluid_face6_8 ;
T_DOF_fluid_faces{6,2} = T_DOF_fluid_face6_20 ;

T_DOF_surface = cell(2,1) ;
T_DOF_surface{1} = T_DOF_beam_faces ;
T_DOF_surface{2} = T_DOF_fluid_faces ;
               
% Assembling the fsimesh structure

nodelist = cell(2,2) ;
nodelist{1,1} = nodes_beam_faces_8 ;
nodelist{1,2} = nodes_beam_faces_20 ;
nodelist{2,1} = nodes_fluid_faces_8 ;
nodelist{2,2} = nodes_fluid_faces_20 ;

T_X = cell(2,1) ;
T_X{1} = T_X_8 ;
T_X{2} = T_X_20 ;

fsimesh = struct ;

fsimesh.T_X = T_X ;
fsimesh.nodelist = nodelist ;
fsimesh.T_E_volume = T_E_volume ;
fsimesh.T_E_surface = T_E_surface ;
fsimesh.T_E_fsi = T_E_fsi ;
fsimesh.T_DOF_volume = T_DOF_volume ;
fsimesh.T_DOF_surface = T_DOF_surface ;
fsimesh.T_DOF_fsi = T_DOF_fsi ;

% Saving fsimesh to a .mat file

if isave == 1
    save('testcase_fsimesh.mat','fsimesh') ;
end