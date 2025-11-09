function [group_rp,group_re,group_spad,canopy_idx_all]=canopy_model_spad(group_se,group_sp,plant_spad,facet_idx,plant_num,row_num,plant_dis,row_dis,all_idx,direction)
%% =========================================================================
% Function: canopy_model_spad
% Purpose : Generate a 3D canopy model from multiple plants and compute SPAD values
%  Author       : Fusang Liu
%  Affiliation  : Shanghai Institute of Plant Physiology and Ecology, Chinese Academy of Sciences 
%                 Wageningen University & Research (WUR)
% Inputs  :
%   group_se    - cell array of plant mesh faces (triangles)
%   group_sp    - cell array of plant point coordinates (Nx3)
%   plant_spad  - cell array of SPAD values for each plant facet
%   facet_idx   - cell array of facet indices for each plant
%   plant_num   - number of plants per row
%   row_num     - number of rows in the canopy
%   plant_dis   - distance between plants in a row (X-axis)
%   row_dis     - distance between rows (Y-axis)
%   all_idx     - matrix of plant indices (plant arrangement)
%   direction   - matrix of rotation angles (0-1, scaled to 360°)
% Outputs :
%   group_rp        - all rotated and translated plant point coordinates in the canopy
%   group_re        - all rotated and translated plant mesh faces in the canopy
%   group_spad      - SPAD values for all canopy facets
%   canopy_idx_all  - mapping of plant index to facet index in the canopy
%% =========================================================================
facet_idx_all=[];
plant_num_all=[];
for i=1:plant_num
    for j=1:row_num
plant_idx=all_idx(i,j);        
facet_idx_one=facet_idx{plant_idx};
plant_num_one=ones(length(facet_idx{plant_idx}),1).*plant_idx;
facet_idx_all=[facet_idx_all;facet_idx_one];
plant_num_all=[plant_num_all;plant_num_one];
    end
end
canopy_idx_all=[plant_num_all,facet_idx_all];

%% actual canopy simulation
for i=1:plant_num
    for j=1:row_num
plant_idx=all_idx(i,j);
group_tr=coordinate_rotate(group_sp{plant_idx},direction(i,j)*360,[0,0,0],3);       
group_cell{i,j}(:,1)=group_tr(:,1)+i*plant_dis;%%
group_cell{i,j}(:,2)=group_tr(:,2)+j*row_dis;%%加上行距
group_cell{i,j}(:,3)=group_tr(:,3);
    end
end
group_rp=[];

for i=1:plant_num
    for j=1:row_num
group_rp1=group_cell{i,j};
group_rp=[group_rp;group_rp1];
    end
end
group_re=[];
n=0;
group_spad=[];

for i=1:plant_num
    for j=1:row_num
plant_idx=all_idx(i,j); 
l=length(group_sp{plant_idx});
group_e1(:,1)=group_se{plant_idx}(:,1)+n;
group_e1(:,2)=group_se{plant_idx}(:,2)+n;
group_e1(:,3)=group_se{plant_idx}(:,3)+n;
n=n+l;
group_re=[group_re;group_e1]; 
group_e1=[];
group_spad=[group_spad;plant_spad{plant_idx}];
    end
end

end