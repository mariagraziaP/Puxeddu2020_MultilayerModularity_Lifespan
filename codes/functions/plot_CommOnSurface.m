function plot_CommOnSurface(path,part,lines,colors)

% INPUT: path -> directory of data and annotation
%        part -> vector of communitiy membership (N°nodes x 1)
%        lines -> specify if you want to draw lines separating ROIs
%                   'yes' = you want line
%                   'no' = no lines
%       colors -> matrix with colors used to fill communities (N°Clust x 3) 
%                   (optional, if not specified parula is used)

% OUTPUT: communities sructure plotted on the brain surface 

if nargin<4
    colors = parula;
end

% load data 
LH_surf_string = fullfile(path,'lh.inflated');
RH_surf_string = fullfile(path,'rh.inflated');

% load annotation
switch lines
    case 'yes'
        LH_annot_string = fullfile(path,'lh.YeoUpsample.annot');
        RH_annot_string = fullfile(path,'rh.YeoUpsample.annot');
    case 'no'
        LH_annot_string = fullfile(path,'lh.yeo17dil.annot');
        RH_annot_string = fullfile(path,'rh.yeo17dil.annot');
end

% read data
surfData = struct();

[surfData.LH.verts,surfData.LH.faces] = read_surf(LH_surf_string) ;
[surfData.RH.verts,surfData.RH.faces] = read_surf(RH_surf_string) ;

surfData.LH.faces = surfData.LH.faces(:,1:3) + 1 ;
surfData.RH.faces = surfData.RH.faces(:,1:3) + 1 ;

% get annotation
annotData = struct();

[annotData.LH.verts,annotData.LH.labs,annotData.LH.ct] = read_annotation(LH_annot_string) ;
[annotData.RH.verts,annotData.RH.labs,annotData.RH.ct] = read_annotation(RH_annot_string) ;

% check consistency between annotation and partition
nodes = length(part);
nrois = size(annotData.LH.ct.table,1);
if nodes~=nrois-1
    display('Error: different number of nodes in annotation and partition')
end

% get weights
weights_unknown = -1 ;

weights_LH = ones(size(annotData.LH.verts,1),1) * weights_unknown;
weights_RH = ones(size(annotData.RH.verts,1),1) * weights_unknown;

% new weights according to communities
new_col = [length(unique(part))+1; part];
annotData.LH.ct.table(:,6) = new_col;   
annotData.RH.ct.table(:,6) = new_col;
for idx = 1:nrois 
    weights_LH(annotData.LH.labs == annotData.LH.ct.table(idx,5)) = annotData.LH.ct.table(idx,6);
    weights_RH(annotData.RH.labs == annotData.RH.ct.table(idx,5)) = annotData.RH.ct.table(idx,6);
end

% viz
figure;

subplot(2,2,1) % right hemisphere lateral view
tmpViz = trisurf(surfData.RH.faces,...
    surfData.RH.verts(:,1),...
    surfData.RH.verts(:,2),...
    surfData.RH.verts(:,3),...
    weights_RH);
set(tmpViz,'EdgeColor','none');
axis equal; axis off
tmpViz.CDataMapping = 'direct' ;
view(90,0)
camlight headlight; material dull; lighting gouraud
title('RH, lat')

subplot(2,2,3) % right hemisphere medial view
tmpViz = trisurf(surfData.RH.faces,...
    surfData.RH.verts(:,1),...
    surfData.RH.verts(:,2),...
    surfData.RH.verts(:,3),...
    weights_RH);
set(tmpViz,'EdgeColor','none');
axis equal; axis off
tmpViz.CDataMapping = 'direct' ;
view(-90,0)
camlight headlight; material dull; lighting gouraud
title('RH, med')

subplot(2,2,2) % left hemisphere lateral view
tmpViz = trisurf(surfData.LH.faces,...
    surfData.LH.verts(:,1),...
    surfData.LH.verts(:,2),...
    surfData.LH.verts(:,3),...
    weights_LH);
set(tmpViz,'EdgeColor','none');
axis equal; axis off
tmpViz.CDataMapping = 'direct' ;
view(-90,0)
camlight headlight; material dull; lighting gouraud
title('LH, lat')

subplot(2,2,4) % left hemisphere medial view
tmpViz = trisurf(surfData.LH.faces,...
    surfData.LH.verts(:,1),...
    surfData.LH.verts(:,2),...
    surfData.LH.verts(:,3),...
    weights_LH);
set(tmpViz,'EdgeColor','none');
axis equal; axis off
tmpViz.CDataMapping = 'direct' ;
view(90,0)
camlight headlight; material dull; lighting gouraud
title('LH, med')


if nargin<4
    colormap([parula(length(unique(part))); 0 0 0])
else
    colormap([colors; 0 0 0])
end