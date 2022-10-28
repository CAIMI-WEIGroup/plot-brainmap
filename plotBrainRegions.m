function yPlotRegionMap(values, regions, atlas, varargin)
% ======================== Function description ===========================
% This function plot brain map for DK atlas regions or economo
% Input
%       values: N by 1 vectors of values to be ploted. N: number of regions 
%       regions: N by 1 cell of region descriptions
%       atlas: string of the atlas name, e.g., 'lausanne120'
%
% Optional input
%       limits: [minimum value, maximum value] to be plotted
%           Default: [min(values), max(values)]
%       subject: string of the subject name (folder name of the FS output)
%           Default: fsaverage
%       subjects_dir: path to the directory of the subject folder.
%           Default: SUBJECTS_DIR in bash
%       type: colormap type of cbrewer. e.g., 'div', 'seq', 'qual'
%           Default: 'seq'
%       color: colormap color of cbrewer. e.g., 'RdBu', 'Blues'
%           Default: 'Blues'
%       background: color of the background. e.g., 'k', 'w'
%           Default: 'w'
%       surf: string of the surface type. e.g., 'pial', 'white', 'inflated'
%           Default: 'pial'
%       out_dir: string of the outputfile.
%           default: fig.png under the current folder
% by Yongbin Wei, Jul 2020
% =========================================================================

% ======================== Function starts here ===========================

while ~isempty(varargin)
    if numel(varargin) == 1
        error('lscatter:missing_option', ...
        'Optional arguments must come in pairs.');
    end 
    
    switch varargin{1}
        case 'limits'
            limits = varargin{2};
        case 'color'
            color = varargin{2};
        case 'type'
            type = varargin{2};
        case 'subjects_dir'
            subjects_dir = varargin{2};
        case 'subject'
            subject = varargin{2};
        case 'surf'
            surf = varargin{2};
        case 'background'
            background = varargin{2};
        case 'out_dir'
            out_dir = varargin{2};
    end
    varargin(1:2) = [];  
end

if ~exist('limits','var')
   limits = [min(values), max(values)]; 
else
   if limits(1) > min(values) || limits(2) < max(values)
       error('## Wrong arguments: limits are in the range of values');
       return
   end
end

if ~exist('color', 'var')
    color = 'Blues';
end

if ~exist('type', 'var')
    type = 'seq';
end

if ~exist('surf', 'var')
    surf = 'pial';
end

if ~exist('background', 'var')
    background = 'w';
end

if ~exist('subject', 'var')
    subject = 'fsaverage';
end

if ~exist('subjects_dir', 'var')
    subjects_dir = getenv('SUBJECTS_DIR');
    disp('## SUBJECT_DIR');
    disp(subjects_dir);
end

if ~exist('out_dir', 'var')
    out_dir = './fig.png';
end

datapath = fullfile(subjects_dir, subject);
annotPath = {fullfile(datapath, 'label', ['lh.', atlas, '.annot']); ...
    fullfile(datapath, 'label', ['rh.', atlas, '.annot'])};
disp('## Annotation files: ')
disp(annotPath);
if ~exist(annotPath{1}, 'file') || ~exist(annotPath{2}, 'file')
    error('## Annotation files do not exist');
    return
end

surfPath = {fullfile(datapath, 'surf', ['lh.', surf]); ...
    fullfile(datapath, 'surf', ['rh.', surf])};
disp('## Surface files: ')
disp(surfPath);

% val2rgb
if size(values, 1) == 1
    warning('## Transpose VALUES as ths number of rows is 1');
    values = values';    
end
values = [values; limits'];
colors = cbrewer(type, color, 100);
colors = colors(10:90, :);
CB = squeeze(real2rgb(values, colors));
CB = CB(1:end-2, :);

h = plotSurface4(surfPath, annotPath, regions, ...
    CB, background, colors, limits);

% save figure
saveas(h, out_dir, 'png');
disp('## Save figure as:')
disp(out_dir);
end


function h = plotSurface4(freesurferSurface, freesurferAnnotation, ...
    regions, colorMatrix, background, colors, limits)

    brain = 0.85; % 0~1, brain background
    figure('units', 'normalized');
    set(gcf,'Color', background)
    set(gcf, 'InvertHardcopy', 'off');
    tiledlayout(5, 4, 'TileSpacing', 'Compact', 'Padding', 'compact');

    % Left hemishpere
    [~, label, colortable] = read_annotation(...
        freesurferAnnotation{1});

    II = contains(regions, 'ctx-lh-') | ...
        contains(regions, '_LH') | contains(regions, '_L');

    if nnz(II) ~= 0
        regions_lh = regions(II);
        colorMatrix_lh = [colorMatrix(II, :); brain,brain,brain];
        colorMatrix_lh = floor(colorMatrix_lh*255);

        % adjust region names
        if contains(regions_lh{1}, 'ctx-lh')
            for ii = 1:numel(regions_lh)
                tmp = regions_lh{ii};
                regions_lh{ii} = tmp(numel('ctx-lh-')+1 : end);
            end
        end

        % adjust color table
        [~, J] = ismember(colortable.struct_names, regions_lh);
        J(J == 0) = size(colorMatrix_lh, 1);
        colortable.table(:, 1:3) = colorMatrix_lh(J, :);
        oldLabels = colortable.table(:, 5);
        newLabels = colortable.table(:, 1) + ...
            colortable.table(:, 2)*2^8 + ...
            colortable.table(:, 3)*2^16 + ...
            colortable.table(:, 4)*2^24;
        colortable.table(:, 5) = newLabels;
        [I, J] = ismember(label, oldLabels);
        label(I) = newLabels(J(I));
    else
        colortable.table(:, 1:3) = ones(size(colortable.table, 1), 3)...
            * brain * 255;
    end
    
    % read surface file
    [vertices, faces] = read_surf(freesurferSurface{1});
    v(:,3) = vertices(:,1) - nanmean(vertices(:,1));
    v(:,1) = vertices(:,2) - nanmean(vertices(:,2));
    v(:,2) = vertices(:,3) - nanmean(vertices(:,3));
    v(:,3) = -v(:,3);    
    v(:,1) = -v(:,1);
        
    [~, colorIndices] = ismember(label, colortable.table(:, 5));
    
    % lateral view
    nexttile([2,2])
    h1 = patch('Vertices', v, 'Faces', faces + 1, ...
        'FaceVertexCData', colorIndices, 'CDataMapping', 'direct', ...
        'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
        'EdgeColor', 'none');
    colormap(gca, colortable.table(:, 1:3)/255);
    axis off;
    axis equal;
    view([0 90]);
    lightangle(0, 90)
    material dull;

    % medial view
    nexttile([2,2])
    v(:,3) = -v(:,3);    
    v(:,1) = -v(:,1);
    h2 = patch('Vertices', v, 'Faces', faces + 1, ...
        'FaceVertexCData', colorIndices, 'CDataMapping', 'direct', ...
        'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
        'EdgeColor', 'none');
    colormap(gca, colortable.table(:, 1:3)/255);
    axis off;
    view([0 90]);
    lightangle(0, 90)
    axis equal;
    material dull;    
    
    % Right hemishpere
    [~, label, colortable] = read_annotation(freesurferAnnotation{2});

    II = contains(regions, 'ctx-rh-') | ...
        contains(regions, '_RH') | contains(regions, '_R');

    if nnz(II) ~= 0
        regions_rh = regions(II);
        colorMatrix_rh = [colorMatrix(II, :); brain,brain,brain];
        colorMatrix_rh = floor(colorMatrix_rh*255);

        % adjust region names
        if contains(regions_rh{1}, 'ctx-rh')
            for ii = 1:numel(regions_rh)
                tmp = regions_rh{ii};
                regions_rh{ii} = tmp(numel('ctx-rh-')+1 : end);
            end
        end

        % adjust color table
        [~, J] = ismember(colortable.struct_names, regions_rh);
        J(J == 0) = size(colorMatrix_rh, 1);
        colortable.table(:, 1:3) = colorMatrix_rh(J, :);
        oldLabels = colortable.table(:, 5);
        newLabels = colortable.table(:, 1) + colortable.table(:, 2)*2^8 + ...
            colortable.table(:, 3)*2^16 + colortable.table(:, 4)*2^24;
        colortable.table(:, 5) = newLabels;
        [I, J] = ismember(label, oldLabels);
        label(I) = newLabels(J(I));
    else
        colortable.table(:, 1:3) = ones(size(colortable.table, 1), 3) ...
            * brain * 255;
    end
    
    % read surface file
    [vertices_rh, faces_rh] = read_surf(freesurferSurface{2});
    v(:,3) = vertices_rh(:,1) - nanmean(vertices_rh(:,1));
    v(:,1) = vertices_rh(:,2) - nanmean(vertices_rh(:,2));
    v(:,2) = vertices_rh(:,3) - nanmean(vertices_rh(:,3));

    [~, colorIndices] = ismember(label, colortable.table(:, 5));

    % lateral
    nexttile([2,2])
    h3 = patch('Vertices', v, 'Faces', faces_rh + 1, ...
        'FaceVertexCData', colorIndices, 'CDataMapping', 'direct', ...
        'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
        'EdgeColor', 'none');
    colormap(gca, colortable.table(:, 1:3)/255);
    axis equal; 
    axis off;
    view([0 90]);
    lightangle(0, 90)
    material dull;

    % medial
    v(:,3) = -v(:,3);    
    v(:,1) = -v(:,1);
    nexttile([2,2])
    h4 = patch('Vertices', v, 'Faces', faces_rh + 1, ...
        'FaceVertexCData', colorIndices, 'CDataMapping', 'direct', ...
        'FaceColor', 'flat', 'FaceLighting', 'gouraud', ...
        'EdgeColor', 'none');
    colormap(gca, colortable.table(:, 1:3)/255);
    axis equal; 
    axis off;
    view([0 90]);
    lightangle(0, 90)
    material dull;    
 
    % colorbar
    nexttile
    axis off;

    nexttile([1,2])
    M = nan(3, 100);
    M(2, 20:80) = 20:80;
    imagesc(M);
    hold on
    colors = [1, 1, 1; colors];
    colormap(gca, colors);
    axis off;
    text(10, 2, num2str(round(limits(1), 3)));
    text(90, 2, num2str(round(limits(2), 3)));
    
    h = gcf;
    return
end

