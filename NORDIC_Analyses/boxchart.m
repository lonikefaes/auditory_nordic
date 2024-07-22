function h = boxchart(varargin)
% BOXCHART Creates box and whisker plot.
%
%   BOXCHART(YDATA) creates a box chart, or box plot, for each column of
%   the matrix YDATA. If ydata is a vector, then boxchart returns a single
%   box chart. Each box chart displays the following information: the
%   median, first and third quartiles, outliers (computed using the
%   interquartile range), and minimum and maximum values that are not
%   outliers.
%
%   BOXCHART(XGROUPDATA,YDATA) groups the data in vector ydata according to
%   the unique values in XGROUPDATA and plots each group of data as a
%   separate box chart. XGROUPDATA determines the position of each box
%   chart along the x-axis. XGROUPDATA must be a vector of the same length
%   as YDATA.
%
%   BOXCHART(___,Name,Value) specifies additional chart options using one
%   or more name-value pair arguments. Specify the name-value pair
%   arguments after all other input arguments.
%
%   BOXCHART(AX,___) plots into the axes specified by AX instead of into
%   the current axes (gca). The option AX can precede any of the input
%   argument combinations in the previous syntaxes.
% 
%   H = BOXCHART(___) returns a BoxChart object. Use H to set properties of
%   the box charts after creating them.
%
%   H = BOXCHART(___,'GroupByColor',cgroupdata) uses the data specified by
%   cgroupdata for grouping ydata. The function returns a vector of
%   BoxChart objects, one for each distinct category in cgroupdata.
%
%   Example: 
%
%   tiledlayout(2,1);
%   nexttile
%   boxchart(rand(10,5))
%   nexttile
%   boxchart(randi(2,50,1),rand(50,1))
%
%   See also HISTOGRAM, PLOT, BAR, BARH, BAR3, BAR3H.
    
%   Copyright 2019-2020 The MathWorks, Inc. 

% Pass inputs to axescheck. axescheck is not guaranteed to return an axes
% handle but may return any valid graphics handle.
[parent,args] = axescheck(varargin{:});
if isempty(parent) && ~isempty(args) && isa(args{1},'matlab.graphics.Graphics')
    parent = args{1};
    args = args(2:end);
end

% Parent may not be an axes handle. If it is empty or is an axes, pass to
% newplot. Otherwise, find the axes in the graphics tree via ancestor.
if isempty(parent) || isgraphics(parent,'axes')
    parent = newplot(parent);
    cax = parent;
else
    cax = ancestor(parent,'axes');
end

% Ensure that supplied parent is valid
if isempty(cax)
    shortname = regexprep(class(parent), '.*\.', '');
    error(message('MATLAB:handle_graphics:exceptions:HandleGraphicsException',...
        getString(message("MATLAB:HandleGraphics:hgError", 'BoxChart', shortname))));
end

% Obtain y from the input arguments
if numel(args)<1
    error(message('MATLAB:minrhs'));
end

% Obtain xgroupdata and ydata from args
xdatamode = 'auto';
ydata = args{1};
if numel(args)>1 && ~(ischar(args{2}) || isstring(args{2}))
    xgroupdata = ydata;
    ydata = args{2};
    
    % xgroupdata may only be supplied when ydata is a vector
    if ~isvector(ydata)
        error(message('MATLAB:graphics:boxchart:NoXGDataWhenY2D'));
    end
    xdatamode = 'manual';
    args = args(3:end);
else % args == 1
    % If xgroupdata is empty create a vector of ones
    if isvector(ydata)
        xgroupdata = ones(size(ydata));
    else
        xgroupdata = 1:size(ydata,2);
    end
    
    % The default x-ruler is categorical
    xgroupdata = categorical(xgroupdata);
    
    args = args(2:end);
end
validateattributes(ydata,{'numeric'},{'2d','real'},mfilename,'ydata');
validateattributes(xgroupdata,{'numeric','categorical'},{'vector','real'},mfilename,'xgroupdata');

% Check that the sizes of xgroupdata and ydata are consistent
if isvector(ydata) && (numel(xgroupdata) ~= numel(ydata))
    error(message('MATLAB:graphics:boxchart:BadXVectorY'));
end

% Configure x-axis and convert xgroupdata to numeric
matlab.graphics.internal.configureAxes(cax,xgroupdata,ydata);

% Get the automatic color
autoColor = true;
colorPropF = 'BoxFaceColor_I';
for i = 1:2:numel(args)
    if startsWith('BoxFaceColor',args{i},'IgnoreCase',true) && ~strcmp(args{i+1},'flat')
        autoColor = false;
        colorPropF = 'BoxFaceColor';
    end
end
[~,fc] = matlab.graphics.chart.internal.nextstyle(cax,autoColor,false,true);

autoColor = true;
colorPropM = 'MarkerColor_I';
for i = 1:2:numel(args)
    if startsWith('MarkerColor',args{i},'IgnoreCase',true) && ~strcmp(args{i+1},'flat')
        autoColor = false;
        colorPropM = 'MarkerColor';
    end
end
[~,mc] = matlab.graphics.chart.internal.nextstyle(cax,autoColor,false,true);

ngrp = 1;
gnum = ones(size(xgroupdata));
colGrpIdx = [];
for i = 1:2:numel(args)
    if startsWith('GroupByColor',args{i},'IgnoreCase',true) && (i+1 <= numel(args))
        % Obtain color group data
        grp = args{i+1};
        colGrpIdx = [colGrpIdx,i];
        % Validate GroupByColor
        validateattributes(grp,{'numeric','categorical','logical',...
                'char','string','cell'},{'real','nonsparse','vector'},'','GroupByColor');
        
        % GroupByColor is only supported when ydata is a vector
        if ~isvector(ydata)
            error(message('MATLAB:graphics:boxchart:NoColGDataWhenY2D'));
        end
        
        % Check that the sizes of grp and ydata are consistent
        if isvector(ydata) && (numel(grp) ~= numel(ydata))
            error(message('MATLAB:graphics:boxchart:BadColGroupVectorY'));
        end
        
        % Obtain group indices and names
        [gnum,gnames] = findgroups(grp);
        gnames = gnames(:);
        ngrp = numel(gnames);
    end
end
args([colGrpIdx,colGrpIdx+1]) = [];

if isvector(ydata)
    grpargs = {};
    H = gobjects(ngrp,1);
    for idx = 1:ngrp
        ind = gnum == idx;
        x = xgroupdata(ind);
        y = ydata(ind);
        if ~isempty(colGrpIdx)
            % Display name is used to populate the legend's text
            grpargs = {'DisplayName',string(gnames(idx,:)),...
                'NumColorGroups', ngrp, 'GroupByColorMode','manual'};
        end

        % Call the class constructor
        H(idx) = matlab.graphics.chart.primitive.BoxChart('Parent', cax,...
            'XData', x, 'YData', y, colorPropF, fc, colorPropM, mc, ...
            'XDataMode', xdatamode, 'PeerID', idx, args{:}, grpargs{:});
        H(idx).assignSeriesIndex();
    end
    
    % If the user passed colgroupdata, ensure that each boxchart stores
    % handles to its peers
    if ngrp > 1
        for idx = 1:ngrp
            H(idx).BoxPeers = H(idx ~= 1:ngrp);
        end
    end
else
    H = matlab.graphics.chart.primitive.BoxChart('Parent', cax,...
            'XData', xgroupdata, 'YData', ydata, colorPropF, fc, ...
            colorPropM, mc, 'XDataMode', xdatamode, 'PeerID', 1, args{:});
    H.assignSeriesIndex();
end

% Return handle only if the user asks for it
if nargout > 0
    h = H;
end
end