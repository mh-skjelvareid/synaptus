function polFocusedImage(im,phi,rr,phiGrid,rGrid,varargin)
% polFocusedImage - create polar plot of focused image
%
%   Usage:
%   polMultiLayerImage(im,phi,rr,phiGrid,rGrid)
%
%   Input parameters:
%   im      -   NxM image
%   phi     -   Mx1 vector of angles [rad]
%   rr      -   Nx1 vector of ranges [m]
%   phiGrid -   vector of angle coordinates for vertical grid lines
%   rGrid   -   vector of range coordinates for lateral grid lines
%
%   The function will create a polar plot where the standard rectangular
%   axes have been replaced with axes following the contour of the data.
%   The coordinates of end points are marked, and additional "tick marks"
%   can be made by specifying grid coordinates rGrid and phiGrid.
%
%   2012-05-10  Martin H. Skjelvareid

%% Parse optional input arguments
p = inputParser;                                % Create input parser
p.addParamValue('rInterface',[]);               % Ranges of interfaces
p.addParamValue('fontSize',14);
p.parse(varargin{:});                           % Parse param.-value pairs
param = p.Results;                              % Transfer to "param" structure
clear p

%% Constants (make into input parameters?)
gridStyle = ':';
gridColor = 0.5*ones(1,3);

%% Get size of input
nRGrid = length(rGrid);
nPhiGrid = length(phiGrid);

rMin = rr(1);
rMax = rr(end);

%% Plot image
[PHI,RR] = meshgrid(phi,rr);
[XX,YY] = pol2cart(PHI,RR);

% Plot as surface
surf(XX,YY,im,'edgecolor','none')
% hpc = pcolor(XX,YY,im);
% set(hpc,'edgecolor','none')
view(90,90)
hold on

%% Find range of x-y coordinates
[xbi,ybi] = pol2cart(phi,rr(1));            % x-y coordinates of "inner circle"
[xbo,ybo] = pol2cart(phi,rr(end));          % x-y coordinates of "outer circle"

xRange = max([xbi(:);xbo(:)])-min([xbi(:);xbo(:)]);     % Test x extreme values
yRange = max([ybi(:);ybo(:)])-min([ybi(:);ybo(:)]);     % Test y extreme values

%% Add borders
[xLeft yLeft] = pol2cart(phi(1),[rMin rMax]);
[xRight yRight] = pol2cart(phi(end),[rMin rMax]);
[xTop yTop] = pol2cart(phi,rMin);
[xBottom yBottom] = pol2cart(phi,rMax);

plot([xBottom(:) xTop(:)],[yBottom(:) yTop(:)],'k')
plot([xLeft(:) xRight(:)],[yLeft(:) yRight(:)],'k')

%% Add interfaces (optional)
if ~isempty(param.rInterface)
    for ii = 1:length(param.rInterface)
        [xTick yTick] = pol2cart(phi,param.rInterface(ii));
        plot(xTick,yTick,'k--');
    end
end

%% Add grid lines
for ii = 1:nPhiGrid
    [xTick yTick] = pol2cart(phiGrid(ii)*[1 1],[rMin rMax]);
    plot(xTick,yTick,gridStyle,'Color',gridColor,'LineWidth',1.5);
end

for ii = 1:nRGrid
    [xTick yTick] = pol2cart(phi,rGrid(ii));
    plot(xTick,yTick,gridStyle,'Color',gridColor,'LineWidth',1.5);
end

%% Add grid text labels (TODO: make placement adjustable?)
[xTextR,yTextR] = pol2cart(phi(1),rGrid);
text(xTextR,yTextR-0.08*yRange,num2str(rGrid(:)*1e3),...
    'FontSize',param.fontSize)

[xTextPhi,yTextPhi] = pol2cart(phiGrid,rMax);
text(xTextPhi+0.05*xRange,yTextPhi-0.03*yRange,num2str(phiGrid(:)*(180/pi)),...
    'FontSize',param.fontSize)

%% Add range label (TODO: make placement adjustable?)
[xRLabel, yRLabel] = pol2cart(phi(1),(rMin+rMax)/2);
text(xRLabel,yRLabel-0.14*yRange,'r [mm]','rotation',90+phi(1)*(180/pi),...
    'FontSize',param.fontSize)

%% Add angle label (TODO: make placement adjustable?)
[xPhiLabel, yPhiLabel] = pol2cart((phi(1)+phi(end))/2,rMax);
text(xPhiLabel+0.1*xRange, yPhiLabel,'\phi [deg]','FontSize',param.fontSize)

%% Set axes
axis tight
% axis image
set(gca,'Visible','off')
