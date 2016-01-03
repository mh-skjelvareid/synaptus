function polRawDataImage(ptp,tt,phi,tGrid,phiGrid,R1,c1,varargin)
% polRawDataImage - polar plot of ultrasonic raw data acquired on circular arc
%
%   Usage:
%   polRawDataImage(ptp,tt,phi,tGrid,phiGrid,R1,c1)
%
%   Input parameters:
%   ptp     -   2D matrix of ultrasonic data, time along vertical axis,
%               angle along lateral axis
%   tt      -   time axis vector [s]
%   phi     -   angle axis vector [rad]
%   tGrid   -   vector of time coordinates for lateral grid lines
%   phiGrid -   vector of angle coordinates for vertical grid lines
%   R1      -   scan radius [m]
%   c1      -   sound speed in first layer (if several layers)
%
%   Optional parameter-value pairs:
%   'tInterface'    -   t coordinate for interface, marked with dashed line
%   'tLabelOffset'  -   offset of time axis label, normalized with width of
%                       image. Default: 0.15.
%
%   The function will create a polar plot where the standard rectangular
%   axes have been replaced with axes following the contour of the data.
%   The coordinates of end points are marked, and additional "tick marks"
%   can be made by specifying grid coordinates tGrid and phiGrid.

%% Parse optional input arguments
p = inputParser;                                    % Create input parser
p.addParamValue('tInterface',[]);                   % T coordinate layer interf.
p.addParamValue('tLabelOffset',0.15);               % time axis label offset
p.addParamValue('fontSize',14);
p.parse(varargin{:});                               % Parse param.-value pairs
param = p.Results;                                  % Transfer to "param"
clear p

%% Constants (make into input parameters?)
gridStyle = ':';
gridColor = 0.5*ones(1,3);

%% Get size of input
[nT,nPhi] = size(ptp);
nTGrid = length(tGrid);
nPhiGrid = length(phiGrid);

%% % Create time axis referenced to origo
ttRef0 = tt + 2*(R1/c1);
tGridRef0 = tGrid + 2*(R1/c1);

%% Create x and y coordinates for polar image
[PHI,TT] = meshgrid(phi,ttRef0);
[XX,YY] = pol2cart(PHI,TT);

xRange = max(XX(end,:))-min(XX(1,:));
yRange = max(YY(:,end))-min(YY(:,1));

%% Plot as surface
surf(XX,YY,zeros(size(XX)),ptp,'edgecolor','none')
view(90,90)
hold on

%% Add borders
[xLeft yLeft] = pol2cart(phi(1),[ttRef0(1) ttRef0(end)]);
[xRight yRight] = pol2cart(phi(end),[ttRef0(1) ttRef0(end)]);

[xTop yTop] = pol2cart(phi,ttRef0(1));
[xBottom yBottom] = pol2cart(phi,ttRef0(end));

plot([xBottom(:) xTop(:)],[yBottom(:) yTop(:)],'k')
plot([xLeft(:) xRight(:)],[yLeft(:) yRight(:)],'k')

%% Add interfaces (optional)
if ~isempty(param.tInterface)
    for ii = 1:length(param.tInterface)
        [xTick yTick] = pol2cart(phi,param.tInterface(ii)+2*(R1/c1));
        plot(xTick,yTick,'k--');
    end
end

%% Add grid lines
for ii = 1:nPhiGrid
    [xTick yTick] = pol2cart(phiGrid(ii)*[1 1],ttRef0([1 end]));
    plot(xTick,yTick,gridStyle,'Color',gridColor);
end

for ii = 1:nTGrid
    [xTick yTick] = pol2cart(phi,tGrid(ii)+2*(R1/c1));
    plot(xTick,yTick,gridStyle,'Color',gridColor);
end

%% Add grid text labels
[xTextT,yTextT] = pol2cart(phi(1),tGridRef0);
% text(xTextT,yTextT-0.1*yRange,num2str(tGrid(:)*1e6))
text(xTextT,yTextT-0.6*param.tLabelOffset*yRange,num2str(tGrid(:)*1e6),...
    'FontSize',param.fontSize)

[xTextPhi,yTextPhi] = pol2cart(phiGrid,ttRef0(end));
text(xTextPhi+0.05*xRange,yTextPhi-0.03*yRange,num2str(phiGrid(:)*(180/pi)),...
    'FontSize',param.fontSize)

%% Add time label
[xTLabel, yTLabel] = pol2cart(phi(1),ttRef0(round(nT/2)));
text(xTLabel,yTLabel-param.tLabelOffset*yRange,'t [us]',...
    'rotation',90+phi(1)*(180/pi),'FontSize',param.fontSize)

%% Add angle label
[xPhiLabel, yPhiLabel] = pol2cart(phi(round(nPhi/2)),ttRef0(end));
text(xPhiLabel+0.1*xRange, yPhiLabel,'\phi [deg]','FontSize',param.fontSize)

%% Set axes
axis normal
set(gca,'Visible','off')
