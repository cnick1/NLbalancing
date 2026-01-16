function runExample32_LTI_example(reduction)
%runExample32_LTI_example Runs the 3D example to visualize the linear balancing transformations.
%
%   Usage:  runExample32_LTI_example()
%
%   Inputs:
%       reduction - boolean, whether or not to apply reduction
%
%   Description: This simple 3D examples capture the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated. Consider the
%   linear system from [1]:
%           ẋ₁ = −x₁ + 100 x₃ + u,
%           ẋ₂ = −2 x₂ + 100 x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%            y = x₁ + x₂ + x₃,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to approximate this model
%   with a 2D model. However, x₃ strongly drives the states x₁ and x₂. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states).
%
%   In this script, we use linear balancing to illustrate the presence of a
%   2D model that approximates the input-output behavior of the system.
%
%   References: [1] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 23\n')
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultLineLineWidth',1.5)


if nargin < 1
    reduction = false;
end

%% Get system dynamics
[f, g, h] = getSystem32();
A = full(f{1}); B = full(g{1}); C = full(h{1});

%%  Compute the Gramians
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing Gramians (square-roots):  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
L = lyapchol(A,B); R = lyapchol(A.',C.');

%% Compute the balancing transformation
fprintf(" ~~~~~~~~~~~ Computing transformation:  ~~~~~~~~~~~~ \n")
[U,S,V] = svd(L*R.');

U1 = U(:,1:2); S1 = S(1:2,1:2); V1 = V(:,1:2);
T = S1.^(1/2) \ V1.' * R; Tinv = L.' * U1 / S1.^(1/2);

%% Construct ROM
Ar = T*A*Tinv; Br = T*B; Cr = C*Tinv;

%% Create ss models for FOM and ROM
FOM = ss(A,B,C,0)
ROM = ss(Ar,Br,Cr,0);

%% Plot results
t = 0:.05:10;
[y1,t1] = impulse(FOM,t);
[y2,t2] = impulse(ROM,t);

figure; subplot(2,1,1); hold on;
plot(t1,y1)
plot(t2,y2,'r--')
title('Impulse response'); xlabel('t'); ylabel('y(t)')
legend('FOM','ROM')


[y1,t1] = step(FOM,t);
[y2,t2] = step(ROM,t);

subplot(2,1,2); hold on;
plot(t1,y1)
plot(t2,y2,'r--')
title('Step response'); xlabel('t'); ylabel('y(t)')
legend('FOM','ROM')


figure; hold on;
bode(FOM)
bode(ROM,'r--')
% title('Transfer Function Bode Plot'); xlabel('t'); ylabel('y(t)')

set(findall(gcf, 'type', 'line'), 'LineWidth', 1.5);
legend('FOM','ROM')

%% Plot balanced subspace
[X,Y,Z] = meshgrid(-2:.05:2,-2:.1:2,-2:.1:2);
% v = x.*exp(-x.^2-y.^2-z.^2) + 1;

v = zeros(size(X)); %V2 = 0.1*eye(3); V2(end) = 100; V2(2,2) = 10; V2 = S;

for i=1:size(v,1)
    for j=1:size(v,2)
        for k=1:size(v,3)
            x = [X(i,j,k),Y(i,j,k),Z(i,j,k)].';
            z = (S.^(1/2) \ V.' * R)*x;
            v(i,j,k) = 1/2 * z.' * inv(S) * z;
        end
    end
end

figure;
plot3([0 Tinv(1,1)]/2,[0 Tinv(2,1)]/2,[0 Tinv(3,1)]/2, 'r', 'LineWidth', 2);hold on;
plot3([0 Tinv(1,2)]/2,[0 Tinv(2,2)]/2,[0 Tinv(3,2)]/2, 'r', 'LineWidth', 2);

% X-Y projection (Z = 0)
plot3([0 Tinv(1,1)]/2, [0 0], [0 0], 'k--', 'LineWidth', 1);
plot3([0 Tinv(1,2)]/2, [0 0], [0 0], 'k--', 'LineWidth', 1);

% Y-Z projection (X = 0)
plot3([Tinv(1,1) Tinv(1,1)]/2, [0 Tinv(2,1)]/2, [0 0], 'k--', 'LineWidth', 1);
plot3([Tinv(1,2) Tinv(1,2)]/2, [0 Tinv(2,2)]/2, [0 0], 'k--', 'LineWidth', 1);

% X-Z projection (Y = 0)
plot3([Tinv(1,1) Tinv(1,1)]/2, [Tinv(2,1) Tinv(2,1)]/2, [0 Tinv(3,1)]/2, 'k--', 'LineWidth', 1);
plot3([Tinv(1,2) Tinv(1,2)]/2, [Tinv(2,2) Tinv(2,2)]/2, [0 Tinv(3,2)]/2, 'k--', 'LineWidth', 1);

xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

figure;
pcolor3(X,Y,Z,v)
hold on;
rgbmap('blue','white')
caxis([0.8 1.2])

plot3([0 Tinv(1,1)]/2,[0 Tinv(2,1)]/2,[0 Tinv(3,1)]/2, 'r', 'LineWidth', 2);hold on;
plot3([0 Tinv(1,2)]/2,[0 Tinv(2,2)]/2,[0 Tinv(3,2)]/2, 'r', 'LineWidth', 2);

%% Simulate random noise, add to plot
figure;
pcolor3(X,Y,Z,v)
hold on;
rgbmap('blue','white')
caxis([0.8 1.2])

plot3([0 Tinv(1,1)]/2,[0 Tinv(2,1)]/2,[0 Tinv(3,1)]/2, 'r', 'LineWidth', 2);hold on;
plot3([0 Tinv(1,2)]/2,[0 Tinv(2,2)]/2,[0 Tinv(3,2)]/2, 'r', 'LineWidth', 2);

t = 0:.05:1000;
u = .2*randn(size(t));  % White noise input
[~, ~, x3] = lsim(FOM, u, t);
% Plot singular vectors
plot3(x3(:,1),x3(:,2),x3(:,3),'k:','LineWidth',.5)
% legend('Balanced (controllable) subspace','Response to white noise')
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

%% Simulate random noise, add to plot
figure;
pcolor3(X,Y,Z,v)
hold on;
rgbmap('blue','white')
caxis([0.8 1.2])

plot3([0 Tinv(1,1)]/2,[0 Tinv(2,1)]/2,[0 Tinv(3,1)]/2, 'r', 'LineWidth', 2);hold on;
plot3([0 Tinv(1,2)]/2,[0 Tinv(2,2)]/2,[0 Tinv(3,2)]/2, 'r', 'LineWidth', 2);

% Plot singular vectors
plot3(x3(:,1),x3(:,2),x3(:,3),'k:','LineWidth',.5)
% legend('Balanced (controllable) subspace','Response to white noise')
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
view([1 0 0])

%% Simulate random noise, add to plot
figure;
pcolor3(X,Y,Z,v)
hold on;
rgbmap('blue','white')
caxis([0.8 1.2])

plot3([0 Tinv(1,1)]/2,[0 Tinv(2,1)]/2,[0 Tinv(3,1)]/2, 'r', 'LineWidth', 2);hold on;
plot3([0 Tinv(1,2)]/2,[0 Tinv(2,2)]/2,[0 Tinv(3,2)]/2, 'r', 'LineWidth', 2);

% Plot singular vectors
plot3(x3(:,1),x3(:,2),x3(:,3),'k:','LineWidth',.5)
% legend('Balanced (controllable) subspace','Response to white noise')
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
view(2)

end

function h = pcolor3(varargin)
% pcolor3 plots a 3D data volume as 100 color-scaled semitransparent
% surface planes in each dimension.
%
%% Syntax
%
%  pcolor3(V)
%  pcolor3(X,Y,Z,V)
%  pcolor3(...,'alpha',AlphaValue)
%  pcolor3(...,'edgealpha',EdgeAlphaValue)
%  pcolor3(...,'alphalim',AlphaLimits)
%  pcolor3(...,InterpolationMethod)
%  pcolor3(...,'N',NumberOfSlices)
%  pcolor3(...,'Nx',NumberOfXSlices)
%  pcolor3(...,'Ny',NumberOfYSlices)
%  pcolor3(...,'Nz',NumberOfZSlices)
%  h = pcolor3(...)
%
%% Description
%
% pcolor3(V) plots a field of 3D volume V.
%
% pcolor3(X,Y,Z,V) plots 3D volume V at locations given by X,Y,Z. X, Y, and
% Z can be 3D matrices matching the dimensions of V, or 1D arrays.
%
% pcolor3(...,'alpha',AlphaValue) specifies a volume transparency value between 0
% (completely transparent) and 1 (completely opaque). Default AlphaValue is
% 0.01. This value may seem surprisingly low, but remember that you'll be
% looking through 100 slices--they add up.
%
% pcolor3(...,'edgealpha',EdgeAlphaValue) specifies transparency of sides of
% the volume faces of the volume. An EdgeAlphaValue greater than the volume
% AlphaValue helps define corners and edges, especially in the presence of
% lighting objects. Default EdgeAlphaValue is 0.05.
%
% pcolor3(...,'alphalim',AlphaLimits) scales transparency values with
% values of V. This can help highlight a variable of interest by making
% low V values invisible. AlphaLimits is a two-element array
% corresponding of values in V. If AlphaLimits is 'auto',
% AlphaLimits is taken as [min(V(:)) max(V(:))].
%
%     Tip: If interesting values diverge about an uninteresting mean (e.g.,
%     temperature of 25 is not interesting whereas T = 10 is interesting and T = 40 is also
%     interesting), use 'alphalim',[25 40] and select a colormap that
%     diverges from 25. Although T = 10 is well below the minimum
%     AlphaLimits, 10 and 40 are equidistant from 25 and are therefore given
%     equal opacity.
%
% pcolor3(...,InterpolationMethod) specifies an interpolation method as
%   'linear'  trilinear slice interpolation (default),
%   'cubic'   tricubic slice interpolation,
%   'nearest' nearest-neighbor slice interpolation, or
%   'direct'  plots data directly instead of interpolated slices.
%
% pcolor3(...,'N',NumberOfSlices) specifies a number of slices in each
% direction. Default value is 100. Increasing number of slices can make a
% smoother, higher quality graphic, but may slow performance.
%
% pcolor3(...,'Nx',NumberOfXSlices) specifies a number of slices in the x
% direction. Default value is 100.
%
% pcolor3(...,'Ny',NumberOfYSlices) specifies a number of slices in the y
% direction. Default value is 100.
%
% pcolor3(...,'Nz',NumberOfZSlices) specifies a number of slices in the z
% direction. Default value is 100.
%
% h = pcolor3(...) returns a vector of handles to surface graphics objects.
%
%% Examples (Type showdemo pcolor3_documentation for more examples. )
%
% % Using this sample data:
%   [x,y,z] = meshgrid(-1:.2:3,-2:.25:2,-2:.16:2);
%   v = x.*exp(-x.^2-y.^2-z.^2);
%
% % Plot a simple field:
%
%   pcolor3(v)
%
% % Or specify x,y,z values and set alpha limits:
%
%   pcolor3(x,y,z,v,'alphalim',[0 0.2],'cubic','edgealpha',.1)
%   camlight
%   view(-34,56)
%
%% Author Info
% This function was written by Chad A. Greene of the University of
% Texas at Austin's Institute for Geophysics (UTIG) March 2015.
% http://www.chadagreene.com
%
% See also slice, surf, alpha.
%% Initial input checks:
assert(nargin>0,'pcolor3 requires at least one input.')
assert(isnumeric(varargin{1})==1,'First argument of pcolor3 must be numeric.')

%% Set defaults:
Alpha = 0.01;
EdgeAlpha = 0.05;
nx = 100;
ny = 100;
nz = 100;
InterpolationMethod = 'linear';
setAlphaLim = false;
%% Parse inputs:

% Is input format pcolor3(X,Y,Z,V,...) or simply pcolor3(V,...)?
if nargin>1 && isnumeric(varargin{2})
    assert(nargin>3,'Input error. If the second input to pcolor3 is numeric, inputs are assumed to be in the form pcolor3(X,Y,Z,V,...). You have either entered too few or too many inputs.')
    X = varargin{1};
    Y = varargin{2};
    Z = varargin{3};
    V = varargin{4};
    assert(isnumeric(Z)==1,'The pcolor3 function has interpreted inputs in the form pcolor3(X,Y,Z,V,...), but your third input here is not numeric. I am confused.')
    assert(isnumeric(V)==1,'The pcolor3 function has interpreted inputs in the form pcolor3(X,Y,Z,V,...), but your fourth input here is not numeric. I am confused.')
else
    V = varargin{1};
    [X,Y,Z] = meshgrid(1:size(V,2),1:size(V,1),1:size(V,3));
end
% Set user-defined volume (body) transparency:
tmp = strcmpi(varargin,'alpha');
if any(tmp)
    Alpha = varargin{find(tmp)+1};
    assert(Alpha>=0,'Alpha value must be between zero and one.')
    assert(Alpha<=1,'Alpha value must be between zero and one.')
end
% Set user-defined edge (sides, top, and bottom) transparency:
tmp = strcmpi(varargin,'edgealpha');
if any(tmp)
    EdgeAlpha = varargin{find(tmp)+1};
    assert(EdgeAlpha>=0,'EdgeAlpha value must be between zero and one.')
    assert(EdgeAlpha<=1,'EdgeAlpha value must be between zero and one.')
end
% Set user-defined volume (body) transparency:
tmp = strcmpi(varargin,'alphalim');
if any(tmp)
    AlphaLim = varargin{find(tmp)+1};
    setAlphaLim = true;
end
% Number of slices:
tmp = strcmpi(varargin,'n');
if any(tmp)
    nx = varargin{find(tmp)+1};
    ny = nx;
    nz = nx;
    assert(isscalar(nx)==1,'Invalid input after N declaration. Must be a scalar.')
    assert(nx>=0,'Number of slices N must be greater than zero.')
end
tmp = strcmpi(varargin,'nx');
if any(tmp)
    nx = varargin{find(tmp)+1};
    assert(isscalar(nx)==1,'Invalid input after Nx declaration. Must be a scalar.')
    assert(nx>=0,'Number of slices Nx must be greater than zero.')
end
tmp = strcmpi(varargin,'ny');
if any(tmp)
    ny = varargin{find(tmp)+1};
    assert(isscalar(ny)==1,'Invalid input after Ny declaration. Must be a scalar.')
    assert(ny>=0,'Number of slices Ny must be greater than zero.')
end
tmp = strcmpi(varargin,'nz');
if any(tmp)
    nz = varargin{find(tmp)+1};
    assert(isscalar(nz)==1,'Invalid input after Nz declaration. Must be a scalar.')
    assert(nz>=0,'Number of slices Nz must be greater than zero.')
end
% Interpolation method:
if any(strncmpi(varargin,'cubic',3))
    InterpolationMethod = 'cubic';
end
if any(strncmpi(varargin,'nearest',4))
    InterpolationMethod = 'nearest';
end
if any(strncmpi(varargin,'direct',3))
    InterpolationMethod = 'direct';
end
%% Some more checks now that all inputs are parsed:
assert(ndims(V)==3,'Input volume matrix V must be 3 dimensional.')
% Allow inputs as vectors:
if isvector(X)
    assert(isvector(Y)==1,'If X is a vector, Y must be a vector. Check your input X,Y,Z values.')
    assert(isvector(Z)==1,'If X is a vector, Z must be a vector. Check your input X,Y,Z values.')
    [X,Y,Z] = meshgrid(X,Y,Z);
end
% Make sure no 2D X grid slipped in there:
assert(ndims(X)==3,'Currently, X must be 1D or 3D with dimensions corresponding to V. This might change in the future, but until then, use meshgrid.')
%%
switch InterpolationMethod
    case 'direct'
        nx = size(V,2);
        ny = size(V,1);
        nz = size(V,3);
        
        % Make direct alpha roughly the same total value as when 300 slices are interpolated:
        Alpha = Alpha*(300/(nx+ny+nz));
        hold on
        
        % Set 3D view if not already 3D:
        [az,el] = view;
        if az==0 && el==90
            view(3)
        end
        
        % Plot x slices:
        for k = 1:nx
            h(k) = surface(squeeze(X(:,k,:)),squeeze(Y(:,k,:)),squeeze(Z(:,k,:)),squeeze(V(:,k,:)));
        end
        
        % Plot y slices:
        for k2 = 1:ny
            k = k+1;
            h(k) = surface(squeeze(X(k2,:,:)),squeeze(Y(k2,:,:)),squeeze(Z(k2,:,:)),squeeze(V(k2,:,:)));
        end
        
        % Plot z slices:
        for k3 = 1:nz
            k = k+1;
            h(k) = surface(squeeze(X(:,:,k3)),squeeze(Y(:,:,k3)),squeeze(Z(:,:,k3)),squeeze(V(:,:,k3)));
        end
        grid on
        
    otherwise
        
        % Generate slices:
        xslice = linspace(min(X(:)),max(X(:)),nx);
        yslice = linspace(min(Y(:)),max(Y(:)),ny);
        zslice = linspace(min(Z(:)),max(Z(:)),nz);
        % Plot slices:
        h = slice(X,Y,Z,V,xslice,yslice,zslice,InterpolationMethod);
end
% Set formatting:
shading interp
if setAlphaLim
    
    if strcmpi(AlphaLim,'auto')
        AlphaLim = [min(V(:)) max(V(:))];
    end
    assert(numel(AlphaLim)==2,'AlphaLim can only be a two-element array or ''auto''.')
    assert(AlphaLim(2)>AlphaLim(1),'AlphaLim values must be in the order [minAlphaLim maxAlphaLim].')
    
    switch InterpolationMethod
        case 'direct'
            
            % Plot x slices:
            for k = 1:nx
                set(h(k),'alphadata',Alpha*abs((squeeze(V(:,k,:))-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
            % Plot y slices:
            for k2 = 1:ny
                k = k+1;
                set(h(k),'alphadata',Alpha*abs((squeeze(V(k2,:,:))-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
            % Plot z slices:
            for k3 = 1:nz
                k = k+1;
                set(h(k),'alphadata',Alpha*abs((squeeze(V(:,:,k3))-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
            
        otherwise
            
            % Plot x slices:
            for k = 1:nx
                Vi = get(h(k),'Cdata');
                set(h(k),'alphadata',Alpha*abs((Vi-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
            % Plot y slices:
            for k2 = 1:ny
                k = k+1;
                Vi = get(h(k),'Cdata');
                set(h(k),'alphadata',Alpha*abs((Vi-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
            % Plot z slices:
            for k3 = 1:nz
                k = k+1;
                Vi = get(h(k),'Cdata');
                set(h(k),'alphadata',Alpha*abs((Vi-AlphaLim(1)))/(AlphaLim(2)-AlphaLim(1)),...
                    'AlphaDataMapping','none','facealpha','flat','edgecolor','none')
            end
    end
    
else
    
    set(h,'edgecolor','none','facealpha',Alpha)
    
end
% Set different transparency (typically slightly more opaque) for sides, top and bottom:
set(h([1 nx nx+1 nx+ny nx+ny+1 nx+ny+nz]),'facealpha',EdgeAlpha)
%% Clean up:
if nargout==0
    clear h
end
end
function [cmap] = rgbmap(varargin)
% RGBMAP creates color maps from any list of colors given as their common
% names.  Include a scalar anywhere in the list to specify the total number
% of levels in the color map.
%
% This function requires the rgb function found on the Mathworks File
% Exchange site here: http://www.mathworks.com/matlabcentral/fileexchange/46872
%
%
%% Syntax
%
%  cmap = rgbmap('first color name','second color name')
%  cmap = rgbmap('first color name','second color name',...,'nth color name')
%  cmap = rgbmap(...,levels)
%  rgbmap(...)
%
%% Description
%
% |cmap = rgbmap('first color name','second color name')| creates an RGB
% color map |cmap| from some first color to a second color.
%
% |cmap = rgbmap('first color name','second color name',...,'nth color
% name')| creates a color map linearly scaled between any number of colors.
%
% |cmap = rgbmap(...,M)| specifies the approximate number of levels |M| of the M x 3
% output colormap. Default value is 256.
%
% |rgbmap(...)| sets the color map without creating an array in the
% workspace.
%
%
%% Example 1: Create a nice color map for some scattered data
%
% x = 1:50;
% y = cos(x*pi/50);
% colors = rgbmap('blue','taupe','silver',50) ;
% scatter(x,y,30,colors,'filled')
%
%% Example 2: Get a 12-color map matrix from red to blue
%
% cmap = rgbmap('red','blue',12)
%
% h = surf(peaks);
% colorbar
% colormap(cmap)
% shading interp
% set(h,'edgecolor','k','edgealpha',.2)
% axis tight
%
%% Example 3: Plot blue to white to red
%
% h = surf(peaks);
% colorbar
% rgbmap('blue','white','red')
% shading interp
% set(h,'edgecolor','k','edgealpha',.2)
% caxis([-5 5])
% axis tight
%
%% Author Info
% This function was written by Chad A. Greene on June 5, 2014.
%
% See also rgb and colormap.
%
%% Input check:
assert(nargin>1,'Must have at least two inputs as strings.')
assert(exist('rgb','file')==2,'Cannot find the rgb function. Make sure Matlab knows where it is. If you don''t know where it is, you can download it here: http://www.mathworks.com/matlabcentral/fileexchange/46872')
%% Create color map:
levels = 256; % default
for k = 1:length(varargin)
    if isnumeric(varargin{k})
        levels = varargin{k};
        varargin(k)=[];
    end
end
for k = 1:length(varargin)
    try
        cp(k,:) = rgb(varargin{k});
    catch err
        error('MATLAB:rgbmap:rgb',['Cannot find an rgb value for the color ',varargin{k},'.'])
    end
end
numColors = length(varargin);
levPerCol = floor(levels/(numColors - 1));
cmap=[];
for k = 1:numColors-1
    cmap(length(cmap)+1:length(cmap)+levPerCol,:) = [linspace(cp(k,1),cp(k+1,1),levPerCol)' linspace(cp(k,2),cp(k+1,2),levPerCol)' linspace(cp(k,3),cp(k+1,3),levPerCol)'];
end
% The next two while loops are a somewhat clunky fix in case the number of
% colors does not equal the number of levels specified by the user. If that
% happens, the last line(s) will be deleted or repeated until the correct
% number of colors is obtained.
while length(cmap(:,1))>levels
    cmap(end,:)=[];
end
while length(cmap(:,1))<levels
    cmap(end+1,:)=cmap(end);
end
if nargout==0
    colormap(cmap)
    clear cmap
end
end
function [RGB] = rgb(ColorNames,varargin)
% RGB returns the RGB triple of the 949 most common names for colors,
% according to the results of the XKCD survey described here:
% http://blog.xkcd.com/2010/05/03/color-survey-results/
%
% To see a chart of color options, check out this internet website:
% http://xkcd.com/color/rgb/
%
% SYNTAX
% RGB = rgb('Color Name')
% RGB = rgb('Color Name 1','Color Name 2',...,'Color Name N')
% RGB = rgb({'Color Name 1','Color Name 2',...,'Color Name N'})
%
%
% DESCRIPTION
% RGB = rgb('Color Name') returns the rgb triple of the color described by the
% string 'Color Name'. Try any color name you can think of, it'll probably
% work.
%
% RGB = rgb('Color Name 1','Color Name 2',...,'Color Name N') returns an
% N by 3 matrix containing RGB triplets for each color name.
%
% RGB = rgb({'Color Name 1','Color Name 2',...,'Color Name N'}) accepts
% list of color names as a character array.
%
%
% EXAMPLE 1
% rgb('baby blue')
%
%
% EXAMPLE 2
% rgb('wintergreen','sunflower yellow','sapphire','radioactive green')
%
%
% EXAMPLE 3
% x = -pi:.1:pi;
% y = cos(x);
% plot(x,y,'linewidth',4,'color',rgb('cornflower'))
% hold on
% plot(x,y-1,'*','color',rgb('plum'))
% plot(x,y-2,'linewidth',4,'color',rgb('puke green'))
% legend('cornflower','plum','puke green')
% set(gca,'color',rgb('azure'))
% text(0,-2,'This text is burnt orange.','fontweight','bold','color',rgb('burnt orange'));
%
%
% AUTHOR INFO
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics.  I do not claim credit for the data
% from the color survey. http://www.chadagreene.com.
%
% Updated July 2015 to fix an installation error. Thanks for the tip,
% Florian Klimm!
%
% See also ColorSpec, hex2rgb, rgb2hex.
%% Make sure the function can find the data:
if exist('xkcd_rgb_data.mat','file')~=2
    disp 'Cannot find xkcd_rgb_data.mat. I will try to install it from rgb.txt now.'
    rgb_install
    disp 'Installation complete.'
end
%% Check inputs:
if iscellstr(ColorNames)==0 && iscellstr({ColorNames})==1
    ColorNames = {ColorNames};
end
if ~isempty(varargin)
    ColorNames = [{ColorNames} varargin];
end
assert(isnumeric(ColorNames)==0,'Input must be color name(s) as string(s).')
%% Search for color, return rgb value:
% Load data created by rgb_install.m:
load xkcd_rgb_data.mat
% Number of input color names:
numcolors = length(ColorNames);
% Preallocate a matrix to fill with RGB triplets:
RGB = NaN(numcolors,3);
% Find rgb triplet for each input color string:
for k = 1:numcolors
    ColorName = ColorNames{k};
    ColorName = strrep(ColorName,'gray','grey'); % because many users spell it 'gray'.
    % If a color name is not found in the database, display error message
    % and look for near matches:
    if sum(strcmpi(colorlist,ColorName))==0
        disp(['Color ''',ColorName,''' not found. Consider one of these options:'])
        
        % Special thanks to Cedric Wannaz for writing this bit of code. He came up with a
        % quite clever solution wherein the spectrum of the input color
        % name is compared to the spectra of available options.  So cool.
        spec = @(name) accumarray(upper(name.')-31, ones(size(name)), [60 1]) ;
        spec_mycol = spec(ColorName); % spectrum of input color name
        spec_dist = cellfun(@(name) norm(spec(name)-spec_mycol), colorlist);
        [sds,idx]   = sort(spec_dist) ;
        nearbyNames = colorlist(idx(sds<=1.5));
        if isempty(nearbyNames)
            nearbyNames = colorlist(idx(1:3));
        end
        disp(nearbyNames);
        
        clear RGB
        return % gives up and makes the user try again.
    end
    RGB(k,:) = rgblist(strcmpi(colorlist,ColorName),:);
end
end
%% Installation subfunction:
function rgb_install
if ~exist('rgb.txt','file')
    disp 'Cannot find rgb.txt file. I will try to download it from the internet now.'
    try
        urlwrite('http://xkcd.com/color/rgb.txt','rgb.txt');
    catch
        disp('Having trouble downloading the data file. You''ll need to download it manually')
        disp('from http://xkcd.com/color/rgb.txt and place it in your current folder.')
        return
    end
end

fid = fopen('rgb.txt');
RGB = textscan(fid,'%s %s','delimiter','\t');
fclose(fid);
colorlist = RGB{1};
hex = RGB{2};
rgblist = newhex2rgb(char(hex));
save('xkcd_rgb_data.mat','colorlist','rgblist')
end
%% newhex2rgb subfunction:
function [ rgb ] = newhex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1.
%
%
% * * * * * * * * * * * * * * * * * * * *
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default.
%
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255.
%
%
% * * * * * * * * * * * * * * * * * * * *
% EXAMPLES:
%
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
%
%
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional
%    = 0.2000    0.3020    0.4000
%
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
%
%
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
%
%
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
%
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'};
% rgbvals = hex2rgb(HexValsAsACharacterArray)
%
% * * * * * * * * * * * * * * * * * * * *
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input.
%
% * * * * * * * * * * * * * * * * * * * *
% See also rgb2hex, dec2hex, hex2num, and ColorSpec.
%
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.')
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary:
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex';
    end
    
    % If input is cell, convert to matrix:
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
r = hex2dec(hex(:,2:3));
g = hex2dec(hex(:,4:5));
b = hex2dec(hex(:,6:7));
rgb = [r g b]/255;
end