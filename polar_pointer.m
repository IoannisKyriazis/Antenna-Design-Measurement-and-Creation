function hpol = polar_pointer(varargin)
%POLAR  Polar coordinate plot.
%   POLAR(THETA, RHO) makes a plot using polar coordinates of
%   the angle THETA, in radians, versus the radius RHO.
%   POLAR(THETA,RHO,R) uses the radial limits specified by the two element
%   vector R. If R is a 4-element vector, the second two elements are
%   angular limits.
%   POLAR(THETA,RHO,S) uses the linestyle specified in string S.
%   See PLOT for a description of legal linestyles. 
%   POLAR(THETA,RHO,R,S) uses the linestyle specified in string S and the
%   radial limits in R.
%
%   POLAR(AX,...) plots into AX instead of GCA.
%
%   H = POLAR(...) returns a handle to the plotted object in H.
%
%   Example:
%      t = 0:.01:2*pi;
%      polar(t,sin(2*t).*cos(2*t),'--r')
%
%   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
%
%   Revised version by Daniel Armyr, 2009. Based on Mathworks original.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 5.22.4.9 $  $Date: 2007/08/27 17:06:52 $

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,4,nargs,'struct'));

%global size_rho;
%global theta_fg;
%global theta_f;

if nargs < 1 || nargs > 4
    error('MATLAB:polar:InvalidInput', 'Requires 2 to 4 data arguments.')
elseif nargs == 2 
    theta = args{1};
    rho = args{2};
    
    %size_rho=size(rho);
    
    if ischar(rho)
        line_style = rho;
        rho = theta;
        [mr,nr] = size(rho);
        if mr == 1
            theta = 1:nr;
        else
            th = (1:mr)';
            theta = th(:,ones(1,nr));
        end
    else
        line_style = 'auto';
    end
    radial_limits = [];
elseif nargs == 1
    theta = args{1};
    line_style = 'auto';
    rho = theta;
    [mr,nr] = size(rho);
    if mr == 1
        theta = 1:nr;
    else
        th = (1:mr)';
        theta = th(:,ones(1,nr));
    end
    radial_limits = [];
elseif nargs == 3
    if ( ischar(args{3}) )
        [theta,rho,line_style] = deal(args{1:3});
        radial_limits = [];
    else
        [theta,rho,radial_limits] = deal(args{1:3});
        line_style = 'auto';
        if ( numel(radial_limits) == 2 )
			radial_limits = [radial_limits(:); 0; 2*pi];
		elseif ( ~(numel(radial_limits) == 4) )
            error ( 'R must be a 2 element vector' );
        end
    end
else %nargs == 4
	[theta,rho,radial_limits,line_style] = deal(args{1:4});
	if ( numel(radial_limits) == 2 )
		radial_limits = [radial_limits(:); 0; 2*pi];
	elseif ( ~(numel(radial_limits) == 4) )
		error ( 'R must be a 2 element vector' );
	end
end

if ischar(theta) || ischar(rho)
    error('MATLAB:polar:InvalidInputType', 'Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    error('MATLAB:polar:InvalidInput', 'THETA and RHO must be the same size.');
end

% get hold state
cax = newplot(cax);
set(gcf,'Position',[100 100 800 500]);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

title('Polar Plot (Pointer)')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold(cax,'on');
    set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    
    % ensure that Inf values don't enter into the limit calculation.
    arho = abs(rho(:));    
    if ( isempty(radial_limits) )
        maxrho = max(arho(arho ~= Inf));
        minrho = 0;        
        hhh=line([minrho minrho maxrho maxrho],[minrho maxrho maxrho minrho],'parent',cax);        
        v = [get(cax,'xlim') get(cax,'ylim')];
        ticks = numel(get(cax,'ytick'));
        delete(hhh);
        % check radial limits and ticks
        rmin = v(1); rmax = v(4); rticks = max(ticks-1,2);
        if rticks > 5   % see if we can reduce the number
            if rem(rticks,2) == 0
                rticks = rticks/2;
            elseif rem(rticks,3) == 0
                rticks = rticks/3;
            end
        end
        rinc = (rmax-rmin)/rticks;
        
	else
		
		thmax = radial_limits(4);
        thmin = radial_limits(3);
		rmax = radial_limits(2);
        rmin = radial_limits(1);
		
        
        
        
        order = (10^floor(log10(rmax-rmin)));
        firstDigit = floor((rmax-rmin)/order);
        if ( firstDigit <= 1 )
            step = 0.2*order;
        elseif ( firstDigit <= 3 )
            step = 0.5*order;
        elseif ( firstDigit <= 7 )
            step = order;
        else
            step = 2*order;
        end
        rinc = step;
		
		%Also, make sure to drop any values outside the limits.
		if ( rmin == 0 )
			subset = rho >= -rmax & rho <= rmax & theta >= thmin & theta <= thmax;
		else
			subset = rho >= rmin & rho <= rmax & theta >= thmin & theta <= thmax;
		end
		theta(~subset) = nan;
		rho(~subset) = nan;		
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~ischar(get(cax,'color')),
       patch('xdata',xunit*(rmax-rmin),'ydata',yunit*(rmax-rmin), ...
             'edgecolor',tc,'facecolor',get(cax,'color'),...
             'handlevisibility','off','parent',cax);
    end

% draw radial circles
    c0 = cos(0*pi/180);
    s0 = sin(0*pi/180);    
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*(i-rmin),yunit*(i-rmin),'linestyle',ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off','parent',cax);
        text((i-rmin-2)*c0,(i-rmin-2)*s0, ...
             num2str(i),'verticalalignment','top','horizontalalignment','center',...
             'fontsize',7,'fontname','Times New Roman','color',[0.2 0 0.8],...
            'handlevisibility','off','parent',cax)
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line((rmax-rmin)*cs,(rmax-rmin)*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax)

% annotate spokes in degrees
    rt = 1.1*(rmax-rmin);
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),...
             'horizontalalignment','center',...
             'handlevisibility','off','parent',cax);
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','parent',cax)
    end

% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,(rmax-rmin)*[-1 1 -1.15 1.15]);
    
    setappdata( cax, 'rMin', rmin );
    
else
    %Try to find the inner radius of the current axis.
    if ( isappdata ( cax, 'rMin' ) )
        rmin = getappdata( cax, 'rMin' );
    else
        rmin = 0;
    end
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits);

% transform data to Cartesian coordinates.
xx = (rho - rmin).*cos(theta);
yy = (rho - rmin).*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
    q = plot(xx,yy,'parent',cax,'linewidth',2);
    hold on
%q_m=plot([0 abs(rmin)],[0 0],'r','parent',cax,'linewidth',2);
q_m=plot([0 abs(rmin)],[0 0],'color',[0 1 1],'parent',cax,'linewidth',2);
hold off
else
    q = plot(xx,yy,line_style,'parent',cax);
    hold on
%q_m=plot([0 abs(rmin)],[0 0],line_style,'r','parent',cax,'linewidth',2);
q_m=plot([0 abs(rmin)],[0 0],'color',[0 1 1],'parent',cax,'linewidth',2);
hold off
end



theta_f=0;


set(q_m,'ButtonDownFcn',@startDragFcn);

set (gcf, 'WindowButtonUpFcn', @stopDragFcn);



if nargout == 1
    hpol = q;
end

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')

if ~isempty(q) && ~isdeployed
    makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
end


size_rho=size(rho);
v2=[get(cax,'xlim') get(cax,'ylim')];

function startDragFcn(varargin)
    set (cax, 'HitTest','off');
set(gcf,'WindowButtonMotionFcn',{@mouseMove,q_m,size_rho,v2})

end



function mouseMove(object, eventdata,q_m,size_rho,v2)
C = get (gca, 'CurrentPoint');
[theta2 rho_m] = cart2pol(C(1,1), C(1,2));
%set(q_m,'r')

th_m=theta2;

if theta2<0
    theta2=2*pi+theta2;
end





theta_f=round(rad2deg(theta2));



theta_r=abs(round(theta_f*(size_rho/360)));

if theta_r==0
    theta_r=1;
end

rho_text=rho(abs(theta_r(1,1)),1);





    
text('string',{['Angle=',num2str(theta_f),'deg'],['Power= ',num2str(rho_text),'dB         ']},...
    'position',[-2.5*v2(1)/2,0],'backgroundcolor',[1 1 1],'color',[1 0 0]);





theta_m=[0 th_m];
        r_m=[0 abs(rmin)];


        [x y]=pol2cart(theta_m,r_m);
        
set(q_m,'xData',x,'yData',y);

end



   function stopDragFcn(varargin)
        set(gcf,'WindowButtonMotionFcn','');
   end




end

