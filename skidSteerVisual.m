function skidSteerVisual()
% SKIDSTEERVISUAL  Interactive 4-wheel skid-steer model.
%
%   >> skidSteerVisual
%
% You choose:
%   - the Instantaneous Center of Rotation (ICR) of the chassis
%     (DRAG the red star in the plot, or use the sliders)
%   - yaw rate, friction coefficient, chassis geometry, CG offset, mass
% It shows, for each wheel (FL, FR, RL, RR):
%   - wheel rotation speed (rad/s, rpm)  [animated tread marks]
%   - rolling velocity, lateral skid velocity, normal load, friction force
%
% Body frame: x forward, y left, yaw rate w positive CCW.
% ICR at (xc, yc) in the body frame  ->  body velocity at origin:
%       vx0 =  w*yc,   vy0 = -w*xc
% Wheel i at (xi, yi) moves with
%       vix = vx0 - w*yi    (along the wheel's rolling direction)
%       viy = vy0 + w*xi    (sideways -> the wheel must skid)
% Wheel rotation:  omega_i = vix / r      (+ = rolling forward)
%
% Friction model: wheels roll freely along x (motors supply that force);
% the lateral skid is resisted by Coulomb friction
%       Fy_i = -mu * N_i * tanh(viy_i / vs)
% with normal loads N_i set by mass and CG position.

% ----------------------------- state -----------------------------------
S = struct('xc',0.10,'yc',1.00,'w',1.00,'mu',0.60, ...
           'L',0.60,'b',0.50,'r',0.10,'xcg',0.00,'m',20);
S0 = S;
lim   = 2.5;                    % plot half-size [m]
ww    = 0.07;                   % wheel width (drawing) [m]
names = {'FL','FR','RL','RR'};
g     = 9.81;
dragging = false;
sliders  = struct();
phase    = zeros(1,4);          % tread animation phase
cur      = compute();           %#ok<NASGU>

% ----------------------------- figure ----------------------------------
fig = figure('Name','Skid-steer visual model','NumberTitle','off', ...
             'Position',[40 40 1300 760],'Color','w');
ax = axes('Parent',fig,'Position',[0.05 0.32 0.54 0.64]);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
xlabel(ax,'x forward [m]'); ylabel(ax,'y left [m]');
title(ax,'Drag in the plot to move the ICR (red star)');

o  = cur;
th = linspace(0,2*pi,240);
hCirc = plot(ax,nan,nan,'--','Color',[0.6 0.6 0.6]);
for j = 1:4
    hRay(j) = line('Parent',ax,'XData',[0 0],'YData',[0 0], ...
                   'Color',[0.75 0.75 0.75],'LineStyle',':'); %#ok<AGROW>
end
hBody = patch('Parent',ax,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
              'FaceColor',[0.85 0.9 1],'EdgeColor',[0.2 0.3 0.6],'LineWidth',1.5);
for j = 1:4
    hWh(j) = patch('Parent',ax,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
                   'FaceColor',[0.25 0.25 0.25]); %#ok<AGROW>
    hTr(j) = line('Parent',ax,'XData',nan,'YData',nan, ...
                  'Color',[0.95 0.95 0.95],'LineWidth',1.5); %#ok<AGROW>
    hTx(j) = text('Parent',ax,'Position',[0 0],'String','', ...
                  'HorizontalAlignment','center','FontSize',9); %#ok<AGROW>
end
hVel  = quiver(ax,o.xi,o.yi,zeros(1,4),zeros(1,4),0,'Color',[0 0 0],   'LineWidth',1.0,'MaxHeadSize',0.6);
hRoll = quiver(ax,o.xi,o.yi,zeros(1,4),zeros(1,4),0,'Color',[0 .35 .85],'LineWidth',2.0,'MaxHeadSize',0.6);
hSkid = quiver(ax,o.xi,o.yi,zeros(1,4),zeros(1,4),0,'Color',[.85 .1 .1],'LineWidth',2.0,'MaxHeadSize',0.6);
hFric = quiver(ax,o.xi,o.yi,zeros(1,4),zeros(1,4),0,'Color',[.75 0 .75],'LineWidth',2.0,'MaxHeadSize',0.6);
hCG   = plot(ax,0,0,'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
hICR  = plot(ax,S.xc,S.yc,'p','MarkerSize',18,'MarkerFaceColor','r','MarkerEdgeColor','k');
legend(ax,[hVel hRoll hSkid hFric],{'total wheel velocity','rolling (wheel spin)', ...
       'lateral skid','friction force'},'Location','northeast','FontSize',8);

% ----------------------------- controls --------------------------------
defs = { 'xc',  'ICR x  [m]  (+ forward)',  -lim, lim, '%+.2f';
         'yc',  'ICR y  [m]  (+ left)',     -lim, lim, '%+.2f';
         'w',   'Yaw rate  w  [rad/s]',     -5,   5,   '%+.2f';
         'mu',  'Friction coefficient  mu',  0,   1.5, '%.2f';
         'L',   'Wheelbase  L  [m]',         0.2, 1.5, '%.2f';
         'b',   'Track width  b  [m]',       0.2, 1.5, '%.2f';
         'r',   'Wheel radius  r  [m]',      0.03,0.3, '%.3f';
         'xcg', 'CG offset x  [m]',         -0.3, 0.3, '%+.2f';
         'm',   'Mass  [kg]',                1,   150, '%.0f' };
for i = 1:size(defs,1)
    mkSlider(defs{i,1},defs{i,2},defs{i,3},defs{i,4},defs{i,5},0.90-(i-1)*0.065);
end
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.62 0.30 0.115 0.045], ...
    'String','Spin in place','Callback',@(~,~) preset(0,0));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.745 0.30 0.115 0.045], ...
    'String','Wide turn','Callback',@(~,~) preset(0,2));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.87 0.30 0.11 0.045], ...
    'String','Reset all','Callback',@(~,~) resetAll());

hTab = uitable('Parent',fig,'Units','normalized','Position',[0.05 0.04 0.54 0.22], ...
    'ColumnName',{'Wheel','w_wheel [rad/s]','rpm','roll v [m/s]','skid v [m/s]','N [N]','Fy fric [N]'}, ...
    'RowName',[],'ColumnWidth',{55,105,70,90,90,70,85},'FontSize',11);
hSum = uicontrol(fig,'Style','text','Units','normalized','Position',[0.62 0.04 0.36 0.24], ...
    'HorizontalAlignment','left','BackgroundColor','w','FontName','FixedWidth','FontSize',10);

set(fig,'WindowButtonDownFcn',  @(~,~) mouseDown(), ...
        'WindowButtonMotionFcn',@(~,~) mouseMove(), ...
        'WindowButtonUpFcn',    @(~,~) mouseUp(), ...
        'CloseRequestFcn',      @(~,~) closeAll());
tm = timer('ExecutionMode','fixedSpacing','Period',0.05,'BusyMode','drop', ...
           'TimerFcn',@(~,~) tick());

redraw();
start(tm);

% ======================= nested helper functions =======================
    function o = compute()
        L = S.L; b = S.b; r = S.r;
        xi = [ L/2,  L/2, -L/2, -L/2];          % FL FR RL RR
        yi = [ b/2, -b/2,  b/2, -b/2];
        o.xi = xi; o.yi = yi;
        o.vx0 =  S.w*S.yc;
        o.vy0 = -S.w*S.xc;
        o.vix = o.vx0 - S.w*yi;
        o.viy = o.vy0 + S.w*xi;
        o.om  = o.vix/r;
        xcg = max(min(S.xcg, L/2), -L/2);
        a   = (L/2 + xcg)/L;                    % front axle load share
        o.N = S.m*g*[a a 1-a 1-a]/2;
        vs  = 0.02;                             % friction smoothing speed
        o.Fy   = -S.mu*o.N.*tanh(o.viy/vs);
        o.Fnet = sum(o.Fy);
        o.Mz   = sum((xi - xcg).*o.Fy);
        o.P    = sum(-o.Fy.*o.viy);             % skid dissipation power
    end

    function redraw()
        o = compute(); cur = o;
        L = S.L; r = S.r;
        set(hBody,'XData',[-L/2-.03 L/2+.03 L/2+.03 -L/2-.03], ...
                  'YData',[-S.b/2+ww/2 -S.b/2+ww/2 S.b/2-ww/2 S.b/2-ww/2]);
        sv = 0.7/max([abs(o.vix) abs(o.viy) 0.2]);   % arrow scale (m per m/s)
        sf = 0.5/(S.m*g/4*1.5);                      % friction arrow scale
        for j = 1:4
            set(hWh(j),'XData',o.xi(j)+[-r r r -r], ...
                       'YData',o.yi(j)+[-ww/2 -ww/2 ww/2 ww/2]);
            set(hTx(j),'Position',[o.xi(j), o.yi(j)+sign(o.yi(j))*(ww/2+0.15)], ...
                       'String',sprintf('%s\n%+.1f rad/s',names{j},o.om(j)));
            set(hRay(j),'XData',[S.xc o.xi(j)],'YData',[S.yc o.yi(j)]);
        end
        set(hVel, 'XData',o.xi,'YData',o.yi,'UData',o.vix*sv,'VData',o.viy*sv);
        set(hRoll,'XData',o.xi,'YData',o.yi,'UData',o.vix*sv,'VData',zeros(1,4));
        set(hSkid,'XData',o.xi,'YData',o.yi,'UData',zeros(1,4),'VData',o.viy*sv);
        set(hFric,'XData',o.xi,'YData',o.yi,'UData',zeros(1,4),'VData',o.Fy*sf);
        R = hypot(S.xc,S.yc);
        set(hCirc,'XData',S.xc+R*cos(th),'YData',S.yc+R*sin(th));
        set(hICR,'XData',S.xc,'YData',S.yc);
        set(hCG,'XData',max(min(S.xcg,L/2),-L/2),'YData',0);

        rows = cell(4,7);
        for j = 1:4
            rows(j,:) = {names{j}, sprintf('%+.2f',o.om(j)), sprintf('%+.1f',o.om(j)*30/pi), ...
                         sprintf('%+.3f',o.vix(j)), sprintf('%+.3f',o.viy(j)), ...
                         sprintf('%.1f',o.N(j)), sprintf('%+.1f',o.Fy(j))};
        end
        set(hTab,'Data',rows);
        set(hSum,'String',{ ...
            sprintf('Body velocity   vx = %+.3f   vy = %+.3f m/s',o.vx0,o.vy0), ...
            sprintf('Body slip angle = %+.1f deg',rad2deg(atan2(o.vy0,o.vx0))), ...
            sprintf('Yaw rate        w  = %+.3f rad/s',S.w), ...
            sprintf('ICR distance    R  = %.3f m',R), ...
            sprintf('Net lateral friction  sum(Fy) = %+.2f N',o.Fnet), ...
            sprintf('Friction yaw moment   Mz      = %+.2f N*m',o.Mz), ...
            sprintf('Skid power dissipated        = %.2f W',o.P), ...
            '', ...
            'Wheel rotation: + = rolling forward (x).', ...
            'Skid v: + = sliding to the left (y).'});
        drawTreads();
    end

    function tick()
        if ~ishandle(fig), return; end
        sp = S.r/2;
        step = max(min(cur.vix*0.05, sp/2), -sp/2);   % cap to avoid aliasing
        phase = mod(phase + step, sp);
        drawTreads();
    end

    function drawTreads()
        r = S.r; sp = r/2;
        for j = 1:4
            xs = -r + mod(phase(j),sp) + (0:3)*sp;
            X = [xs; xs; nan(1,4)];
            Y = repmat([-ww/2; ww/2; NaN],1,4);
            set(hTr(j),'XData',cur.xi(j)+X(:),'YData',cur.yi(j)+Y(:));
        end
    end

    function mkSlider(key,label,lo,hi,fmt,yPos)
        uicontrol(fig,'Style','text','Units','normalized', ...
            'Position',[0.62 yPos+0.027 0.24 0.025],'String',label, ...
            'HorizontalAlignment','left','BackgroundColor','w');
        vt = uicontrol(fig,'Style','text','Units','normalized', ...
            'Position',[0.88 yPos+0.027 0.10 0.025],'String',sprintf(fmt,S.(key)), ...
            'HorizontalAlignment','right','BackgroundColor','w');
        sl = uicontrol(fig,'Style','slider','Units','normalized', ...
            'Position',[0.62 yPos 0.36 0.025],'Min',lo,'Max',hi,'Value',S.(key), ...
            'Callback',@(src,~) onSlide(key,src,vt,fmt));
        sliders.(key) = struct('sl',sl,'vt',vt,'fmt',fmt);
    end

    function onSlide(key,src,vt,fmt)
        S.(key) = get(src,'Value');
        set(vt,'String',sprintf(fmt,S.(key)));
        redraw();
    end

    function syncSliders()
        f = fieldnames(sliders);
        for i = 1:numel(f)
            s = sliders.(f{i});
            set(s.sl,'Value',S.(f{i}));
            set(s.vt,'String',sprintf(s.fmt,S.(f{i})));
        end
    end

    function setICR(p)
        S.xc = p(1); S.yc = p(2);
        syncSliders(); redraw();
    end

    function preset(x,y)
        setICR([x y]);
    end

    function resetAll()
        S = S0; syncSliders(); redraw();
    end

    function mouseDown()
        p = get(ax,'CurrentPoint'); p = p(1,1:2);
        if abs(p(1)) <= lim && abs(p(2)) <= lim
            dragging = true; setICR(p);
        end
    end

    function mouseMove()
        if ~dragging, return; end
        p = get(ax,'CurrentPoint'); p = p(1,1:2);
        setICR(max(min(p,lim),-lim));
    end

    function mouseUp()
        dragging = false;
    end

    function closeAll()
        try, stop(tm); delete(tm); catch, end %#ok<NOCOM>
        delete(fig);
    end
end
