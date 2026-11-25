function skidSteerPro()
% SKIDSTEERPRO  Extended-physics 4-wheel skid-steer simulator.
%
%   >> skidSteerPro
%
% You command wheel speeds (v, w -> naive ideal diff-drive inverse), and the
% rest is physics. The ICR, effective track factor chi, skid, loads, currents
% and tipping risk are all RESULTS.
%
% Physics included
%   Chassis      planar rigid body (vx, vy, yaw) with centripetal terms,
%                CG offset, yaw-inertia factor, aerodynamic drag
%   Load         dynamic wheel loads: longitudinal + lateral load transfer
%   transfer     via CG height h, front/rear split of lateral transfer,
%                first-order suspension lag, wheel lift-off / tip warning
%   Slope        incline angle + uphill direction (gravity components + normal
%                load reduction), static tip limits
%   Friction     anisotropic friction ellipse (mux along wheel, muy across)
%                with a soil-shear "bristle" model: friction builds with slip
%                distance (shear modulus K, Janosi-Shields-like) and gives
%                true static friction (no creep on a slope)
%   Terrain      rolling resistance coeff + Bekker compaction resistance
%                (sinkage from kc, kphi, n) for sand / loam / snow presets
%   Scrub        contact-patch spin (scrub) yaw moment vs wheel width
%   Drivetrain   PI wheel-speed loop -> duty -> DC motor (R, kt, back-EMF),
%                gearbox (ratio, efficiency, reflected inertia), current
%                limit, battery internal resistance (voltage sag)
%
% Frames: x forward, y left, yaw w positive CCW. Wheel order FL FR RL RR.
% Wheel spin positive = rolling forward.
%
% Simplifications that remain: flat-plane planar motion (no pitch/roll
% dynamics, only quasi-static load transfer), lumped soil model, no thermal
% or inductance effects.

% ----------------------------- constants --------------------------------
g = 9.81;  Jw = 0.01;  Jm = 6e-5;     % gravity, wheel inertia, motor rotor inertia
Kc = 0.2;  Ki = 1.5;                  % wheel-speed PI gains (duty per rad/s)
tauS = 0.05;                          % suspension load-transfer lag [s]
rhoAir = 1.2;
dt = 5e-4; nSub = 60;                 % physics step, steps per frame (0.03 s)
limB = 2.0; limW = 2.5;
names = {'FL','FR','RL','RR'};

TER = struct( ...
  'name',{'Hard (asphalt)','Gravel','Dry sand','Sandy loam / grass','Snow'}, ...
  'mx',  {0.80,0.55,0.55,0.50,0.25}, ...
  'my',  {0.70,0.50,0.50,0.45,0.20}, ...
  'K',   {0.004,0.010,0.010,0.025,0.025}, ...
  'crr', {0.02,0.04,0.01,0.01,0.01}, ...
  'bek', {[],[],[1.1 0.99 1528.43],[0.7 5.27 1515.04],[1.6 4.37 196.72]}, ...
  'col', {[0.93 0.93 0.93],[0.86 0.85 0.80],[0.96 0.90 0.70],[0.80 0.90 0.75],[0.95 0.97 1.0]});

S = struct('vcmd',0.8,'wcmd',0.8, ...
           'L',0.6,'b',0.5,'r',0.1,'bw',0.06,'xcg',0,'h',0.12,'m',20,'kz',1.0, ...
           'terr',1,'mx',0.8,'my',0.7,'K',0.004,'crr',0.02,'bek',[], ...
           'slope',0,'psi',0,'CdA',0.1, ...
           'Vbat',24,'Rbat',0.05,'Rm',0.6,'kt',0.05,'Ng',10,'eta',0.9,'Imax',20,'rho',0.5);
S0 = S;

% ----------------------------- state ------------------------------------
Q   = zeros(1,6);          % [X Y theta vx vy w]  (CG pose + body velocity)
om  = zeros(1,4);          % wheel speeds [rad/s]
zi  = zeros(1,4);          % PI integrators
zx  = zeros(1,4);          % soil shear deformation (bristle) x
zy  = zeros(1,4);          % soil shear deformation (bristle) y
Nld = zeros(1,4);          % wheel normal loads [N]
Icur = zeros(1,4);         % motor currents [A]
Ibat = 0; Vterm = S.Vbat; liftFlag = false;
lastFx = zeros(1,4); lastFy = zeros(1,4);
lastSx = zeros(1,4); lastSy = zeros(1,4);
lastTauW = zeros(1,4); lastTrr = zeros(1,4); lastRr = zeros(1,4);
G = struct();
phase = zeros(1,4); frame = 0; running = true;
trail = zeros(0,2);
sliders = struct();
groupH = struct('chassis',[],'terrain',[],'drive',[]);
th = linspace(0,2*pi,240);

updateDerived();
Nld = G.Nstatic;

% ----------------------------- figure -----------------------------------
fig = figure('Name','Skid-steer pro','NumberTitle','off', ...
             'Position',[20 20 1500 840],'Color','w');

% Chassis-frame view
axB = axes('Parent',fig,'Position',[0.04 0.40 0.28 0.55]);
hold(axB,'on'); grid(axB,'on'); axis(axB,'equal');
xlim(axB,[-limB limB]); ylim(axB,[-limB limB]);
xlabel(axB,'x forward [m]'); ylabel(axB,'y left [m]');
title(axB,'Chassis frame: forces and ICR (result)');
hCirc = plot(axB,nan,nan,'--','Color',[0.6 0.6 0.6]);
for pj = 1:4
    hRay(pj) = line('Parent',axB,'XData',nan,'YData',nan, ...
                    'Color',[0.75 0.75 0.75],'LineStyle',':'); %#ok<AGROW>
end
hBody = patch('Parent',axB,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
              'FaceColor',[0.85 0.9 1],'EdgeColor',[0.2 0.3 0.6],'LineWidth',1.5);
for pj = 1:4
    hWh(pj) = patch('Parent',axB,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
                    'FaceColor',[0.25 0.25 0.25]); %#ok<AGROW>
    hTr(pj) = line('Parent',axB,'XData',nan,'YData',nan, ...
                   'Color',[0.95 0.95 0.95],'LineWidth',1.5); %#ok<AGROW>
    hTx(pj) = text('Parent',axB,'Position',[0 0],'String','', ...
                   'HorizontalAlignment','center','FontSize',8); %#ok<AGROW>
end
z4 = zeros(1,4);
hVel  = quiver(axB,G.xg,G.yg,z4,z4,0,'Color',[0 0 0],     'LineWidth',1.0,'MaxHeadSize',0.6);
hRoll = quiver(axB,G.xg,G.yg,z4,z4,0,'Color',[0 .35 .85], 'LineWidth',2.0,'MaxHeadSize',0.6);
hSlip = quiver(axB,G.xg,G.yg,z4,z4,0,'Color',[.85 .1 .1], 'LineWidth',2.0,'MaxHeadSize',0.6);
hFric = quiver(axB,G.xg,G.yg,z4,z4,0,'Color',[.75 0 .75], 'LineWidth',2.0,'MaxHeadSize',0.6);
hCG   = plot(axB,0,0,'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
hICR  = plot(axB,nan,nan,'p','MarkerSize',18,'MarkerFaceColor','r','MarkerEdgeColor','k');
legend(axB,[hVel hRoll hSlip hFric],{'wheel ground velocity','wheel surface speed (r*w)', ...
       'slip velocity','friction force on chassis'},'Location','northeast','FontSize',7);

% World view
axW = axes('Parent',fig,'Position',[0.345 0.40 0.28 0.55]);
hold(axW,'on'); axis(axW,'equal'); box(axW,'on');
[GX,GY] = meshgrid(-60:2:60,-60:2:60);
plot(axW,GX(:),GY(:),'.','Color',[0.6 0.6 0.6],'MarkerSize',4);
xlabel(axW,'X [m]'); ylabel(axW,'Y [m]'); title(axW,'World view (follows robot)');
hTrail = plot(axW,nan,nan,'-','Color',[0.1 0.4 0.8],'LineWidth',1.2);
hWCirc = plot(axW,nan,nan,'--','Color',[0.5 0.5 0.5]);
hWBody = patch('Parent',axW,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
               'FaceColor',[0.85 0.9 1],'EdgeColor',[0.2 0.3 0.6],'LineWidth',1.5);
for pj = 1:4
    hWWh(pj) = patch('Parent',axW,'XData',[0 1 1 0],'YData',[0 0 1 1], ...
                     'FaceColor',[0.25 0.25 0.25]); %#ok<AGROW>
end
hWICR = plot(axW,nan,nan,'p','MarkerSize',14,'MarkerFaceColor','r','MarkerEdgeColor','k');
hUp   = quiver(axW,0,0,0,0,0,'Color',[0.6 0.1 0.1],'LineWidth',2,'MaxHeadSize',1);
hUpTx = text('Parent',axW,'Position',[0 0],'String','uphill','Color',[0.6 0.1 0.1],'Visible','off');
set(axW,'Color',TER(1).col);

% ----------------------------- controls ---------------------------------
xR = 0.68; wR = 0.30;
mkSlider('cmd','vcmd','Command speed  v  [m/s]',      -2,   2,    '%+.2f', 0.93);
mkSlider('cmd','wcmd','Command yaw rate  w  [rad/s]', -3,   3,    '%+.2f', 0.87);

btn = {'Straight',@(~,~) preset(0.8,0); 'Turn',@(~,~) preset(0.8,0.8); ...
       'Spin',@(~,~) preset(0,2);       'Stop',@(~,~) preset(0,0); ...
       'Reset',@(~,~) resetAll();       'Pause',@(~,~) togglePause()};
for pk = 1:size(btn,1)
    hBtn(pk) = uicontrol(fig,'Style','pushbutton','Units','normalized', ...
        'Position',[xR+(pk-1)*0.05 0.80 0.048 0.04],'String',btn{pk,1}, ...
        'FontSize',8,'Callback',btn{pk,2}); %#ok<AGROW>
end
gnames = {'chassis','terrain','drive'}; glabels = {'Chassis','Terrain','Drivetrain'};
for pk = 1:3
    hGrp(pk) = uicontrol(fig,'Style','pushbutton','Units','normalized', ...
        'Position',[xR+(pk-1)*0.1 0.75 0.098 0.04],'String',glabels{pk}, ...
        'FontSize',9,'Callback',@(~,~) showGroup(gnames{pk})); %#ok<AGROW>
end

% Chassis group
mkSlider('chassis','L',  'Wheelbase  L  [m]',               0.2, 1.2, '%.2f', slotY(1));
mkSlider('chassis','b',  'Track width  b  [m]',             0.2, 1.2, '%.2f', slotY(2));
mkSlider('chassis','r',  'Wheel radius  r  [m]',            0.04,0.25,'%.3f', slotY(3));
mkSlider('chassis','bw', 'Wheel width  [m]',                0.02,0.20,'%.3f', slotY(4));
mkSlider('chassis','xcg','CG offset x  [m]',               -0.3, 0.3, '%+.2f',slotY(5));
mkSlider('chassis','h',  'CG height  h  [m]',               0.02,0.5, '%.2f', slotY(6));
mkSlider('chassis','m',  'Mass  [kg]',                      2,   150, '%.0f', slotY(7));
mkSlider('chassis','kz', 'Yaw inertia factor (1 = box)',    0.4, 2.0, '%.2f', slotY(8));

% Terrain group
uTx = uicontrol(fig,'Style','text','Units','normalized','Position',[xR slotY(1)+0.027 0.2 0.025], ...
    'String','Terrain preset','HorizontalAlignment','left','BackgroundColor','w');
hTerr = uicontrol(fig,'Style','popupmenu','Units','normalized','Position',[xR slotY(1) wR 0.03], ...
    'String',{TER.name},'Value',1,'Callback',@(src,~) onTerrain(src));
groupH.terrain = [uTx hTerr];
mkSlider('terrain','mx',   'Friction along wheel  mu_x',     0.05,1.5, '%.2f', slotY(2));
mkSlider('terrain','my',   'Friction across wheel  mu_y',    0.05,1.5, '%.2f', slotY(3));
mkSlider('terrain','K',    'Soil shear modulus K [m]',       0.002,0.04,'%.3f',slotY(4));
mkSlider('terrain','crr',  'Rolling resistance coeff.',      0,   0.15,'%.3f', slotY(5));
mkSlider('terrain','slope','Slope angle  [deg]',             0,   45,  '%.1f', slotY(6));
mkSlider('terrain','psi',  'Uphill direction (world) [deg]', -180,180, '%.0f', slotY(7));
mkSlider('terrain','CdA',  'Aero drag  Cd*A  [m^2]',         0,   1.5, '%.2f', slotY(8));

% Drivetrain group
mkSlider('drive','Vbat','Battery voltage [V]',               6,   60,  '%.0f', slotY(1));
mkSlider('drive','Rbat','Battery internal R [ohm]',          0,   0.5, '%.3f', slotY(2));
mkSlider('drive','Rm',  'Motor terminal R [ohm]',            0.05,3,   '%.2f', slotY(3));
mkSlider('drive','kt',  'Motor torque const kt [N*m/A]',     0.01,0.2, '%.3f', slotY(4));
mkSlider('drive','Ng',  'Gear ratio',                        1,   60,  '%.1f', slotY(5));
mkSlider('drive','eta', 'Gearbox efficiency',                0.5, 1.0, '%.2f', slotY(6));
mkSlider('drive','Imax','Current limit / motor [A]',         1,   80,  '%.0f', slotY(7));
mkSlider('drive','rho', 'Front share of lateral load xfer',  0,   1,   '%.2f', slotY(8));

hTab = uitable('Parent',fig,'Units','normalized','Position',[0.04 0.03 0.60 0.30], ...
    'ColumnName',{'Wheel','w [rad/s]','rpm','I [A]','wheel tau [N*m]','slip x [m/s]', ...
                  'slip y [m/s]','N [N]','Fx [N]','Fy [N]','R_roll [N]'}, ...
    'RowName',[],'ColumnWidth',{50,75,55,55,100,80,80,55,60,60,75},'FontSize',10);
hSum = uicontrol(fig,'Style','text','Units','normalized','Position',[xR 0.02 wR 0.27], ...
    'HorizontalAlignment','left','BackgroundColor','w','FontName','FixedWidth','FontSize',8);

showGroup('chassis');
set(fig,'CloseRequestFcn',@(~,~) closeAll());
tm = timer('ExecutionMode','fixedSpacing','Period',0.03,'BusyMode','drop', ...
           'TimerFcn',@(~,~) tick());
drawAll(); updateTable();
start(tm);

% ======================= nested helper functions =======================
    function y = slotY(k)
        y = 0.67 - (k-1)*0.052;
    end

    function updateDerived()
        L = S.L; b = S.b;
        xcg = max(min(S.xcg, L/2), -L/2);
        G.xcg = xcg;
        G.xg = [ L/2,  L/2, -L/2, -L/2];        % wheel positions, geometric-centre frame
        G.yg = [ b/2, -b/2,  b/2, -b/2];
        G.xr = G.xg - xcg;                      % relative to CG
        G.yr = G.yg;
        G.phi = deg2rad(S.slope);
        G.psi = deg2rad(S.psi);
        G.W   = S.m*g*cos(G.phi);               % load normal to the slope
        a     = (L/2 + xcg)/L;
        G.Nstatic = G.W*[a a 1-a 1-a]/2;
        G.Iz  = S.kz*S.m*(L^2 + b^2)/12 + S.m*xcg^2;
        G.Jeff = Jw + Jm*S.Ng^2;
        wl = (S.vcmd - S.wcmd*b/2)/S.r;         % naive ideal inverse kinematics
        wr = (S.vcmd + S.wcmd*b/2)/S.r;
        G.omc = [wl wr wl wr];
    end

    function stepSim(dt)
        r = S.r; vx = Q(4); vy = Q(5); w = Q(6);
        % --- contact-patch slip and soil-shear friction (anisotropic) ---
        u  = vx - w*G.yr;
        v  = vy + w*G.xr;
        sx = u - r*om;
        sy = v;
        smag = sqrt(sx.^2 + sy.^2);
        zx = zx + dt*(sx - smag.*zx/S.K);
        zy = zy + dt*(sy - smag.*zy/S.K);
        kk = S.mx*Nld/S.K;                       % shear stiffness [N/m]
        mi = max(Nld/g, 0.2);
        cdamp = 2*0.7*sqrt(max(kk,0).*mi);          % damping for stable stick
        fxn = -zx/S.K - cdamp.*sx./max(S.mx*Nld,1e-6);
        fyn = -zy/S.K - cdamp.*sy./max(S.my*Nld,1e-6);
        mag = hypot(fxn,fyn);
        sc  = ones(1,4); idx = mag > 1; sc(idx) = 1./mag(idx);   % friction ellipse limit
        Fx = S.mx*Nld.*fxn.*sc;
        Fy = S.my*Nld.*fyn.*sc;
        % --- rolling resistance (constant coeff + Bekker compaction) ---
        Rr = S.crr*Nld;
        if ~isempty(S.bek)
            n  = S.bek(1);
            kp = (S.bek(2) + S.bw*S.bek(3))*1000;
            zs = (3*Nld/((3-n)*kp*sqrt(2*r))).^(2/(2*n+1));
            Rr = Rr + kp*zs.^(n+1)/(n+1);
        end
        Trr = r*Rr.*tanh(om/0.3);
        % --- drivetrain: PI speed loop, DC motor, gearbox, battery ---
        e    = G.omc - om;
        zi   = max(min(zi + dt*e, 1/Ki), -1/Ki);
        duty = max(min(Kc*e + Ki*zi, 1), -1);
        Vt   = max(min(S.Vbat - S.Rbat*Ibat, 1.3*S.Vbat), 0);
        omm  = S.Ng*om;
        Icur = max(min((duty*Vt - S.kt*omm)/S.Rm, S.Imax), -S.Imax);
        taum = S.kt*Icur;
        tauw = S.Ng*taum;
        mot  = (taum.*omm) >= 0;                 % motoring vs regenerating
        tauw(mot)  = tauw(mot)*S.eta;
        tauw(~mot) = tauw(~mot)/S.eta;
        Ibat  = sum(duty.*Icur);
        Vterm = Vt;
        % --- chassis dynamics ---
        gx = -g*sin(G.phi)*cos(Q(3) - G.psi);
        gy =  g*sin(G.phi)*sin(Q(3) - G.psi);
        V  = hypot(vx,vy);
        Dx = -0.5*rhoAir*S.CdA*V*vx;
        Dy = -0.5*rhoAir*S.CdA*V*vy;
        Msc = -(2/3)*0.5*(S.mx + S.my)*Nld*(S.bw/2).*tanh(w/0.1);   % patch scrub moment
        ax = (sum(Fx) + Dx)/S.m + gx + w*vy;
        ay = (sum(Fy) + Dy)/S.m + gy - w*vx;
        aw = (sum(G.xr.*Fy - G.yr.*Fx) + sum(Msc))/G.Iz;
        % --- quasi-static load transfer (target loads) ---
        SFx = sum(Fx); SFy = sum(Fy);
        Nf  = (-S.h*SFx + G.W*(S.L/2 + G.xcg))/S.L;
        Nr  = G.W - Nf;
        Tl  = -S.h*SFy/S.b;
        Nt  = [Nf/2 + S.rho*Tl, Nf/2 - S.rho*Tl, Nr/2 + (1-S.rho)*Tl, Nr/2 - (1-S.rho)*Tl];
        liftFlag = min(Nt) <= 0;
        Nld = Nld + dt*(max(Nt,0) - Nld)/tauS;
        % --- integrate ---
        c = cos(Q(3)); s = sin(Q(3));
        om   = om + dt*(tauw - r*Fx - Trr)/G.Jeff;
        Q(1) = Q(1) + dt*(vx*c - vy*s);
        Q(2) = Q(2) + dt*(vx*s + vy*c);
        Q(3) = Q(3) + dt*w;
        Q(4) = vx + dt*ax;
        Q(5) = vy + dt*ay;
        Q(6) = w  + dt*aw;
        lastFx = Fx; lastFy = Fy; lastSx = sx; lastSy = sy;
        lastTauW = tauw; lastTrr = Trr; lastRr = Rr;
    end

    function tick()
        if ~ishandle(fig) || ~running, return; end
        for kk2 = 1:nSub, stepSim(dt); end
        if any(~isfinite([Q om Nld])), resetAll(); return; end
        frame = frame + 1;
        sp = S.r/2;                                   % tread animation
        phase = mod(phase + max(min(S.r*om*nSub*dt, sp/2), -sp/2), sp);
        trail = [trail; Q(1) Q(2)]; %#ok<AGROW>
        if size(trail,1) > 4000, trail = trail(end-3999:end,:); end
        drawAll();
        if mod(frame,5) == 0, updateTable(); end
    end

    function [wx,wy] = toWorld(px,py)
        c = cos(Q(3)); s = sin(Q(3));
        dx = px - G.xcg; dy = py;                      % geometric frame -> CG frame
        wx = Q(1) + c*dx - s*dy;
        wy = Q(2) + s*dx + c*dy;
    end

    function [xI,yI,R] = icrGeom()
        w = Q(6);
        if abs(w) > 0.02
            xI = -Q(5)/w + G.xcg;                      % geometric-centre frame
            yI =  Q(4)/w;
            R  = hypot(Q(4),Q(5))/abs(w);
        else
            xI = NaN; yI = NaN; R = Inf;
        end
    end

    function drawAll()
        L = S.L; r = S.r; xi = G.xg; yi = G.yg; wd = S.bw;
        vx = Q(4); w = Q(6);
        u = vx - w*G.yr;  v = Q(5) + w*G.xr;
        [xI,yI,R] = icrGeom();

        % ---- chassis frame ----
        bodyX = [-L/2-.03 L/2+.03 L/2+.03 -L/2-.03];
        bodyY = [-S.b/2+wd/2 -S.b/2+wd/2 S.b/2-wd/2 S.b/2-wd/2];
        set(hBody,'XData',bodyX,'YData',bodyY);
        sp = r/2;
        for j = 1:4
            set(hWh(j),'XData',xi(j)+[-r r r -r],'YData',yi(j)+[-wd/2 -wd/2 wd/2 wd/2]);
            xs = -r + phase(j) + (0:3)*sp;
            X = [xs; xs; nan(1,4)]; Y = repmat([-wd/2; wd/2; NaN],1,4);
            set(hTr(j),'XData',xi(j)+X(:),'YData',yi(j)+Y(:));
            set(hTx(j),'Position',[xi(j), yi(j)+sign(yi(j))*(wd/2+0.2)], ...
                       'String',sprintf('%s\n%+.1f rad/s\nN=%.0f',names{j},om(j),Nld(j)));
            set(hRay(j),'XData',[xI xi(j)],'YData',[yI yi(j)]);
        end
        sv = 0.35; sf = 0.45/(S.m*g/4);
        set(hVel, 'XData',xi,'YData',yi,'UData',u*sv,'VData',v*sv);
        set(hRoll,'XData',xi,'YData',yi,'UData',r*om*sv,'VData',zeros(1,4));
        set(hSlip,'XData',xi,'YData',yi,'UData',lastSx*sv,'VData',lastSy*sv);
        set(hFric,'XData',xi,'YData',yi,'UData',lastFx*sf,'VData',lastFy*sf);
        set(hCG,'XData',G.xcg,'YData',0);
        if isfinite(xI) && abs(xI) < limB && abs(yI) < limB
            set(hICR,'XData',xI,'YData',yI);
        else
            set(hICR,'XData',nan,'YData',nan);
        end
        if isfinite(R) && R < 500
            set(hCirc,'XData',xI + R*cos(th),'YData',yI + R*sin(th));
        else
            set(hCirc,'XData',nan,'YData',nan);
        end

        % ---- world view ----
        [bx,by] = toWorld(bodyX,bodyY);
        set(hWBody,'XData',bx,'YData',by);
        for j = 1:4
            [wxj,wyj] = toWorld(xi(j)+[-r r r -r], yi(j)+[-wd/2 -wd/2 wd/2 wd/2]);
            set(hWWh(j),'XData',wxj,'YData',wyj);
        end
        set(hTrail,'XData',trail(:,1),'YData',trail(:,2));
        if isfinite(xI) && R < 500
            [cx,cy] = toWorld(xI + R*cos(th), yI + R*sin(th));
            set(hWCirc,'XData',cx,'YData',cy);
            [ix,iy] = toWorld(xI,yI);
            set(hWICR,'XData',ix,'YData',iy);
        else
            set(hWCirc,'XData',nan,'YData',nan);
            set(hWICR,'XData',nan,'YData',nan);
        end
        [px,py] = toWorld(0,0);
        set(axW,'XLim',px+[-limW limW],'YLim',py+[-limW limW]);
        if S.slope > 0.1
            bx0 = px - 0.8*limW; by0 = py + 0.8*limW;
            set(hUp,'XData',bx0,'YData',by0,'UData',0.5*cos(G.psi),'VData',0.5*sin(G.psi));
            set(hUpTx,'Position',[bx0 by0+0.15],'Visible','on');
        else
            set(hUp,'UData',0,'VData',0); set(hUpTx,'Visible','off');
        end
        drawnow limitrate;
    end

    function updateTable()
        if ~ishandle(fig), return; end
        rows = cell(4,11);
        for j = 1:4
            rows(j,:) = {names{j}, sprintf('%+.2f',om(j)), sprintf('%+.1f',om(j)*30/pi), ...
                sprintf('%+.1f',Icur(j)), sprintf('%+.2f',lastTauW(j)), ...
                sprintf('%+.3f',lastSx(j)), sprintf('%+.3f',lastSy(j)), ...
                sprintf('%.1f',Nld(j)), sprintf('%+.1f',lastFx(j)), sprintf('%+.1f',lastFy(j)), ...
                sprintf('%.1f',lastRr(j))};
        end
        set(hTab,'Data',rows);

        vx = Q(4); vy = Q(5); w = Q(6); r = S.r; b = S.b;
        [xI,yI,R] = icrGeom();
        ln = { ...
            sprintf('%s | slope %.1f deg', TER(S.terr).name, S.slope), ...
            sprintf('Command  v=%+.2f m/s  w=%+.2f rad/s', S.vcmd, S.wcmd), ...
            sprintf('Actual   vx=%+.3f vy=%+.3f w=%+.3f', vx, vy, w)};
        if hypot(vx,vy) > 0.05
            ln{end+1} = sprintf('Body slip angle   = %+.1f deg', rad2deg(atan2(vy,vx)));
        end
        if abs(w) > 0.02
            vL = r*mean(om([1 3])); vR = r*mean(om([2 4]));
            yL = (vx - vL)/w;  yR = (vx - vR)/w;
            ln{end+1} = sprintf('Vehicle ICR (%+.2f,%+.2f) m  R=%.2f m', xI, yI, R);
            ln{end+1} = sprintf('Track ICRs  yL=%+.3f  yR=%+.3f m', yL, yR);
            ln{end+1} = sprintf('Effective chi = %.2f', (yL-yR)/b);
            if abs(S.wcmd) > 0.02
                ln{end+1} = sprintf('Yaw achieved = %.0f %% of commanded', 100*w/S.wcmd);
            end
        else
            ln{end+1} = 'Vehicle ICR at infinity (straight)';
        end
        tipLat  = (S.b/2)/S.h;
        tipLong = min(S.L/2 - G.xcg, S.L/2 + G.xcg)/S.h;
        ln{end+1} = sprintf('Static tip limit: lat %.2f g  long %.2f g', tipLat, tipLong);
        if liftFlag
            ln{end+1} = '!! WHEEL LIFT-OFF / TIPPING RISK';
        end
        nlim = sum(abs(Icur) >= 0.98*S.Imax);
        if nlim > 0
            ln{end+1} = sprintf('!! CURRENT LIMIT on %d/4 motors', nlim);
        end
        Pel   = Vterm*Ibat;
        Pshaft = sum(lastTauW.*om);
        Pskid = sum(-(lastFx.*lastSx + lastFy.*lastSy));
        Proll = sum(lastTrr.*om);
        Pcu   = S.Rm*sum(Icur.^2) + S.Rbat*Ibat^2;
        ln{end+1} = sprintf('Battery %.1f V  %.1f A  = %.0f W', Vterm, Ibat, Pel);
        ln{end+1} = sprintf('Wheel shaft %.0f W | copper loss %.0f W', Pshaft, Pcu);
        ln{end+1} = sprintf('Skid loss %.0f W | rolling loss %.0f W', Pskid, Proll);
        set(hSum,'String',ln);
    end

    function mkSlider(group,key,label,lo,hi,fmt,y)
        t1 = uicontrol(fig,'Style','text','Units','normalized', ...
            'Position',[xR y+0.027 0.20 0.025],'String',label, ...
            'HorizontalAlignment','left','BackgroundColor','w');
        vt = uicontrol(fig,'Style','text','Units','normalized', ...
            'Position',[xR+0.22 y+0.027 0.08 0.025],'String',sprintf(fmt,S.(key)), ...
            'HorizontalAlignment','right','BackgroundColor','w');
        sl = uicontrol(fig,'Style','slider','Units','normalized', ...
            'Position',[xR y wR 0.025],'Min',lo,'Max',hi,'Value',S.(key), ...
            'Callback',@(src,~) onSlide(key,src,vt,fmt));
        sliders.(key) = struct('sl',sl,'vt',vt,'fmt',fmt);
        if ~strcmp(group,'cmd')
            groupH.(group) = [groupH.(group) t1 vt sl];
        end
    end

    function showGroup(name)
        f = fieldnames(groupH);
        for gi = 1:numel(f)
            if strcmp(f{gi},name), vis = 'on'; else, vis = 'off'; end
            set(groupH.(f{gi}),'Visible',vis);
        end
        for gi = 1:3
            if strcmp(gnames{gi},name), set(hGrp(gi),'FontWeight','bold');
            else, set(hGrp(gi),'FontWeight','normal'); end
        end
    end

    function onSlide(key,src,vt,fmt)
        S.(key) = get(src,'Value');
        set(vt,'String',sprintf(fmt,S.(key)));
        updateDerived();
    end

    function onTerrain(src)
        applyTerrain(get(src,'Value'));
        syncSliders(); updateDerived();
    end

    function applyTerrain(k)
        T = TER(k);
        S.terr = k; S.mx = T.mx; S.my = T.my; S.K = T.K; S.crr = T.crr; S.bek = T.bek;
        set(axW,'Color',T.col);
    end

    function syncSliders()
        f = fieldnames(sliders);
        for si = 1:numel(f)
            s = sliders.(f{si});
            set(s.sl,'Value',S.(f{si}));
            set(s.vt,'String',sprintf(s.fmt,S.(f{si})));
        end
        set(hTerr,'Value',S.terr);
    end

    function preset(v,w)
        S.vcmd = v; S.wcmd = w;
        syncSliders(); updateDerived();
    end

    function resetAll()
        S = S0; applyTerrain(S.terr); updateDerived(); syncSliders();
        Q = zeros(1,6); om = zeros(1,4); zi = zeros(1,4); zx = zeros(1,4); zy = zeros(1,4);
        Nld = G.Nstatic; Icur = zeros(1,4); Ibat = 0; Vterm = S.Vbat; liftFlag = false;
        trail = zeros(0,2); phase = zeros(1,4);
        lastFx = zeros(1,4); lastFy = zeros(1,4); lastSx = zeros(1,4); lastSy = zeros(1,4);
        lastTauW = zeros(1,4); lastTrr = zeros(1,4); lastRr = zeros(1,4);
    end

    function togglePause()
        running = ~running;
        if running, set(hBtn(6),'String','Pause'); else, set(hBtn(6),'String','Resume'); end
        updateTable();
    end

    function closeAll()
        try, stop(tm); delete(tm); catch, end
        delete(fig);
    end
end
