function skidSteerSim()
% SKIDSTEERSIM  Four-wheel skid-steer kinematics playground.
%
% Edit defaultParams() below, then run:  >> skidSteerSim
%
% Models (all driven by the SAME commanded wheel speeds, so you can compare):
%   'ideal' : textbook differential drive, effective track = b
%   'chi'   : ICR model, effective track = chi*b  (chi >= 1 typical)
%   'slip'  : chi model + longitudinal slip ratios on left/right sides
%   'icr'   : full Instantaneous-Center-of-Rotation model (Mandow et al.)
%             yL, yR = lateral ICR coords of left/right tracks,
%             xv = longitudinal ICR offset -> produces lateral skid vy
%
% Body-frame convention: x forward, y left, w = yaw rate (CCW +).
% Track velocities vL, vR (m/s) -> body velocity:
%   w  = (vR - vL)/(yL - yR)
%   vx = (yL*vR - yR*vL)/(yL - yR)
%   vy = -xv*w
% Ideal diff-drive is the case yL = b/2, yR = -b/2, xv = 0.

close all;
P = defaultParams();

t = (0:P.dt:P.T)';
[wl, wr, vCmd, wCmd] = commandProfile(P, t);

nM = numel(P.models);
R = cell(1, nM);
for k = 1:nM
    R{k} = simulate(P.models{k}, P, t, wl, wr);
end

printSummary(P, R);
plotResults(P, t, R, wl, wr, vCmd, wCmd);
if P.animate,     animateRun(P, t, R); end
if P.sweep.enable, runSweep(P, t, wl, wr); end
end

%% ========================================================================
%  PARAMETERS  (edit here)
%  ========================================================================
function P = defaultParams()
% --- Geometry ---
P.r  = 0.10;     % wheel radius [m]
P.b  = 0.40;     % nominal track width (left-right wheel spacing) [m]
P.L  = 0.30;     % wheelbase (front-rear wheel spacing) [m]
P.ww = 0.05;     % wheel width, drawing only [m]

% --- Time ---
P.dt = 0.01;     % step [s]
P.T  = 30;       % duration [s]

% --- Initial pose [x y theta] ---
P.x0 = [0 0 0];

% --- Actuators ---
P.wmax = 20;     % wheel speed saturation [rad/s]
P.tau  = 0.15;   % motor first-order lag time const [s] (0 = off)

% --- Command profile ---
% 'circle' | 'figure8' | 'sine' | 'square' | 'spin' | 'straight'
% | 'custom' (c.fcn(t) -> [v,w]) | 'wheels' (c.wheelFcn(t) -> [wl,wr] rad/s)
P.cmd.profile = 'figure8';
P.cmd.v       = 0.6;      % forward speed [m/s]
P.cmd.w       = 0.5;      % yaw rate [rad/s]
P.cmd.period  = 15;       % period for figure8 / sine [s]
P.cmd.fcn      = @(t) deal(0.5*ones(size(t)), 0.4*sin(0.5*t));
P.cmd.wheelFcn = @(t) deal(4*ones(size(t)), 5*ones(size(t)));

% --- Kinematic models to compare ---
% mkModel(name, type, 'param', value, ...)
% params: chi, iL, iR, yL, yR, xv
P.models = {
    mkModel('Ideal diff-drive', 'ideal')
    mkModel('ICR chi=1.8',      'chi',  'chi', 1.8)
    mkModel('Slip + chi',       'slip', 'chi', 1.5, 'iL', 0.05, 'iR', 0.10)
    mkModel('Full ICR',         'icr',  'yL', 0.28, 'yR', -0.20, 'xv', 0.05, ...
                                        'iL', 0.03, 'iR', 0.05)
    };

% --- Output ---
P.animate  = true;
P.animSkip = 8;       % draw every N-th sample
P.sweep.enable = true;            % sweep chi on the 'chi' model
P.sweep.chi    = [1 1.2 1.5 2 2.5 3];
end

function m = mkModel(name, type, varargin)
m = struct('name', name, 'type', type, 'chi', 1, 'iL', 0, 'iR', 0, ...
           'yL', 0.2, 'yR', -0.2, 'xv', 0);
for i = 1:2:numel(varargin)
    m.(varargin{i}) = varargin{i+1};
end
end

%% ========================================================================
%  COMMAND GENERATION
%  ========================================================================
function [wl, wr, v, w] = commandProfile(P, t)
c = P.cmd; N = numel(t);
v = zeros(N,1); w = zeros(N,1);

switch lower(c.profile)
    case 'circle'
        v(:) = c.v;  w(:) = c.w;
    case 'figure8'
        v(:) = c.v;  w = c.w * tanh(6*sin(2*pi*t/c.period));
    case 'sine'
        v(:) = c.v;  w = c.w * sin(2*pi*t/c.period);
    case 'square'
        Ts = 2.0;                      % straight leg duration [s]
        Tt = (pi/2)/abs(c.w);          % 90 deg turn duration [s]
        ph = mod(t, Ts+Tt);
        straight = ph < Ts;
        v(straight)  = c.v;
        w(~straight) = c.w;
    case 'spin'
        w(:) = c.w;
    case 'straight'
        v(:) = c.v;
    case 'custom'
        [v, w] = c.fcn(t);  v = v(:); w = w(:);
    case 'wheels'
        [wl, wr] = c.wheelFcn(t);  wl = wl(:); wr = wr(:);
        [wl, wr] = actuators(P, wl, wr);
        v = P.r*(wl+wr)/2;  w = P.r*(wr-wl)/P.b;   % ideal equivalents
        return
    otherwise
        error('Unknown profile "%s"', c.profile);
end

% Ideal inverse kinematics (what a naive controller would command)
vL = v - w*P.b/2;
vR = v + w*P.b/2;
[wl, wr] = actuators(P, vL/P.r, vR/P.r);
end

function [wl, wr] = actuators(P, wl, wr)
wl = max(min(wl, P.wmax), -P.wmax);
wr = max(min(wr, P.wmax), -P.wmax);
if P.tau > 0
    a = P.dt/(P.tau + P.dt);
    for k = 2:numel(wl)
        wl(k) = wl(k-1) + a*(wl(k) - wl(k-1));
        wr(k) = wr(k-1) + a*(wr(k) - wr(k-1));
    end
end
end

%% ========================================================================
%  KINEMATIC MODELS
%  ========================================================================
function [vx, vy, w, info] = bodyVel(m, P, wl, wr)
b = P.b;
vLi = P.r*wl;  vRi = P.r*wr;          % ideal (no-slip) track speeds

switch m.type
    case 'ideal'
        vL = vLi; vR = vRi;
        yL = b/2; yR = -b/2; xv = 0;
    case 'chi'
        vL = vLi; vR = vRi;
        yL = m.chi*b/2; yR = -yL; xv = 0;
    case 'slip'
        vL = (1-m.iL)*vLi; vR = (1-m.iR)*vRi;
        yL = m.chi*b/2; yR = -yL; xv = 0;
    case 'icr'
        vL = (1-m.iL)*vLi; vR = (1-m.iR)*vRi;
        yL = m.yL; yR = m.yR; xv = m.xv;
    otherwise
        error('Unknown model type "%s"', m.type);
end

d  = yL - yR;
w  = (vR - vL)/d;
vx = (yL*vR - yR*vL)/d;
vy = -xv*w;
info.slipAngle = atan2(vy, max(abs(vx),1e-9)*sign(vx + (vx==0)));
info.vL = vL;  info.vR = vR;
end

function S = simulate(m, P, t, wl, wr)
N = numel(t);
X = zeros(N,3); V = zeros(N,3); beta = zeros(N,1);
q = P.x0(:)';

for k = 1:N
    [vx, vy, w, info] = bodyVel(m, P, wl(k), wr(k));
    X(k,:) = q;  V(k,:) = [vx vy w];  beta(k) = info.slipAngle;
    if k < N
        dt  = t(k+1) - t(k);
        thm = q(3) + 0.5*w*dt;                     % midpoint heading
        q(1) = q(1) + (vx*cos(thm) - vy*sin(thm))*dt;
        q(2) = q(2) + (vx*sin(thm) + vy*cos(thm))*dt;
        q(3) = q(3) + w*dt;
    end
end
S.name = m.name; S.X = X; S.V = V; S.beta = beta;
end

%% ========================================================================
%  OUTPUT
%  ========================================================================
function printSummary(P, R)
fprintf('\n%-20s %9s %9s %10s %10s\n', 'Model', 'x_end', 'y_end', 'th_end[deg]', 'pos err');
ref = R{1}.X(end,1:2);
for k = 1:numel(R)
    q = R{k}.X(end,:);
    fprintf('%-20s %9.3f %9.3f %10.1f %10.3f\n', R{k}.name, q(1), q(2), ...
            rad2deg(q(3)), norm(q(1:2)-ref));
end
fprintf('(pos err is relative to model 1: %s)\n\n', R{1}.name);
end

function plotResults(P, t, R, wl, wr, vCmd, wCmd)
nM = numel(R); col = lines(nM);

figure('Name','Skid-steer results','Position',[80 80 1250 700]);

subplot(2,3,1); hold on; grid on; axis equal;
for k = 1:nM
    plot(R{k}.X(:,1), R{k}.X(:,2), 'Color', col(k,:), 'LineWidth', 1.6);
    plot(R{k}.X(end,1), R{k}.X(end,2), 'o', 'Color', col(k,:), 'MarkerFaceColor', col(k,:));
end
xlabel('x [m]'); ylabel('y [m]'); title('Trajectories');
legend(cellfun(@(s) s.name, R, 'uni', 0), 'Location', 'best');

subplot(2,3,2); hold on; grid on;
plot(t, wl, 'b', t, wr, 'r');
xlabel('t [s]'); ylabel('rad/s'); title('Commanded wheel speeds'); legend('left','right');

subplot(2,3,3); hold on; grid on;
plot(t, vCmd, 'k--', 'LineWidth', 1);
for k = 1:nM, plot(t, R{k}.V(:,1), 'Color', col(k,:)); end
xlabel('t [s]'); ylabel('v_x [m/s]'); title('Forward velocity (dashed = command)');

subplot(2,3,4); hold on; grid on;
plot(t, wCmd, 'k--', 'LineWidth', 1);
for k = 1:nM, plot(t, R{k}.V(:,3), 'Color', col(k,:)); end
xlabel('t [s]'); ylabel('\omega [rad/s]'); title('Yaw rate (dashed = command)');

subplot(2,3,5); hold on; grid on;
for k = 1:nM, plot(t, rad2deg(R{k}.beta), 'Color', col(k,:)); end
xlabel('t [s]'); ylabel('deg'); title('Body slip angle atan2(v_y, v_x)');

subplot(2,3,6); hold on; grid on;
for k = 2:nM
    e = vecnorm(R{k}.X(:,1:2) - R{1}.X(:,1:2), 2, 2);
    plot(t, e, 'Color', col(k,:), 'LineWidth', 1.4);
end
xlabel('t [s]'); ylabel('m'); title(['Position error vs ' R{1}.name]);
if nM > 1, legend(cellfun(@(s) s.name, R(2:end), 'uni', 0), 'Location', 'northwest'); end
end

function animateRun(P, t, R)
nM = numel(R); col = lines(nM);
allX = cell2mat(cellfun(@(s) s.X(:,1:2), R(:), 'uni', 0));
pad = 1;
figure('Name','Animation','Position',[200 120 800 700]);
axis equal; grid on; hold on;
axis([min(allX(:,1))-pad max(allX(:,1))+pad min(allX(:,2))-pad max(allX(:,2))+pad]);
xlabel('x [m]'); ylabel('y [m]');

hPath = gobjects(nM,1); hBody = gobjects(nM,1); hWh = gobjects(nM,4);
for k = 1:nM
    hPath(k) = plot(nan, nan, 'Color', col(k,:), 'LineWidth', 1.3);
    hBody(k) = patch(nan, nan, col(k,:), 'FaceAlpha', 0.25, 'EdgeColor', col(k,:));
    for j = 1:4
        hWh(k,j) = patch(nan, nan, 'k', 'FaceAlpha', 0.8);
    end
end
legend(hPath, cellfun(@(s) s.name, R, 'uni', 0), 'Location', 'best', 'AutoUpdate', 'off');

for i = 1:P.animSkip:numel(t)
    for k = 1:nM
        q = R{k}.X(i,:);
        set(hPath(k), 'XData', R{k}.X(1:i,1), 'YData', R{k}.X(1:i,2));
        [bx, by, wx, wy] = robotPolys(q, P);
        set(hBody(k), 'XData', bx, 'YData', by);
        for j = 1:4, set(hWh(k,j), 'XData', wx(:,j), 'YData', wy(:,j)); end
    end
    title(sprintf('t = %.2f s', t(i)));
    drawnow limitrate;
end
end

function [bx, by, wx, wy] = robotPolys(q, P)
rect = @(cx, cy, lx, ly) [cx-lx/2 cx+lx/2 cx+lx/2 cx-lx/2; ...
                          cy-ly/2 cy-ly/2 cy+ly/2 cy+ly/2];
Rz = [cos(q(3)) -sin(q(3)); sin(q(3)) cos(q(3))];
tf = @(p) Rz*p + q(1:2)'; 

B  = tf(rect(0, 0, P.L + 0.02, P.b - P.ww));
bx = B(1,:)'; by = B(2,:)';
wx = zeros(4,4); wy = zeros(4,4); j = 0;
for sx = [1 -1]
    for sy = [1 -1]
        j = j + 1;
        Wp = tf(rect(sx*P.L/2, sy*P.b/2, 2*P.r, P.ww));
        wx(:,j) = Wp(1,:)'; wy(:,j) = Wp(2,:)';
    end
end
% patch wants columns of vertices; 4 vertices x 4 wheels
end

function runSweep(P, t, wl, wr)
figure('Name','chi sweep','Position',[260 140 700 600]); hold on; grid on; axis equal;
col = parula(numel(P.sweep.chi));
% reference: ideal
Ri = simulate(mkModel('ideal','ideal'), P, t, wl, wr);
plot(Ri.X(:,1), Ri.X(:,2), 'k--', 'LineWidth', 1.5);
names = {'ideal'};
for i = 1:numel(P.sweep.chi)
    Rs = simulate(mkModel('chi','chi','chi',P.sweep.chi(i)), P, t, wl, wr);
    plot(Rs.X(:,1), Rs.X(:,2), 'Color', col(i,:), 'LineWidth', 1.5);
    names{end+1} = sprintf('\\chi = %.2f', P.sweep.chi(i)); %#ok<AGROW>
end
legend(names, 'Location', 'best');
xlabel('x [m]'); ylabel('y [m]');
title('Effect of effective-track factor \chi on trajectory');
end
