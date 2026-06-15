% Confinement with active control - RK4
clear
clc
close all

% Particle
q = 1.6e-19;
m = 1.67e-27;
vi = 1e4; 
dir = rand*2*pi;
v0 = [vi*cos(dir), vi*sin(dir)];
r0 = [0, 0];

% Reactor
N_cables = 12;
R_toroid = 0.1;
R_limit = 0.05;
I_idle = 100;
I_rescue = 5000;

% Cables position
angls_c = linspace(0, 2*pi, N_cables + 1);
angls_c(end) = [];
x_c = R_toroid * cos(angls_c);
y_c = R_toroid * sin(angls_c);

% Simulation
dt = 1e-9;
steps = 100000;
r = zeros(steps,2);
v = zeros(steps,2);
r(1,:) = r0;
v(1,:) = v0;
I = ones(1,N_cables) * I_idle;

% ─── FIGURE SETUP (Dark Theme) ───────────────────────────────────────
fig = figure('Color', [0.1 0.1 0.1], 'Name', 'RK4 Confinement', 'Position', [100 100 700 700],'WindowState', 'maximized');
pause(1);
ax = axes('Parent', fig, 'Color', [0.1 0.1 0.1], 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8]);
hold(ax, 'on');
grid(ax, 'on');
ax.GridColor = [0.4 0.4 0.4];
axis(ax, 'equal');
xlim(ax, [-0.12 0.12]);
ylim(ax, [-0.12 0.12]);

title('Magnetic Confinement with Active Control', 'Color', 'w', 'FontSize', 14);
xlabel('x (m)', 'Color', 'w', 'FontSize', 12);
ylabel('y (m)', 'Color', 'w', 'FontSize', 12);

% Draw Limit Ring (Soft Red)
theta_circ = linspace(0, 2*pi, 200);
plot(ax, R_limit*cos(theta_circ), R_limit*sin(theta_circ), '--', 'Color', [1 0.3 0.3 0.6], 'LineWidth', 1.5);

% Draw Cables
hCables = zeros(1, N_cables);
for i = 1:N_cables
    hCables(i) = plot(ax, x_c(i), y_c(i), 'o', ...
        'MarkerFaceColor', [0.3 0.3 0.3], ...
        'MarkerEdgeColor', 'w', ...
        'MarkerSize', 10);
end

% Animated Line for Trajectory & Marker for Particle
hTrajectory = animatedline(ax, 'Color', [0 1 1], 'LineWidth', 1.2); % Cyan line
hParticle = plot(ax, r0(1), r0(2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', 'MarkerSize', 6);

% UI Status Text
hStatus = text(ax, -0.115, 0.11, 'Status: IDLE', 'Color', [0 1 0], 'FontSize', 12, 'FontWeight', 'bold');

pasos_rescate = 0;

% ── Función auxiliar: calcula aceleración dada posición y corrientes ──
    function a = calcAccel(pos, vel, q, m, x_c, y_c, I, N_cables)
        mu0 = 4*pi*1e-7;
        Bz = 0;
        for j = 1:N_cables
            dx = pos(1) - x_c(j);
            dy = pos(2) - y_c(j);
            dist = sqrt(dx^2 + dy^2);
            if dist < 1e-5, dist = 1e-5; end
            if mod(j,2)==0
                corriente = I(j);
            else
                corriente = -I(j);
            end
            Bz = Bz + (mu0 * corriente) / (2*pi*dist);
        end
        a = [ (q*vel(2)*Bz)/m, -(q*vel(1)*Bz)/m ];
    end

% ─── MAIN LOOP ───────────────────────────────────────────────────────
for n = 1:steps-1
    dist_center = norm(r(n,:));
    I(:) = I_idle;
    
    % ACTIVE CONTROL
    if dist_center > R_limit
        pasos_rescate = pasos_rescate + 1;
        ang_p = atan2(r(n,2), r(n,1));
        if ang_p < 0, ang_p = ang_p + 2*pi; end
        [~, idx] = min(abs(angls_c - ang_p));
        
        idx_izq = mod(idx - 2, N_cables) + 1;
        idx_der = mod(idx,   N_cables) + 1;
        
        I(idx)     = I_rescue;
        I(idx_izq) = I_rescue * 0.5;
        I(idx_der) = I_rescue * 0.5;
        
        % Visual Update: Active
        set(hCables, 'MarkerFaceColor', [0.3 0.3 0.3]);
        set(hCables(idx),              'MarkerFaceColor', [1 0 0]); % Red
        set(hCables([idx_izq idx_der]),'MarkerFaceColor', [1 0.5 0]); % Orange
        set(hStatus, 'String', 'Status: ACTIVE CONTROL', 'Color', [1 0 0]);
    else
        % Visual Update: Idle
        set(hCables, 'MarkerFaceColor', [0.3 0.3 0.3]);
        set(hStatus, 'String', 'Status: IDLE', 'Color', [0 1 0]);
    end
    
    % ── RK4 ──────────────────────────────────────────────────────────
    rn = r(n,:);
    vn = v(n,:);
    
    k1_v = vn;
    k1_a = calcAccel(rn, vn, q, m, x_c, y_c, I, N_cables);
    
    k2_v = vn + 0.5*dt*k1_a;
    k2_a = calcAccel(rn + 0.5*dt*k1_v, k2_v, q, m, x_c, y_c, I, N_cables);
    
    k3_v = vn + 0.5*dt*k2_a;
    k3_a = calcAccel(rn + 0.5*dt*k2_v, k3_v, q, m, x_c, y_c, I, N_cables);
    
    k4_v = vn + dt*k3_a;
    k4_a = calcAccel(rn + dt*k3_v, k4_v, q, m, x_c, y_c, I, N_cables);
    
    r(n+1,:) = rn + (dt/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
    v(n+1,:) = vn + (dt/6)*(k1_a + 2*k2_a + 2*k3_a + k4_a);
    % ─────────────────────────────────────────────────────────────────
    
    % Actualización eficiente de gráficos
    if mod(n, 100) == 0
        addpoints(hTrajectory, r(n+1,1), r(n+1,2));
        set(hParticle, 'XData', r(n+1,1), 'YData', r(n+1,2));
        drawnow limitrate; % 'limitrate' previene saturación gráfica
    end
end

% Añadir los últimos puntos asegurando que se dibuje completo
addpoints(hTrajectory, r(end,1), r(end,2));
set(hParticle, 'XData', r(end,1), 'YData', r(end,2));
drawnow;

porcentaje_rescate = (pasos_rescate / (steps-1)) * 100;
fprintf('La corriente de rescate se activó el %.2f%% del tiempo total.\n', porcentaje_rescate);