function LEHT(k_m, k_i, frac, field, x_cut, y_cut)
clc;
% =========================================================================
% CALCULATION OF THE EFFECTIVE THERMAL CONDUCTIVITY OF COMPOSITES WITH
% ISOTROPIC PHASES VIA LOCALLY-EXACT HOMOGENIZATION THEORY (LEHT)
% =========================================================================

% 1. DIMENSIONS AND PARAMETERS FOR PLOTS
L = 1;          % CELL LENGTH
H = 1;          % CELL HEIGHT
Hx_plot = 1;    % MACROSCOPIC TEMPERATURE GRADIENT IN X
Hy_plot = 0;    % MACROSCOPIC TEMPERATURE GRADIENT IN Y

% 2. EFFECTIVE THERMAL CONDUCTIVITY MATRIX K*
K_ast_LEHT = zeros(2);
for col = 1:2
    if col == 1
        Hx = 1; Hy = 0;
    else
        Hx = 0; Hy = 1;
    end

    % SOLVES FOR THE UNKNOWN EXPANSION COEFFICIENTS
    coef = LEHT_isotropic_square(L, H, L/2, H/2, sqrt(frac*L*H/pi), k_m, k_i, 30, 250, 400, Hx, Hy);

    % COMPUTES THE AVERAGE HEAT FLUX
    [qxbar, qybar] = mean_flux_LEHT_boundary(coef, L, H, k_m, Hx, Hy, 1000);

    % EFFECTIVE THERMAL CONDUCTIVITY MATRIX
    K_ast_LEHT(:,col) = -[qxbar; qybar];
end

fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRIX (K*) \n');
disp(K_ast_LEHT);

% 3. POST-PROCESSING AND VISUALIZATION
if field == 1 || (x_cut > 0 && x_cut <= 1) || (y_cut > 0 && y_cut <= 1)
    coef_plot = LEHT_isotropic_square(L, H, L/2, H/2, sqrt(frac*L*H/pi), k_m, k_i, 30, 250, 400, Hx_plot, Hy_plot);
    plot_thermal_results(coef_plot, L, H, Hx_plot, Hy_plot, field, x_cut, y_cut);
end
end

function plot_thermal_results(coef, L, H, Hx, Hy, field, x_cut, y_cut)


T_func = @(x, y) evaluate_temperature(x, y, coef, Hx, Hy);

if field == 1
    [Xg, Yg] = meshgrid(linspace(0, L, 150), linspace(0, H, 150));
    Z_LEHT = T_func(Xg, Yg);
    figure();
    surf(Xg, Yg, Z_LEHT, 'EdgeColor', 'none'); 
    view(2); axis equal tight;
    colormap('jet'); colorbar;
    title('Total temperature field - LEHT', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca,'FontSize',12);
end

if y_cut > 0 && y_cut <= 1
    xH = linspace(0, L, 500);
    y_eval = y_cut * ones(size(xH));
    TLEHT_X = T_func(xH, y_eval);
    figure();
    plot(xH, TLEHT_X, 'k-', 'LineWidth', 1.5);
    grid on; xlim([0 L]); box on;
    legend({'LEHT'}, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Temperature field at $x_2 = %.2f$', y_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end

if x_cut > 0 && x_cut <= 1
    yV = linspace(0, H, 500); 
    x_eval = x_cut * ones(size(yV));
    TLEHT_Y = T_func(x_eval, yV);
    figure();
    plot(yV, TLEHT_Y, 'k-', 'LineWidth', 1.5);
    grid on; xlim([0 H]); box on;
    ylim([0 1]);
    legend({'LEHT'}, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Temperature field at $x_1 = %.2f$', x_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end
end

% LEHT FUNCTIONS
function T_total = evaluate_temperature(x, y, coef, Hx, Hy)
rx = x - coef.xc; ry = y - coef.yc;
r = sqrt(rx.^2 + ry.^2); th = atan2(ry, rx);
T_tilde = zeros(size(x));
mask_f = (r <= coef.a); mask_m = ~mask_f;

if any(mask_f, 'all')
    xi_f = r(mask_f) / coef.a; th_f = th(mask_f); T_f = zeros(size(xi_f));
    for n = 1:coef.N
        T_f = T_f + coef.a * (coef.Af(n)*(xi_f.^n).*cos(n*th_f) + coef.Bf(n)*(xi_f.^n).*sin(n*th_f));
    end
    T_tilde(mask_f) = T_f;
end

if any(mask_m, 'all')
    xi_m = max(r(mask_m) / coef.a, 1e-14); th_m = th(mask_m); T_m = zeros(size(xi_m));
    for n = 1:coef.N
        T_m = T_m + coef.a * (coef.Am(n)*(xi_m.^n).*cos(n*th_m) + coef.Bm(n)*(xi_m.^n).*sin(n*th_m) + ...
            coef.Cm(n)*(xi_m.^(-n)).*cos(n*th_m) + coef.Dm(n)*(xi_m.^(-n)).*sin(n*th_m));
    end
    T_tilde(mask_m) = T_m;
end
T_total = Hx.*x + Hy.*y + T_tilde;
end

function coef = LEHT_isotropic_square(L,H,xc,yc,a,km,ki,N,Mbc,Mint,Hx,Hy)
yy = (1:Mbc)'/(Mbc+1)*H; xx = (1:Mbc)'/(Mbc+1)*L;
xL = zeros(size(yy)); yL = yy; xR = L + zeros(size(yy)); yR = yy;
xB = xx; yB = zeros(size(xx)); xT = xx; yT = H + zeros(size(xx));
thI = (0:Mint-1)'*(2*pi/Mint);

A_bcT = [ periodic_Ttilde_rows(xR,yR,xL,yL,xc,yc,a,N); periodic_Ttilde_rows(xT,yT,xB,yB,xc,yc,a,N) ];
A_bcQ = [ periodic_flux_rows('x',xR,yR,xL,yL,xc,yc,a,km,N); periodic_flux_rows('y',xT,yT,xB,yB,xc,yc,a,km,N)];
[A_intT, ~] = interface_T_rows(thI,a,N);
[A_intQ, b_intQ] = interface_qn_rows(thI,km,ki,N,Hx,Hy);

A = [A_bcT; A_bcQ; A_intT; A_intQ];
b = [zeros(size(A_bcT,1)+size(A_bcQ,1)+size(A_intT,1),1); b_intQ];

u = lsqminnorm(A,b);
coef.Af = u(1:N); coef.Bf = u(N+1:2*N);
coef.Am = u(2*N+1:3*N); coef.Bm = u(3*N+1:4*N);
coef.Cm = u(4*N+1:5*N); coef.Dm = u(5*N+1:6*N);
coef.N=N; coef.a=a; coef.km=km; coef.ki=ki; coef.xc=xc; coef.yc=yc;
end

function [qxbar,qybar] = mean_flux_LEHT_boundary(coef,L,H,km,Hx,Hy,M)
yy = (1:M)'/(M+1)*H; xx = (1:M)'/(M+1)*L;
xR = L * ones(size(yy)); [dTdxR, ~] = eval_grad_Ttil_matrix(xR, yy, coef);
qxbar = mean(-km * (Hx + dTdxR));
yT = H * ones(size(xx)); [~, dTdyT] = eval_grad_Ttil_matrix(xx, yT, coef);
qybar = mean(-km * (Hy + dTdyT));
end

function [dTdx, dTdy] = eval_grad_Ttil_matrix(x,y,coef)
rx = x - coef.xc; ry = y - coef.yc; r = sqrt(rx.^2 + ry.^2); th = atan2(ry,rx);
n = 1:coef.N; C = cos(th*n); S = sin(th*n); Cn = C.*n; Sn = S.*n;
xi = r/coef.a; xiP1 = bsxfun(@power, xi, n-1); xiM1 = bsxfun(@power, xi, -n-1);
dTdr = (xiP1.*Cn)*coef.Am + (xiP1.*Sn)*coef.Bm - (xiM1.*Cn)*coef.Cm - (xiM1.*Sn)*coef.Dm;
dTdth = coef.a * ( -(bsxfun(@power,xi,n).*Sn)*coef.Am + (bsxfun(@power,xi,n).*Cn)*coef.Bm - (bsxfun(@power,xi,-n).*Sn)*coef.Cm + (bsxfun(@power,xi,-n).*Cn)*coef.Dm );
ct = cos(th); st = sin(th); invr = 1./r;
dTdx = dTdr.*ct - (dTdth.*st).*invr; dTdy = dTdr.*st + (dTdth.*ct).*invr;
end

function Arows = periodic_Ttilde_rows(x2,y2,x1,y1,xc,yc,a,N)
[PhiAm2,PhiBm2,PhiCm2,PhiDm2] = matrix_Ttilde_basis(x2,y2,xc,yc,a,N);
[PhiAm1,PhiBm1,PhiCm1,PhiDm1] = matrix_Ttilde_basis(x1,y1,xc,yc,a,N);
Arows = [zeros(size(PhiAm2,1),2*N), PhiAm2-PhiAm1, PhiBm2-PhiBm1, PhiCm2-PhiCm1, PhiDm2-PhiDm1];
end

function [PhiAm,PhiBm,PhiCm,PhiDm] = matrix_Ttilde_basis(x,y,xc,yc,a,N)
rx = x - xc; ry = y - yc; r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
xi = max(r/a, 1e-14); n = 1:N;
PhiAm = a*(bsxfun(@power,xi,n).*cos(th*n)); PhiBm = a*(bsxfun(@power,xi,n).*sin(th*n));
PhiCm = a*(bsxfun(@power,xi,-n).*cos(th*n)); PhiDm = a*(bsxfun(@power,xi,-n).*sin(th*n));
end

function Arows = periodic_flux_rows(dir,x2,y2,x1,y1,xc,yc,a,km,N)
[QxAm2,QxBm2,QxCm2,QxDm2, QyAm2,QyBm2,QyCm2,QyDm2] = matrix_flux_basis(x2,y2,xc,yc,a,km,N);
[QxAm1,QxBm1,QxCm1,QxDm1, QyAm1,QyBm1,QyCm1,QyDm1] = matrix_flux_basis(x1,y1,xc,yc,a,km,N);
if dir == 'x'
    d = [QxAm2-QxAm1, QxBm2-QxBm1, QxCm2-QxCm1, QxDm2-QxDm1];
else
    d = [QyAm2-QyAm1, QyBm2-QyBm1, QyCm2-QyCm1, QyDm2-QyDm1];
end
Arows = [zeros(size(d,1),2*N), d];
end

function [QxAm,QxBm,QxCm,QxDm, QyAm,QyBm,QyCm,QyDm] = matrix_flux_basis(x,y,xc,yc,a,km,N)
rx = x - xc; ry = y - yc; r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
xi = max(r/a, 1e-14); n = 1:N; C = cos(th*n); S = sin(th*n);
ct = cos(th); st = sin(th); invr = 1./max(r, 1e-14);

dTdr_Am = bsxfun(@power,xi,n-1).*C.*n; dTdth_Am = -a*bsxfun(@power,xi,n).*S.*n;
dTdr_Bm = bsxfun(@power,xi,n-1).*S.*n; dTdth_Bm =  a*bsxfun(@power,xi,n).*C.*n;
dTdr_Cm = -bsxfun(@power,xi,-n-1).*C.*n; dTdth_Cm = -a*bsxfun(@power,xi,-n).*S.*n;
dTdr_Dm = -bsxfun(@power,xi,-n-1).*S.*n; dTdth_Dm =  a*bsxfun(@power,xi,-n).*C.*n;

QxAm = -km*(dTdr_Am.*ct - dTdth_Am.*st.*invr); QxBm = -km*(dTdr_Bm.*ct - dTdth_Bm.*st.*invr);
QxCm = -km*(dTdr_Cm.*ct - dTdth_Cm.*st.*invr); QxDm = -km*(dTdr_Dm.*ct - dTdth_Dm.*st.*invr);
QyAm = -km*(dTdr_Am.*st + dTdth_Am.*ct.*invr); QyBm = -km*(dTdr_Bm.*st + dTdth_Bm.*ct.*invr);
QyCm = -km*(dTdr_Cm.*st + dTdth_Cm.*ct.*invr); QyDm = -km*(dTdr_Dm.*st + dTdth_Dm.*ct.*invr);
end

function [Arows, brows] = interface_T_rows(thI,a,N)
n = 1:N; C = cos(thI*n); S = sin(thI*n);
Arows = [ -a*C, -a*S,  a*C,  a*S,  a*C,  a*S ];
brows = zeros(size(thI));
end

function [Arows, brows] = interface_qn_rows(thI,km,ki,N,Hx,Hy)
n = 1:N; C = cos(thI*n); S = sin(thI*n);
Arows = [ -ki*C.*n, -ki*S.*n,  km*C.*n,  km*S.*n, -km*C.*n, -km*S.*n ];
brows = -(km-ki)*(Hx*cos(thI) + Hy*sin(thI));
end