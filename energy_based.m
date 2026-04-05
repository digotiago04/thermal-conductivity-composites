function energy_based(nx, ny, k_m, k_i, frac, field, x_cut, y_cut)
% =========================================================================
% EFFECTIVE THERMAL CONDUCTIVITY AND TEMPERATURE FIELDS
% COMPARING LEHT (ANALYTICAL) AND FVT - ENERGY BASED (NUMERICAL)
% =========================================================================
clc;

% 1. DIMENSIONS AND PARAMETERS
L = 1; H = 1;            % Dimensions of the Representative Unit Cell (RUC)      
Hx_plot = 1;             % Macroscopic temperature gradient in x1 for plots
Hy_plot = 0;             % Macroscopic temperature gradient in x2 for plots
l = L/nx; h = H/ny;      % Subvolume dimensions
R = sqrt(frac*L*H/pi);   % Radius of the inclusion

% 2. FINITE-VOLUME THEORY
% Periodic boundary conditions and degrees of fredom
[i, j] = meshgrid(1:nx, 1:ny);
s = (i + (j - 1)*nx)';   % Total number of subvolumes
faces = [s(:), s(:)+(nx*ny)+1, s(:)+nx, s(:)+(nx*ny)];
faces(end-nx+1:end, 3) = faces(1:nx, 1);
faces(nx:nx:end, 2)    = faces(1:nx:end-nx+1, 4);
ndof = max(faces(:));   % Total number of faces
dofIS = unique([faces(1:nx, 1); faces(end-nx+1:end, 3)]);
dofIS = [dofIS(1), dofIS(end)];
dofDE = unique([faces(nx:nx:end, 2); faces(1:nx:end-nx+1, 4)]);
dofDE = [dofDE(1), dofDE(end)];
fixed = [dofIS dofDE];          % Fixed degrees
free = setdiff(1:ndof, fixed);  % Free degrees

% Sparse mapping indices
iK = reshape(kron(faces, ones(4,1))', 16*(nx*ny), 1);
jK = reshape(kron(faces, ones(1,4))', 16*(nx*ny), 1);
iF = repmat(faces', 2, 1);
jF = [ones(4, (nx*ny)); 2*ones(4, (nx*ny))];
k0 = eye(2);

% Auxiliary matrices
a = ones(4,1);
N1 = [0,-1]; N2 = [1,0]; N3 = [0,1]; N4 = [-1,0];
N = [N1, zeros(1,6); zeros(1,2), N2, zeros(1,4); ...
    zeros(1,4), N3, zeros(1,2); zeros(1,6), N4];
A = [0 -h/2 0 h^2/4; l/2 0 l^2/4 0; ...
    0 h/2 0 h^2/4; -l/2 0 l^2/4 0];
E = [0 0 0 0; 0 -1 0 3*h/2; -1 0 -3*l/2 0; 0 0 0 0; ...
    0 0 0 0; 0 -1 0 -3*h/2; -1 0 3*l/2 0; 0 0 0 0];
B = N * eye(8) * E;
ab = (B * (A\a)) \ (B/A);
Ab = A \ (eye(4) - a*ab);
K0 = B * Ab;
K0 = [K0(1,:)*l; K0(2,:)*h; K0(3,:)*l; K0(4,:)*h]; % Local thermal conductivity matrix

% Local load vector
H0 = [N1*l; N2*h; N3*l; N4*h] * k0;

% Material Design
[I, J] = ndgrid(1:nx, 1:ny);
Ic = I(:);
Jc = J(:);
% T0_x: Matrix containing the X1 coordinate of the center of each face.
T0_x = [(0.5 + Ic - 1)*l, (1 + Ic - 1)*l, (0.5 + Ic - 1)*l, (Ic - 1)*l];
% T0_y: Matrix containing the X2 coordinate of the center of each face.
T0_y = [ (Jc - 1)*h, (Jc - 0.5)*h,   (Jc)*h, (Jc - 0.5)*h ];
inclusion = ((I - 0.5)*l - L/2).^2 + ((J - 0.5)*h - H/2).^2 < R^2;
x = k_m * ones(nx, ny);
x(inclusion) = k_i;

% Material Interpolation
sK = K0(:) * x(:)';
sF = H0(:) * x(:)';

% Global thermal conductivity matrix
K = sparse(iK, jK, sK, ndof, ndof);
K = (K + K') / 2;

% Assembly of heat flux vectors corresponding to two unit temperature gradient tests
Q0 = sparse(iF(:), jF(:), sF, ndof, 2);

% Compute fluctuating temperatures for two unit temperature gradient tests
Tf = zeros(ndof, 2);
Tf(free, :) = K(free, free) \ Q0(free, :); 
Tfx = Tf(:, 1);
Tfy = Tf(:, 2);

Tx = T0_x + Tfx(faces);
Ty = T0_y + Tfy(faces);
T = cat(3, Tx, Ty); 

% HOMOGENIZATION BASED ON ENERGY THEORY
C  = zeros(2);
for dir_i = 1:2
    Ti = T(:, :, dir_i);
    for dir_j = 1:2
        Tj = T(:, :, dir_j);
        sumE = sum((Ti * K0) .* Tj, 2) / (L * H);
        C(dir_i, dir_j) = -sum(x(:) .* sumE);
    end
end

% 3. LEHT
K_ast_LEHT = zeros(2);
for col = 1:2
    if col==1, Hx=1; Hy=0; else, Hx=0; Hy=1; end
    coef = LEHT_isotropic_square(L, H, L/2, H/2, R, k_m, k_i, 30, 250, 400, Hx, Hy); 
    [qxbar, qybar] = mean_flux_LEHT_boundary(coef, L, H, k_m, Hx, Hy, 1000);
    K_ast_LEHT(:,col) = -[qxbar; qybar];
end

% PRINT OUTPUT (K*)
fprintf('====================================================\n');
fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRICES (K*)\n');
fprintf('====================================================\n');
fprintf('(FVT - BASED ON ENERGY THEORY) =\n');
disp(C);
fprintf('(LEHT - ANALYTICAL) =\n');
disp(K_ast_LEHT);


% 4. POST-PROCESSING FOR PLOTS (FIELDS & CUTS)

if field == 1 || (x_cut > 0 && x_cut <= 1) || (y_cut > 0 && y_cut <= 1)
coef_plot = LEHT_isotropic_square(L, H, L/2, H/2, R, k_m, k_i, 30, 250, 400, Hx_plot, Hy_plot);
end
% FVT Total Temperature Faces
Ttot_faces = Hx_plot * Tx + Hy_plot * Ty;

% ----- Extract Cuts -----
px = max(1, min(nx, round(x_cut / l) + 1));
py = max(1, min(ny, round(y_cut / h) + 1));
xcourt_val = (px-1)*l;
ycourt_val = (py-1)*h;

% Vertical Cut (x_cut)
svX = px + (0:ny-1)'*nx;
yV  = ((1:ny)' - 0.5)*h;
TFVT_Y = Ttot_faces(svX, 4);
TLEHT_Y = evaluate_temperature(xcourt_val*ones(size(yV)), yV, coef_plot, Hx_plot, Hy_plot);

% Horizontal Cut (y_cut)
svY = (1:nx)' + (py-1)*nx;
xH  = ((1:nx)' - 0.5)*l;
TFVT_X = Ttot_faces(svY, 1);
TLEHT_X = evaluate_temperature(xH, ycourt_val*ones(size(xH)), coef_plot, Hx_plot, Hy_plot);

% ----- Nodal Temperature -----
Tcell_tot = Ttot_faces;   
invA = A \ eye(4);                        
T00  = (ab*Tcell_tot')';                      
Taux = (Tcell_tot' - a*T00');                 
TijM = (invA*Taux)';                          
vert = [-l/2,-h/2; l/2,-h/2; l/2, h/2; -l/2, h/2];
phi = zeros(4,4); 
for p = 1:4
    x1 = vert(p,1); x2 = vert(p,2);
    phi(:,p) = [x1; x2; 0.5*(3*(x1^2) - (l^2)/4); 0.5*(3*(x2^2) - (h^2)/4)];
end
Tnodes = T00(:,ones(1,4)) + TijM*phi;    

nNx = nx+1;  nNy = ny+1;                     
nidBL = Ic + (Jc-1)*nNx;
nidBR = (Ic+1) + (Jc-1)*nNx;
nidTR = (Ic+1) + (Jc)*nNx;
nidTL = Ic + (Jc)*nNx;
nodeIDs = [nidBL, nidBR, nidTR, nidTL];       
sumT = accumarray(nodeIDs(:), Tnodes(:), [nNx*nNy,1], @sum, 0);
cntT = accumarray(nodeIDs(:), 1,              [nNx*nNy,1], @sum, 0);
Tnod = sumT ./ max(cntT,1);
Z_FVT = reshape(Tnod, [nNx, nNy])';           
[Xg,Yg] = meshgrid((0:nx)*l, (0:ny)*h);

% 5. PLOTTING

% Total Temperature Field - FVT
if field == 1
figure();
surf(Xg, Yg, Z_FVT);
view(2); shading flat; axis equal tight;
colormap('jet'); colorbar;
xlabel('$x_1$','Interpreter','latex','FontSize',16,'FontWeight','bold');
ylabel('$x_2$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('Total temperature field - FVT', 'Interpreter','latex','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',12);
end

% Micro-temperature fields (Cuts)
if y_cut > 0 && y_cut <= 1
figure();
plot(xH, TLEHT_X, 'k-', 'LineWidth', 1); hold on;
plot(xH, TFVT_X, 'bo', 'MarkerIndices',...
    1:3:length(xH)); % FVT adjustment points, for example: 1:2, 1:3, 1:4
grid on; xlim([0 L]); box on;
legend({'LEHT', 'FVT'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $x_2 = %.2f$', ycourt_val), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end 

if x_cut > 0 && x_cut <= 1
figure();
plot(yV, TLEHT_Y, 'k-', 'LineWidth', 1); hold on;
plot(yV, TFVT_Y, 'bo', 'MarkerIndices',...
    1:3:length(yV)); % FVT adjustment points, for example: 1:2, 1:3, 1:4
grid on; xlim([0 H]); box on;
ylim([0 1]);
legend({'LEHT', 'FVT'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $x_1 = %.2f$', xcourt_val), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end
end

% Auxiliary functions - LEHT

function coef = LEHT_isotropic_square(L,H,xc,yc,a,km,ki,N,Mbc,Mint,Hx,Hy)
    yy = (1:Mbc)'/(Mbc+1)*H;
    xx = (1:Mbc)'/(Mbc+1)*L;
    xL = zeros(size(yy));   yL = yy;
    xR = L + zeros(size(yy)); yR = yy;
    xB = xx;                yB = zeros(size(xx));
    xT = xx;                yT = H + zeros(size(xx));
    
    thI = (0:Mint-1)'*(2*pi/Mint);
    
    A_bcT = [ periodic_Ttilde_rows(xR,yR,xL,yL,xc,yc,a,N); ...
              periodic_Ttilde_rows(xT,yT,xB,yB,xc,yc,a,N) ];
    b_bcT = zeros(size(A_bcT,1),1);
    
    A_bcQ = [ periodic_flux_rows('x',xR,yR,xL,yL,xc,yc,a,km,N); ...
              periodic_flux_rows('y',xT,yT,xB,yB,xc,yc,a,km,N) ];
    b_bcQ = zeros(size(A_bcQ,1),1);
    
    [A_intT, b_intT] = interface_T_rows(thI,a,N);
    [A_intQ, b_intQ] = interface_qn_rows(thI,km,ki,N,Hx,Hy);
    
    A_mat = [A_bcT; A_bcQ; A_intT; A_intQ];
    b_mat = [b_bcT; b_bcQ; b_intT; b_intQ];
    
    normA = sqrt(sum(A_mat.^2, 1)); 
    normA(normA == 0) = 1; 
    A_scaled = A_mat ./ normA; 
    y = lsqminnorm(A_scaled, b_mat);
    u = y ./ normA';  
    coef.Af = u(1:N);       coef.Bf = u(N+1:2*N);
    coef.Am = u(2*N+1:3*N); coef.Bm = u(3*N+1:4*N);
    coef.Cm = u(4*N+1:5*N); coef.Dm = u(5*N+1:6*N);
    coef.N=N; coef.a=a; coef.km=km; coef.ki=ki; coef.xc=xc; coef.yc=yc;
end

function T_total = evaluate_temperature(x, y, coef, Hx, Hy)
    rx = x - coef.xc;
    ry = y - coef.yc;
    r = sqrt(rx.^2 + ry.^2);
    th = atan2(ry, rx);
    T_tilde = zeros(size(x)); 
    mask_f = (r <= coef.a);
    mask_m = ~mask_f;
    
    if any(mask_f, 'all')
        xi_f = r(mask_f) / coef.a; 
        th_f = th(mask_f);
        T_f = zeros(size(xi_f));
        for n = 1:coef.N
            T_f = T_f + coef.a * (coef.Af(n)*(xi_f.^n).*cos(n*th_f) + coef.Bf(n)*(xi_f.^n).*sin(n*th_f));
        end
        T_tilde(mask_f) = T_f;
    end
    
    if any(mask_m, 'all')
        xi_m = max(r(mask_m) / coef.a, 1e-14); 
        th_m = th(mask_m);
        T_m = zeros(size(xi_m));
        for n = 1:coef.N
            T_m = T_m + coef.a * (...
                coef.Am(n)*(xi_m.^n).*cos(n*th_m) + coef.Bm(n)*(xi_m.^n).*sin(n*th_m) + ...
                coef.Cm(n)*(xi_m.^(-n)).*cos(n*th_m) + coef.Dm(n)*(xi_m.^(-n)).*sin(n*th_m));
        end
        T_tilde(mask_m) = T_m;
    end
    T_total = Hx.*x + Hy.*y + T_tilde;
end

function [qxbar,qybar] = mean_flux_LEHT_boundary(coef,L,H,km,Hx,Hy,M)
    yy = (1:M)'/(M+1)*H;              
    xx = (1:M)'/(M+1)*L;              
    xR = L * ones(size(yy));  
    [dTdxR, ~] = eval_grad_Ttil_matrix(xR, yy, coef);
    qxbar = mean(-km * (Hx + dTdxR)); 
    yT = H * ones(size(xx));
    [~, dTdyT] = eval_grad_Ttil_matrix(xx, yT, coef);
    qybar = mean(-km * (Hy + dTdyT)); 
end

function [dTdx, dTdy] = eval_grad_Ttil_matrix(x,y,coef)  
    rx = x - coef.xc;
    ry = y - coef.yc;
    r  = sqrt(rx.^2 + ry.^2);
    th = atan2(ry,rx);
    n = 1:coef.N;
    C = cos(th*n);
    S = sin(th*n);
    Cn = C.*n;
    Sn = S.*n;
    xi = r/coef.a;
    xiP  = bsxfun(@power, xi,  n);
    xiM  = bsxfun(@power, xi, -n);
    xiP1 = bsxfun(@power, xi,  n-1);
    xiM1 = bsxfun(@power, xi, -n-1);
    dTdr  = (xiP1.*Cn)*coef.Am + (xiP1.*Sn)*coef.Bm - (xiM1.*Cn)*coef.Cm - (xiM1.*Sn)*coef.Dm;
    dTdth = coef.a * ( -(xiP.*Sn)*coef.Am + (xiP.*Cn)*coef.Bm - (xiM.*Sn)*coef.Cm + (xiM.*Cn)*coef.Dm );
    ct = cos(th);
    st = sin(th);
    invr = 1./r; 
    dTdx = dTdr.*ct - (dTdth.*st).*invr;
    dTdy = dTdr.*st + (dTdth.*ct).*invr;
end

function Arows = periodic_Ttilde_rows(x2,y2,x1,y1,xc,yc,a,N)
    [PhiAm2,PhiBm2,PhiCm2,PhiDm2] = matrix_Ttilde_basis(x2,y2,xc,yc,a,N);
    [PhiAm1,PhiBm1,PhiCm1,PhiDm1] = matrix_Ttilde_basis(x1,y1,xc,yc,a,N);
    dAm = PhiAm2 - PhiAm1;
    dBm = PhiBm2 - PhiBm1;
    dCm = PhiCm2 - PhiCm1;
    dDm = PhiDm2 - PhiDm1;
    Z = zeros(size(dAm));
    Arows = [Z Z dAm dBm dCm dDm];
end

function [PhiAm,PhiBm,PhiCm,PhiDm] = matrix_Ttilde_basis(x,y,xc,yc,a,N)
    rx = x - xc; ry = y - yc;
    r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
    xi = max(r/a, 1e-14);
    n = 1:N;
    C = cos(th*n); S = sin(th*n);
    xiP = bsxfun(@power,xi,n);
    xiM = bsxfun(@power,xi,-n);
    PhiAm = a*(xiP.*C);
    PhiBm = a*(xiP.*S);
    PhiCm = a*(xiM.*C);
    PhiDm = a*(xiM.*S);
end

function Arows = periodic_flux_rows(dir,x2,y2,x1,y1,xc,yc,a,km,N)
    [QxAm2,QxBm2,QxCm2,QxDm2, QyAm2,QyBm2,QyCm2,QyDm2] = matrix_flux_basis(x2,y2,xc,yc,a,km,N);
    [QxAm1,QxBm1,QxCm1,QxDm1, QyAm1,QyBm1,QyCm1,QyDm1] = matrix_flux_basis(x1,y1,xc,yc,a,km,N);
    
    if dir == 'x'
        dAm = QxAm2 - QxAm1; dBm = QxBm2 - QxBm1;
        dCm = QxCm2 - QxCm1; dDm = QxDm2 - QxDm1;
    else
        dAm = QyAm2 - QyAm1; dBm = QyBm2 - QyBm1;
        dCm = QyCm2 - QyCm1; dDm = QyDm2 - QyDm1;
    end
    
    Z = zeros(size(dAm));
    Arows = [Z Z dAm dBm dCm dDm];
end

function [QxAm,QxBm,QxCm,QxDm, QyAm,QyBm,QyCm,QyDm] = matrix_flux_basis(x,y,xc,yc,a,km,N)
    rx = x - xc; ry = y - yc;
    r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
    xi = max(r/a, 1e-14);
    n = 1:N;
    C = cos(th*n); S = sin(th*n);
    Cn = C.*n; Sn = S.*n;
    xiP  = bsxfun(@power,xi,n);
    xiM  = bsxfun(@power,xi,-n);
    xiP1 = bsxfun(@power,xi,n-1);
    xiM1 = bsxfun(@power,xi,-n-1);
    ct = cos(th); st = sin(th);
    invr = 1./max(r, 1e-14);
   
    dTdr_Am = (xiP1).*Cn;     dTdth_Am = -a*(xiP).*Sn;
    dTdr_Cm = -(xiM1).*Cn;    dTdth_Cm = -a*(xiM).*Sn;
    dTdr_Bm = (xiP1).*Sn;     dTdth_Bm =  a*(xiP).*Cn;
    dTdr_Dm = -(xiM1).*Sn;    dTdth_Dm =  a*(xiM).*Cn;
    
    dTdx_Am = dTdr_Am.*ct - (dTdth_Am.*st).*invr;
    dTdx_Bm = dTdr_Bm.*ct - (dTdth_Bm.*st).*invr;
    dTdx_Cm = dTdr_Cm.*ct - (dTdth_Cm.*st).*invr;
    dTdx_Dm = dTdr_Dm.*ct - (dTdth_Dm.*st).*invr;
    
    dTdy_Am = dTdr_Am.*st + (dTdth_Am.*ct).*invr;
    dTdy_Bm = dTdr_Bm.*st + (dTdth_Bm.*ct).*invr;
    dTdy_Cm = dTdr_Cm.*st + (dTdth_Cm.*ct).*invr;
    dTdy_Dm = dTdr_Dm.*st + (dTdth_Dm.*ct).*invr;
    
    QxAm = -km*dTdx_Am; QxBm = -km*dTdx_Bm; QxCm = -km*dTdx_Cm; QxDm = -km*dTdx_Dm;
    QyAm = -km*dTdy_Am; QyBm = -km*dTdy_Bm; QyCm = -km*dTdy_Cm; QyDm = -km*dTdy_Dm;
end

function [Arows, brows] = interface_T_rows(thI,a,N)
    n = 1:N;
    C = cos(thI*n);
    S = sin(thI*n);
    Arows = [ -a*C, -a*S,  a*C,  a*S,  a*C,  a*S ];
    brows = zeros(size(thI));
end

function [Arows, brows] = interface_qn_rows(thI,km,ki,N,Hx,Hy)
    n = 1:N;
    C = cos(thI*n);
    S = sin(thI*n);
    Cn = C.*n;
    Sn = S.*n;
    Arows = [ -ki*Cn, -ki*Sn,  km*Cn,  km*Sn, -km*Cn, -km*Sn ];
    dT0dr = Hx*cos(thI) + Hy*sin(thI);
    brows = -(km-ki)*dT0dr;
end