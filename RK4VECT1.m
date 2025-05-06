% author: BELOBO BELOBO DIDIER ACAS SEPT 08 2024
% RK4 method for disrete polaron Eq.(25) of Eur. Phys. J. B 90, 155 (2017) A generalized Davydov-Scott model for polarons in linear peptide chains J. Luo and B. M. A. G. Piette
% periodic BCs for bright solitons
% absorbing BCs for kinks and dark solitons

clear variables ;
close all ;
clc ;
tic 
% profile on
global N_0 R Q M_ Omega rho K J_1 chi beta Gamma gamma_ dt Mmax f_p0 omega_0
global f_p sigma delta sigdelta M epsi_ q E Fn c_0 a_0 fn N_00 phi0_p alph 
set(0,'defaultaxesfontsize',45, 'defaultaxesfontWeight','bold', 'defaultaxesLinewidth', 1);
format long g
% constants 
N_00 = 100 ; alph = 1000 ; 
hbar = 1.054e-34; %in Joule seconds, reduced Planck's constant
M_ = 1.774e-25 ; % peptide unit average mass in kg
Omega = 5.5e12 ; % in s^-1, natural angular frequency of slow phonon in alpha-helix
R = (1880*hbar/(M_*Omega))^0.5 ; % lattice spacing at equilibrium in meters
rho = 2.1 ;
K = Omega^2*M_ ;
J_1 = rho*hbar*Omega ;  % nearest-neighbour dipole interation energy
chi = 3.5e-11 ; % adiabaity parameters
beta = 0.6 ; % asymmetry parameter of the model
gamma_ = 0.005 ; %
q = 1.06e-19 ; % eletri charge of eletron
E = 1 ; %
Fn = 0.1 ; %

sigma = R*chi/(2*hbar*Omega) ;
delta = chi/(2*M_*R*Omega^2) ;
Gamma = gamma_*M_*Omega ;
epsi_ = q*E*R/(hbar*Omega) ;
ratios = sigma/delta ; 
neta1 = 2*pi/3 ;
q_ = 0.3; % (neta1 )/R  ; % 3*pi/R ; % good value  0.3 ;
% neta = 5*pi/4 ; m = 0 ; Q = (neta + 2*m*pi)/R ;
% Vg = 2*rho*R*sin(Q*R) ; M = rho*R^2*cos(Q*R) ;
% e_tcrit = - 8*sin(neta)/( 3*(1-beta^2) ) ; 
% e_t = e_tcrit - 2 ;
Q = q_ ; 
q_L = Q ; omega_L = 4*rho*sin(q_L*R/2)^2 ;
sigdelta = sigma*delta ;
M = rho*R^2*cos(Q*R) ; 
fn = Fn/(M_*R*Omega^2) ;
c_0 = 12*R^2*sigma^2*delta^2*(1-beta^2) ;
a_0 = (2*pi)^0.5*M ;
m11 = 4 ; dt = 10^(-m11) ; % in s, time step
% reduced constants for the right hand side
alpha_1 = 1i*sigma*delta*(1-beta^2)*dt ;
alpha_2 = 2*1i*sigma*delta*(1+beta^2)*dt ;
alpha_3 = 1i*rho*dt ;
alpha_4 = -2*1i*rho*dt ;
% right-hand side a1: psi_n+1, b1: psi_n-1, c1: psi_n, d1: time 
RHS_ = @(a1, b1, c1, d1) ( alpha_1*(conj(a1)*a1 + conj(b1)*b1) + ...
       alpha_2*conj(c1)*c1 )*c1 + alpha_3*(a1 + b1) + alpha_4*c1  ;
  % creation of grid points
N_ = 100 ; % number of discrete points
a  = R ; % lattice spacing
n1  = - N_ : N_ ;  % lattice grid
Nmax = length(n1) ; % number spatial grid points
t0 = 0 ; % initial time, in s
tf = 400 ; % final time, s 
omega_0 = 2*pi*5/tf ;
Mmax = round( ( tf-t0 )/dt ) + 1  ; %
t = t0:dt:tf ;
% tl = t(250)
% v=1;
% initial condition
  % x0: center of mass position, k: linear phase shift, 
  % b: nonlinear phase-shift or chirp
  % case 1 stationary bright soliton k=b=0 with homegeneous phase

N_0 = 9e-9 ;
% N_0 = - 2*M*(2*pi)^0.5*( - Vg*R*e_t )^0.5 /( R*sigma*delta*( 3*R*e_t*(1-beta^2) + 4*Vg ) ) ; % 9e-9 ; % good values 8e-9, 9e-9

% f_p0 = ( - Vg*R/e_t )^0.5 ; % (1/(4*sigdelta*N_0))*(a_0 + (c_0*N_0^2 + a_0^2).^0.5 ) ;
% f_p    = f_p0 ; % (1/4)*(M*sqrt(2)*sqrt(pi)+sqrt(-12*N_0^2*a^2*beta^2*delta^2*sigma^2+12*N_0^2*a^2*delta^2*sigma^2+2*M^2*pi))/(N_0*delta*sigma) ;
% alpha = 0 ;
% x01 = - 0.5*( 6*Vg*R*e_t*alpha*( beta^2 - 1 ) + M*R*e_t*( 1 - beta^2 ) + 4*Vg*( 3*M - 2*Vg*alpha ) ) ;
% x02 = Vg*( 3*R*e_t*( 1 - beta^2 ) + 4*Vg ) ;
% x0 = x01/x02 ;
% % phi0_p = ( sigma*delta*N_0/(2*pi)^0.5 )*( a^2*(beta^2-1)/f_p^3 + 4/f_p )...
% %        - M/( 2*f_p^2 )  ; 
% phi0_p = 0 ;
% phi_0p = phi0_p*t0 ; 
% Phi_npt0 = (N_0/(f_p*pi^0.5))^0.5.*exp( -(n1*a-x0).^2/(2*f_p^2) + 1i*phi_0p );
% % 
% % omega = 0.005 ;
% psi_np0 = (Phi_npt0.*exp(1i*(q_L*n1*a - omega_L*t0 ) ))';
% Amplp_max = (N_0/(f_p*pi^0.5))^0.5 ; 
% psint_p = psi_np0 + 0.001*Amplp_max*rand(Nmax,1) ;

% new stationary state : x0 = k = 0

% N_0 = M*sqrt(2*pi)*sqrt(3*(1-beta^2))/( 8*sigma*delta*R*(beta^2-1) ) ; % good values 8e-9, 9e-9
% f_p0 = R*sqrt( 3*( 1-beta^2 ) )/6  ;
% f_p    = f_p0 ;
% phi0_p = -( 8*sqrt(pi)*f_p0^3 )^(-1)*( N_0*sqrt(2)*sigma*delta*( 7*R^2*( 1 - beta^2 ) - 20*f_p0^2 ) + 8*M*f_p0*sqrt(pi) ) ; 
% phi_0p = phi0_p*t0 ; 
% Phi_npt0 = (N_0/(f_p*pi^0.5))^0.5.*exp( -(n1*a).^2/(2*f_p^2) + 1i*phi_0p );
% 
% omega = 0.005 ;
% psi_np0 = (Phi_npt0.*exp(1i*(q_L*n1*a - omega_L*t0 ) ))';
% Amplp_max = (N_0/(f_p*pi^0.5))^0.5 ; 
% psint_p = psi_np0 + 0.001*Amplp_max*rand(Nmax,1) ;

% MI code quest of solitary waves
Psi0 = 2 ; b01 = 0.01 ; omega01 = 0.01 ; q1 = 0.5 ;
psint_p = ( ( Psi0 + 2*b01*cos(q1*n1*R) ).*exp(1i*(q_L*R*n1 + 2*omega01*cos(q1*n1*a) ) ) )' ;
% psinp_t = zeros(Mmax, Nmax ) ;
psinp_t(1, 1:Nmax ) = psint_p ;  % psi_np0 ; 
umax( 1 ) = max( conj(psint_p).*psint_p ) ; 
% Pow_t(1, 1:Nmax) = conj(psint_p).*psint_p ;
% Phase_(1, 1:Nmax)  = angle( psint_p ) ; 
s = 1/6 ;
% construction of sol via variational approximation
% % stationary solution
% VA = VA_para_VECT('hom_p',t) ;
% VA_Q = ones(1,Mmax)'*Q ; 
% VA_Om = ones(1,Mmax)'*omega_L ;
% VA = [VA VA_Q VA_Om] ;
% Nt(1) = VA(1,1) ; % initial value of N  
% width_p(1) = VA(1,2) ; % initial value of f_p 
% width_m(1) = VA(1,2) ; % initial value of f_m numerical  
% phit(1) = VA(1,3) ;  % initial value of phi 
% Qt(1) = VA(1,4) ; % initial value of Q 
% omega_lt(1) = VA(1,5) ; % initial value of omega_l

% linear phase solution

% VA = VA_para_VECT('chirp_sta',t) ;
% VA_Q = ones(1,Mmax)'*Q ; 
% VA_Om = ones(1,Mmax)'*omega_L ;
% VA = [VA VA_Q VA_Om] ;

figure(1)
plot(n1, conj(psint_p).*psint_p, 'r', 'LineStyle', '-', 'LineWidth', 12)
xlabel('n')
% xticks([-600 -300 0 300 600])
yticks([4.0800 4.0804 4.0810])
ylabel('|\Psi_{n}(0)|^2')
% axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(a)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3A.fig' )

% figure
% plot(n1, conj(psi_np0).*psi_np0, 'r', 'LineStyle', '-.', 'LineWidth', 8)
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(0)|^2')
% % axis([-600 600 0 .00065]);
% % legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(a)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1A.fig' )

% clmn = length(VA(1,:)) ; % VA parameters
% Err = zeros(Mmax,clmn) ; % initial value of error matrix between anal ans numerics
 
% K1 = zeros(1,Nmax) ; K2 = zeros(1,Nmax) ; 
% K3 = zeros(1,Nmax) ; K4 = zeros(1,Nmax) ; 

% RK4 LOOP WITH VECTORIZE RHS
e22 = ones(Nmax,1); 
AA = spdiags([e22 -2*e22 e22], [-1 0 1], Nmax, Nmax) ; % vect for centered finite differences
AA(1,Nmax) = 1 ;  AA(Nmax,1) = 1 ;  % periodic BCs
BB = spdiags([e22*0 -2*e22*0 e22], [-1 0 1], Nmax, Nmax) ; %  upper digaonal 
BB(Nmax,1) = 1 ;  % periodic BCs
QQ = spdiags([e22 -2*e22*0 e22*0], [-1 0 1], Nmax, Nmax) ; %  lower diagonal
QQ(1,Nmax) = 1 ;    % periodic BCs

m1 = 1 ; t1(m1) = 0 ; m0 = (Mmax-1)/1000 ; 
Egy = zeros(1,Mmax) ; 
s00 = 0 ; s001 = 0 ; 
for n = 1:Nmax-1
        s00 = s00 -rho*( conj(psint_p(n+1))*psint_p(n) + ...
              conj(psint_p(n))*psint_p(n+1) ) ;
        s001 = s001 + 0.5*chi^2*(conj(psint_p(n)).*psint_p(n))^2 ; 
end
    
    s001 = s001 + 0.5*chi^2*(conj(psint_p(Nmax)).*psint_p(Nmax))^2 ; 
    
    Egy(1) = s00 + s001 ;
for m = 2: Mmax
    % 1st RK function
%     P11 = conj(psint_p).*psint_p ;
    K1 = alpha_1*( BB*( conj(psint_p).*psint_p ) + QQ*( conj(psint_p).*psint_p ) ... % alpha_1*( QQ*P11 )*psint_p ...
         + alpha_2*conj(psint_p).*psint_p ).*psint_p + alpha_3*AA*psint_p ;
    % actualize the vaiable \psi_n
    psint_p1 = psint_p + 0.5*K1 ;
    % 2nd RK function
    K2 = alpha_1*( BB*( conj(psint_p1).*psint_p1 ) + QQ*( conj(psint_p1).*psint_p1 ) ... % alpha_1*( QQ*P11 )*psint_p ...
         + alpha_2*conj(psint_p1).*psint_p1 ).*psint_p1 + alpha_3*AA*psint_p1 ;
    % actualize the vaiable \psi_n
    psint_p2 = psint_p + 0.5*K2 ;
    % 3rd RK function
    K3 = alpha_1*( BB*( conj(psint_p2).*psint_p2 ) + QQ*( conj(psint_p2).*psint_p2 ) ... % alpha_1*( QQ*P11 )*psint_p ...
         + alpha_2*conj(psint_p2).*psint_p2 ).*psint_p2 + alpha_3*AA*psint_p2 ;
    % actualize the vaiable \psi_n
    psint_p3 = psint_p + K3 ;
    % 4th RK function
    K4 = alpha_1*( BB*( conj(psint_p3).*psint_p3 ) + QQ*( conj(psint_p3).*psint_p3 ) ... % alpha_1*( QQ*P11 )*psint_p ...
         + alpha_2*conj(psint_p3).*psint_p3 ).*psint_p3 + alpha_3*AA*psint_p3 ;
    % adance solution in time
    psint_p = psint_p + s*(K1 +2*K2 + 2*K3 + K4 ) ; 
    % energy computation
%     E1a = 0 ; E1b = 0 ;
%     for n = 1:Nmax-1
%         E1a = E1a -rho*( conj(psint_p(n+1))*psint_p(n) + ...
%               conj(psint_p(n))*psint_p(n+1) ) ;
%         E1b = E1b + 0.5*chi^2*(conj(psint_p(n)).*psint_p(n))^2 ; 
%     end
%     
%     E1b = E1b + 0.5*chi^2*(conj(psint_p(Nmax)).*psint_p(Nmax))^2 ; 
%     Egy(m) = E1a + E1b ;
        
%     psinp_t(m, 1:Nmax ) = psint_p ;
    % saving sol in time-space matrix
    if mod(m,m0) == 0
        m1 = m1 + 1 ;
        t1(m1) = m*dt ;
        psinp_t(m1, 1:Nmax ) = psint_p ;  c45 = max(abs(psinp_t).^2) ;
        umax( m1 ) = max( conj(psint_p).*psint_p ) ;
    end
    
    % LEAST SQUARE FITTING 
% X0 = VA(1, :) ; % abritrary initial guess
%       % case 1 stationary bright soliton x0=0, k=b=0 with homegeneous
% func_model = @(X,u1) ( X(1)/( X(2)*pi^0.5 ) )^0.5*exp( -u1.^2*a^2/(2*X(2)^2) + ...
%      1i*( X(3) + X(4)*u1*a - X(5)*t(m) )  ) ;
% % re22 = size(func_model(X0,n1));
% % re33 = size(psint_p);
% opts = optimset('Algorithm','levenberg-marquardt','Display','off');% 'Algorithm','levenberg-marquardt' since complex values
% % [VAestimated] = lsqnonlin(Fittinfunc1,VA0,[],[],opts);    
% [Xestimated] = lsqcurvefit(func_model,X0, n1, psint_p', [],[],opts); 
% Err(m,:) = abs(VA(m,:)-Xestimated) ;
    
end

% uu = psinp_t(1:250:Mmax, : ) ; t1 = t(1:250:Mmax) ;
% figure(99)
% plot(t, Egy, 'b', 'LineStyle', '-', 'LineWidth', 8 )                                                                                                                                         
% xlabel('t')
% ylabel('energy')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(a)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG2.fig' )
figure(2)
plot(t1, umax, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('t')
% xticks([-600 -300 0 300 600])
ylabel('Max|\Psi_{n}(t)|^2')
% axis([-600 600 0 .00065]);
%legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(a)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG9B.fig')

figure(11)
mesh(n1, t1, abs(psinp_t).^2)
% mesh(n1, t1, abs(uu).^2)
grid off
xlabel('n')
% xticks([-4 -2 0 2 4])
% xticks([-40 -20 0 20 40]) 
ylabel('t')
% yticks([0 20 40 60 80 100])
% yticks([0 7 14 18])
zlabel('|\Psi_{n}(x,t)|^2')
% axis([-4 4 0 Z(end) 0 6.5]);
% axis([-40 40 0 t1(end) 0 1]);
% title('|\Psi_{1}(n,t)|^2')
rotate_labels(gca);
call_when_rotates = @(~,~,hax)(rotate_labels(hax));
hrot = rotate3d;
set(hrot,'ActionPreCallBack',...
@(hfig,event_data)(set(hfig,'WindowButtonMotionFcn',{call_when_rotates event_data.Axes})));
set(hrot,'ActionPostCallBack',@(hfig,~)(set(hfig,'WindowButtonMotionFcn','')));
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(f)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3F.fig' )

% toc

n_00= m0 ; n_2 = 2*n_00 ; n_3 = 3*n_00 ; n_4 = 4*n_00 ;
al_25 = 2.5*10^(m11) ; m_al = al_25/m0 ; t_25 = 1 + m_al ;
% figure(1000)
% % plot(n1, abs(psint_anal_n_00).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psinp_t(t_25,:)).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(2.5)|^2')
% % axis([-600 600 0 .00065]);
% %legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(b)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1B.fig' )

% case 1: stationary analytical sol at a specified time

% psint_anal_n_00 = ( VA(n_00+1,1)/(VA(n_00+1,2)*pi^0.5) )^0.5*...
%                  exp( -(n1*R).^2/( 2*VA(n_00+1,2)^2 ) + ... 
%                  1i*( VA(n_00+1,3) + VA(n_00+1,4)*n1*a - ...
%                  VA(n_00+1,5)*t(n_00+1) ) ) ;
psint_p_n_00 = psinp_t(450, : ) ;

figure(100)
% plot(n1, abs(psint_anal_n_00).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
plot(n1, abs(psint_p_n_00).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
% xticks([-600 -300 0 300 600])
ylabel('|\Psi_{n}(45)|^2')
% axis([-600 600 0 .00065]);
%legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(b)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3B.fig' )

% psint_anal_200 = ( VA(n_2+1,1)/(VA(n_2+1,2)*pi^0.5) )^0.5*...
%                  exp( -(n1*R).^2/( 2*VA(n_2+1,2)^2 ) + ... 
%                  1i*( VA(n_2+1,3) + VA(n_2+1,4)*n1*a - ...
%                  VA(n_2+1,5)*t(n_2+1) ) ) ;
psint_p_200 = psinp_t(600, : ) ;

figure(200)
% plot(n1, abs(psint_anal_200).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
plot(n1, abs(psint_p_200).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
% xticks([-600 -300 0 300 600])
ylabel('|\Psi_{n}(60)|^2')
% axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(c)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3C.fig' )

% psint_anal_300 = ( VA(n_3+1,1)/(VA(n_3+1,2)*pi^0.5) )^0.5*...
%                  exp( -(n1*R).^2/( 2*VA(n_3+1,2)^2 ) + ... 
%                  1i*( VA(n_3+1,3) + VA(n_3+1,4)*n1*a - ...
%                  VA(n_3+1,5)*t(n_3+1) ) ) ;
psint_p_300 = psinp_t(800, : ) ;

figure(300)
% plot(n1, abs(psint_anal_300).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
plot(n1, abs(psint_p_300).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
% xticks([-600 -300 0 300 600])
ylabel('|\Psi_{n}(80)|^2')
% axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(d)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3D.fig' )

% psint_anal_400 = ( VA(n_4+1,1)/(VA(n_4+1,2)*pi^0.5) )^0.5*...
%                  exp( -(n1*R).^2/( 2*VA(n_4+1,2)^2 ) + ... 
%                  1i*( VA(n_4+1,3) + VA(n_4+1,4)*n1*a - ...
%                  VA(n_4+1,5)*t(n_4+1) ) ) ;
psint_p_400 = psinp_t(900, : ) ;

figure(400)
% plot(n1, abs(psint_anal_400).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
plot(n1, abs(psint_p_400).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
% xticks([-600 -300 0 300 600])
ylabel('|\Psi_{n}(90)|^2')
% axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
dim = [0.2 0.5 0.3 0.3] ;
Str15 = '(e)' ;
annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
saveas(gcf, 'FIG3E.fig' )

% Err_N = Err(:,1)' ; % abs error anal vs num N parameter
% Err_fp = Err(:,2)' ; % abs error anal vs num width parameter
% Err_phi = Err(:,3)' ; % abs error anal vs num phi parameter
% Err_Q = Err(:,4)' ; % abs error anal vs num Q parameter
% Err_omega_L = Err(:,5)' ; % abs error anal vs num omega parameter
% 
% figure(101)
% plot(t,Err_N, 'b', 'LineStyle', '-', 'LineWidth', 4) ;
% hold on
% plot(t,Err_fp, 'r', 'LineStyle', ':', 'LineWidth', 4 )
% % hold on
% % plot(t,Err_phi, 'g', 'LineStyle', '-.', 'LineWidth', 4 )
% hold on
% plot(t,Err_Q, 'c', 'LineStyle', '--', 'LineWidth', 4 )
% hold on
% plot(t,Err_omega_L, 'm+', 'LineStyle', ':', 'LineWidth', 4 )
% xlabel('t')
% % xticks([-600 -300 0 300 600])
% ylabel('errors')
% % axis([-600 600 0 .00065]);
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(e)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1E.fig' )
%legend({'ErrN(t)', 'Errfp(t)', 'ErrQ(t)', 'Err\omega(t)'}, 'Box', 'off')

% CASE 2 LINEAR PHASE
% psint_anal_n_00 = ( VA(n_00+1,1)/( VA(n_00+1,2)*pi^0.5 ) )^0.5*exp( -( n1*a-VA(n_00+1,3) ).^2/(2*VA(n_00+1,2)^2) + ...
%      1i*( VA(n_00+1,4)*( n1*a-VA(n_00+1,3) ) + VA(n_00+1,5) + n1*a*VA(n_00+1,6) - VA(n_00+1,7)*t(n_00+1) ) )  ;
% 
% psint_p_n_00 = psinp_t(2, : ) ;
% 
% figure(100)
% plot(n1, abs(psint_anal_n_00).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psint_p_n_00).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(100\tau)|^2')
% % axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(a)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1A.fig' )
%  
%  psint_anal_200 = ( VA(n_2+1,1)/( VA(n_2+1,2)*pi^0.5 ) )^0.5*exp( -( n1*a-VA(n_2+1,3) ).^2/(2*VA(n_2+1,2)^2) + ...
%      1i*( VA(n_2+1,4)*( n1*a-VA(n_2+1,3) ) + VA(n_2+1,5) + n1*a*VA(n_2+1,6) - VA(n_2+1,7)*t(n_2+1) ) )  ;
% 
% psint_p_200 = psinp_t(3, : ) ;
% 
% figure(200)
% plot(n1, abs(psint_anal_200).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psint_p_200).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(200\tau)|^2')
% % axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(b)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1B.fig' ) 
% 
% psint_anal_300 = ( VA(n_3+1,1)/( VA(n_3+1,2)*pi^0.5 ) )^0.5*exp( -( n1*a-VA(n_3+1,3) ).^2/(2*VA(n_3+1,2)^2) + ...
%      1i*( VA(n_3+1,4)*( n1*a-VA(n_3+1,3) ) + VA(n_3+1,5) + n1*a*VA(n_3+1,6) - VA(n_3+1,7)*t(n_3+1) ) )  ;
% 
% psint_p_300 = psinp_t(4, : ) ;
% 
% figure(300)
% plot(n1, abs(psint_anal_300).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psint_p_300).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(300\tau)|^2')
% % axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(c)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1C.fig' ) 
% 
% psint_anal_400 = ( VA(n_4+1,1)/( VA(n_4+1,2)*pi^0.5 ) )^0.5*exp( -( n1*a-VA(n_4+1,3) ).^2/(2*VA(n_4+1,2)^2) + ...
%      1i*( VA(n_4+1,4)*( n1*a-VA(n_4+1,3) ) + VA(n_4+1,5) + n1*a*VA(n_4+1,6) - VA(n_4+1,7)*t(n_4+1) ) )  ;
% 
% psint_p_400 = psinp_t(5, : ) ;
% 
% figure(400)
% plot(n1, abs(psint_anal_400).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psint_p_400).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;
% xlabel('n')
% % xticks([-600 -300 0 300 600])
% ylabel('|\Psi_{n}(400\tau)|^2')
% % axis([-600 600 0 .00065]);
% legend({'ANAL', 'NUM'}, 'Box', 'off')
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(d)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1D.fig' ) 

% Err_N = Err(:,1)' ; % abs error anal vs num N parameter
% Err_fp = Err(:,2)' ; % abs error anal vs num width parameter
% Err_x0 = Err(:,3)' ; % abs error anal vs num x0 parameter
% Err_k = Err(:,4)' ; % abs error anal vs num k parameter
% Err_phi = Err(:,5)' ; % abs error anal vs num phi parameter
% Err_Q = Err(:,6)' ; % abs error anal vs num Q parameter
% Err_omega_L = Err(:,7)' ; % abs error anal vs num omega parameter

% figure(101)
% % plot(t,Err(:,1)', 'b', 'LineStyle', '-', 'LineWidth', 4) ;
% plot(t,Err_N, 'b', 'LineStyle', '-', 'LineWidth', 4) ;
% hold on
% % plot(t,Err(:,2)', 'r', 'LineStyle', ':', 'LineWidth', 4 )
% plot(t,Err_fp, 'r', 'LineStyle', ':', 'LineWidth', 4 )
% hold on
% % plot(t,Err_phi, 'g', 'LineStyle', '--', 'LineWidth', 4 )
% % plot(t,Err(:,4)', 'c', 'LineStyle', '-.', 'LineWidth', 4 )
% plot(t,Err_x0, 'g', 'LineStyle', '--', 'LineWidth', 4 )
% hold on
% % plot(t,Err_k, 'y', 'LineStyle', '-.', 'LineWidth', 4 )
% % hold on
% plot(t,Err_Q, 'yo', 'LineStyle', ':', 'LineWidth', 4 )
% hold on
% % plot(t,Err(:,5)', 'm', 'LineStyle', '-', 'LineWidth', 4 )
% plot(t,Err_omega_L, 'm+', 'LineStyle', '--', 'LineWidth', 4 )
% xlabel('t')
% % xticks([-600 -300 0 300 600])
% ylabel('errors')
% % axis([-600 600 0 .00065]);
% dim = [0.2 0.5 0.3 0.3] ;
% Str15 = '(e)' ;
% annotation('textbox', dim, 'String', Str15, 'FitBoxToText', 'on') ;
% saveas(gcf, 'FIG1E.fig' )
% % legend({'Err_N(t)', 'Err_f_p(t)', 'Err_x0(t)', 'Err_k(t)', 'Err_phi(t)', 'Err_Q(t)', 'Err_\omega(t)'}, 'Box', 'off')
% legend({'ErrN(t)', 'Errfp(t)', 'Errx0(t)', 'ErrQ(t)', 'Err\omega(t)'}, 'Box', 'off')
 
% profile viewer


% psint_anal_1000 = ( VA(1000,1)/(VA(1000,2)*pi^0.5) )^0.5*...
%                  exp( -(n1*R).^2/( 2*VA(1000,2)^2 ) + ... 
%                  1i*( VA(1000,3) + VA(1000,4)*n1*a - ...
%                  VA(1000,5)*t(1000) ) ) ;
% psint_p_1000 = psinp_t(1000, : ) ;
% figure(1000)
% plot(n1, abs(psint_anal_1000).^2, 'b', 'LineStyle', '-', 'LineWidth', 8) ; 
% hold on
% plot(n1, abs(psint_p_1000).^2, 'r', 'LineStyle', ':', 'LineWidth', 6) ;




toc   
  
    
    
 