% author: BELOBO BELOBO DIDIER ACAS SEPT 08 2024
% RK4 method for disrete polaron Eq.(25) of Eur. Phys. J. B 90, 155 (2017) A generalized Davydov-Scott model for polarons in linear peptide chains J. Luo and B. M. A. G. Piette
% periodic BCs for bright solitons


clear variables ;
close all ;
clc ;
tic 

set(0,'defaultaxesfontsize',45, 'defaultaxesfontWeight','bold', 'defaultaxesLinewidth', 1);
format long g
% constants 

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
q_ = 0.3; 
Q = q_ ; 
q_L = Q ; omega_L = 4*rho*sin(q_L*R/2)^2 ;
sigdelta = sigma*delta ;
M = rho*R^2*cos(Q*R) ; 
fn = Fn/(M_*R*Omega^2) ;
c_0 = 12*R^2*sigma^2*delta^2*(1-beta^2) ;
a_0 = (2*pi)^0.5*M ;
m11 = 4 ; dt = 10^(-m11) ; %  time step

% reduced constants for the right hand side

alpha_1 = 1i*sigma*delta*(1-beta^2)*dt ;
alpha_2 = 2*1i*sigma*delta*(1+beta^2)*dt ;
alpha_3 = 1i*rho*dt ;
alpha_4 = -2*1i*rho*dt ;

% creation of grid points

N_ = 40 ; % number of discrete points
a  = R ; % lattice spacing
n1  = - N_ : N_ ;  % lattice grid
Nmax = length(n1) ; % number spatial grid points
t0 = 0 ; % initial time, in s
tf = 100 ; % final time, s 
omega_0 = 2*pi*5/tf ;
Mmax = round( ( tf-t0 )/dt ) + 1  ; %
t = t0:dt:tf ;

% Modulation instability: quest of solitary waves

Psi0 = 2 ; b01 = 0.01 ; omega01 = 0.01 ; q1 = 0.5 ;
psint_p = ( ( Psi0 + 2*b01*cos(q1*n1*R) ).*exp(1i*(q_L*R*n1 + 2*omega01*cos(q1*n1*a) ) ) )' ;
psinp_t(1, 1:Nmax ) = psint_p ;  
umax( 1 ) = max( conj(psint_p).*psint_p ) ; 
s = 1/6 ;


figure(1) % densty at initial time
plot(n1, conj(psint_p).*psint_p, 'r', 'LineStyle', '-', 'LineWidth', 12)
xlabel('n')
yticks([4.0800 4.0804 4.0810])
ylabel('|\Psi_{n}(0)|^2')

% RK4 LOOP WITH VECTORIZE RHS
e22 = ones(Nmax,1); 
AA = spdiags([e22 -2*e22 e22], [-1 0 1], Nmax, Nmax) ; % vect for centered finite differences
AA(1,Nmax) = 1 ;  AA(Nmax,1) = 1 ;  % periodic BCs
BB = spdiags([e22*0 -2*e22*0 e22], [-1 0 1], Nmax, Nmax) ; %  upper digaonal 
BB(Nmax,1) = 1 ;  % periodic BCs
QQ = spdiags([e22 -2*e22*0 e22*0], [-1 0 1], Nmax, Nmax) ; %  lower diagonal
QQ(1,Nmax) = 1 ;    % periodic BCs

m1 = 1 ; t1(m1) = 0 ; m0 = (Mmax-1)/1000 ; 

for m = 2: Mmax
    % 1st RK function
    K1 = alpha_1*( BB*( conj(psint_p).*psint_p ) + QQ*( conj(psint_p).*psint_p ) ... 
         + alpha_2*conj(psint_p).*psint_p ).*psint_p + alpha_3*AA*psint_p ;
     
    % actualize the vaiable \psi_n
    
    psint_p1 = psint_p + 0.5*K1 ;
    
    % 2nd RK function
    
    K2 = alpha_1*( BB*( conj(psint_p1).*psint_p1 ) + QQ*( conj(psint_p1).*psint_p1 ) ... 
         + alpha_2*conj(psint_p1).*psint_p1 ).*psint_p1 + alpha_3*AA*psint_p1 ;
     
    % actualize the vaiable \psi_n
    
    psint_p2 = psint_p + 0.5*K2 ;
    
    % 3rd RK function
    
    K3 = alpha_1*( BB*( conj(psint_p2).*psint_p2 ) + QQ*( conj(psint_p2).*psint_p2 ) ... 
         + alpha_2*conj(psint_p2).*psint_p2 ).*psint_p2 + alpha_3*AA*psint_p2 ;
     
    % actualize the vaiable \psi_n
    
    psint_p3 = psint_p + K3 ;
    
    % 4th RK function
    
    K4 = alpha_1*( BB*( conj(psint_p3).*psint_p3 ) + QQ*( conj(psint_p3).*psint_p3 ) ... 
         + alpha_2*conj(psint_p3).*psint_p3 ).*psint_p3 + alpha_3*AA*psint_p3 ;
     
    % advance solution in time
    
    psint_p = psint_p + s*(K1 +2*K2 + 2*K3 + K4 ) ; 
    
    
    % saving solution in time-space matrix
    if mod(m,m0) == 0
        m1 = m1 + 1 ;
        t1(m1) = m*dt ;
        psinp_t(m1, 1:Nmax ) = psint_p ;  
        umax( m1 ) = max( conj(psint_p).*psint_p ) ;
    end
    
end


figure(2)  % maximum of square density vs time 
plot(t1, umax, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('t')
ylabel('Max|\Psi_{n}(t)|^2')

figure(11)  % spatiotemporal evolution square density 
mesh(n1, t1, abs(psinp_t).^2)
grid off
xlabel('n')
ylabel('t')


psint_p_n_00 = psinp_t(450, : ) ;

figure(100)

plot(n1, abs(psint_p_n_00).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
ylabel('|\Psi_{n}(45)|^2')

psint_p_200 = psinp_t(600, : ) ;

figure(200)
plot(n1, abs(psint_p_200).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
ylabel('|\Psi_{n}(60)|^2')

psint_p_300 = psinp_t(800, : ) ;

figure(300)
plot(n1, abs(psint_p_300).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
ylabel('|\Psi_{n}(80)|^2')

psint_p_400 = psinp_t(900, : ) ;

figure(400)
plot(n1, abs(psint_p_400).^2, 'r', 'LineStyle', '-', 'LineWidth', 12) ;
xlabel('n')
ylabel('|\Psi_{n}(90)|^2')

toc   
  
    
    
 
