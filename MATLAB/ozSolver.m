% Code for solving Onrstein-Zernike equation using the 
% Picard interation method

clc, clear all, close all

% defining the number of points in the function
N = 800;

% generating coordinates r and k
Rmax=15.0;

dr=Rmax/N;
dk=pi/Rmax;

imax=N-1; jmax=imax;

for i=1:imax
    k(i)=i*dk;
    r(i)=i*dr;
end

% defining thermodynamic state: density and temperature

rho=0.9; % rho can be an array of values
Tr=2.0;

% defining initial guess for gamma function

gammai(:,1)=r.*0;

% calculating the interaction potential energy

Ur=(4/Tr)*((1./r').^12 - (1./r').^6);

% iterating

for w=1:size(rho,2) % loop over the densities defined previously
    
erro=1000;
t=0;
clear Ch=0 Ch2=0 gammah=0 gamma2=0 gamma3=0 C h gammaf

if w<2 % for the first density it is considered the initial guess for gamma function
    gammaf(:,1)=gammai(:,1);
else % from the second density the previous result is considered as initial guess
    gammaf(:,1)=gammaf2(:,w-1);
end

%% calculating the function C(r)
% Considering Percus-Yevick approximation
%F=exp(Ur);
%C(:,1)=(gammaf(:,1)./r'+1).*((1./F)-1);
%C(:,1)=C(:,1).*r';

% Considering the Hypernetted-Chain approximation
C(:,1)=exp(-Ur + gammaf(:,1)./r') - gammaf(:,1)./r' - 1;
C(:,1)=C(:,1).*r';

%% calculating the function h(r)
h(:,1)=C(:,1).*r'+gammaf(:,1);

tic
while erro > 1e-4 % iterating to solve OZ equation
t=t+1;

% defining mixing parameter
alfa=1;

if t>40000
    break
end

% calculating Fourier transform of C(r)
for j=1:jmax
    for i=1:imax
        Ch(j,i)=(4*pi)*dr*C(i,t)*sin(k(j)*r(i));
    end
end
Ch2=sum(Ch,2);

% calculating Fourier transform of gamma function using OZ equation
gammah(:,t)=(rho(w)*Ch2(:,1).^2)./(k'-rho(w)*Ch2(:,1));

% calculating inverse Fourier transform of gamma(k)
for i=1:imax
    for j=1:jmax
        gamma3(i,j)=(1/(2*pi^2))*dk*gammah(j,t)*sin(k(j)*r(i));
    end
end
gamma2=sum(gamma3,2);

gammaf(:,t+1)=gamma2;

% mixing the results
gammaf(:,t+1)=(1-alfa)*gammaf(:,t)+alfa*gammaf(:,t+1);

%% calculating new C(r)
% PY
%F=exp(Ur);
%C(:,t+1)=(gammaf(:,t+1)./r'+1).*((1./F)-1);   
%C(:,t+1)=C(:,t+1).*r';

% HNC
C(:,t+1)=exp(log(1./exp(Ur)) + gammaf(:,t+1)./r') - gammaf(:,t+1)./r' - 1;
C(:,t+1)=C(:,t+1).*r';

%% calculating new h(r)
h(:,t+1)=gammaf(:,t+1)+C(:,t+1);

% calculating error
erro=sqrt(sum((gammaf(:,t+1)-gammaf(:,t)).^2));

if isnan(erro)== 1 || norm(C(:,t+1)) == inf || erro == inf
    break
end

values = ['Iteration: ', num2str(t), '. Error: ', num2str(erro)];
disp(values)

end
toc

% storing final values of correlation function for all densities defined

C2(:,w)=C(:,t+1);
h2(:,w)=h(:,t+1);
gammaf2(:,w)=gammaf(:,t+1);

end

% plotting data

figure(1)
plot(r,h2./r','-k')
axis([0 6 -1.1 1.8])
xlabel(['r / \sigma'],'FontSize',18)
ylabel('h(r)','FontSize',18)
set(gca,'FontSize',18)
%print -depsc -r600 hrHNCgl

figure(2)
plot(r,C2./r','-k')
axis([0 6 -22 2])
xlabel(['r / \sigma'],'FontSize',18)
ylabel('C(r)','FontSize',18)
set(gca,'FontSize',18)
%print -depsc -r600 CrHNCgl

figure(3)
plot(r,gammaf2(:,1)./r','-k'),hold on,
axis([0 6 -2 22])
xlabel(['r / \sigma'],'FontSize',18)
ylabel('\gamma(r)','FontSize',18)
set(gca,'FontSize',18)
%print -depsc -r600 gammarHNCgl
