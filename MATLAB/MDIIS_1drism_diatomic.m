%% Code for solving 1D-RISM equations using the MDIIS method
%---------------------------------------------------------------------
% - Solves for diatomic molecules;
% - Examples included are: two models for HCl and one for Cl2;
% - It not seems to be difficult generalize it for larger molecules;
%---------------------------------------------------------------------
%%%%%%%%%%%%%%%% By Felipe Silva Carvalho

clc, clear all, close all

% defining image parameters
f = figure;
f.Position = [112 300 1300 250];

% defining the number of points in the function
N = 1500;

% generating coordinates r and k
Rmax=20.0;

dr=Rmax/N;
dk=pi/Rmax;

imax=N-1; jmax=imax;

for i=1:imax
    k(i)=i*dk;
    r(i)=i*dr;
end

% calculating the interaction potential energy and defining thermodynamic 
% state: density and temperature

%{
%%%%%%%%%%%%%%%% HCl model 1 %%%%%%%%%%%%%%%%%%
% J. Chem. Phys. 103, 3680 (1995)
rho = 0.836; % g/cm^3
MM = 36.458; % g/mol
NA = 6.022e23; % 1/mol

rho = ((rho/MM)*NA)/(1e8)^3; % 1/Angs^3
T= 297;

L=1.30; 

% 1 = H, 2 = Cl
sigma(1,1) = 2.81;
sigma(2,2) = 3.47;
sigma(1,2) = 0.5*(sigma(1,1)+sigma(2,2));
sigma(2,1) = sigma(1,2);

eps_kb(1,1) = 14.7;
eps_kb(2,2) = 190.9;
eps_kb(1,2) = sqrt(eps_kb(1,1)*eps_kb(2,2));
eps_kb(2,1) = eps_kb(1,2);

eta=0.4; 
maxdim=5;
fat=1;

% load data
data11=dlmread("dataHH.txt");
data12=dlmread("dataClH.txt");
data22=dlmread("dataClCl.txt");

title1='H-H';
title2='H-Cl';
title3='Cl-Cl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
%%%%%%%%%%%%%%%% Cl2 %%%%%%%%%%%%%%%%%%
%  J. Chem. Phys. 102, 2092 (1995)
rho = 1.660; % g/cm^3
MM = 70.906; % g/mol
NA = 6.022e23; % 1/mol

rho = ((rho/MM)*NA)/(1e8)^3; % 1/Angs^3
%rho = rho*3.332^3;
T= 200;

L=2.10;

% 1 = H, 2 = Cl
sigma(1,1) = 3.332;
sigma(2,2) = 3.332;
sigma(1,2) = 3.332;
sigma(2,1) = 3.332;

eps_kb(1,1) = 178.3;
eps_kb(2,2) = 178.3;
eps_kb(1,2) = 178.3;
eps_kb(2,1) = 178.3;

eta=0.5; 
maxdim=5;
fat=1;

% reading data to compare final results
data11=dlmread("dataCl2.txt");
data11(:,1) = data11(:,1)*3.332;
data12=dlmread("dataCl2.txt");
data12(:,1) = data12(:,1)*3.332;
data22=dlmread("dataCl2.txt");
data22(:,1) = data22(:,1)*3.332;

title1='Cl-Cl';
title2='Cl-Cl';
title3='Cl-Cl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%%%%%%%%%%%%%%%% HCl model 2 %%%%%%%%%%%%%%%%%%
%  J. Chem. Phys. 77, 509 (1982)
rho = 0.018;
T= 210; %297;

L=1.30; 

% 1 = H, 2 = Cl
sigma(1,1) = 0.400; %2.81;
sigma(2,2) = 3.353; %3.47;
sigma(1,2) = 0.5*(sigma(1,1)+sigma(2,2));
sigma(2,1) = sigma(1,2);

eps_kb(1,1) = 20; %14.7;
eps_kb(2,2) = 259; %190.9;
eps_kb(1,2) = sqrt(eps_kb(1,1)*eps_kb(2,2));
eps_kb(2,1) = eps_kb(1,2);

eta=0.4; 
maxdim=5;
fat=1;

% load data
data11=dlmread("dataHH-2.txt");
data12=dlmread("dataClH-2.txt");
data22=dlmread("dataClCl-2.txt");

title1='H-H';
title2='H-Cl';
title3='Cl-Cl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% i and j stands for the particles and third dimension for the coordinates

for i=1:2
    for j=1:2
        Ur(i,j,:) = (4*eps_kb(i,j)/T)*((sigma(i,j)./r').^12 - (sigma(i,j)./r').^6);
    end
end

% defining initial guess for c_{ij}(r) function. 

for i=1:2
    for j=1:2
        c(i,j,:)=(exp(-Ur(i,j,:))-1)*fat;
    end
end

% iterating up to n and defining m to generate the DIIS mxm matrix
m=1;
n=1;

I=eye(2,2);

% calculating intramolecular correlation
for i=1:2
    for j=1:2
        if i==j
            w(i,j,:) = k*0+1;
        else
            w(i,j,:) = sin(k*L)./(k*L);
        end
    end
end

while n <= 1000
    %% getting the Fourier transform of c_{ij}(r), C(k)
    for i=1:2
        for j=1:2
            for jj=1:jmax
                for ii=1:imax
                    Ch(i,j,jj,ii)=(4*pi/k(jj))*dr*c(i,j,ii,m)*sin(k(jj)*r(ii))*r(ii);
                end
            end
        end
    end
    Ch2=sum(Ch,4);

    %% calculating H(k) using the Fourier trnasformed 1D-RISM equations
    for pp=1:size(k,2)
        hh(:,:,pp) = w(:,:,pp)*Ch2(:,:,pp)*((I-rho*w(:,:,pp)*Ch2(:,:,pp))^-1)*w(:,:,pp);
    end

    %% calculating inverse Fourier transform of H(k), h(r)
    for i=1:2
        for j=1:2
            for ii=1:imax
                for jj=1:jmax
                    gamma3(i,j,ii,jj)=(1/(r(ii)*2*pi^2))*dk*hh(i,j,jj)*sin(k(jj)*r(ii))*k(jj);
                end
            end
        end
    end
    h(:,:,:,m)=sum(gamma3,4);

    %% Calculating h(r) using the closure relation
    for i=1:2
        for j=1:2
            h2(i,j,:,m)=exp(-Ur(i,j,:) + h(i,j,:,m) - c(i,j,:,m)) - 1;
        end
    end
    
    %% calculating residue

    for i=1:2
        for j=1:2
            R(i,j,:,m) = h2(i,j,:,m) - h(i,j,:,m);
        end
    end

    %% defining stop criterion
    mean_norm = (norm(squeeze(R(1,1,:,m)))+norm(squeeze(R(1,2,:,m)))+norm(squeeze(R(2,2,:,m))))/3;
    values = ['Iteration: ', num2str(n), '. Residual norm: ', num2str(mean_norm)];
    disp(values)
    if mean_norm < 1d-4
        break
    end
    
    %% getting S matrix
    for p=1:2
        for o=1:2
            for i=1:m
                for j=1:m
                    s(i,j,p,o) = dot(R(p,o,:,i),R(p,o,:,j),3);
                end
            end
        end
    end
    
    for p=1:2
        for o=1:2
            for i=1:m+1
                s(i,m+1,p,o) = -1; s(m+1,i,p,o) = -1;
            end
        end
    end
       
    s(m+1,m+1,:,:) = 0;

    %% calculating b vector
    for i=1:m
        b(i,1) = 0;
    end

    b(m+1,1) = -1;

    %% getting coefficients

    % solving with singular value decomposition
    
    for p=1:2
        for o=1:2
            [U,S,V] = svd(s(:,:,p,o));
            sol = V*(S^-1)*U'*b;cs(1:m,p,o) = sol(1:m,1);
        end
    end
    
    %% getting new c(r)

    for p=1:2
        for o=1:2
            r_star(:,p,o) = squeeze(R(p,o,:,1:m))*squeeze(cs(:,p,o));
        end
    end

    for p=1:2
        for o=1:2
            c(p,o,:,m+1) = squeeze(c(p,o,:,1:m))*squeeze(cs(:,p,o)) +...
                eta*squeeze(r_star(:,p,o));
        end
    end

    m=m+1;

    %% defining maximum dimenstion for matrix S
    if m > maxdim;
        m=maxdim;
        for i = 1:m
            c(:,:,:,i) = c(:,:,:,i+1);
        end

        for i=1:m-1
            R(:,:,:,i) = R(:,:,:,i+1);
            h(:,:,:,i) = h(:,:,:,i+1);
            h2(:,:,:,i) = h2(:,:,:,i+1);
        end

    end
    
    n = n+1;

    %% ploting
    h11=squeeze(h(1,1,:,m-1));
    h12=squeeze(h(1,2,:,m-1));
    h22=squeeze(h(2,2,:,m-1));

    t = tiledlayout(2,9);

    nexttile([2 3])
    plot(r,h11+1,'k-',data11(:,1),data11(:,2),'.k','MarkerSize',10)
    title(title1)
    xlim([0 10])

    nexttile([2 3])
    plot(r,h12+1,'k-',data12(:,1),data12(:,2),'.k','MarkerSize',10)
    title(title2)
    xlim([0 10])

    nexttile([2 3])
    plot(r,h22+1,'k-',data22(:,1),data22(:,2),'.k','MarkerSize',10)
    title(title3)
    xlim([0 10])
    
    lgd = legend('1D-RISM','Monte Carlo');
    lgd.Layout.Tile = 'east';

    pause(0.001)
    
end
