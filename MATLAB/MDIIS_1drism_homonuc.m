%% Code for solving 1D-RISM equations using the MDIIS method
%{
- Solves for diatomic molecules;
- All examples are homonuclear, but it can be easily modified for
heteronuclear molecules;
- It needs to change from S matrix part, so all linear equations can be 
stored in a single tensor;
- It not seems to be difficult generalize it for larger molecules;

By Felipe Silva Carvalho
%}

clc, clear all, close all

% defining the number of points in the function
N = 1000;

% generating coordinates r and k
Rmax=15.0;

dr=Rmax/N;
dk=pi/Rmax;

imax=N-1; jmax=imax;

for i=1:imax
    k(i)=i*dk;
    r(i)=i*dr;
end

% calculating the interaction potential energy and defining thermodynamic 
% state: density and temperature


% comment or uncomment this block to change system.
%%%%%%%%%%%%%%%% Cl2 %%%%%%%%%%%%%%%%%%
%  J. Chem. Phys. 102, 2092 (1995)
rho=0.5216;
Tr=1.122;
L=2.10; 
sigma1=3.332;

eta=0.5; 
maxdim=10;

% reading data to compare final results
data=dlmread("dataCl2.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
% comment or uncomment this block to change system.
%%%%%%%%%%%%%%%%% LJ %%%%%%%%%%%%%%%%%%
% J. Phys.: Condens. Matter 28 (2016) 414007
rho=0.4;
Tr=1.7;
L=1; 
sigma1=1;

eta=0.9;
maxdim=10;

% reading data to compare final results
data=dlmread("dataT1p7Dens0p4.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
% comment or uncomment this block to change system.
%%%%%%%%%%%%%%%%% LJ %%%%%%%%%%%%%%%%%%
% not getting solution for this example...
% J. Phys.: Condens. Matter 28 (2016) 414007
rho=0.01;
Tr=1.7;
L=1; 
sigma1=1;

eta=0.5;
maxdim=3;

% reading data to compare final results
data=dlmread("datat1p7Dens0p01.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

sigma = 1; % for reduced coordinates

% i and j stands for the particles and third dimension for the coordinates

for i=1:2
    for j=1:2
        Ur(i,j,:) = (4/Tr)*((sigma./r').^12 - (sigma./r').^6);
    end
end

% defining initial guess for c_{ij}(r) function. 

for i=1:2
    for j=1:2
        c(i,j,:)=(exp(-Ur(i,j,:))-1)*0;
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
            w(i,j,:) = sin(k*L/sigma1)./(k*L/sigma1);
        end
    end
end

while n <= 500
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
    for i=1:m
        for j=1:m
            s11(i,j) = dot(R(1,1,:,i),R(1,1,:,j),3);
            s12(i,j) = dot(R(1,2,:,i),R(1,2,:,j),3);
            s21(i,j) = s12(i,j);
            s22(i,j) = dot(R(2,2,:,i),R(2,2,:,j),3);
        end
    end

    for i=1:m
        s11(i,m+1) = -1;s11(m+1,i) = -1;
        s12(i,m+1) = -1;s12(m+1,i) = -1;
        s21(i,m+1) = -1;s21(m+1,i) = -1;
        s22(i,m+1) = -1;s22(m+1,i) = -1;
    end

    s11(m+1,m+1) = 0;s12(m+1,m+1) = 0;s21(m+1,m+1) = 0;s22(m+1,m+1) = 0;

    %% calculating b vector
    for i=1:m
        b(i,1) = 0;
    end

    b(m+1,1) = -1;

    %% getting coefficients


    % solving with mldivide: https://www.mathworks.com/help/matlab/ref/mldivide.html
    %{
    sol = s11\b;cs11 = sol(1:m,1);
    sol = s12\b;cs12 = sol(1:m,1);
    sol = s22\b;cs22 = sol(1:m,1);
    %}

    % solving with singular value decomposition
    
    [U,S,V] = svd(s11);
    sol = V*(S^-1)*U'*b;cs11 = sol(1:m,1);

    [U,S,V] = svd(s12);
    sol = V*(S^-1)*U'*b;cs12 = sol(1:m,1);

    [U,S,V] = svd(s22);
    sol = V*(S^-1)*U'*b;cs22 = sol(1:m,1);
    %}

    %% getting new c(r)

    r_star11 = squeeze(R(1,1,:,1:m))*cs11;
    r_star12 = squeeze(R(1,2,:,1:m))*cs12;
    r_star21 = r_star12;
    r_star22 = squeeze(R(2,2,:,1:m))*cs22;

    c(1,1,:,m+1) = squeeze(c(1,1,:,1:m))*cs11 + eta*r_star11;
    c(1,2,:,m+1) = squeeze(c(1,2,:,1:m))*cs12 + eta*r_star12;
    c(2,1,:,m+1) = c(1,2,:,m+1);
    c(2,2,:,m+1) = squeeze(c(2,2,:,1:m))*cs22 + eta*r_star22;

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

    figure(2)
    plot(data(:,1),data(:,2),'ob',r,h11+1,'k-',r,h12+1,'r--',r,h22+1,'m-.')
    axis([0 5 -0.25 2.5])
    pause(0.001)

end
