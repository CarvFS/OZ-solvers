%% Code for solving Onrstein-Zernike equation using the MDIIS method

clc, clear all, close all

% reading data to compare final results
data=dlmread("hrHNC.txt");

%% falta verificar o resultado para h(r) e deixar o código mais limpo...

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

% defining thermodynamic state: density and temperature
rho=0.8;
Tr=2.0;

% defining initial guess for c(r) function

c(:,1)=zeros(length(r),1);

% calculating the interaction potential energy

eps = 1;
sigma = 1;
Ur=(4*eps/Tr)*((sigma./r').^12 - (sigma./r').^6);

% iterating up to n and defining m to generate the DIIS mxm matrix
m=1;
n=1;

while n <= 200
    %% Calculating h(r) using OZ equation in the reciprocal space
    % getting the Fourier transform of c(r), C(k)
    for j=1:jmax
        for i=1:imax
            Ch(j,i)=(4*pi/k(j))*dr*c(i,m)*sin(k(j)*r(i))*r(i);
        end
    end
    Ch2=sum(Ch,2);
    
    % acquiring H(k)
    hh = Ch2./(1 - rho.*Ch2);
    %plot(k,hh,'-')
    %pause()

    % calculating inverse Fourier transform of H(k), h(r)
    for i=1:imax
        for j=1:jmax
            gamma3(i,j)=(1/(r(i)*2*pi^2))*dk*hh(j)*sin(k(j)*r(i))*k(j);
        end
    end
    h(:,m)=sum(gamma3,2);

    %% Calculating h(r) using the closure relation
    h2(:,m)=exp(-Ur + h(:,m) - c(:,m)) - 1;

    %% Calculating residue
    R(:,m) = h2(:,m) - h(:,m);
    
    %% setting stop criterium
    if norm(R(:,m)) < 1e-4
        break
    end

    %% getting the S matrix

    for i=1:m
        for j=1:m
            s(i,j) = dot(R(:,i),R(:,j));
        end
    end

    for i=1:m
        s(i,m+1) = -1;
        s(m+1,i) = -1;
    end

    s(m+1,m+1) = 0;

    for i=1:m
        b(i,1) = 0;
    end

    b(m+1,1) = -1;

    %% getting coefficients

    sol = s\b;
    cs = sol(1:m,1);
    
    %% getting new c(r)

    r_star = R(:,1:m)*cs;

    c(:,m+1) = c(:,1:m)*cs + 0.05*r_star;
    m=m+1;

    %% defining maximum dimenstion for matrix S
    if m > 10;
        m=10;
        for i = 1:m
            c(:,i) = c(:,i+1);
        end

        for i=1:m-1
            R(:,i) = R(:,i+1);
            h(:,i) = h(:,i+1);
            h2(:,i) = h2(:,i+1);
        end

    end
    
    %% printing numerical values for the last residue norm
    values = ['Iteration: ', num2str(n), '. Residual norm: ', num2str(norm(R(:,m-1)))];
    disp(values)

    n=n+1;

    %% printing
    figure(2)
    plot(data(:,1),data(:,2),'k-',r,h(:,m-1),'r--')
    axis([0 10 -1.5 1.5])
    pause(0.001)

end