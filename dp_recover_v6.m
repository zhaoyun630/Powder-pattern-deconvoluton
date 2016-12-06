% This piece of code is used to recover diffraction pattern from correlation functions.
% last modified by Yun Zhao. 2013/4/26
% sample m over 4. In the research report, I used m=2. I choose 4 here as I found that only 4*m are non-zero values practically.


clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% loading data  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_crystal=input('number of crystals per shot is ');
N_dp=input('number of diffraction pattern is ');

b2=[int2str(N_dp) 'N_dp_' int2str(N_crystal) '_parameter.mat'];
b3=[int2str(N_dp) 'dp_SAXS_R.mat'];
b4=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_ac.mat'];
b5=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_tc.mat'];
b6=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_pc.mat'];
parameter=importdata(b2);
r_th=parameter(3);  % The reference ring
SAXS_R=importdata(b3);
C1=importdata(b4);   % load auto correlation function
C2_r=importdata(b6); % load pair correlation function
C3_r=importdata(b5); % load triple correlation function

N_rot=360;
Detector_size=401;
N_r=(Detector_size-1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Fourier transform of pair and triple correlation  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_max=19;        % set the number of non-zeros terms in Fourier transform
m_max=M_max*4+1; % set cut off value for M.

% fourier transform of angular pair correlation function
B_M=zeros(N_r,m_max);
for i=1:N_r
    for j=1:m_max
        for k=1:N_rot
            B_M(i,j)=B_M(i,j)+C2_r(i,k)*cos(j*k*pi/180);  % note that m starts from 1 rather than 0.
        end
    end 
end
B_M=B_M/360;

% set cut off value for M as 38.
B_M_ref_obs=zeros(1,M_max);
for i=1:M_max
    B_M_ref_obs(i)=B_M(r_th,4*i);
end
I_M_ref_mag=sqrt(abs(B_M_ref_obs));


% fourier transform of angular triple correlation function for one ring
T_M=zeros(N_r,m_max);
for i=1:N_r
    for j=1:m_max
        for k=1:N_rot
            T_M(i,j)=T_M(i,j)+C3_r(i,k)*cos(j*k*pi/180);
        end
    end
end
T_M=T_M/360;

% set cut off value for M as 38.
T_M_ref_obs=zeros(1,M_max);
for i=1:M_max
    T_M_ref_obs(i)=T_M(r_th,4*i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Searching over all possible combination of signs of I_M %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get all 2^19 possible phases for I_M, which are different +/- sign combinations for M_max I_M values
Phi_q=[1;-1];
for i=1:M_max-1
    b=size(Phi_q);
    a1=ones(b(1),1);
    a2=[[Phi_q a1];[Phi_q -a1]];
    Phi_q=a2;
end

mat_size=size(Phi_q);
T_M_calc_1=zeros(mat_size);
T_M_calc_2=zeros(mat_size);

% Calculate the T_M with all possible phase combinations
for i=1:mat_size(1)
    % for those m >0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
    for j=1:M_max
        % for those M >0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
        for k=1:M_max
            if j~=k
               T_M_calc_1(i,j)=T_M_calc_1(i,j)+I_M_ref_mag(k)*I_M_ref_mag(abs(k-j))*Phi_q(i,k)*Phi_q(i,abs(k-j));
            end
        end
         % for those M <0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
        for k=1:M_max
            if j+k < M_max+1
               T_M_calc_1(i,j)=T_M_calc_1(i,j)+I_M_ref_mag(k)*I_M_ref_mag(k+j)*Phi_q(i,k)*Phi_q(i,k+j);
            end
        end
        T_M_calc_1(i,j)=T_M_calc_1(i,j)*(I_M_ref_mag(j)*Phi_q(i,j));
    end
    
    % for those m <0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
    for j=1:M_max
        % for those M >0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
        for k=1:M_max
            if j+k < M_max+1
               T_M_calc_2(i,j)=T_M_calc_2(i,j)+I_M_ref_mag(k)*I_M_ref_mag(k+j)*Phi_q(i,k)*Phi_q(i,k+j);
            end
        end
         % for those M <0. refer to eqn(10) in (Saldin, 2010, New Journal of Physics
        for k=1:M_max
            if j~=k
               T_M_calc_2(i,j)=T_M_calc_2(i,j)+I_M_ref_mag(k)*I_M_ref_mag(abs(k-j))*Phi_q(i,k)*Phi_q(i,abs(k-j));
            end
        end
        T_M_calc_2(i,j)=T_M_calc_2(i,j)*(I_M_ref_mag(j)*Phi_q(i,j));
    end
end

T_error=zeros(mat_size(1),1);
% Searching the combination of signs which could minimize the difference between T_cacl and T_obs.

% calculate the error between T_M_calc and T_M_obs
for i=1:mat_size(1)
    for j=1:M_max
        T_error(i)=T_error(i)+abs(T_M_calc_1(i,j)-T_M_ref_obs(j))+abs(T_M_calc_2(i,j)-T_M_ref_obs(j));
    end
end

% searching smallest error and save the corresponding sign as Phi.
for i=2:mat_size(1)
    if T_error(i)<T_error(i-1)
        Phi=Phi_q(i,:);
        n_th=i;
    end
end

T_obs=0;
for j=1:M_max
        T_obs=T_obs+2*abs(T_M_ref_obs(j));
end

% Reconstruct the intensity circular expansion coefficients for the reference ring
I_ref=zeros(N_rot,1);
I_M=Phi.*I_M_ref_mag;
for i=1:N_rot
    for j=1:M_max
        I_ref(i)=I_ref(i)+2*I_M(j)*cos(4*j*i*pi/180);
    end
end

% Reconstruct the intensity circular expansion coefficients for other rings with respect to the reference ring
I_M_r=zeros(N_r,M_max);
I_r=zeros(N_r,N_rot);
for r=1:N_r
    for j=1:M_max
        I_M_r(r,j)=B_M(r,4*j)/I_M(j);
    end
end

% Reconstruct the diffraction pattern in polar coordinates
for r=1:N_r
    for i=1:N_rot
        for j=1:M_max
            I_r(r,i)=I_r(r,i)+2*I_M_r(r,j)*cos(4*j*i*pi/180);
        end
           % I_r(r,i)=I_r(r,i)+SAXS_R(r);
    end
end


% transform from polar coordinates to cartesian coordinates
I_c=zeros(Detector_size);
for l=1:N_r
    for m=1:N_rot
         x1=round(l*cos(m*pi/180)+201);
         y1=round(l*sin(m*pi/180)+201); 
         if x1>0 & x1<402 & y1>0 & y1<402
            I_c(x1,y1)=I_r(l,m);
         end
    end
end

imagesc(abs(I_c))
axis image
title('reconstructed single crystal diffraction pattern')

