% Used to simulate diffraction pattern, calculate pair correlation and triple correlation
% Xstals are considered as finite. Flat Ewald surface is assumed. 2D demonstration
% Last edited in 2013/5/10.  By Yun Zhao.

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Set input parameters and global variables %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_crystal=input('set number of crystals per shot as ');
N_dp=input('set number of diffraction pattern as ');

% coordinates of all atoms in one unit cell
x=3/8;
C_Mg=[0 0 0;1/4 1/4 1/4];
C_Al=[5/8 5/8 5/8;5/8 7/8 7/8;7/8 5/8 7/8;7/8 7/8 5/8];
C_O=[x x x;x -x -x;-x x -x;-x -x x;1/4-x 1/4-x 1/4-x;1/4-x 1/4+x 1/4+x;1/4+x 1/4-x 1/4+x;1/4+x 1/4+x 1/4-x];
C=[C_Mg;C_Al;C_O];
C_rot=zeros(14,3);

r_th=50;  % set the second ring as reference ring
Detector_size=401; % with 401*401 pixels
k_max=4; % maximum h,k,l value
N_r=(Detector_size-1)/2;
N_rot=360;  % possible of orientation of each crystal
N_unit_cell=15;   % number of unit cells in one crystal.
angle_all=zeros(1,N_rot);
delta_angle=2*pi/N_rot;

for i=1:N_rot
    angle_all(i)=(i-1)*delta_angle;
end
Rot_z=zeros(3,3,360);   % creat rotation matrix
for i=1:N_rot
    Rot_z(:,:,i)=[cos(angle_all(i)) sin(angle_all(i)) 0;-sin(angle_all(i)) cos(angle_all(i)) 0; 0 0 1];
end


% set mask the central beam
mask=ones(Detector_size,Detector_size,'single');   
for i=190:210
    for j=190:210
         mask(i,j)=0;
    end
end

% transform from cartesian coordinates to spherical coordinates for Detector.
k_x=linspace(-k_max,k_max,Detector_size);
k_y=linspace(-k_max,k_max,Detector_size);
[k_x,k_y]=meshgrid(k_x,k_y);

% transform from cartesian coordinates to polar
p_i=1:1:N_r;
p_j=0:2*pi/N_rot:2*pi-2*pi/N_rot;
[p_j,p_i]=meshgrid(p_j,p_i);
[c_x,c_y] = pol2cart(p_j,p_i); % c_i is the x component in cartesian coordinate; c_j is the y component in cartesian coordinate
c_x=round(c_x)+(Detector_size+1)/2*ones(N_r,N_rot);   % map the coordinates into matrix indice.
c_y=round(c_y)+(Detector_size+1)/2*ones(N_r,N_rot);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%        Simulate powder diffraction        %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assume flat Ewald sphere and small angle scattering.
f_O=8*ones(Detector_size);     % atomic scattering factor for oxygen
f_Mg=12*ones(Detector_size);
f_Al=13*ones(Detector_size);

cart_powder_diff=zeros(Detector_size,Detector_size,N_dp,'single');  % cartesian coordinates of intensity distribution on detector
for iter=1:N_dp                               % N_dp is the total number of diffraction patterns in simulation.
    diff_lattice_1=zeros(Detector_size);
    for n=1:N_crystal                         % N_crystals per shot
        Rot=squeeze(Rot_z(:,:,randi(N_rot))); % assign a random rotation operation on crystal along z axis.
        C_rot=C*Rot; % X'=Rot_z*X, X' are the basis after rotation                
        % structure factor for one unit cell
        F_cell=0;
        for i=1:2
            F_cell=f_Mg.*exp(1i*2*pi*(C_rot(i,1).*k_x+C_rot(i,2).*k_y))+F_cell;
        end
        for i=3:6
             F_cell=f_Al.*exp(1i*2*pi*(C_rot(i,1).*k_x+C_rot(i,2).*k_y))+F_cell;
        end
        for i=7:14
            F_cell=f_O.*exp(1i*2*pi*(C_rot(i,1).*k_x+C_rot(i,2).*k_y))+F_cell;
        end 
       
        % coordinates of unit cell after rotation by theta
        k_x1=k_x.*Rot(1,1)+k_y.*Rot(2,1);
        k_y1=k_x.*Rot(1,2)+k_y.*Rot(2,2);
       
        % structure factor for lattice
        F_lattice=((sin(pi*N_unit_cell*(k_x1+0.000001)).*sin(pi*N_unit_cell*(k_y1+0.000001))).^2./(sin(pi*(k_x1+0.000001)).*sin(pi*(k_y1+0.000001))).^2)/(N_unit_cell^4); 
        % structure factor for the whole crystal
        F_xtal=F_cell.*F_lattice;
        % total intensity on detector is the sum of scattering intensity from different crystals
        diff_lattice_1=diff_lattice_1+abs(F_xtal).^2.*mask;   % sum over the intensity from each crystals as they scatter x-ray incoherently.
    end
    cart_powder_diff(:,:,iter)=diff_lattice_1;   % save the diffraction pattern from many crystals in this matrix.
end

parameter=[N_crystal N_dp r_th]; % save basic parameters

b1=[int2str(N_dp) 'N_dp_' int2str(N_crystal) '_parameter.mat'];
b2=[int2str(N_dp) 'N_dp_' int2str(N_crystal) '_powder_diff.mat'];
save(b1,'parameter')   % save basic parameters
save(b2,'cart_powder_diff')     % save all simulated diffraction patterns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Calculate pair and triple autocorrelation   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_phi=2*pi/N_rot;
%N_r=round(sqrt(2)*800);
c_phi=zeros(N_r,N_rot); % used to calculate autocorrelation for one ring
t_phi=zeros(N_r,N_rot);  % used to calculate triple autocorrelation for one ring
p_phi=zeros(N_r,N_rot);  % used to calculate pair correlation between a given ring and other ringss.

% Calculate SAXS term
SAXS_C=mean(cart_powder_diff,3);   % SAXS in Cartesian coordinate
SAXS_P=zeros(N_r,N_rot);   % SAXS in polar coordinate
for i=1:N_r
    for j=1:N_rot
        SAXS_P(i,j)=SAXS_C(c_x(i,j),c_y(i,j));
    end
end
SAXS_R=mean(SAXS_P,2);   % average SAXS value for one ring

% Convert coordinates of powder diffraction pattern in polar coordinates
pol_powder_diff=zeros(N_r,N_rot,N_dp);
for k=1:N_dp
    for i=1:N_r
        for j=1:N_rot
            pol_powder_diff(i,j,k)=cart_powder_diff(c_x(i,j),c_y(i,j),k);
        end
    end
end

% calculate experimental angular autocorrelation, pair-correlation and triple autocorrelation for one ring
auto_corr_exp=zeros(N_r,N_rot);
pair_corr_exp=zeros(N_r,N_rot);
tri_corr_exp=zeros(N_r,N_rot);
SAXS_1=repmat(SAXS_R,[1 N_rot N_dp]);
I_ref=repmat(pol_powder_diff(r_th,:,:),[N_r 1 1]);
for phi=1:N_rot
    auto_corr_exp(:,phi)=mean(mean((pol_powder_diff-SAXS_1).*(circshift(pol_powder_diff,[0,phi,0])-SAXS_1),3),2);
    tri_corr_exp(:,phi)=mean(mean((pol_powder_diff-SAXS_1).^2.*(circshift(pol_powder_diff,[0,phi,0])-SAXS_1),3),2); 
    pair_corr_exp(:,phi)=mean(mean((pol_powder_diff-SAXS_1).*(circshift(I_ref,[0,phi,0])-SAXS_1),3),2);
end

% calculate angular autocorrelation, pair-correlation and triple autocorrelation for one ring
SAXS_2=repmat(SAXS_R,[1 N_rot]);
SAXS_rth=repmat(SAXS_R(r_th),[N_r,N_rot]);
auto_corr_the=auto_corr_exp/N_crystal+(SAXS_2/N_crystal).^2;
pair_corr_the=pair_corr_exp/N_crystal+SAXS_2.*SAXS_rth./(N_crystal.^2);
tri_corr_the=tri_corr_exp/N_crystal+2*SAXS_2.*tri_corr_exp/N_crystal-(SAXS_2/N_crystal).^3;

b3=[int2str(N_dp) 'dp_SAXS_R.mat'];
b4=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_ac.mat'];
b5=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_tc.mat'];
b6=[int2str(N_dp) '_dp_N_xstal' int2str(N_crystal) '_pc.mat'];
save(b3,'SAXS_R')             % save SAXS term
save(b4,'auto_corr_the')      % save auto correlation function of 1 xtal which is calculated from experimental data
save(b5,'tri_corr_the')       % save triple correlation function of 1 xtal which is calculated from experimental data
save(b6,'pair_corr_the')      % save pair correlation function of 1 xtal which is calculated from experimental data
