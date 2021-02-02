function [A_l,b_l,Omega_l]=SLR_ungm_dynamic(meank,Pk,alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,k)


Nz=Nx;

%We obtain the sigma points
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];

sigma_points=repmat(meank,1,length(weights))+sigma_points;


%Transformed sigma points
sigma_points_z=alfa_mod*sigma_points+beta_mod*sigma_points./(1+sigma_points.^2)+gamma_mod*cos(1.2*k);

z_pred=sigma_points_z*weights';


var_pred=zeros(Nz);
var_xz=zeros(Nx,Nz);

for j=1:length(weights)
    sub_z_j=sigma_points_z(:,j)-z_pred;
    var_pred=var_pred+weights(j)*(sub_z_j*sub_z_j');
    var_xz=var_xz+weights(j)*(sigma_points(:,j)-meank)*sub_z_j';
end

%Statistical linearisaltion

A_l=var_xz'/Pk;
b_l=z_pred-A_l*meank;
Omega_l=var_pred-A_l*Pk*A_l';