function [A_l,b_l,Omega_l]=SLR_measurement_ax3(meank,Pk,a,weights,W0,Nx,Nz)

%Statistical linear regression for the function h(x)=a*x^3

chol_var_mult=chol((Nx/(1-W0)*Pk));

sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];

sigma_points=repmat(meank,1,length(weights))+sigma_points;

%Transformation of sigma points
sigma_points_z=a*sigma_points.^3;

%Required moments
z_pred=sigma_points_z*weights';
var_pred=zeros(Nz);
var_xz=zeros(Nx,Nz);

for j=1:length(weights)
    sub_z_j=sigma_points_z(:,j)-z_pred;
    var_pred=var_pred+weights(j)*(sub_z_j*sub_z_j');
    var_xz=var_xz+weights(j)*(sigma_points(:,j)-meank)*sub_z_j';
end

%Statistical linear regression parameters

A_l=var_xz'/Pk;
b_l=z_pred-A_l*meank;
Omega_l=var_pred-A_l*Pk*A_l';