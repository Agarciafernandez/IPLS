function [mean_u,P_u]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R,z_real)

%Kalman filter update 

z_pred_j=A_l*meank+b_l;
Psi_j=Pk*A_l';

S_j=A_l*Pk*A_l'+R+Omega_l;

resta_z_ukf=z_real-z_pred_j;

Gain=Psi_j/S_j;

mean_u=meank+Gain*resta_z_ukf;
P_u=Pk-Gain*Psi_j';

P_u=(P_u+P_u')/2;