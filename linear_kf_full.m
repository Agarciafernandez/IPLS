function [meank_t,Pk_t]=linear_kf_full(mean_ini,P_ini,A_m,b_m,Omega_m,A_dyn,b_dyn,Omega_dyn,R,Q,z_real_t)

%Kalman filter on the whole sequence

Nsteps=size(b_m,2);
Nx=size(mean_ini,1);

meank_t=zeros(Nx,Nsteps);
Pk_t=zeros(Nx,Nx,Nsteps);

meank=mean_ini;
Pk=P_ini;

for k=1:Nsteps
   
    z_real=z_real_t(:,k);
    %Update
    [mean_ukf_act,var_ukf_act]=linear_kf_update(meank,Pk,A_m(:,:,k),b_m(:,k),Omega_m(:,:,k),R,z_real);

    meank_t(:,k)=mean_ukf_act;
    Pk_t(:,:,k)=var_ukf_act;
    
    %Prediction
    A_dyn_k=A_dyn(:,:,k);
    b_dyn_k=b_dyn(:,k);
    Omega_dyn_k=Omega_dyn(:,:,k);
    
    
    meank=A_dyn_k*mean_ukf_act+b_dyn_k;
    Pk=A_dyn_k*var_ukf_act*A_dyn_k'+Q+Omega_dyn_k;
    Pk=(Pk+Pk')/2;
    
    
    
end
