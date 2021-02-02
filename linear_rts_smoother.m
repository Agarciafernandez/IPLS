function   [meank_smoothed_t,Pk_smoothed_t]=linear_rts_smoother(meank_t,Pk_t,A_dyn,b_dyn,Omega_dyn,Q)

%RTS smoother

Nsteps=size(meank_t,2);
Nx=size(b_dyn,1);
 
meank_smoothed_t=zeros(Nx,Nsteps);
Pk_smoothed_t=zeros(Nx,Nx,Nsteps);
 
meank_smoothed=meank_t(:,Nsteps);
Pk_smoothed=Pk_t(:,:,Nsteps);

meank_smoothed_t(:,Nsteps)=meank_smoothed;
Pk_smoothed_t(:,:,Nsteps)=Pk_smoothed;


for k=Nsteps-1:-1:1
    
    meank=meank_t(:,k);
    Pk=Pk_t(:,:,k);
    
    A_dyn_k=A_dyn(:,:,k);
    b_dyn_k=b_dyn(:,k);
    Omega_dyn_k=Omega_dyn(:,:,k);
    
    
    meank_pred=A_dyn_k*meank+b_dyn_k;
    cov_pred=A_dyn_k*Pk*A_dyn_k'+Omega_dyn_k+Q;
    Gain_smoother=(Pk*A_dyn_k')/cov_pred;
    
   

    meank_1_smoothed=meank+Gain_smoother*(meank_smoothed-meank_pred);
    Pk_1_smoothed=Pk+Gain_smoother*(Pk_smoothed-cov_pred)*Gain_smoother';
    
    Pk_1_smoothed=(Pk_1_smoothed+Pk_1_smoothed')/2;
    
    meank_smoothed=meank_1_smoothed;
    Pk_smoothed=Pk_1_smoothed;
    meank_smoothed_t(:,k)=meank_smoothed;
    Pk_smoothed_t(:,:,k)=Pk_smoothed;
    
    
end
