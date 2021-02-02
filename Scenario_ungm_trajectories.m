

Nsteps=50; %Number of time steps
Nmc=1000; %Number of Monte Carlo runs
Nmc_trajectory=50; %Number of Monte Carlo runs per trajectory
N_trajectories=Nmc/Nmc_trajectory; %Number of simulated trajectories


x0=5;
Q=1;
a=1/20;%Parameter multiplying x^3
R=1;
P_ini=4;

alfa_mod=0.9;
beta_mod=10;
gamma_mod=8;

chol_ini=chol(P_ini)';
chol_Q=chol(Q)';
chol_R=chol(R)';


%Ground truth 
X_multi_series=zeros(N_trajectories,Nsteps);

for i=1:N_trajectories
    X_multi_i=zeros(1,Nsteps);
    xk=x0+chol_ini*randn(1);
    X_multi_i(1)=xk;
    for k=1:Nsteps-1
        xk_pred=alfa_mod*xk+beta_mod*xk/(1+xk^2)+gamma_mod*cos(1.2*k)+chol_Q*randn(1);
        X_multi_i(k+1)=xk_pred;
        xk=xk_pred;
    end
    X_multi_series(i,:)=X_multi_i;
end

noise_z=randn(1,Nsteps*Nmc); %Measurement noise