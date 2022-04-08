%Demo of the iterated posterior linearisation smoother with the same scenario with
%a univariate non-stationary growth model as in

%Á. F. García-Fernández, L. Svensson and S. Särkkä, "Iterated Posterior Linearization Smoother,"
%in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 2056-2063, April 2017

%First, statistical linear regression (SLR) using IPLF, then smoother and then iterated
%SLRs, with filtering and smoothing

clear
randn('seed',9)
rand('seed',9)

Scenario_ungm_trajectories;

%Parameters of the sigma-point method (unscented transform)
Nx=1; %Dimension of the state
W0=1/3; %Weight of sigma-point at the mean
Nz=1; %Dimension of the measurement

Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)]; %Sigma-points weights (according to the unscented transform)
square_error_t_tot=zeros(1,Nsteps);


nees_t_tot=zeros(1,Nsteps);
rms_t_series=zeros(1,Nmc);
square_error_t_series=zeros(Nmc,Nsteps);

square_error_t_tot_smoothing=zeros(1,Nsteps);


%Nit_s-> Number of iterations in smoothing
Nit_s=10;
Nit_iplf=1; %Number of iterations in filtering (Nit_iplf=1 is similar to the UKF)

square_error_t_tot_smoothing_series=zeros(Nmc,Nsteps,Nit_s+1); %For filtering and the N_it_s iterations
NLL_smoothing_series=zeros(Nmc,Nsteps,Nit_s+1); %Negative log-likelihood (See Marc Deisenroth PhD thesis)
Nees_smoothing_series=zeros(Nmc,Nsteps,Nit_s+1); %NEES
cte_NLL=Nx/2*log(2*pi); %Constant for NLL

mean_ini=x0;


randn('seed',9)
rand('seed',9)

%Number of Monte Carlo runs

for i=1:Nmc
    tic
    
    
    n_trajectory=fix((i-1)/Nmc_trajectory)+1;
    X_multi=X_multi_series(n_trajectory,:);
    
    meank=mean_ini;
    Pk=P_ini;
    
    
    square_error_t=zeros(1,Nsteps);
    square_error_t_smoothing=zeros(1,Nsteps);
    
    nees_t=zeros(1,Nsteps);
    
    meank_t=zeros(Nx,Nsteps);
    Pk_t=zeros(Nx,Nx,Nsteps);
    
    
    %SLR parameters for dynamics
    A_dyn=zeros(Nx,Nx,Nsteps);
    b_dyn=zeros(Nx,Nsteps);
    Omega_dyn=zeros(Nx,Nx,Nsteps);
    
    %SLR parameters for measurements
    A_m=zeros(Nz,Nx,Nsteps);
    b_m=zeros(Nz,Nsteps);
    Omega_m=zeros(Nz,Nz,Nsteps);
    
    
    
    %Generation of measurements
    z_real_t=zeros(Nz,Nsteps);
    for k=1:Nsteps
        state_k=X_multi(:,k);
        noise=chol_R*noise_z(:,k+Nsteps*(i-1));
        z_real=a*state_k^3+noise;
        z_real_t(:,k)=z_real;
        
    end
    
    %Filtering
    
    for k=1:Nsteps
        pos_x=X_multi(1,k);
        
        z_real=z_real_t(:,k);
        
        %Calculate iterated SLR
        
        meank_j=meank;
        Pk_j=Pk;
        
        for p=1:Nit_iplf
            
            %SLR of function a*x^3 w.r.t. current approximation of the
            %posterior, with mean meank_j and covariance Pk_j
            [A_l,b_l,Omega_l]=SLR_measurement_ax3(meank_j,Pk_j,a,weights,W0,Nx,Nz);
            
            
            %KF update for the linearised model
            [mean_ukf_act,var_ukf_act]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R,z_real);
            
            meank_j=mean_ukf_act;
            Pk_j=var_ukf_act;
            
        end
        
        
        A_m(:,:,k)=A_l;
        b_m(:,k)=b_l;
        Omega_m(:,:,k)=Omega_l;
        
        meank_t(:,k)=mean_ukf_act;
        Pk_t(:,:,k)=var_ukf_act;
        
        
        square_error_t(k)=(mean_ukf_act-pos_x)^2;
        pos_error=mean_ukf_act-pos_x;
        var_pos_act=var_ukf_act(1);
        %nees_t(k)=state_error'*inv(var_mp_tukf_act)*state_error;
        nees_t(k)=pos_error'/var_pos_act*pos_error;

        square_error_t_tot_smoothing_series(i,k,1)=square_error_t(k);
        NLL_smoothing_series(i,k,1)=1/2*log(var_ukf_act(1))+1/2*square_error_t(k)/var_ukf_act(1)+cte_NLL;
        Nees_smoothing_series(i,k,1)=square_error_t(k)/var_ukf_act(1);
        
        
        %Prediction
        [A_dyn_k,b_dyn_k,Omega_dyn_k]=SLR_ungm_dynamic(mean_ukf_act,var_ukf_act,alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,k);
        

        meank=A_dyn_k*mean_ukf_act+b_dyn_k;
        Pk=A_dyn_k*var_ukf_act*A_dyn_k'+Omega_dyn_k+Q;
        Pk=(Pk+Pk')/2;

        
        A_dyn(:,:,k)=A_dyn_k;
        b_dyn(:,k)=b_dyn_k;
        Omega_dyn(:,:,k)=Omega_dyn_k;
        
        
    end
    
    
    %Smoothing
    [meank_smoothed_t,Pk_smoothed_t]=linear_rts_smoother(meank_t,Pk_t,A_dyn,b_dyn,Omega_dyn,Q);
    
    for k=1:Nsteps
        pos_x=X_multi(1,k);
        square_error_t_tot_smoothing_series(i,k,2)=(meank_smoothed_t(1,k)-pos_x)^2;
        NLL_smoothing_series(i,k,2)=1/2*log(Pk_smoothed_t(1,1,k))+1/2*square_error_t_tot_smoothing_series(i,k,2)/Pk_smoothed_t(1,1,k)+cte_NLL;
        Nees_smoothing_series(i,k,2)=square_error_t_tot_smoothing_series(i,k,2)/Pk_smoothed_t(1,1,k);
    end
    
    for p=1:Nit_s-1
        
        %Iterated SLR using the current posterior
        
        for k=1:Nsteps
            %Generation of sigma points
            meank=meank_smoothed_t(:,k);
            Pk=Pk_smoothed_t(:,:,k);
            
            %SLR for the measurement
            [A_l,b_l,Omega_l]=SLR_measurement_ax3(meank,Pk,a,weights,W0,Nx,Nz);
            
            
            A_m(:,:,k)=A_l;
            b_m(:,k)=b_l;
            Omega_m(:,:,k)=Omega_l;
            
            %SLR for the dynamics
            [A_dyn_k,b_dyn_k,Omega_dyn_k]=SLR_ungm_dynamic(meank,Pk,alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,k);
            
            A_dyn(:,:,k)=A_dyn_k;
            b_dyn(:,k)=b_dyn_k;
            Omega_dyn(:,:,k)=Omega_dyn_k;
            
            
        end
        
        %Filter with the linearised model
        
        [meank_t,Pk_t]=linear_kf_full(mean_ini,P_ini,A_m,b_m,Omega_m,A_dyn,b_dyn,Omega_dyn,R,Q,z_real_t);
        
        
        
        %Smoother
        
        [meank_smoothed_t,Pk_smoothed_t]=linear_rts_smoother(meank_t,Pk_t,A_dyn,b_dyn,Omega_dyn,Q);
        
        for k=1:Nsteps
            pos_x=X_multi(1,k);
            square_error_t_tot_smoothing_series(i,k,p+2)=(meank_smoothed_t(1,k)-pos_x)^2;
            NLL_smoothing_series(i,k,p+2)=1/2*log(Pk_smoothed_t(1,1,k))+1/2*square_error_t_tot_smoothing_series(i,k,p+2)/Pk_smoothed_t(1,1,k)+cte_NLL;
            Nees_smoothing_series(i,k,p+2)=square_error_t_tot_smoothing_series(i,k,p+2)/Pk_smoothed_t(1,1,k);
            
        end
        
    end
    
    
    
    %Square error calculation
    
    for k=1:Nsteps
        pos_x=X_multi(1,k);
        square_error_t_smoothing(k)=(meank_smoothed_t(1,k)-pos_x)^2;  
    end
    
    
    
    square_error_t_tot=square_error_t_tot+square_error_t;
    square_error_t_series(i,:)=square_error_t; 
    square_error_t_tot_smoothing=square_error_t_tot_smoothing+square_error_t_smoothing;
    
    
    nees_t_tot=nees_t_tot+nees_t;
    rms_t_series(i)=sqrt(sum(square_error_t)/(Nsteps));
    
    t=toc;
    display(['Completed iteration nº ', num2str(i),' time ', num2str(t), ' seconds'])
end

square_error_t_tot=square_error_t_tot/Nmc;
rmse_filtering_tot=sqrt(sum(square_error_t_tot)/(Nsteps))

square_error_t_tot_smoothing=square_error_t_tot_smoothing/Nmc;

rmse_tot_smoothing=sqrt(sum(square_error_t_tot_smoothing)/(Nsteps))

rmse_t_tot_smoothing_series=sqrt(sum(squeeze(sum(square_error_t_tot_smoothing_series,1)),1)/(Nmc*Nsteps));


nees_t_tot=nees_t_tot/Nmc;



%Smoothing error for different J
rmse_smoothing_1=sqrt(sum(square_error_t_tot_smoothing_series(:,:,2),1)/Nmc);
rmse_smoothing_5=sqrt(sum(square_error_t_tot_smoothing_series(:,:,6),1)/Nmc);
rmse_smoothing_10=sqrt(sum(square_error_t_tot_smoothing_series(:,:,11),1)/Nmc);


%Output figure (IPLS(i)-J denotes a IPLS with i SLR iterations for the IPLF and
%J SLR smoothing iterations)


figure(1)
plot(1:Nsteps,sqrt(square_error_t_tot),'b',...
    1:Nsteps,rmse_smoothing_1,'--r',...
    1:Nsteps,rmse_smoothing_5,'-.black',...
    1:Nsteps,rmse_smoothing_10,'-*g','Linewidth',1.3)
grid on
legend('IPLS(1)-0','IPLS(1)-1','IPLS(1)-5','IPLS(1)-10')
xlabel('Time step')
ylabel('RMS error')



NLL_average_list=zeros(1,Nit_s);

for i=1:Nit_s+1
    NLL=NLL_smoothing_series(:,:,i);
    
    NLL_average=sum(sum(NLL))/(Nsteps*Nmc);
    NLL_average_list(i)=NLL_average;
    
end


