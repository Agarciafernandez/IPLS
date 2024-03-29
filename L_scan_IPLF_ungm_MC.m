%Demo of the L-scan iterated posterior linearisation filter (IPLF) with the same scenario with
%a univariate non-stationary growth model as in 

%�. F. Garc�a-Fern�ndez, L. Svensson and S. S�rkk�, "Iterated Posterior Linearization Smoother,"
%in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 2056-2063, April 2017

%The L-scan IPLF is described in Section III.D in the above paper. The idea
%is to re-linearise dynamic and measurement functions in a sliding window of
%the past L time steps to improve filtering performance.

%Author: �ngel F. Garc�a-Fern�ndez


clear
randn('seed',9)
rand('seed',9)


Scenario_ungm_trajectories;

%Parametros sigma points
Nx=1;
W0=1/3;
Nz=1;

Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];
square_error_t_tot=zeros(1,Nsteps);


rms_t_series=zeros(1,Nmc);
square_error_t_series=zeros(Nmc,Nsteps);



%Nit_s-> Number of iterations in the L-scan IPLF
Nit_s=5;
L=5; %Length of the L-scan window

NLL_series=zeros(Nmc,Nsteps); %Negative log-likelihood (See Marc Deisenroth PhD thesis)
Nees_series=zeros(Nmc,Nsteps); %NEES
cte_NLL=Nx/2*log(2*pi); %Constant for NLL

mean_ini=x0;


randn('seed',9)
rand('seed',9)


%Number of Monte Carlo runs

for i=1:Nmc
    tic
    
    
    n_trajectory=fix((i-1)/Nmc_trajectory)+1;
    X_multi=X_multi_series(n_trajectory,:);
    
  
    %Means and covariances in the L-scan sliding window

    meank_L=zeros(Nx,L);
    Pk_L=zeros(Nx,Nx,L);
    
    meank_L(:,L)=mean_ini;
    Pk_L(:,:,L)=P_ini;
    
    
    square_error_t=zeros(1,Nsteps);
    
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
        
        kmin=max(1,k-L+1);
        N_w=k-kmin+1; %N_w is the length of the window to relinearise
        
        if(k<=L)
            mean_prior=mean_ini; %Prior mean at the beginning of the L-scan window
            P_prior=P_ini; %Prior covariance at the beginning of the L-scan window
      
        end
        
        
        %We first use a sigma-point Kalman filter update step

        meank=meank_L(:,L);
        Pk=Pk_L(:,:,L);
        
        [A_l,b_l,Omega_l]=SLR_measurement_ax3(meank,Pk,a,weights,W0,Nx,Nz);
        
        %KF update
        [mean_ukf_act,var_ukf_act]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R,z_real);
        
        meank_L(:,L)=mean_ukf_act;
        Pk_L(:,:,L)=var_ukf_act;
        
        
        
        
        if(Nit_s>1)

            %We perform RTS smoothing for the current linearisation in the
            %L-scan window (the window has length N_w, as at the beginning
            %the length of the window is smaller than L)
            
            [meank_smoothed_t,Pk_smoothed_t]=linear_rts_smoother(meank_L(:,L-N_w+1:end),Pk_L(:,:,L-N_w+1:end),A_dyn(:,:,kmin:k),b_dyn(:,kmin:k),Omega_dyn(:,:,kmin:k),Q);
            
            %We perform the iterated SLR in the area of the posterior in
            %the L-scan window
            for p=1:Nit_s-1
                
                for kk=1:N_w
                    %Generation of sigma points
                    meank=meank_smoothed_t(:,kk);
                    Pk=Pk_smoothed_t(:,:,kk);
                    
                    %SLR of measurement function
                    [A_l,b_l,Omega_l]=SLR_measurement_ax3(meank,Pk,a,weights,W0,Nx,Nz);
                    
                    
                    A_m(:,:,kmin+kk-1)=A_l;
                    b_m(:,kmin+kk-1)=b_l;
                    Omega_m(:,:,kmin+kk-1)=Omega_l;
                    
                    %SLR of dynamic function
                    [A_dyn_k,b_dyn_k,Omega_dyn_k]=SLR_ungm_dynamic(meank,Pk,alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,kmin+kk-1);
                    
                    A_dyn(:,:,kmin+kk-1)=A_dyn_k;
                    b_dyn(:,kmin+kk-1)=b_dyn_k;
                    Omega_dyn(:,:,kmin+kk-1)=Omega_dyn_k;
                    
                    
                end
                
                %Filter with the linearised model in the L-scan window
                [meank_t,Pk_t]=linear_kf_full(mean_prior,P_prior,A_m(:,:,kmin:k),b_m(:,kmin:k),...
                    Omega_m(:,:,kmin:k),A_dyn(:,:,kmin:k),b_dyn(:,kmin:k),Omega_dyn(:,:,kmin:k),R,Q,z_real_t(:,kmin:k));
                         
                %Smoother with the linearised model in the L-scan window
                [meank_smoothed_t,Pk_smoothed_t]=linear_rts_smoother(meank_t,Pk_t,A_dyn(:,:,kmin:k),b_dyn(:,kmin:k),Omega_dyn(:,:,kmin:k),Q);
                          
            end
            
            meank_L(:,L-N_w+1:end)=meank_t; %We store the values of the filtering output
            Pk_L(:,:,L-N_w+1:end)=Pk_t;
            
            mean_ukf_act=meank_t(:,end);
            var_ukf_act=Pk_t(:,:,end);
            
            %If k>=L, we set mean_prior and P_prior using the prediction of the
            %oldest state in the L-scan window (whose information
            %disappears from the L-scan window)

             if(k>=L)              
                 [A_dyn_k,b_dyn_k,Omega_dyn_k]=SLR_ungm_dynamic(meank_t(:,1),Pk_t(:,:,1),alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,kmin);
                 
                 mean_prior= A_dyn_k*meank_t(:,1)+b_dyn_k;
                 P_prior=A_dyn_k*Pk_t(:,:,1)*A_dyn_k'+Omega_dyn_k+Q;
                 P_prior=(P_prior+P_prior')/2;
               
             end
            
            
        else
            meank_L(:,L-N_w+1:end)=mean_ukf_act; %We store the values of the filtering output
            Pk_L(:,:,L-N_w+1:end)=var_ukf_act;

            %NOTE: mean_prior and P_prior are not used in this case as
            %there is not re-linearisation of past dynamic and measurement
            %functions
        end
        
        
        square_error_t(k)=(mean_ukf_act-pos_x)^2;
        pos_error=mean_ukf_act-pos_x;
        var_pos_act=var_ukf_act(1);
        nees_t(k)=pos_error'/var_pos_act*pos_error;
        square_error_t_series(i,k)=square_error_t(k);
        NLL_series(i,k)=1/2*log(var_ukf_act(1))+1/2*square_error_t(k)/var_ukf_act(1)+cte_NLL;
        Nees_series(i,k)=square_error_t(k)/var_ukf_act(1);
           
        
        %Prediction
        
        [A_dyn_k,b_dyn_k,Omega_dyn_k]=SLR_ungm_dynamic(mean_ukf_act,var_ukf_act,alfa_mod,beta_mod,gamma_mod,weights,W0,Nx,k);
        
        
        meank=A_dyn_k*mean_ukf_act+b_dyn_k;
        Pk=A_dyn_k*var_ukf_act*A_dyn_k'+Omega_dyn_k+Q;
        Pk=(Pk+Pk')/2;
        
        
              
        A_dyn(:,:,k)=A_dyn_k;
        b_dyn(:,k)=b_dyn_k;
        Omega_dyn(:,:,k)=Omega_dyn_k;
        
        %We store the prediction with the rest of the values in the L-scan
        %window
        
        meank_L_pred=zeros(Nx,L);
        Pk_L_pred=zeros(Nx,Nx,L);
        
        meank_L_pred(:,1:L-1)=meank_L(:,2:L);
        meank_L_pred(:,L)=meank;
        
        Pk_L_pred(:,:,1:L-1)=Pk_L(:,:,2:L);
        Pk_L_pred(:,:,L)=Pk;
        
        meank_L=meank_L_pred;
        Pk_L=Pk_L_pred;
        
        
    end
    
    square_error_t_tot=square_error_t_tot+square_error_t;
    
    t=toc;
    display(['Completed iteration n� ', num2str(i),' time ', num2str(t), ' seconds'])
end

square_error_t_tot=square_error_t_tot/Nmc;
rmse_filtering_tot=sqrt(sum(square_error_t_tot)/(Nsteps))






figure(1)
plot(1:Nsteps,sqrt(square_error_t_tot),'b','Linewidth',1.3)
ylabel ('RMS error (m)')
xlabel('Time step')
grid on
legend('L-scan IPLF')




