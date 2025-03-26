% ensemble stoc-stat filter for the triad system
clear
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 20)

% Monte Carlo - integration parameters
MC=1*10^4;
Dt=1e-3;
tstep=10;
T=10;
Nt=round(T/Dt);

% struct for useful model parameters
params = struct('MC',MC,'Dt',Dt,'T',T,'Nt',Nt,'tstep',tstep);

% statistical regime
params.L=zeros(3);
params.G=diag([0.02,0.01,0.01]);
params.B(1)=1; params.B(2)=-0.6; params.B(3)=-0.4;
params.req=2.5;
params.SS=params.req*sqrt(2*params.G);
params.SSt=zeros(3,3);
for j=1:params.Nt
    params.tmS(j)=0;
    params.tmM(j,:)=0*[1;1;1];
end

% Initial conditions
params.m0 = [3,.1,-.1];
params.var0 = [.5,.01,.01];

%% run direct Monte-Carlo simulation 
% direct MC
params.MC=1*10^5;
[u,TT_MC, samps,Ts_MC, mean_MC, cov_MC, M3_MC_norm, M3_MC, Ene_MC, Ene_dyn] = MC_triad_direct(params, params.m0, params.var0);

% draw results
% stats 
figure;
subplot(2,2,1);
 plot(TT_MC,mean_MC(1,:),'LineWidth',1); hold on;
 plot(TT_MC,mean_MC(2,:),'LineWidth',1);
 plot(TT_MC,mean_MC(3,:),'LineWidth',1);
 legend('<u_1>','<u_2>','<u_3>');
 xlabel('time'); ylabel('mean')
subplot(2,2,2);
 plot(TT_MC,reshape(cov_MC(1,1,:),1,length(TT_MC)), 'LineWidth',1); hold on;
 plot(TT_MC,reshape(cov_MC(2,2,:),1,length(TT_MC)), 'LineWidth',1);
 plot(TT_MC,reshape(cov_MC(3,3,:),1,length(TT_MC)), 'LineWidth',1);
 legend('var u_1','var u_2','var u_3'); %,'total var');
 xlabel('time'); ylabel('variance')
subplot(2,2,3);
 plot(TT_MC,reshape(cov_MC(1,2,:),1,length(TT_MC)), 'LineWidth',1); hold on;
 plot(TT_MC,reshape(cov_MC(1,3,:),1,length(TT_MC)), 'LineWidth',1);
 plot(TT_MC,reshape(cov_MC(2,3,:),1,length(TT_MC)), 'LineWidth',1);
 legend('cov(u_1,u_2)','cov(u_1,u_3)','cov(u_2,u_3)');
 xlabel('time'); ylabel('cross-covariance')
subplot(2,2,4);
 plot(TT_MC,reshape(M3_MC(1,2,3,:),1,length(TT_MC)),'LineWidth',2); hold on;
 plot(TT_MC,reshape(M3_MC(1,2,2,:),1,length(TT_MC)),'LineWidth',1);
 plot(TT_MC,reshape(M3_MC(1,3,3,:),1,length(TT_MC)),'LineWidth',1);
 plot(TT_MC,reshape(M3_MC(2,2,3,:),1,length(TT_MC)),'LineWidth',1);
 legend('<M_1_2_3>','<M_1_2_2>','<M_1_3_3>','<M_2_2_3>');
 xlabel('time'); ylabel('3rd order central moments')

 figure
 subplot(1,3,1)
 scatter(u(:,1),u(:,2),'.'); xlabel('u_1'); ylabel('u_2');
 title('scatter plots for the joint distribution')
 subplot(1,3,2)
 scatter(u(:,2),u(:,3),'.'); xlabel('u_2'); ylabel('u_3');
 subplot(1,3,3)
 scatter(u(:,3),u(:,1),'.'); xlabel('u_3'); ylabel('u_1');

 % PDFs
 figure
 for ii=1:3
     subplot(1,3,ii)
     pdf_u = histogram(u(:,ii),200,'Normalization','pdf'); 
     xx= pdf_u.BinEdges(2:end);
     plot(xx,pdf_u.Values,'LineWidth',2); set(gca,'yscale','log'); hold on
     xlabel(['u_',num2str(ii)]); ylabel(['p(u_',num2str(ii),')']);
    
     % gaussian fit
     sm = mean(u(:,ii));
     sv = var(u(:,ii));
     pdf_norm = normpdf(xx,sm,sqrt(sv));
     plot(xx,pdf_norm,'k--','LineWidth',1);
 end

 % linear instability
 Bsis=eye(3);
 for j1=1:length(TT_MC)
    U=mean_MC(:,j1);
    for i=1:3
        for ii=1:3
            Lvs(i,ii)=(Buu2(params.B,U,Bsis(:,ii))+Buu2(params.B,Bsis(:,ii),U)+(params.L-params.G)*Bsis(:,ii))'*Bsis(:,i);
        end
    end
    EIGLR(:,j1)=sort(real(eig(Lvs)));
    EIGLI(:,j1)=sort(imag(eig(Lvs)));

 end
 
figure; 
subplot(1,2,1);
 plot(TT_MC,(EIGLR),'k'); 
 xlabel('time'); title('Real part of eigenvalues of L_v(<u>)')
subplot(1,2,2);
 plot(TT_MC,(EIGLI),'k'); 
 xlabel('time'); title('Imaginary part of eigenvalues of L_v(<u>)')


%% generate observation data
params.MC=1*10^5;
params.MC1=100;
params.eps=.1;
params.tstep=1;
m0=params.m0;
var0=params.var0;
[zt,TT,z_samps,Ts, mean_u,cov_u, samp_U,samp_R,mean_samps,cov_samps, mom_u] = triad_stocstat_samps(params, m0, var0);

figure;
subplot(3,1,1);
 plot(TT_MC,mean_MC(1,:),'LineWidth',1); hold on; plot(TT,mean_u(1,:),'--');
 xlabel('time'); ylabel('u_1')
title('comparison for the mean')
subplot(3,1,2);
 plot(TT_MC,mean_MC(2,:),'LineWidth',1); hold on; plot(TT,mean_u(2,:),'--');
 xlabel('time'); ylabel('u_2')
subplot(3,1,3);
 plot(TT_MC,mean_MC(3,:),'LineWidth',1); hold on; plot(TT,mean_u(3,:),'--');
 xlabel('time'); ylabel('u_3')
figure
subplot(3,1,1);
 plot(TT_MC,reshape(cov_MC(1,1,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(1,1,:),1,length(TT)), '--'); 
 xlabel('time'); ylabel('var u_1')
title('comparison for the variance')
subplot(3,1,2);
 plot(TT_MC,reshape(cov_MC(2,2,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(2,2,:),1,length(TT)), '--'); 
 xlabel('time'); ylabel('var u_2')
subplot(3,1,3);
 plot(TT_MC,reshape(cov_MC(3,3,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(3,3,:),1,length(TT)), '--'); 
 xlabel('time'); ylabel('var u_3')
figure
subplot(3,1,1);
 plot(TT_MC,reshape(cov_MC(1,2,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(1,2,:),1,length(TT)), '--');
 xlabel('time'); ylabel('cov(u_1,u_2)')
title('comparison for the covariance')
subplot(3,1,2);
 plot(TT_MC,reshape(cov_MC(1,3,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(1,3,:),1,length(TT)), '--');
 xlabel('time'); ylabel('cov(u_1,u_3)')
subplot(3,1,3);
 plot(TT_MC,reshape(cov_MC(2,3,:),1,length(TT_MC)), 'LineWidth',1); hold on; plot(TT,reshape(cov_u(2,3,:),1,length(TT)), '--');
 xlabel('time'); ylabel('cov(u_2,u_3)')


nsamp=100;
figure
subplot(3,1,1)
plot(TT(1:end),squeeze(mean_samps(1,1:nsamp,:)),'LineWidth',.5); hold on; plot(TT,mean_u(1,:),'k','LineWidth',3);
title('stoc-stat model forecast in leading moments, regime II'); ylabel('mean');
subplot(3,1,2)
plot(TT(1:end),squeeze(cov_samps(1,1,1:nsamp,:)),'LineWidth',.5); hold on; plot(TT,reshape(cov_u(1,1,:),1,length(TT)), 'k','LineWidth',3);
ylabel('variance')
subplot(3,1,3);
 plot(TT_MC,(EIGLR),'k'); 
 ylabel('eig[L(u)]')
xlabel('time'); 

%% filtering models
params.MC=100;   % small ensemble size
params.eps=.1;   % relaxation parameter

% observation noises
params.sig_m=diag([0.2853,0.1188,0.0967]);
params.sig_v=[0.7188,0.5202,0.4370;
              0.5202,0.4126,0.3093;
              0.4370,0.3093,0.2752];

params.skip=1;
obs_mean = mean_u(:,1:end);
obs_cov = cov_u(:,:,1:end);
m0=params.m0;
var0=params.var0;
rng(177);

% high-order filter
[zt_fil,TT_fil,z_fil,Ts_fil, fil_mean,fil_cov, fil_mom] = triad_filtering_model(params, m0,var0, obs_mean,obs_cov);
% EnKF
[zt_fil1,TT_fil1,z_fil1,Ts_fil1, fil_mean1,fil_cov1, fil1_mom] = triad_filtering_enkf(params, m0,var0, obs_mean,obs_cov);
% no filter
[zt_nofil,TT_nofil,z_nofil,Ts_nofil, nofil_mean,nofil_cov, nofil_mom] = triad_stocstat_model(params, m0, var0);


% filter result
figure;
subplot(3,1,1);
 plot(TT,mean_u(1,:),'k','LineWidth',3); hold on; 
 plot(TT_fil,fil_mean(1,:),'--','LineWidth',2);
 plot(TT_fil1,fil_mean1(1,:),'-.','LineWidth',1.5);
 plot(TT_nofil,nofil_mean(1,:),':','LineWidth',2);
 legend('truth','optimal filter','EnKF','no filter')
 ylabel('mean(u_1)'); 
title('prediction of mean state')
subplot(3,1,2);
 plot(TT,mean_u(2,:),'k','LineWidth',3); hold on; 
 plot(TT_fil,fil_mean(2,:),'--','LineWidth',2); 
 plot(TT_fil1,fil_mean1(2,:),'-.','LineWidth',1.5);
 plot(TT_nofil,nofil_mean(2,:),':','LineWidth',2);
 ylabel('mean(u_2)');
subplot(3,1,3);
 plot(TT,mean_u(3,:),'k','LineWidth',3); hold on; 
 plot(TT_fil,fil_mean(3,:),'--','LineWidth',2); 
 plot(TT_fil1,fil_mean1(3,:),'-.','LineWidth',1.5);
 plot(TT_nofil,nofil_mean(3,:),':','LineWidth',2);
 ylabel('mean(u_3)'); xlabel('time');
figure
subplot(3,1,1);
 plot(TT,reshape(cov_u(1,1,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(1,1,:),1,length(TT_fil)), '--','LineWidth',2); 
 plot(TT_fil1,reshape(fil_cov1(1,1,:),1,length(TT_fil1)), '-.','LineWidth',1.5); 
 plot(TT_nofil,reshape(nofil_cov(1,1,:),1,length(TT_nofil)), ':','LineWidth',2); 
 legend('truth','optimal filter','EnKF','no filter')
 ylabel('var(u_1)');
title('prediction of variance')
subplot(3,1,2);
 plot(TT,reshape(cov_u(2,2,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(2,2,:),1,length(TT_fil)), '--','LineWidth',2); 
 plot(TT_fil1,reshape(fil_cov1(2,2,:),1,length(TT_fil1)), '-.','LineWidth',1.5);
 plot(TT_nofil,reshape(nofil_cov(2,2,:),1,length(TT_nofil)), ':','LineWidth',2);
 ylabel('var(u_2)');
subplot(3,1,3);
 plot(TT,reshape(cov_u(3,3,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(3,3,:),1,length(TT_fil)), '--','LineWidth',2); 
 plot(TT_fil1,reshape(fil_cov1(3,3,:),1,length(TT_fil1)), '-.','LineWidth',1.5);
 plot(TT_nofil,reshape(nofil_cov(3,3,:),1,length(TT_nofil)), ':','LineWidth',2);
 ylabel('var(u_3)'); xlabel('time'); 
figure
subplot(3,1,1);
 plot(TT,reshape(cov_u(1,2,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(1,2,:),1,length(TT_fil)), '--','LineWidth',2); 
 plot(TT_fil1,reshape(fil_cov1(1,2,:),1,length(TT_fil1)), '-.','LineWidth',1.5); 
 plot(TT_nofil,reshape(nofil_cov(1,2,:),1,length(TT_nofil)), ':','LineWidth',2);
 legend('truth','optimal filter','EnKF','no filter')
 ylabel('cov(u_1,u_2)');
title('prediction of covariance')
subplot(3,1,2);
 plot(TT,reshape(cov_u(1,3,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(1,3,:),1,length(TT_fil)), '--','LineWidth',2); 
 plot(TT_fil1,reshape(fil_cov1(1,3,:),1,length(TT_fil1)), '-.','LineWidth',1.5); 
 plot(TT_nofil,reshape(nofil_cov(1,3,:),1,length(TT_nofil)), ':','LineWidth',2);
 ylabel('cov(u_1,u_3)');
subplot(3,1,3);
 plot(TT,reshape(cov_u(2,3,:),1,length(TT)),'k','LineWidth',3); hold on; 
 plot(TT_fil,reshape(fil_cov(2,3,:),1,length(TT_fil)), '--','LineWidth',2);
 plot(TT_fil1,reshape(fil_cov1(2,3,:),1,length(TT_fil1)), '-.','LineWidth',1.5); 
 plot(TT_nofil,reshape(nofil_cov(2,3,:),1,length(TT_nofil)), ':','LineWidth',2);
 ylabel('cov(u_2,u_3)'); xlabel('time'); 


figure
 plot(TT,squeeze(mom_u(3,3,:)),'LineWidth',2); hold on
 plot(TT,squeeze(mom_u(1,2,:)+mom_u(2,1,:)),'LineWidth',2);
 plot(TT,squeeze(mom_u(1,3,:)+mom_u(3,1,:)),'LineWidth',2);
 plot(TT,squeeze(mom_u(2,3,:)+mom_u(3,2,:)),'LineWidth',2);
 xlabel('time');
 title('ensemble prediction of the high-order moments')
 plot(TT,squeeze(fil_mom(3,3,:)),'--','LineWidth',1); hold on
 plot(TT,squeeze(fil_mom(1,2,:)+fil_mom(2,1,:)),'--','LineWidth',1);
 plot(TT,squeeze(fil_mom(1,3,:)+fil_mom(3,1,:)),'--','LineWidth',1);
 plot(TT,squeeze(fil_mom(2,3,:)+fil_mom(3,2,:)),'--','LineWidth',1);