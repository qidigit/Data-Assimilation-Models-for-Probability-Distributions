function [u,TT,samps,Ts, mean_MC, cov_MC, M3_MC_norm, M3_MC, Ene_MC, Ene_dyn] = MC_triad_direct(params, m0, var0)
MC=params.MC; Dt=params.Dt; tstep=params.tstep; Nt=params.Nt;
L=params.L; G=params.G; B=params.B; 

% Initial conditions
u(:,1)=sqrt(var0(1))*randn(MC,1)+m0(1);
u(:,2)=sqrt(var0(2))*randn(MC,1)+m0(2);
u(:,3)=sqrt(var0(3))*randn(MC,1)+m0(3);


% Integration scheme
samp=u;
ts=1;
TT(ts)=0;
mean_MC(:,ts)=mean(u,1)';
cov_MC(:,:,ts)=cov(u);
for k1=1:3
    for k2=1:3
        for k3=1:3
            M3_MC(k1,k2,k3,ts) = mean((u(:,k1)-mean(u(:,k1))).*(u(:,k2)-mean(u(:,k2))).*(u(:,k3)-mean(u(:,k3))));
            M3_MC_norm(k1,k2,k3,ts)    = M3_MC(k1,k2,k3,ts)/sqrt(cov_MC(k1,k1,ts)*cov_MC(k2,k2,ts)*cov_MC(k3,k3,ts));
        end
    end
end
% energy
Ene_MC(1,ts)=0.5*(sum(mean_MC(:,ts).^2)+sum(diag(cov_MC(:,:,ts))));
Ene_MC(2,ts)=0.5*(mean_MC(1,ts)^2+cov_MC(1,1,ts));
Ene_MC(3,ts)=0.5*(mean_MC(2,ts)^2+cov_MC(2,2,ts));
Ene_MC(4,ts)=0.5*(mean_MC(3,ts)^2+cov_MC(3,3,ts));

Ene_dyn(1,ts)=0.5*(sum(mean_MC(:,ts).^2)+sum(diag(cov_MC(:,:,ts)))); % total energy
Ene_dyn(2,ts)=Ene_dyn(1,ts);                                         % lower bound
Ene_dyn(3,ts)=Ene_dyn(1,ts);                                         % upper bound
Ene_dyn(4,ts)=Ene_dyn(1,ts);                                         % mean
Ene_update=Ene_dyn;
% save samples
ic=1; ns=floor(1/Dt);
for ii=1:Nt
    t=Dt*(ii);     % time

    % FE
    % samp=samp+Dt*((L*samp')'-samp*G+Buu(B,samp,samp)+ones(MC,1)*params.tmM(ii,:));
    % RK4
    k1 = (L*samp')'-samp*G+Buu(B,samp,samp)+ones(MC,1)*params.tmM(ii,:);
    samp1 = samp+.5*Dt*k1;
    k2 = (L*samp1')'-samp1*G+Buu(B,samp1,samp1)+ones(MC,1)*params.tmM(ii,:);
    samp2 = samp+.5*Dt*k2;
    k3 = (L*samp2')'-samp2*G+Buu(B,samp2,samp2)+ones(MC,1)*params.tmM(ii,:);
    samp3 = samp+Dt*k3;
    k4 = (L*samp3')'-samp3*G+Buu(B,samp3,samp3)+ones(MC,1)*params.tmM(ii,:);
    samp=samp+Dt*(k1/6+k2/3+k3/3+k4/6);

    % noise term
    SS=params.SS+params.tmS(ii)^2*(params.SSt-params.SS);

    Winc(:,1) = randn(MC,1);
    Winc(:,2) = randn(MC,1);
    Winc(:,3) = randn(MC,1);
    samp=samp+sqrt(Dt)*Winc*SS;
    
    Ene_update(1)=Ene_update(1)+Dt * (- sum(diag(G).*((mean(samp,1)').^2+diag(cov(samp)))) ...
                                      + sum(params.tmM(ii,:).*mean(samp,1)) + 0.5*sum(diag(SS^2)));
    damp1=max(diag(G)); damp2=max(min(diag(G)),0); damp3=mean(diag(G));
    Ene_update(2)=Ene_update(2)+Dt * (-2*damp1*Ene_update(2) + sum(params.tmM(ii,:).*mean(samp,1)) + 0.5*sum(diag(SS^2)));
    Ene_update(3)=Ene_update(3)+Dt * (-2*damp2*Ene_update(3) + sum(params.tmM(ii,:).*mean(samp,1)) + 0.5*sum(diag(SS^2)));
    Ene_update(4)=Ene_update(4)+Dt * (-2*damp3*Ene_update(4) + sum(params.tmM(ii,:).*mean(samp,1)) + 0.5*sum(diag(SS^2)));
    
    if ~mod(ii,tstep)
        ts=ts+1;
        TT(ts)=t;

        u(:,1)=samp(:,1);
        u(:,2)=samp(:,2);
        u(:,3)=samp(:,3);
        
        mean_MC(:,ts)=mean(u(:,:),1)';
        cov_MC(:,:,ts)=cov(u(:,:));
        for k1=1:3
            for k2=1:3
                for k3=1:3
                    M3_MC(k1,k2,k3,ts)         = mean((u(:,k1)-mean(u(:,k1))).*(u(:,k2)-mean(u(:,k2))).*(u(:,k3)-mean(u(:,k3))));
                    M3_MC_norm(k1,k2,k3,ts)    = M3_MC(k1,k2,k3,ts)/sqrt(cov_MC(k1,k1,ts)*cov_MC(k2,k2,ts)*cov_MC(k3,k3,ts));

                end
            end
        end
        
        % truth energy from the mean & var
        Ene_MC(1,ts)=0.5*(sum(mean_MC(:,ts).^2)+sum(diag(cov_MC(:,:,ts))));
        Ene_MC(2,ts)=0.5*(mean_MC(1,ts)^2+cov_MC(1,1,ts));
        Ene_MC(3,ts)=0.5*(mean_MC(2,ts)^2+cov_MC(2,2,ts));
        Ene_MC(4,ts)=0.5*(mean_MC(3,ts)^2+cov_MC(3,3,ts));
        
        Ene_dyn(1,ts)=Ene_update(1);
        Ene_dyn(2,ts)=Ene_update(2);
        Ene_dyn(3,ts)=Ene_update(3);
        Ene_dyn(4,ts)=Ene_update(4);
    end
    if mod(ii,ns)==0
        samps(:,:,ic)=samp;
        Ts(ic)=t;
        ic=ic+1;

        display(['MC iter = ',num2str(ii),': ene = ',num2str(Ene_MC(1,ts)), ', ene_dyn = ',num2str(Ene_dyn(1,ts))]);
    end
    
end


