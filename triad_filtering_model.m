function [zt,TT,z_samps,Ts, mean_u,cov_u, mom_u] = triad_filtering_model(params, m0,var0, obs_m,obs_v)
MC=params.MC; Dt=params.Dt; tstep=params.tstep; Nt=params.Nt;
L=params.L; G=params.G; B=params.B; 
skip=params.skip;

% Initial conditions
zt(:,1)=sqrt(var0(1))*randn(MC,1);
zt(:,2)=sqrt(var0(2))*randn(MC,1);
zt(:,3)=sqrt(var0(3))*randn(MC,1);

% model coeffs
Bsis=eye(3);
for i=1:3
    for j=1:3
        coupling_strch(:,i+(j-1)*size(Bsis,2)) =  Buu2(B,Bsis(:,i),Bsis(:,j));
    end
end
params.coupling_strch=coupling_strch;


ts=1;
Ut=m0';
Rt=diag(var0);
zsamp=zt-mean(zt);
TT(ts)=0;
mean_u(:,ts)=Ut;
cov_u(:,:,ts)=Rt;
% save samples
ic=1; ns=floor(.5/Dt);
for ii=1:Nt
    
    t=Dt*(ii);     % time

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % forcing terms
    FF=params.tmM(ii,:);
    SS=params.SS+params.tmS(ii)^2*(params.SSt-params.SS);

    % observations
    dU_obs=(obs_m(:,ii+1)-obs_m(:,ii)) / Dt; 
    dR_obs=(obs_v(:,:,ii+1)-obs_v(:,:,ii)) / Dt;
    % RK4
    % step 1
    [dz1,dU1,dR1] = stocstat_update_step(zsamp,Ut,Rt, FF,SS, params);

    z1 = zsamp+.5*Dt*dz1;
    U1 = Ut   +.5*Dt*dU1;
    R1 = Rt   +.5*Dt*dR1;
    % step 2
    [dz2,dU2,dR2] = stocstat_update_step(z1,U1,R1, FF,SS, params);

    z2 = zsamp+.5*Dt*dz2;
    U2 = Ut   +.5*Dt*dU2;
    R2 = Rt   +.5*Dt*dR2;
    % step 3
    [dz3,dU3,dR3] = stocstat_update_step(z2,U2,R2, FF,SS, params);

    z3 = zsamp+1*Dt*dz3;
    U3 = Ut   +1*Dt*dU3;
    R3 = Rt   +1*Dt*dR3;
    % step 4
    [dz4,dU4,dR4] = stocstat_update_step(z3,U3,R3, FF,SS, params);

    Ut_pre=Ut; Rt_pre=Rt; 
    % final update 
    Ut = Ut + Dt*(dU1/6+dU2/3+dU3/3+dU4/6);
    Rt = Rt + Dt*(dR1/6+dR2/3+dR3/3+dR4/6);
    zsamp = zsamp+Dt*(dz1/6+dz2/3+dz3/3+dz4/6)+sqrt(Dt)*randn(MC,3)*SS;

    z_pre=zsamp;
    % observation update
    if ~mod(ii,skip)
        params.scale_m=1.; params.scale_v=1.;

        dzf1 = filter_update_step(z1,U1,R1, dU_obs,dR_obs,FF,SS, params);
        dzf2 = filter_update_step(z2,U2,R2, dU_obs,dR_obs,FF,SS, params);
        dzf3 = filter_update_step(z3,U3,R3, dU_obs,dR_obs,FF,SS, params);
        dzf4 = filter_update_step(zsamp,Ut,Rt, dU_obs,dR_obs,FF,SS, params);
        zsamp = zsamp+Dt*(dzf1/6+dzf2/3+dzf3/3+dzf4/6);

    end
    

    % normalize samples & correction of the mean
    zsamp=zsamp-mean(zsamp);

    
    % save solution
    if ~mod(ii,tstep)
        ts=ts+1;
        TT(ts)=t;

        zt=zsamp;
        mean_u(:,ts)=Ut;
        cov_u(:,:,ts)=Rt;
        % cov_z(:,:,ts)=cov(zsamp);

        % get higher mom
        zp=zsamp-mean(zsamp);
        Hm=Buu(B,zp,zp);
        Hv=repmat(zp,[1 1 3]).*reshape(Hm,[MC 1 3]);
        mom_u(:,:,ts) = mean(Hv);
    end
    if mod(ii,ns)==0
        z_samps(:,:,ic)=zsamp;
        Ts(ic)=t;
        ic=ic+1;

        display(['High-order Filter iter = ',num2str(ii),': mean = ',num2str(sum(mean_u(:,ts))), ', var = ',num2str(trace(cov_u(:,:,ts)))]);
    end
end