%% Non linear scattering for single pillars in AlGaAs

clear all
tic
format short e
addpath(genpath(pwd));

model = Model; % It initialize the object model with Matlab class Model
model.name='NanodiskScattering.mph'; % COMSOL mph file name
model = model.load(); % It loads COMSOL with Matlab LiveLink and mph file

model.linear_study_name='std1';
model.geometry_name='geom1';
model.linearmesh_name='mesh1';
 
model.dsetFF='dset1';
model.res_domain = [9];
%model.logfile = 'Log_SweepWavelength.txt';

params = struct('ra', 200e-9, 'rb', 200e-9, 'h0', 400e-9,...
                'lambdaFF', 1550e-9,...
                'I0', 1e9,...
                'n_c', 1.0, 'n_s', 1.0,...
                'theta_pil', 0,...% Crystalline axis and ellipsitcal base major sexiaxis orientation wrt x-axis respectively 
                'theta_Eb', 0, 'phi_Eb', 0, 'pol_Eb', 90); %Respectively k tilt wrt to z axis, k tilt wrt x axis, E tilt wrt p polarization direction


model = model.reset_params(params);

l0 = 1550e-9;
l1 = 1555e-9;
dl = 5e-9;
l = l0:dl:l1;

fprintf(['Loaded Model: ' model.name '\n'...
        'Starting model parameters:\n'...
        'radius: ' num2str(model.params.ra) '[m]\n'...
        'height: ' num2str(model.params.h0) '[m]\n'...
        'lambda: ' num2str(model.params.lambdaFF) '[m]\n'...
        'I0: ' num2str(model.params.I0) '[W/cm^2]\n'...
        'Incident Polarization to 100 axis: ' num2str(model.params.pol_Eb) '[deg]\n'])

fprintf('Choose folder destination for log file:\n')
selpath = uigetdir(path,'Choose folder destination for log file');
logID = fopen([selpath,'\',model.logfile],'a+');
model.init_logfile(logID);

lmax = 3;

Csca = zeros(length(l),1);
Cse = zeros(length(l),lmax);
Csm = zeros(length(l),lmax);

figure()
for il=1:length(l)
    
    params.lambdaFF = l(il);
    model = model.reset_params(params);
    
    fprintf(['\r lambda= ',num2str(l(il)*1e9)])

    model.solve();
    [Csca(il)] = sca_efficiency(model);
    fprintf(logID,['\r',num2str(l(il)*1e9),'\t',num2str(Csca(il)),'\t']);
    [Csca2(il),Cse(il,:),Csm(il,:)] = multipolesFF(model,lmax);
    clf('reset')
    plot(l,Csca,'ks'),hold on, plot(l,Cse,'o'), plot(l,Csm,'^');
    legend('Total','ED','EQ','EO','MD','MQ','MO');
    drawnow
    
end


%% NearField plot -- EXAMPLE
figure()
model.plot_nearfieldXZ('real(emw.normE)',model.dsetFF);

%% FarField Computation -- EXAMPLE
[angles_up,angles_dn] = model.compute_farfieldRETOP (model.params.lambdaFF/2,'emw.',[0,pi/2,40],[0,pi,2],model.dsetFF);
figure()
polarplot(angles_dn.theta.*(1-2*sin((angles_dn.delta)/2))-pi/2,abs(angles_dn.F),'b*')
hold on,
polarplot(angles_up.theta.*(1-2*sin((angles_up.delta)/2))+pi/2,abs(angles_up.F),'b*')
set(gca,'FontName','Times','FontSize',16)


[angles_up,angles_dn] = model.compute_farfieldRETOP (model.params.lambdaFF/2,'emw.',[0,pi/2,40],[-pi,pi,50],model.dsetFF);
figure()
Fh=reshape(angles_up.F,size(angles_up.ug));              % the flux at each solid angle, upper space ("h"-->upper)
Fb=reshape(angles_dn.F,size(angles_dn.ug));              % the flux at each solid angle, lower space ("b"-->lower)
surf(angles_up.ug.*Fh,angles_up.vg.*Fh,angles_up.wg.*Fh,Fh,'LineStyle','none');
hold on
surf(angles_dn.ug.*Fb,angles_dn.vg.*Fb,-angles_dn.wg.*Fb,Fb,'LineStyle','none');
set(gca,'projection','perspective');
axis tight;axis equal;xlabel('x');ylabel('y');
shading interp
colormap(plasma())
colorbar
%caxis([0 24e12])
lighting PHONG;
% light('position',[-1,0,1],'color','r');
light('position',[5,1,1],'color','w');
light('position',[5,1,1],'color','y');
view(60,20)
axis off


time = toc;
disp(['Execution time: ' num2str(floor(time/3600)) 'h ' num2str(floor(time/60-floor(time/3600)*60)) 'm ' num2str(floor(time-floor(time/60)*60)) 's'])

figure()
plot_Csca(l,0,Csca,model.params,'lambda')

% I = logspace(7,10,10);
% [Csca,CscaSH] = model.sweep_intensity(I);
