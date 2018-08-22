
%% Multipole expansion linear scattering for single pillars in AlGaAs
% The following script receives from COMSOL the Incident and Scattered
% electric field from a Cylindrical Nanoantenna (centered in (0,0,0) in COMSOL)
% and perform the multipole expansion of
% scattered cross section and extinction cross section. 
% Input from COMSOL:
% - Ex,Ey,Ez interpolated
% - k0 for all considered wavelength (wavevector in hosting medium)

% - Geometry of the Nanoantenna

clear all

c=2.9979*10^8;      %[m/s]
eps0=8.85*10^-12;   %[F/m]
mu0=4*pi*10^(-7);   %[H/m]
E0=1;
nh=1;
eta=sqrt(mu0/(eps0*nh^2));

lmax=3;

model_name='LinearMultipolesSubstrate.mph';

import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar

model = mphload(model_name); %we reload the model to refresh it

linear_study_name='std1';
geometry_name='geom1';
mesh_name='mesh1';
linear_solver_name='sol1';

multipoles_study_name='std2';
multipoles_solver_name='sol2';
% Set Initial Geometry
r_cyl=245:2:295;              %Pillar Radius [nm]
h_cyl=400;              %Pillar Height   [nm]
lambda0=1000;   %Fundamental Frequency Wavelength [nm]

k=2*pi./(lambda0.*10^(-9))*nh;
N=size(lambda0,2); 
omega=k*c;

resolution=2;% (0.5nm) pixel size
    
Nr=length(r_cyl);
Nl=length(lambda0);

QFF=zeros(Nr,Nl);

Csca=zeros(Nr,Nl);
Cse=zeros(Nr,Nl,lmax);
Csm=zeros(Nr,Nl,lmax);

for ir=1:Nr
    
    model.param.set('ra',[num2str(r_cyl(ir)) '[nm]']);   
    model.param.set('rb',[num2str(r_cyl(ir)) '[nm]']);   
    model.param.set('h0',[num2str(h_cyl) '[nm]']);
    
    model.geom(geometry_name).run();
    
    model.mesh(mesh_name).run();
    
    for il=1:Nl
        
    model.param.set('lambda0',[num2str(lambda0(il)) '[nm]']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do Not Change %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Solving linear problem (%d/%d)...\n',Nl*(ir-1)+il,Nl*Nr);
    
 %% Run the COMSOL simulation - Linear
    model.study(linear_study_name).run(); 

 %% Extract Results    
    fprintf('Multipole decomposition (%d)...\n', il)

    QFF(ir,il)=mphglobal(model,'QFF'); % Computing Scattering Efficiency FF
    model.study(multipoles_study_name).run();    
    %ae=mpheval(model,'ae');
    
    model.result.numerical('gev1').setResult;
    tabl=mphtable(model,'tbl9');
    %[ae] = mphglobal(model,'ae','Dataset','dset2','outersolnum','','complexout','on');
    Ce(il,:,:)=reshape(tabl.data(:,3)',[7,3])';
    clear tabl
    
    model.result.numerical('gev2').setResult;
    tabl=mphtable(model,'tbl10');
    %[ae] = mphglobal(model,'ae','Dataset','dset2','outersolnum','','complexout','on');
    Cm(il,:,:)=reshape(tabl.data(:,3)',[7,3])';
    clear tabl
    
    end %End for cycle on wavelength
      
    sigma_geom=pi*r_cyl(ir).^2*10^(-18); 
    
    for i=1:3 %% Cycle on l index
        Cse(ir,:,i)=sum(Ce(:,i,:),3)./sigma_geom;
        Csm(ir,:,i)=sum(Cm(:,i,:),3)./sigma_geom;
        Csca(ir,:)=Csca(ir,:)+Cse(ir,:,i)+Csm(ir,:,i);
    end
end %End for cycle on radius

if Nr==1 && Nl>1 
    str=sprintf('Nanopillar Scattering r:%3.0f nm h:%3.0f nm', r_cyl(1), h_cyl);
    hold on
    grid on
    plot(lambda0',Cse(1,:,1),'r',lambda0',Csm(1,:,1),'b')
    plot(lambda0',Cse(1,:,2),'g',lambda0',Csm(1,:,2),'y')
    plot(lambda0',Cse(1,:,3),'c',lambda0',Csm(1,:,3),'m')
    plot(lambda0',Csca(1,:),'k')
    title(str)
    xlabel('Wavelength [nm]')
    ylabel('Scattering efficiency')
    legend('ED','MD','EQ','MQ','EO','MO','Total')
    hold off
elseif Nr>1 && Nl==1
    str=sprintf('Nanopillar Scattering h:%3.0f nm', h_cyl);
        hold on
    grid on
    plot(r_cyl,Cse(:,1,1),'r',r_cyl,Csm(:,1,1),'b')
    plot(r_cyl,Cse(:,1,2),'g',r_cyl,Csm(:,1,2),'y')
    plot(r_cyl,Cse(:,1,3),'c',r_cyl,Csm(:,1,3),'m')
    plot(r_cyl,Csca(:,1),'k')
    title(str)
    xlabel('Radius [nm]')
    ylabel('Scattering efficiency')
    legend('ED','MD','EQ','MQ','EO','MO','Total')
    hold off
end    
