%% Non linear scattering for single pillars in AlGaAs

clear all

model_name='NonlinearStable.mph';

import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar

model = mphload(model_name); %we reload the model to refresh it

linear_study_name='std1';
nonlinear_study_name='std2';
geometry_name='geom1';
linearmesh_name='mesh1';
nonlinearmesh_name='mesh2';
linear_solver_name='sol1';

dsetFF='dset1';
dsetSH='dset2';

pltflag=0; %(0 no plot - 1 plot)

%% Set Initial Geometry
ra=245:2:295;        %Pillar major Radius [nm]
rb=245:2:295;        %Pillar minor Radius [nm]
h_pil=400;     %Pillar Height   [nm]
lambda0=2000;   %Fundamental Frequency Wavelength [nm]

Nr=length(ra);
Nl=length(lambda0);

QFF=zeros(Nr,Nl);
QSH=zeros(Nr,Nl);
eta=zeros(Nr,Nl);
eff=zeros(Nr,Nl);

for ir=1:Nr
    
    for il=1:Nl
        
    model.param.set('ra',[num2str(ra(ir)) '[nm]']);
    model.param.set('rb',[num2str(rb(ir)) '[nm]']);
    model.param.set('h0',[num2str(h_pil) '[nm]']);
    model.param.set('lambda',[num2str(lambda0(il)) '[nm]']);

    % component of the fields (Ex, Ey, Ez, Hx, Hy, Hz) used to find the pole (Ez in OExpress 2013)
    tested_field_comp='normE';
    % window & resolution to visualize the field on xOz plane
    window=[-1.5*ra 1.5*ra -0.5*h_pil 1.5*h_pil].*10^(-9);% [xmin xmax zmin zmax] (m)
    resolution=0.5*10^(-9);% (0.5nm) pixel size

    model.geom(geometry_name).run();
    model.mesh(linearmesh_name).run();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do Not Change %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xp=window(1):resolution:window(2);% create the x axis 
    zp=window(3):resolution:window(4);% create the z axis
    [Xp,Zp]=meshgrid(xp,zp);
    [Xxp,Yyp]=meshgrid(xp,xp);
    
    if pltflag==1
        XZ=[Xp(:),0*Xp(:),Zp(:)];% create coordinates to which we evaluate the field map
        YZ=[0*Xp(:),Xp(:),Zp(:)];% create coordinates to which we evaluate the field map
        XY=[Xxp(:),Yyp(:),0*Xxp(:)];% create coordinates to which we evaluate the field map
        
        xpq=linspace(window(1),window(2),20);% create the x axis 
        zpq=linspace(window(3),window(4),20);% create the z axis
        [Xpq,Zpq]=meshgrid(xpq,zpq);
        [Xxpq,Yypq]=meshgrid(xpq,xpq);
        XZq=[Xpq(:),0*Xpq(:),Zpq(:)];% create coordinates to which we evaluate the field map
        YZq=[0*Xpq(:),Xpq(:),Zpq(:)];% create coordinates to which we evaluate the field map
        XYq=[Xxpq(:),Yypq(:),0*Xxpq(:)];% create coordinates to which we evaluate the field map
        
    end
    
    fprintf('Solving linear problem (%d/%d)...\n',Nl*(ir-1)+il,Nl*Nr);

    r=0:resolution:ra(ir)*10^(-9);
    theta=0:pi/20:2*pi;
    zeta=0:20*resolution:h_pil*10^(-9);

    [R,T,Zr]=meshgrid(r,theta,zeta);
    Xr=R.*cos(T);
    Yr=R.*sin(T);
    XYZ=[Xr(:),Yr(:),Zr(:)];

 %% Run the COMSOL simulation - Linear
    model.study(linear_study_name).run(); 

 %% Extract Results
    fprintf('Interpolating Results Fundamental Frequency...\n');

    QFF(ir,il)=mphglobal(model,'QFF','Dataset',dsetFF); % Computing Scattering Efficiency FF

    [ExFF]=mphinterp(model,'ewfd.Ex','coord',XYZ','complexout','on');
    [EyFF]=mphinterp(model,'ewfd.Ey','coord',XYZ','complexout','on');
    [EzFF]=mphinterp(model,'ewfd.Ez','coord',XYZ','complexout','on');
    
    if pltflag==1
        
        [total_fieldXZ]=mphinterp(model,['ewfd.' tested_field_comp],'coord',XZ');
        [total_fieldYZ]=mphinterp(model,['ewfd.' tested_field_comp],'coord',YZ');
        [total_fieldXY]=mphinterp(model,['ewfd.' tested_field_comp],'coord',XY');

        % XZ plane interpolation
        
        [Ex]=mphinterp(model,'ewfd.Ex','coord',XZq');
        [Ez]=mphinterp(model,'ewfd.Ez','coord',XZq');

        ExFF_XZ=reshape(Ex,size(zpq,2),size(xpq,2));
        EzFF_XZ=reshape(Ez,size(zpq,2),size(xpq,2));

        clear Ex Ez

        ExFF_XZn=ExFF_XZ./(sqrt(ExFF_XZ.^2+EzFF_XZ.^2));
        EzFF_XZn=EzFF_XZ./(sqrt(ExFF_XZ.^2+EzFF_XZ.^2));

        % YZ plane interpolation
        
        [Ey]=mphinterp(model,'ewfd.Ey','coord',YZq');
        [Ez]=mphinterp(model,'ewfd.Ez','coord',YZq');

        EyFF_YZ=reshape(Ey,size(zpq,2),size(xpq,2));
        EzFF_YZ=reshape(Ez,size(zpq,2),size(xpq,2));

        clear Ey Ez

        EyFF_YZn=EyFF_YZ./(sqrt(EyFF_YZ.^2+EzFF_YZ.^2));
        EzFF_YZn=EzFF_YZ./(sqrt(EyFF_YZ.^2+EzFF_YZ.^2));

        % XY plane interpolation
        
        [Ex]=mphinterp(model,'ewfd.Ex','coord',XYq');
        [Ey]=mphinterp(model,'ewfd.Ey','coord',XYq');

        ExFF_XY=reshape(Ex,size(xpq,2),size(xpq,2));
        EyFF_XY=reshape(Ey,size(xpq,2),size(xpq,2));

        clear Ex Ey

        ExFF_XYn=ExFF_XY./(sqrt(ExFF_XY.^2+EyFF_XY.^2));
        EyFF_XYn=EyFF_XY./(sqrt(ExFF_XY.^2+EyFF_XY.^2));
        
      %% Plot Field Distribution
        figure (1)
        field_normXZ=reshape(total_fieldXZ,size(zp,2),size(xp,2))./max(abs(total_fieldXZ));

        hold on

        quiver(Xpq,Zpq,ExFF_XZn,EzFF_XZn)
        axis equal
        shading interp
        rectangle('Position',[-ra*10^(-9) 0 2*ra*10^(-9) h_pil*10^(-9)]','LineStyle','--')
        h=surf(Xp,Zp,field_normXZ,'LineStyle','none');
        colorbar
        view(2)
        hz = get(h,'ZData');
        set(h,'ZData',hz-max(max(field_normXZ))-1)  

        xlim([window(1) window(2)])
        ylim([window(3) window(4)])
        title('Electric Field - Norm distribution - XZ plane (y=0)')
        xlabel(['X position [nm] - radius ' num2str(ra) ' nm'])
        ylabel(['Z position [nm] - height ' num2str(h_pil) ' nm'])
        hold off

        figure (2)
        field_normYZ=reshape(total_fieldYZ,size(zp,2),size(xp,2))./max(abs(total_fieldYZ));
        
        hold on
        quiver(Xpq,Zpq,EyFF_YZn,EzFF_YZn)
        axis equal
        rectangle('Position',[-rb*10^(-9) 0 2*rb*10^(-9) h_pil*10^(-9)]','EdgeColor','k')
        g=surf(Xp,Zp,field_normYZ,'LineStyle','none');
        gz = get(g,'ZData');
        set(g,'ZData',gz-max(max(field_normYZ))-1) 
        xlim([window(1) window(2)])
        ylim([window(3) window(4)])
        colorbar
        view(2)
        shading interp


        title('Electric Field - Norm distribution - YZ plane (x=0)')
        xlabel(['Y position [nm] - radius ' num2str(rb) ' nm'])
        ylabel(['Z position [nm] - height ' num2str(h_pil) ' nm'])
        hold off
        
        figure (3)
        field_normXY=reshape(total_fieldXY,size(xp,2),size(xp,2))./max(abs(total_fieldXY));
        
        hold on
        quiver(Xxpq,Yypq,ExFF_XYn,EyFF_XYn)
        axis equal
        rectangle('Position',[-ra*10^(-9) -rb*10^(-9) 2*ra*10^(-9) 2*rb*10^(-9)]','Curvature',[1,1],'EdgeColor','k')
        i=surf(Xxp,Yyp,field_normXY,'LineStyle','none');
        iz = get(i,'ZData');
        set(i,'ZData',iz-max(max(field_normXY))-1) 
        xlim([window(1) window(2)])
        ylim([window(1) window(2)])
        colorbar
        view(2)
        shading interp

        title('Electric Field - Norm distribution - XY plane (z=0)')
        xlabel(['X position [nm] - radius ' num2str(ra) ' nm'])
        ylabel(['Y position [nm] - height ' num2str(h_pil) ' nm'])
        hold off
    
    end % End if plot cycle
    
 %% Run the COMSOL simulation - NonLinear
    fprintf('Solving nonlinear problem...\n');
    model.mesh(nonlinearmesh_name).run();
    model.study(nonlinear_study_name).run();

 %% Extract Results
    fprintf('Interpolating Results Second Harmonic...\n');
 
    QSH(ir,il)=mphglobal(model,'QSH','Dataset',dsetSH); %Computing SH Scattering Efficiency
    
    if pltflag==1
        
        [SH_fieldXZ]=mphinterp(model,['ewfd2.' tested_field_comp],'coord',XZ','Dataset',dsetSH);
        [SH_fieldYZ]=mphinterp(model,['ewfd2.' tested_field_comp],'coord',YZ','Dataset',dsetSH);
        [SH_fieldXY]=mphinterp(model,['ewfd2.' tested_field_comp],'coord',XY','Dataset',dsetSH);
        
        [Ex]=mphinterp(model,'ewfd2.Ex','coord',XZq','Dataset',dsetSH);
        [Ez]=mphinterp(model,'ewfd2.Ez','coord',XZq','Dataset',dsetSH);

        ExSH_XZ=reshape(Ex,size(zpq,2),size(xpq,2));
        EzSH_XZ=reshape(Ez,size(zpq,2),size(xpq,2));

        clear Ex Ez

        ExSH_XZn=ExSH_XZ./(sqrt(ExSH_XZ.^2+EzSH_XZ.^2));
        EzSH_XZn=EzSH_XZ./(sqrt(ExSH_XZ.^2+EzSH_XZ.^2));

        [Ey]=mphinterp(model,'ewfd2.Ey','coord',YZq','Dataset',dsetSH);
        [Ez]=mphinterp(model,'ewfd2.Ez','coord',YZq','Dataset',dsetSH);

        EySH_YZ=reshape(Ey,size(zpq,2),size(xpq,2));
        EzSH_YZ=reshape(Ez,size(zpq,2),size(xpq,2));

        clear Ey Ez

        EySH_YZn=EySH_YZ./(sqrt(EySH_YZ.^2+EzSH_YZ.^2));
        EzSH_YZn=EzSH_YZ./(sqrt(EySH_YZ.^2+EzSH_YZ.^2));
  
        % XY plane interpolation
        
        [Ex]=mphinterp(model,'ewfd2.Ex','coord',XYq','Dataset',dsetSH);
        [Ey]=mphinterp(model,'ewfd2.Ey','coord',XYq','Dataset',dsetSH);

        ExSH_XY=reshape(Ex,size(xpq,2),size(xpq,2));
        EySH_XY=reshape(Ey,size(xpq,2),size(xpq,2));
        
        clear Ex Ey

        ExSH_XYn=ExSH_XY./(sqrt(ExSH_XY.^2+EySH_XY.^2));
        EySH_XYn=EySH_XY./(sqrt(ExSH_XY.^2+EySH_XY.^2));
        
    end
    
 %   [ExSH]=mphinterp(model,'ewfd2.Ex','coord',XYZ','complexout','on','Dataset',dsetSH);
 %   [EySH]=mphinterp(model,'ewfd2.Ey','coord',XYZ','complexout','on','Dataset',dsetSH);
 %   [EzSH]=mphinterp(model,'ewfd2.Ez','coord',XYZ','complexout','on','Dataset',dsetSH);

  %  eta(ir,il)=abs(sum(sum(sum(ExFF.*EzFF.*conj(EySH)))))/(sqrt(sum(sum(sum((abs(ExFF.*EzFF)).^2))))+sqrt(sum(sum(sum(abs((EySH).^2))))));

   % eff(ir,il)=eta(ir,il).*QFF(ir,il).^2.*QSH(ir,il);

    if pltflag==1
 %% Plot Field Distribution
        figure (4)
        fieldSH_normXZ=reshape(SH_fieldXZ,size(zp,2),size(xp,2))./max(abs(SH_fieldXZ));

        hold on
        quiver(Xpq,Zpq,ExSH_XZn,EzSH_XZn)
        axis equal
        shading interp
        rectangle('Position',[-ra*10^(-9) 0 2*ra*10^(-9) h_pil*10^(-9)]','LineStyle','--')
        h=surf(Xp,Zp,fieldSH_normXZ,'LineStyle','none');
        colorbar
        view(2)
        hz = get(h,'ZData');
        set(h,'ZData',hz-10)  

        xlim([window(1) window(2)])
        ylim([window(3) window(4)])
        title('Electric Field SH - Norm distribution - XZ plane (y=0)')
        xlabel(['X position [nm] - radius ' num2str(ra) ' nm'])
        ylabel(['Z position [nm] - height ' num2str(h_pil) ' nm'])
        hold off

        figure (5)
        fieldSH_normYZ=reshape(SH_fieldYZ,size(zp,2),size(xp,2))./max(abs(SH_fieldYZ));

        hold on
        quiver(Xpq,Zpq,EySH_YZn,EzSH_YZn)
        axis equal
        rectangle('Position',[-rb*10^(-9) 0 2*rb*10^(-9) h_pil*10^(-9)]','EdgeColor','k')
        g=surf(Xp,Zp,fieldSH_normYZ,'LineStyle','none');
        gz = get(g,'ZData');
        set(g,'ZData',gz-10) 
        xlim([window(1) window(2)])
        ylim([window(3) window(4)])
        colorbar
        view(2)
        shading interp

        title('Electric Field SH - Norm distribution - YZ plane (x=0)')
        xlabel(['Y position [nm] - radius ' num2str(rb) ' nm'])
        ylabel(['Z position [nm] - height ' num2str(h_pil) ' nm'])
        hold off
        
        figure (6)
        fieldSH_normXY=reshape(SH_fieldXY,size(xp,2),size(xp,2))./max(abs(SH_fieldXY));
        
        hold on
        quiver(Xxpq,Yypq,ExSH_XYn,EySH_XYn)
        axis equal
        rectangle('Position',[-ra*10^(-9) -rb*10^(-9) 2*ra*10^(-9) 2*rb*10^(-9)]','Curvature',[1,1],'EdgeColor','k')
        i=surf(Xxp,Yyp,fieldSH_normXY,'LineStyle','none');
        iz = get(i,'ZData');
        set(i,'ZData',iz-max(max(fieldSH_normXY))-1) 
        xlim([window(1) window(2)])
        ylim([window(1) window(2)])
        colorbar
        view(2)
        shading interp

        title('Electric Field - Norm distribution - XY plane (z=0)')
        xlabel(['X position [nm] - radius ' num2str(ra) ' nm'])
        ylabel(['Y position [nm] - height ' num2str(h_pil) ' nm'])
        hold off
    end
    
    end %End for cycle on wavelength
end % End for cycle on radius 
