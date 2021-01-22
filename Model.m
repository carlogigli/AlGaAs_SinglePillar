classdef Model
    % COMSOL simulation object
    % This class sets all the parameters of the model
    % and executes functions to run simulations
    
    properties
        
    name
    linear_study_name = 'std1'
    geometry_name = 'geom1'
    linearmesh_name = 'mesh1'
    mph
    dsetFF = 'dset1'
    res_domain
    logfile = ['Log_', datestr(date),'.txt']
    params = struct('ra', 0, 'rb', 0, 'h0', 0,...
                    'lambdaFF', 0,...
                    'I0', 0,...
                    'n_c', 0, 'n_s', 0,...
                    'theta_pil', 0,... % Crystalline axis and ellipsitcal base major sexiaxis orientation wrt x-axis respectively 
                    'theta_Eb', 0, 'phi_Eb', 0, 'pol_Eb', 0); %Resctively k tilt wrt to z axis, k tilt wrt x axis, E tilt wrt p polarization direction
    end
    
    methods
        
        function self = Model(name)
            if nargin > 0
                self.name = name;
            end
        end
        
        function self = load(self)
            import com.comsol.model.* 		% Load APIs
            import com.comsol.model.util.*	% idem
            ModelUtil.showProgress(true) ;  % Progress bar
            self.mph = mphload(self.name);
            
            %%% Read all the initial parameters in the mph file
            
            self.params.ra = self.mph.param.evaluate('ra');
            self.params.rb = self.mph.param.evaluate('rb');
            self.params.h0 = self.mph.param.evaluate('h0');
            self.params.theta_pil = self.mph.param.evaluate('theta_pil');
            self.params.I0 = self.mph.param.evaluate('I0');
            self.params.lambdaFF = self.mph.param.evaluate('lambdaFF');
            self.params.n_c = self.mph.param.evaluate('n_c');
            self.params.n_s = self.mph.param.evaluate('n_s');
            self.params.theta_Eb = self.mph.param.evaluate('theta_Eb');
            self.params.phi_Eb = self.mph.param.evaluate('phi_Eb');
            self.params.pol_Eb = self.mph.param.evaluate('pol_Eb');
            
        end

        function self = init_logfile(self,logID)
            
             fprintf(logID,['New simulation from COMSOL file:', self.name,...
             '\rDate: ',datestr(datetime),...
             '\rParameters:\r',...
             'ra = ',num2str(self.params.ra*1e9), ' [nm]\r',...
             'rb = ',num2str(self.params.rb*1e9), ' [nm]\r',...
             'h0 = ',num2str(self.params.h0*1e9), ' [nm]\r',...
             'lambdaFF = ',num2str(self.params.lambdaFF*1e9), ' [nm]\r',...
             'I0 = ',num2str(self.params.I0), ' [W/cm^2]\r',...
             'n_c = ',num2str(self.params.n_c),'\r',...
             'n_s = ',num2str(self.params.n_s),'\r',...
             'theta_Eb = ',num2str(self.params.theta_Eb),'\r',...
             'phi_Eb = ',num2str(self.params.phi_Eb),'\r',...
             'pol_Eb = ',num2str(self.params.pol_Eb),'\r',...
             'theta_pil = ',num2str(self.params.theta_pil),'\r']);
            
        end
        
        function self = reset_params(self, params)
            
            self.params.ra = params.ra;
            self.params.rb = params.rb;
            self.params.h0 = params.h0;
            self.params.I0 = params.I0;
            self.params.lambdaFF = params.lambdaFF;
            self.params.n_c = params.n_c;
            self.params.n_s = params.n_s;
            self.params.theta_pil = params.theta_pil;
            self.params.theta_Eb = params.theta_Eb;
            self.params.pol_Eb = params.pol_Eb;
            self.params.phi_Eb = params.phi_Eb;
            
            self.mph.param.set('ra',[num2str(params.ra) '[m]']);
            self.mph.param.set('rb',[num2str(params.rb) '[m]']);
            self.mph.param.set('h0',[num2str(params.h0) '[m]']);
            self.mph.param.set('I0',[num2str(params.I0) '[W/cm^2]']);
            self.mph.param.set('lambdaFF',[num2str(params.lambdaFF) '[m]']);
            self.mph.param.set('n_c',num2str(params.n_c));
            self.mph.param.set('n_s',num2str(params.n_s));
            self.mph.param.set('theta_pil',[num2str(params.theta_pil) '[deg]']);
            self.mph.param.set('phi_Eb',[num2str(params.phi_Eb) '[deg]']);
            self.mph.param.set('theta_Eb',[num2str(params.theta_Eb) '[deg]']);
            self.mph.param.set('pol_Eb',[num2str(params.pol_Eb) '[deg]']);
        
        end
            
        function [Csca] = sweep_radius(self, ra, rb)
            
            if nargin ~= 3
                error('Error input arguments calling function sweep_radius!!');
            else
                Csca=zeros(1,length(ra));
                                
                logID = fopen(self.logfile,'a+');
                fprintf(logID,['Run: sweep_radius ', datestr(datetime), ' \r']);
                self.init_logfile(logID)
                
                for ir=1:length(ra)
                    
                    self.params.ra = ra(ir);
                    self.params.rb = rb(ir);

                    self = self.reset_params(self.params);

                    fprintf('Solving radius sweep r: %.3e[m] (%d/%d)...\n',self.mph.param.evaluate('ra'),ir,length(ra));
                       
                    solve(self);
                    [Csca(ir)] = sca_efficiency(self);
                end
                
            end
        end
            
        function [Csca] = sweep_lambda(self, l)
            
            if nargin ~= 2
                error('Error input arguments calling function sweep_lambda!!');
            else              
                Csca=zeros(1,length(l));
                
                for il=1:length(l)
                    
                    self.params.lambdaFF = l(il);
                    self = self.reset_params(self.params);

                    fprintf('Solving wavelength sweep lambda: %.3e[m] (%d/%d)...\n',self.mph.param.evaluate('lambdaFF'),il,length(l));
                    solve(self);
                    [Csca(il)] = sca_efficiency(self);
                end
            end            
        end
        
          
        function  self = solve(self)
                    
                    self.mph.geom(self.geometry_name).run();
                    self.mph.mesh(self.linearmesh_name).run();
                        
                    self.mph.study(self.linear_study_name).run();
                     
        end
        
        function [Csca] = sca_efficiency(self)
            
                    Csca=mphglobal(self.mph,'C_sca','Dataset',self.dsetFF); % Computing Scattering Efficiency FF
                                        
        end
        
        function [Xp,Zp,field_XZ,hAx] = plot_nearfieldXZ(self,field_component,dataset)
        
            z0 = self.mph.param.evaluate('z0');
            window=[-1.5*self.params.rb 1.5*self.params.rb (self.params.h0/2+z0)-self.params.h0 (self.params.h0/2+z0)+self.params.h0];% [xmin xmax zmin zmax] (m)
            resolution=5*10^(-9);% (0.5nm) pixel size
            xp=window(1):resolution:window(2);% create the x axis 
            zp=window(3):resolution:window(4);% create the z axis
            [Xp,Zp]=meshgrid(xp,zp);
            XZ=[Xp(:),0*Xp(:),Zp(:)];% create coordinates to which we evaluate the field map
            [total_fieldXZ]=mphinterp(self.mph,field_component,'coord',XZ','dataset',dataset);
            field_XZ=reshape(total_fieldXZ,size(zp,2),size(xp,2));
            hAx = pcolor(Xp,Zp,field_XZ);
            hAx.LineStyle='none';
            rectangle('Position',[-self.params.ra z0 2*self.params.ra self.params.h0]','LineStyle','--')
            %hAx = surf(Xp,Zp,field_XZ,'LineStyle','none');
            axis equal
            shading interp
            colormap(RWB())
            colorbar
            view(2)
            xlim([window(1) window(2)])
            ylim([window(3) window(4)])
            title([field_component,' - XZ plane (y=0)'])
            xlabel(['X position [nm] - radius ' num2str(self.params.ra*1e9) ' nm'])
            ylabel(['Z position [nm] - height ' num2str(self.params.h0*1e9) ' nm'])
            set(gca,'FontName','Times','FontSize',14);
            box on
            
        end
        
        function [Yp,Zp,field_YZ,hAx] = plot_nearfieldYZ(self,field_component,dataset)
        
            z0 = self.mph.param.evaluate('z0');
            window=[-1.5*self.params.rb 1.5*self.params.rb (self.params.h0/2+z0)-self.params.h0 (self.params.h0/2+z0)+self.params.h0];% [xmin xmax zmin zmax] (m)
            resolution=5*10^(-9);% (0.5nm) pixel size
            yp=window(1):resolution:window(2);% create the x axis 
            zp=window(3):resolution:window(4);% create the z axis
            [Yp,Zp]=meshgrid(yp,zp);
            YZ=[0*Yp(:),Yp(:),Zp(:)];% create coordinates to which we evaluate the field map
            [total_fieldYZ]=mphinterp(self.mph,field_component,'coord',YZ','dataset',dataset);
            field_YZ=reshape(total_fieldYZ,size(zp,2),size(yp,2));
            %hAx = surf(Yp,Zp,field_YZ,'LineStyle','none');
            hAx = pcolor(Yp,Zp,field_YZ);
            hAx.LineStyle='none';
            rectangle('Position',[-self.params.rb z0 2*self.params.rb self.params.h0]','LineStyle','--')
            axis equal
            shading interp
            colormap(RWB())
            colorbar
            view(2)
            xlim([window(1) window(2)])
            ylim([window(3) window(4)])
            title([field_component,' - YZ plane (x=0)'])
            xlabel(['Y position [nm] - radius ' num2str(self.params.rb*1e9) ' nm'])
            ylabel(['Z position [nm] - height ' num2str(self.params.h0*1e9) ' nm'])
            set(gca,'FontName','Times','FontSize',14);
            box on
            
        end
              
        function EH=extract_comsol_field(self,x,y,z,physics,dataset)
        %%%%% extract_comsol_field.m - RETOP subroutine %%%%%
        % this function extracts the (scattered or total) field from COMSOL multiphysics
        % 'model' : the COMSOL model that provides the field on the box
        % x,y,z are vectors that define the x, y and z coordinates of the points for which we extract the COMSOL field
        % 'physics': COMSOL physics name to extract near field. It can
        % assume generally three values:    'ewfd.' -> total field at FF
        %                                   'ewfd.rel' -> scatteredfield at FF
        %                                   'ewfd2.' -> total field at SH
        % Dataset has to be passed coherently with physics choice to this
        % function.
        
            coord=[x(:),y(:),z(:)];
            [Ex,Ey,Ez,Hx,Hy,Hz]=mphinterp(self.mph,{[physics,'Ex'], [physics,'Ey'],[physics,'Ez'],[physics,'Hx'],[physics,'Hy'],[physics,'Hz']},'coord',coord.','Complexout','on','Dataset',dataset);

            EH=([Ex;Ey;Ez;Hx;Hy;Hz].');
        end
        
        function [angles_up,angles_dn] = compute_farfieldRETOP (self,wavelength,physics,theta_arr,phi_arr,dataset)
            %% RETOP Subroutine to compute farfield
            % It implements RETOP package to compute the near to far field
            % transformation from COMSOL solution 
            % Inputs:   - wavelength working wavelength of COMSOL solution 
            %           - physics COMSOL physics name to extract near field.
            %           - theta_arr: 3 elements array containing
            %             [min_theta,max_theta,Npoints] expressing minimum
            %             maximum and number of points for the polar angle.
            %             theta can span from 0 to pi/2
            %           - phi_arr: 3 elements array containing
            %             [min_phi,max_phi,Npoints] expressing minimum
            %             maximum and number of points for the azimuthal angle
            %             phi can span from -pi to pi
            %           - dataset: COMSOL dataset for nearfield extraction
            
            wavenumber=2*pi/wavelength;
            refractive_indices=[self.params.n_c,self.params.n_s];       % the refractive indices of the stratified medium
            z_layers=self.mph.param.evaluate('z0');                           %  z coordinate of each discontinuity (unit: meter)

            center_x=0; center_y=0; center_z=self.mph.param.evaluate('zbox');     % center of the box (in RETOP unit : meter)
            box_center=[center_x,center_y,center_z];        
            Lx=self.mph.param.evaluate('Lbox');Ly=Lx;Lz=Lx;
            box_size=[Lx,Ly,Lz];                    % Box Size (in RETOP unit : meter)--> the box has the same size as the one defined in COMSOL model
            option_i=1;                             % 'option_i=1' : input field has convention exp(i*omega*t)
                                        % 'option_i=-1': input field has convention exp(-i*omega*t)
            N_gauss=[8,8;8,8;8,8];                  % A parameter for generating the sampling locations on the "box" (normally, the present setting is fine)
            % define the (theta --> polar angle, phi --> azimuthal angle) (for details, see the user guide)

            [teta,wteta]=retgauss(theta_arr(1),theta_arr(2),theta_arr(3),1); % polar angle
            % [phi,wphi]=retgauss(0,pi,2,1);   % azimuthal angle
            if phi_arr(3)==2; phi = [phi_arr(1), phi_arr(2)]; else [phi,wphi]=retgauss(phi_arr(1), phi_arr(2),phi_arr(3),1);end % [0 pi]: 2D XZ plane; [pi/2 3*pi/2]: 2D YZ plane;
            [Teta,Phi]=ndgrid(teta,phi);        % generate the grid
            u=sin(Teta).*cos(Phi);              % weight in x-coordinate
            v=sin(Teta).*sin(Phi);              % weight in y-coordinate
            w=cos(Teta);                        % weight in z-coordinate

            cal_champ=@(x,y,z)self.extract_comsol_field(x,y,z,physics,dataset); 

            % Initialize RETOP by inputting material and geometrical parameters
            init=retop(cal_champ,wavenumber,refractive_indices,z_layers,box_center,box_size,N_gauss,struct('option_cal_champ',2,'option_i',option_i));
            
            %%% Planewave decomposition in the upper free-space (if the corresponding refractive index is real)
            if imag(refractive_indices(1))==0
                uh=wavenumber*refractive_indices(1)*u;vh=wavenumber*refractive_indices(1)*v;  %  wavenumber kx=uh, ky=vh, refractive_indices(1): refractive index of the upper free-space
                direction=1;                                           %  Attention: direction=1 (upper space, +z)
                angles_up=retop(init,uh,vh,direction); %  planewave decomposition for the upper space
            % the option "struct('test',1)" allows displaying the field on the box in order to check the quality of sampling
            end
            angles_up.ug = u;angles_up.vg = v;angles_up.wg = w;
            %%% Planewave decomposition in the lower free-space (if the corresponding refractive index is real)
            if imag(refractive_indices(end))==0
                ub=wavenumber*refractive_indices(end)*u;vb=wavenumber*refractive_indices(end)*v;  %  wavenumber kx=ub, ky=vb, refractive_indices(end): refractive index of the lower free-space
                direction= -1;                                   %  Attention: direction= -1 (lower space, -z)
                angles_dn=retop(init,ub,vb,direction);            %  planewave decomposition for the lower space
            end
            angles_dn.ug = u;angles_dn.vg = v;angles_dn.wg = w;
            
        end
              
        function [Csca,Cse,Csm] = multipolesFF(self,lmax)
            %% Sub routine to automatically compute multipoles at FF
            % It call the sub routines "extract_field" to evaluate the total
            % field @FF and "compute_multipoles" to compute the multipoles
            % in spherical coordinates. lmax define the maximum angular
            % moment for the decomposition (typical lmax=3 to stop at octupoles)
            [field,mesh] = extract_field(self,'emw.',self.dsetFF);
            [Csca,Cse,Csm] = compute_multipoles(self,field,mesh,self.params.lambdaFF,lmax);
            
        end
               
        function [Csca,Cse,Csm] = compute_multipoles(self,field,mesh,lmd,lmax)
            %% Subroutine to compute multipoles in spherical coordinates
            % It implements the paper by Grahn et al. (2012)
            % Inputs:   field: total electric field (Ex,Ey,Ez) - it can be
            %           computed with the subroutine "extract_field"
            %           mesh: struct containing evaluation points and mesh
            %           elements volume for integration
            %           lmd: working vacuum-wavelength
            %           lmax: maximum angular moment for the decomposition
            if self.params.n_c~=self.params.n_s,error('Decomposition is not valid in an inhomogeneous surrounding');else nh = self.params.n_c;end
            [ae,am] = compute_multipolecoeff(self,field,mesh,lmd,lmax);
            k0 = 2*pi/lmd;
            sigma_geom = pi*self.params.ra*self.params.rb;
            for il=1:lmax
                Cse(il)=pi./((k0*nh)'.^2)*(2*il+1).*sum(abs(ae(il,:)).^2,2)/sigma_geom;
                Csm(il)=pi./((k0*nh)'.^2)*(2*il+1).*sum(abs(am(il,:)).^2,2)/sigma_geom;                
            end
            Csca=sum(Cse)+sum(Csm);
        end
        
        function [ae,am] = compute_multipolecoeff(self,field,mesh,lmd,lmax)
            %% Subroutine to compute multipoles in spherical coordinates
            % It implements the paper by Grahn et al. (2012)
            % Inputs:   field: total electric field (Ex,Ey,Ez) - it can be
            %           computed with the subroutine "extract_field"
            %           mesh: struct containing evaluation points and mesh
            %           elements volume for integration
            %           lmd: working vacuum-wavelength
            %           lmax: maximum angular moment for the decomposition
            
            if self.params.n_c~=self.params.n_s,error('Decomposition is not valid in an inhomogeneous surrounding');else nh = self.params.n_c;end
            r=sqrt(mesh.cord(1,:).^2+mesh.cord(2,:).^2+mesh.cord(3,:).^2);
            theta=acos(mesh.cord(3,:)./r);
            phi=atan2(mesh.cord(2,:),mesh.cord(1,:));
        
            Er=(sin(theta).*cos(phi)).*(field.Ex)+(sin(theta).*sin(phi)).*(field.Ey)+(cos(theta)).*(field.Ez);
            Et=(cos(theta).*cos(phi)).*(field.Ex)+(cos(theta).*sin(phi)).*(field.Ey)-(sin(theta)).*(field.Ez);
            Ep=(-sin(phi)).*(field.Ex)+cos(phi).*(field.Ey);

            n_res = mphglobal(self.mph,['nAlGaAs(',num2str(lmd),')'],'Dataset',self.dsetFF);
            k_res = mphglobal(self.mph,['kAlGaAs(',num2str(lmd),')'],'Dataset',self.dsetFF);
            eps0 = 8.85e-12;
            eps_res=(n_res+1j*k_res).^2;
            k0 = 2*pi/lmd;
            omega = k0*3e8;
            
%             Jr=(sin(theta).*cos(phi)).*(field.Jx)+(sin(theta).*sin(phi)).*(field.Jy)+(cos(theta)).*(field.Jz);
%             Jt=(cos(theta).*cos(phi)).*(field.Jx)+(cos(theta).*sin(phi)).*(field.Jy)-(sin(theta)).*(field.Jz);
%             Jp=(-sin(phi)).*(field.Jx)+cos(phi).*(field.Jy);

            Jr=1j*eps0*Er*omega*(eps_res-nh^2);
            Jt=1j*eps0*Et*omega*(eps_res-nh^2);
            Jp=1j*eps0*Ep*omega*(eps_res-nh^2);

            ae=zeros(lmax,lmax*2+1);
            am=zeros(lmax,lmax*2+1);

            kr=r*k0*nh;
            eta = 377/nh;
            E0 = sqrt(self.params.I0*2*eta*1e4);
            for l=1:1:lmax
                psi=Psi(l,kr);
                psi1=Psi1(l,kr);
                psi2=Psi2(l,kr);
                jl=sphericalBessel(l,kr);
                for m=-l:1:l
                    plm=Plm(cos(theta),l,m);
                    pilm=m./sin(theta).*plm;
                    taulm=Taulm(theta,l,m);        
                    Olm=1/(l*(l+1))^(1/2)*((2*l+1)/(4*pi)*gamma(l-m+1)/gamma(l+m+1))^(1/2);
                    intE=(mesh.vol.*exp(1j*m*phi)).*(plm.*(psi+psi2).*Jr+(psi1./kr).*(taulm.*Jt+1j*pilm.*Jp));
                    intM=(mesh.vol.*exp(1j*m*phi)).*jl.*(-1j*pilm.*Jt+taulm.*Jp);
                    mi=lmax+m+1;
                    ae(l,mi)=(1j)^(l-1).*(k0*nh).^2*Olm*eta/(E0*sqrt(pi*(2*l+1))).*(sum(intE,2));
                    am(l,mi)=(1j)^(l+1).*(k0*nh).^2*Olm*eta/(E0*sqrt(pi*(2*l+1))).*(sum(intM,2));
                    %clear plm pilm taulm Olm intE intM der
                end
                %   clear psi psi1 psi2 jl  
            end
            ae=-conj(ae); % COMSOL works with iwt convention
            am=-conj(am); % this is the convertion to -jwt convention
        end
        
        function [field,mesh] = extract_field(self,physics,dataset)
            
            temp=mpheval(self.mph, 'meshvol','pattern','gauss','dataset',dataset,'selection', self.res_domain);
            mesh.vol=temp.d1;   % mesh volume   
            mesh.cord=temp.p;        % mesh coordinates   
            z0 = self.mph.param.evaluate('z0');
            mesh.cord(3,:) = mesh.cord(3,:)-(z0+self.params.h0/2);
            
            temp=mpheval(self.mph,[physics,'Ex'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Ex=temp.d1;         % Electric field x component 
            temp=mpheval(self.mph,[physics,'Ey'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Ey=temp.d1;         % Electric field x component 
            temp=mpheval(self.mph,[physics,'Ez'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Ez=temp.d1;         % Electric field x component 
            
            temp=mpheval(self.mph,[physics,'Jx'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Jx=temp.d1;         % Electric field x component 
            temp=mpheval(self.mph,[physics,'Jy'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Jy=temp.d1;         % Electric field x component 
            temp=mpheval(self.mph,[physics,'Jz'],'pattern','gauss','dataset',dataset,'selection', self.res_domain,'Complexout','on');
            field.Jz=temp.d1;         % Electric field x component 
            clear temp
        end
         
    end
end
