function plotHandle = plot_Csca(x,y,Csca,params,kind)

 if strcmp(kind,'radius')
    r = x.*1e9; 
    str=sprintf('Linear Scattering \\lambda:%1.3e m - h:%1.3e m - I_0:%1.2e GW/cm^2', params.lambdaFF, params.h0, params.I0*1e-9);
    hold on
    grid on
    box on
    hAx = plot(r,Csca);
    plotHandle = hAx;
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1)
    title(str)
    xlabel('Radius [nm]')
    
    ylabel('\sigma_{sca}/\sigma_{geom}') % left y-axis 
    hold off
    
 elseif strcmp(kind,'radii')
    r1 = x.*1e9;
    r2 = y.*1e9;
    str=sprintf('Linear Scattering \\lambda:%1.3e m - h:%1.3e m - I_0:%1.2e GW/cm^2', params.lambdaFF, params.h0, params.I0*1e-9);
   
    [X,Y] = meshgrid(r2,r1);
    
    figure(1)
   
    mesh(X,Y,Csca)
    ylabel('Ra [nm]')
    xlabel('Rb [nm]')
    %zlabel('\sigma_{sca}/\sigma_{geom}')
    colorbar
    
    view(2)
    title(str)
    
 elseif strcmp(kind,'lambda')
    lmd = x.*1e9; 
    str=sprintf('Linear Scattering r:%1.3e m - h:%1.3e m - I_0:%1.2e GW/cm^2', params.ra, params.h0, params.I0*1e-9);
    hold on
    grid on
    box on
    hAx = plot(lmd,Csca,'linewidth',1.5);
    plotHandle = hAx;
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1)
    
    title(str)
    xlabel('Wavelength FF \lambda_{FF} [nm]')
    
    ylabel('\sigma_{sca}/\sigma_{geom}') % left y-axis 


    hold off
    
 else
    error('Plot type not accepted: use radius/lambda')
 end
end