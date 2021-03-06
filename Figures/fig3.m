function fig3(prn)

%% Local Variables
% Impostazione Figura
nx = 1;			% Numero di figure in orizzontale
ny = 4;			% Numero di figure in verticale
bx = 0.48*1.10;	    % Margini sinistro  e destro    [cm]
by = 0.30;		% Margini superiore e inferiore [cm]
fx = 7.5*0+7.5;		% Larghezza primo riquadro interno [cm]
fy = 1.8;		    % Altezza riquadri interni   [cm]
dx = 0.0;		% Distanza orizzontale tra i riquadri [cm]
dy = 0.32;		% Distanza verticale tra i riquadri   [cm]
ox = 0.30;		% Offset orizzontale [cm]
oy = 0.03;		% Offset verticale   [cm]
ospl = 0.0;
fnt  = 8;		% Dimensione dei font per gli assi
fnl  = 8;		% Dimensione dei font per le labels
fnl2 = 7;		% Dimensione dei font per le labels

fwidth  = nx*fx+(nx-1)*dx+2*bx;
fheigth = ny*fy+(ny-1)*dy+2*by+ospl;
%%
lincol1 = 'b';
lincol2 = [50 205 50]/256;
lincol3 = 'c';
lincol7 = 'm';
lincol8 = 'g';
lincol5=[0 0 0];          
linwid1=0.75;              
linwid2=1.55;              
linwid0=1.5;
%% Creazione Figura
figure('Renderer', 'Painters',...
    'units','centimeters',...
    'position',[bx, by, fwidth, fheigth], ...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPosition',[bx, by, fwidth, fheigth], ...
    'PaperType','A4', ...
    'Color','white' ...
    );
%%
% Creazione Riquadri e Grafici
% #format table ## [WaveView Analyzer] saved 16:33:34 Sun Mar 17 2019

GRID_ON    = 1;
GRID_MINOR = 0;
CONST      = 180/pi;

DCBUSINF.f    = struct2array(load('dcInfAc.mat','freq'));
DCBUSINF.Y_AA = -struct2array(load('dcInfAc.mat','Yaa_0'));
DCBUSINF.Y_BA = -struct2array(load('dcInfAc.mat','Yba_0'));

DCLINE.f    = struct2array(load('dcLineAc.mat','freq'));
DCLINE.Y_AA = -struct2array(load('dcLineAc.mat','Yaa_0'));
DCLINE.Y_BA = -struct2array(load('dcLineAc.mat','Yba_0'));

DCCOMPLETE.f    = struct2array(load('Dcs1Ac.mat','freq'));
DCCOMPLETE.Y_AA = -struct2array(load('Dcs1Ac.mat','Yaa_0'));
DCCOMPLETE.Y_BA = -struct2array(load('Dcs1Ac.mat','Yba_0'));

min([angle(DCBUSINF.Y_AA)*CONST; angle(DCLINE.Y_AA)*CONST; angle(DCCOMPLETE.Y_AA)*CONST])
max([angle(DCBUSINF.Y_AA)*CONST; angle(DCLINE.Y_AA)*CONST; angle(DCCOMPLETE.Y_AA)*CONST])


%%


for i=1:nx
    for j=1:ny
        h(i,j)=axes(...
            'units','centimeters',...
            'position',[bx+(i-1)*(fx+dx)+ox, by+(ny-j)*(fy+dy)+oy, fx fy],...
            'fontsize',fnt,...
            'fontname','helvetica',...
            'fontangle','normal',...
            'visible','on', ...
            'Color','none', ...
            'box','on', ...
            'XTickLabel', [], ...
            'YTickLabel', [], ...
            'XTick', [], ...
            'Ytick', []);
        
        OFF = 0;
        
        
        if j==1
            
            line([23 23], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');

            hold on
            
            line([61 61], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            line([152 152], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            line(DCBUSINF.f, abs(DCBUSINF.Y_AA), 'Color', lincol5, ...
                'LineWidth',  linwid2,...
                'Linestyle','-');
            
            line(DCLINE.f, abs(DCLINE.Y_AA), 'Color', lincol7, ...
                'LineWidth',  linwid2,...
                'Linestyle','-');
            
            line(DCCOMPLETE.f, abs(DCCOMPLETE.Y_AA), 'Color', lincol2, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            h(i,j).YScale = 'log';
            
            axis([0 350 0.001 1])
            
            if GRID_ON
                grid on;
            end
            
            if GRID_MINOR
                grid minor;
            end

        elseif j == 2
            line(DCBUSINF.f, angle(DCBUSINF.Y_AA)*CONST, 'Color', lincol5, ...
                'LineWidth',  linwid0,...
                'Linestyle','-');
            
            hold on
            
            line(DCLINE.f, angle(DCLINE.Y_AA)*CONST, 'Color', lincol7, ...
                'LineWidth',  linwid2,...
                'Linestyle','-');
            
            line(DCCOMPLETE.f, angle(DCCOMPLETE.Y_AA)*CONST, 'Color', lincol2, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
        elseif j == 3
            
            line([23 23], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');

            hold on
            
            line([61 61], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            line([152 152], [0.001 1], 'Color', lincol8, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            line(DCBUSINF.f, abs(DCBUSINF.Y_BA), 'Color', lincol5, ...
                'LineWidth',  linwid0,...
                'Linestyle','-');
            
            line(DCLINE.f, abs(DCLINE.Y_BA), 'Color', lincol7, ...
                'LineWidth',  linwid2,...
                'Linestyle','-');
            
            line(DCCOMPLETE.f, abs(DCCOMPLETE.Y_BA), 'Color', lincol2, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
            h(i,j).YScale = 'log';
            
            axis([0 350 0.001 1])
                       
        elseif j == 4
            line(DCBUSINF.f, angle(DCBUSINF.Y_BA)*CONST, 'Color', lincol5, ...
                'LineWidth',  linwid0,...
                'Linestyle','-');
            
            hold on
            
            line(DCLINE.f, angle(DCLINE.Y_BA)*CONST, 'Color', lincol7, ...
                'LineWidth',  linwid2,...
                'Linestyle','-');
            
            line(DCCOMPLETE.f, angle(DCCOMPLETE.Y_BA)*CONST, 'Color', lincol2, ...
                'LineWidth',  linwid1,...
                'Linestyle','--');
            
        end
        
        hold on
        axis xy
        
        
        if GRID_ON
            grid on;
        end
        
        if GRID_MINOR
            grid minor;
        end
        %
        if j==1
            format_ticks(h(i,j),{'$ $','$23$','$ $','$61$','$ $','$ $','$152$','$ $','$ $','$ $','$ $'},...
                {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'},...
                sort([0:50:350 23 61 152]), logspace(-3,0,4),...
                [],[],[0.000020 0.020]);
            axlbl={'$$ $$','$$ |Y_\mathrm{aa}^{(0)}| [\mathrm{S}]$$'};
            
            lt=[0.005, 0.77];
            text(lt(1),lt(2),axlbl{2}, ...
                'fontname','times',...
                'fontangle','italic',...
                'fontsize',fnl, ...
                'fontweight', 'normal', ...
                'interpreter','latex',...
                'HorizontalAlignment','l',...
                'VerticalAlignment','bottom', ...
                'Unit', 'normalized' ...
                );
        end
        
        if j==2
            format_ticks(h(i,j),{'$ $','$ $','$ $','$ $','$ $','$ $','$ $','$ $'},...
                {'$-200$','$-100$','$0$','$100$','$200$','$300$'},...
                ([0:50:350]), ([-200:100:300]),...
                [],[],[0.020 0.020]);
            axlbl={'$$ $$','$$ \angle Y_\mathrm{aa}^{(0)} [\mathrm{deg}]$$'};
            
            lt=[0.005, 0.77];
            text(lt(1),lt(2),axlbl{2}, ...
                'fontname','times',...
                'fontangle','italic',...
                'fontsize',fnl, ...
                'fontweight', 'normal', ...
                'interpreter','latex',...
                'HorizontalAlignment','l',...
                'VerticalAlignment','bottom', ...
                'Unit', 'normalized' ...
                );
        end
        
        if j==3
            format_ticks(h(i,j),{'$ $','$23$','$ $','$61$','$ $','$ $','$152$','$ $','$ $','$ $','$ $'},...
                {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'},...
                sort([0:50:350 23 61 152]), logspace(-3,0,4),...
                [],[],[0.000020 0.020]);
            axlbl={'$$ $$','$$ |Y_\mathrm{ba}^{(0)}| [\mathrm{S}]$$'};
            
            lt=[0.005, 0.77];
            text(lt(1),lt(2),axlbl{2}, ...
                'fontname','times',...
                'fontangle','italic',...
                'fontsize',fnl, ...
                'fontweight', 'normal', ...
                'interpreter','latex',...
                'HorizontalAlignment','l',...
                'VerticalAlignment','bottom', ...
                'Unit', 'normalized' ...
                );
        end
        
        
        
        if j==4
            format_ticks(h(i,j),{'$0$','$50$','$100$','$150$','$200$','$250$','$300$','$350$'},...
                {'$-200$','$-100$','$0$','$100$','$200$','$300$'},...
                ([0:50:350]), ([-200:100:300]),...
                 [],[],[0.020 0.020]);
                axlbl={'$$ f [\mathrm{Hz}] $$','$$ \angle Y_\mathrm{ba}^{(0)} [\mathrm{deg}] $$'};
            
            lt=[0.999, 0.003];
            text(lt(1),lt(2),axlbl{1}, ...
                'fontname','times',...
                'fontangle','italic',...
                'fontsize',fnl, ...
                'fontweight', 'normal', ...
                'interpreter','latex',...
                'HorizontalAlignment','right',...
                'VerticalAlignment','bottom', ...
                'Color', 'black', ...
                'Unit', 'normalized' ...
                );
            
            lt=[0.005, 0.77];
            text(lt(1),lt(2),axlbl{2}, ...
                'fontname','times',...
                'fontangle','italic',...
                'fontsize',fnl, ...
                'fontweight', 'normal', ...
                'interpreter','latex',...
                'HorizontalAlignment','l',...
                'VerticalAlignment','bottom', ...
                'Unit', 'normalized' ...
                );
        end
        
    end
    
    if(prn)
        print -loose -painters -depsc2 fig3.eps
    end
end