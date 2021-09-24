function fig7(prn)

%% Local Variables
% Impostazione Figura
nx = 2;			% Numero di figure in orizzontale
ny = 2;			% Numero di figure in verticale
bx = 0.5;	    % Margini sinistro  e destro    [cm]
by = 0.01;		% Margini superiore e inferiore [cm]
fx = 6.4*[0.5 0.5];	% Larghezza primo riquadro interno [cm]
fy = [1 1.5];	% Altezza riquadri interni   [cm]
dx = 0.7;		% Distanza orizzontale tra i riquadri [cm]
dy = 1;% Distanza verticale tra i riquadri   [cm]
ox = 0.25;		% Offset orizzontale [cm]
oy = 0.3;		% Offset verticale   [cm]
ospl = 0.0;
fnt  = 8;		% Dimensione dei font per gli assi
fnl  = 8;		% Dimensione dei font per le labels
fnl2 = 7;		% Dimensione dei font per le labels

fwidth  = sum(fx)+(nx-1)*dx+2*bx+ ospl;
fheigth = sum(fy)+(ny-1)*dy+2*by+ ospl;

if fwidth > 8.5
    disp('Error. Plot wider than single column')
    return
end

%%
lincol5=[0 0 0];           
lincol3=0.75*[1 1 1];       
lincol7=0.5*[1 1 1];        
lincol8=0.25*[1 1 1];
linwid1=0.5;              
linwid2=1.5;              
MrkSize1 = 3.0;
MrkSize2 = 3.0;

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

load('Dcs1_800M.mat'); 

Ypp_freq = freq.';
Ypp  = -Ypp_0.' ;
Ypp_w =   2*pi*Ypp_freq;
EQUAL_WEIGHT = 1;
Npoles = 10;                                                                %In the previous version of the paper it was set to 20
Niter  = 100;
Ka   = 2;                                                                   %1.-Strictly proper, 2.-Proper, 3.-Improper

if EQUAL_WEIGHT                                                             %Vector of 1 is used as custom_weights
    weight = ones(size(Ypp_freq));
else                                                                        %Previously computed custom_weights are loaded
    weight = 1./(abs(Ypp));
end

[Poles,Residues,Constant,Proportional] = ...
    fitting(Ypp,Ypp_freq,Npoles,Niter,Ka,weight);

Ypp_fit             = sum(Residues./(1i*Ypp_w - Poles)) + ...
    Constant + Proportional*1i*Ypp_w;

%%

Ypp_tf = 0;
k = 1;
while (k <= length(Residues))
    if imag(Residues(k)) ~= 0
        a = real(Residues(k));
        b = imag(Residues(k));
        c = real(Poles(k));
        d = imag(Poles(k));
        
        Ypp_tf = Ypp_tf + ...
            tf([2*a -2*(a*c + b*d)],[1 -2*c (c^2+d^2)]);
        k = k + 2;
        
    else
        Ypp_tf = Ypp_tf + tf([Residues(k)],[1 -Poles(k)]);
        k = k + 1;
    end
end
Ypp_tf = Ypp_tf + Constant + tf([Proportional 0],[1]);

%%

z = tf([1 0],[1]);

L_dc = logspace(-3,5,1e5);
poles = rlocus(Ypp_tf*z, L_dc);

idx(1) = find(L_dc >= 100e-3,1);
idx(2) = find(L_dc >= 150e-3,1);
idx(3) = find(L_dc >= 200e-3,1);
idx(4) = find(L_dc >= 250e-3,1);

% keyboard

%%

[re,im] = nyquist(Ypp_tf*z*200e-3);
re = squeeze(re);
im = squeeze(im);


for i=1:nx
    for j=1:ny
        [bx+(i-1)*dx+sum(fx(1:i))-fx(i)+ox  by+(ny-j)*dy+sum(fy(j:end-1))+oy fx(i) fy(j) i j];
        h(i,j)=axes(...
            'units','centimeters',...
            'position', [bx+(i-1)*dx+sum(fx(1:i))-fx(i)+ox  by+(ny-j)*dy+sum(fy(j:end-1))+oy fx(i) fy(j)],...
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
        
        
        if j == 1 && i == 1
            for k = 1: size(poles,1)
                line(real(poles(k,2:end-1)),imag(poles(k,2:end-1)),'Color', lincol5)
                hold on
                line(real(poles(k,1)),imag(poles(k,1)),'Marker','x','MarkerSize',MrkSize1,'Color', lincol5)
                line(real(poles(k,end)),imag(poles(k,end)),'Marker','o','MarkerSize',MrkSize1,'Color', lincol5)
                line([0 0],[-2e3 2e3],'Linestyle','-.','Color', lincol5)
            end
            
            h(i,j).XLim = [-130 20];
            h(i,j).YLim = [-240.5 240.5];
            h(i,j).XGrid = 'on';
            h(i,j).YGrid = 'on';
            h(i,j).XMinorGrid = 'off';
            h(i,j).YMinorGrid = 'off';
        end
        
        if j == 1 && i == 2
            for k = 1: size(poles,1)
                line(real(poles(k,2:end-1)),imag(poles(k,2:end-1)),'Color', lincol5)
                hold on
                line(real(poles(k,1)),imag(poles(k,1)),'Marker','x','MarkerSize',MrkSize1,'Color', lincol5)
                line(real(poles(k,end)),imag(poles(k,end)),'Marker','o','MarkerSize',MrkSize1,'Color', lincol5)
                line([0 0],[-2e3 2e3],'Linestyle','-.','Color', lincol5)
                
                line(real(poles(k,idx(1))),imag(poles(k,idx(1))),'Marker','*','MarkerSize',MrkSize2,'Color', lincol7)
                line(real(poles(k,idx(2))),imag(poles(k,idx(2))),'Marker','s','MarkerSize',MrkSize2,'Color', lincol7)
                line(real(poles(k,idx(3))),imag(poles(k,idx(3))),'Marker','d','MarkerSize',MrkSize2,'Color', lincol7)
                line(real(poles(k,idx(4))),imag(poles(k,idx(4))),'Marker','^','MarkerSize',MrkSize2,'Color', lincol7)
            end
            
            h(i,j).XLim = [-5 10];
            h(i,j).YLim = [-100 100];
            h(i,j).XGrid = 'on';
            h(i,j).YGrid = 'on';
            h(i,j).XMinorGrid = 'off';
            h(i,j).YMinorGrid = 'off';
            
        end
        
        if j == 2 && i == 1
            line(re,im,'Color', lincol5)
            hold on
            line(re,-im,'Color', lincol5)
            line(-1,0,'Marker','+','MarkerSize',MrkSize1,'Color', 'r')
            
            h(i,j).XLim = [-15 10];
            h(i,j).YLim = [-30 30];
            h(i,j).XGrid = 'on';
            h(i,j).YGrid = 'on';
            h(i,j).XMinorGrid = 'off';
            h(i,j).YMinorGrid = 'off';
        end
        
        if j == 2 && i == 2
            line(re,im,'Color', lincol5)
            hold on
            line(re,-im,'Color', lincol5)
            line(-1,0,'Marker','+','MarkerSize',MrkSize1,'Color', 'r')
            
            h(i,j).XLim = [-1.5 0.15];
            h(i,j).YLim = [-2 2];
            h(i,j).XGrid = 'on';
            h(i,j).YGrid = 'on';
            h(i,j).XMinorGrid = 'off';
            h(i,j).YMinorGrid = 'off';
        end
        
        
        if j == 1 && i == 1
            format_ticks(h(i,j),{'$-130$','$-80$','$-30$','$20$'},...
                {'$-240$','$-120$','$0$','$120$','$240$'},...
                ([-130:50:20]), ([-240:120:240]),...
                [],[],[0.020 0.020]);
            axlbl={'$$ $$','$$ $$'};
            
            lt=[0.005, 0.75];
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
        
        if j == 1 && i == 2
            format_ticks(h(i,j),{'$-5$','$-2$','$1$','$4$','$7$','$10$'},...
                {'$-100$','$-50$','$0$','$50$','$100$'},...
                ([-5:3:10]), ([-100:50:100]),...
                [],[],[0.020 0.020]);
            axlbl={'$$ $$','$$ $$'};
            
            lt=[0.005, 0.75];
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
        
        if j == 2 && i == 1
            format_ticks(h(i,j),{'$-15$','$-10$','$-5$','$0$','$5$','$10$'},...
                {'$-30$','$ $','$-10$','$ $','$10$','$ $','$30$'},...
                ([-15:5:10]), ([-30:10:30]),...
                [],[],[0.020 0.020]);
            axlbl={'$$ $$','$$ $$'};
            
            lt=[0.005, 0.75];
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
        
                if j == 2 && i == 2
            format_ticks(h(i,j),{'$-1.5$','$-1$','$-0.5$','$0.5$'},...
                {'$-2$','$-1$','$0$','$1$','$2$'},...
                ([-1.5:0.5:0]), ([-2:1:2]),...
                [],[],[0.020 0.020]);
            axlbl={'$$ $$','$$ $$'};

            lt=[0.005, 0.75];
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
    
end

if(prn)
    print -loose -painters -depsc2 ..\figs\fig7.eps
end