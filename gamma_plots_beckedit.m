clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% mu O to PO2 %%%%%%%%%%%%%%%%%%%%%%%%%%
aa = 201;
T = linspace(0,2000,aa); %Torr
T = transpose(T);
%{
cc = 161;
P_torr = logspace(-50,0,cc);
P_torr = transpose(P_torr);
for bb = 1:aa
    H(bb) = -8E-12*T(bb)^3+5E-08*T(bb)^2+0.0003*T(bb)-0.0889;
    TxS(bb)=3E-07*T(bb)^2+0.0023*T(bb)-0.0596;
end
H = transpose(H);
TxS = transpose(TxS);


hO = 1.36; %eV
P_MPa = P_torr*0.000133322; %MPa
P0_MPa = 0.1; %MPa
P0_torr = 750.062; %torr
T_r = 298.15; %K
E_O2_VASP = -9.6887; %eV
k_B = 8.617E-5; %eV/K


for ii = 1:aa
    for kk = 1:cc
        plotmu_O(ii,kk) = (E_O2_VASP+hO+H(ii)-TxS(ii)...
            +k_B*T(ii)*log(P_MPa(kk)/P0_MPa))/2;
    end
end


figure
contourf(P_torr,T,plotmu_O,100,'LineColor','none')
colorbar
xlabel({'P_{O2} (Torr)'},'fontsize', 32);
ylabel({'Temperature (K)'},'fontsize', 32);
set(gca,'XScale','log','FontSize',28);
ax = gca;
ax.LineWidth = 3;
box on
hold on
%}
%%%%%%%%%%%%%%%%%%%%%%%%%% surface energies %%%%%%%%%%%%%%%%%%%%%%%%%%

load mu_O.mat
load mu_BaO.mat
load mu_metallicBa.mat
load mu_metW.mat
load mu_WO3.mat
load Ba4O8W001_F.mat
load Ba2O8W110_F.mat
load Ba2O4W112_F.mat
load W112_O4_trough_Ba2_F.mat
load O8W001_F.mat
load O8W110_F.mat
load O4W112_F.mat
load Wcutoff.mat
load Bacutoff.mat
load gamma_bareW001.mat
load gamma_W110.mat
load gamma_bareW112.mat
load Ba1Sc1O4W112_F.mat
load Sccutoff.mat
load mu_Sc2O3.mat
load mu_metSc.mat
load Ba2Sc2O4W112_F.mat
load Ba0p5_O_tri_F.mat
load Ba2O8W001_F.mat
load Ba2Sc2O8W001_F.mat
load W112_Ba2Sc2O6_F.mat

Ba4O8W001_area = 80.5; %this is the 2A value
Ba4O8W001_Watoms = 52;
Ba4O8W001_Baatoms = 4;
Ba4O8W001_Oatoms = 8;
Ba4O8W001_config = 5.808E-6;

Ba2O8W001_area = 80.5;
Ba2O8W001_Watoms = 52;
Ba2O8W001_Baatoms = 2;
Ba2O8W001_Oatoms = 8;
Ba2O8W001_config = 4.712E-6;

Ba2Sc2O8W001_area = 80.5;
Ba2Sc2O8W001_Watoms = 52;
Ba2Sc2O8W001_Baatoms = 2;
Ba2Sc2O8W001_Scatoms = 2;
Ba2Sc2O8W001_Oatoms = 8;
Ba2Sc2O8W001_config = 8.712E-6;

Ba2O8W110_area = 56.88; 
Ba2O8W110_Watoms = 36;
Ba2O8W110_Baatoms = 2;
Ba2O8W110_Oatoms = 8;
Ba2O8W110_config = 6.67E-6;

Ba2O4W112_area = 49.26; 
Ba2O4W112_Watoms = 38;
Ba2O4W112_Baatoms = 2;
Ba2O4W112_Oatoms = 4;
Ba2O4W112_config = 4.745E-6;

Ba1Sc1O4W112_area = 98.52;
Ba1Sc1O4W112_Watoms = 76;
Ba1Sc1O4W112_Baatoms = 2;
Ba1Sc1O4W112_Oatoms = 8;
Ba1Sc1O4W112_Scatoms = 2;
Ba1Sc1O4W112_config = 7.118E-6;

Ba2Sc2O4W112_area = 98.52;
Ba2Sc2O4W112_Watoms = 76;
Ba2Sc2O4W112_Baatoms = 4;
Ba2Sc2O4W112_Oatoms = 8;
Ba2Sc2O4W112_Scatoms = 4;
Ba2Sc2O4W112_config = 4.745E-6;

Ba2Sc2O6W112_area = 147.7998;
Ba2Sc2O6W112_Watoms = 114;
Ba2Sc2O6W112_Baatoms = 4;
Ba2Sc2O6W112_Oatoms = 12;
Ba2Sc2O6W112_Scatoms = 4;
Ba2Sc2O6W112_config = 7.52E-6;

%{
k = 1;
q(k,1,:) = zeros(1,4);
q(k,2,:) = zeros(1,4);
goodshape=false;
Sc_O = true;
Ba_O= true;
W_O = true;
%}
for i = 201
    %{
    mu_O = plotmu_O(i,:);
    goodshape = false;
    Sc_O = true;
    Ba_O= true;
    W_O = true;
    %}
        for j = 1:length(mu_O)
            
            if mu_O(j)< Bacutoff(i)
                mu_Ba(j) = mu_metallicBa(i);
            else 
                mu_Ba(j) = mu_BaO(i)- mu_O(j);
            end
            
            if mu_O(j) < Wcutoff(i)
                mu_W(j) = mu_metW(i);
            else 
                mu_W(j) = mu_WO3(i)- 3*mu_O(j);
            end
            
            if mu_O(j) < Sccutoff(i)
                mu_Sc(j) = mu_metSc(i);
            else 
                mu_Sc(j) = (mu_Sc2O3(i)- 3*mu_O(j))/2;
            end
%%%%001%%%%
            Ba4O8W001_gamma(j) = ((Ba4O8W001_F(i)-(Ba4O8W001_Watoms*mu_W(j)...
                +Ba4O8W001_Baatoms*mu_Ba(j)+Ba4O8W001_Oatoms*mu_O(j)))/Ba4O8W001_area)+(T(i)*Ba4O8W001_config);
            
            Ba2O8W001_gamma(j) = ((Ba2O8W001_F(i)-(Ba2O8W001_Watoms*mu_W(j)...
                +Ba2O8W001_Baatoms*mu_Ba(j)+Ba2O8W001_Oatoms*mu_O(j)))/Ba2O8W001_area+(T(i)*Ba2O8W001_config));
            
            Ba2Sc2O8W001_gamma(j) = ((Ba2Sc2O8W001_F(i)-(Ba2Sc2O8W001_Watoms*mu_W(j)...
                +Ba2Sc2O8W001_Baatoms*mu_Ba(j)+Ba2Sc2O8W001_Oatoms*mu_O(j)+Ba2Sc2O8W001_Scatoms*mu_Sc(j)))/...
                Ba2O8W001_area+(T(i)*Ba2Sc2O8W001_config));
            
            O8W001_gamma(j) = (O8W001_F(i)-(Ba4O8W001_Watoms*mu_W(j)...
                +Ba4O8W001_Oatoms*mu_O(j)))/Ba4O8W001_area;                     
%%%%110%%%%

            Ba2O8W110_gamma(j) = ((Ba2O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Baatoms*mu_Ba(j)+Ba2O8W110_Oatoms*mu_O(j)))/Ba2O8W110_area)+(T(i)*Ba2O8W110_config);
                
            O8W110_gamma(j) = (O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Oatoms*mu_O(j)))/Ba2O8W110_area;
%%%%112%%%%
            
            Ba2O4W112_gamma(j) = ((Ba2O4W112_F(i)-(Ba2O4W112_Watoms*mu_W(j)...
                +Ba2O4W112_Baatoms*mu_Ba(j)+Ba2O4W112_Oatoms*mu_O(j)))/Ba2O4W112_area)+(T(i)*Ba2O4W112_config);

            %W112_O4_trough_Ba2_gamma(n) = (W112_O4_trough_Ba2_F(i)-(Ba2O4W112_Watoms*mu_W(n)...
            %    +Ba2O4W112_Baatoms*mu_Ba(n)+Ba2O4W112_Oatoms*mu_O(n)))/Ba2O4W112area;
            
            O4W112_gamma(j) = (O4W112_F(i)-(Ba2O4W112_Watoms*mu_W(j)...
                +Ba2O4W112_Oatoms*mu_O(j)))/Ba2O4W112_area;
            
            Ba1Sc1O4W112_gamma(j) = ((Ba1Sc1O4W112_F(i)-(Ba1Sc1O4W112_Watoms*mu_W(j)...
                +Ba1Sc1O4W112_Baatoms*mu_Ba(j)+Ba1Sc1O4W112_Oatoms*mu_O(j)+...
                Ba1Sc1O4W112_Scatoms*mu_Sc(j)))/Ba1Sc1O4W112_area)+(T(i)*Ba1Sc1O4W112_config);
            
            Ba2Sc2O4W112_gamma(j) = ((Ba2Sc2O4W112_F(i)-(Ba2Sc2O4W112_Watoms*mu_W(j)...
                +Ba2Sc2O4W112_Baatoms*mu_Ba(j)+Ba2Sc2O4W112_Oatoms*mu_O(j)+...
                Ba2Sc2O4W112_Scatoms*mu_Sc(j)))/Ba2Sc2O4W112_area)+(T(i)*Ba2Sc2O4W112_config);
            
            %Ba0p5_O_tri_gamma(j) = (Ba0p5_O_tri_F(i)-(Ba2O4W112_Watoms*mu_W(j)...
            %    +Ba2O4W112_Baatoms*mu_Ba(j)+Ba2O4W112_Oatoms*mu_O(j)))/Ba2O4W112area;
            
            W112_Ba2Sc2O6_gamma(j) = ((W112_Ba2Sc2O6_F(i)-(Ba2Sc2O6W112_Watoms*mu_W(j)...
                +Ba2Sc2O6W112_Baatoms*mu_Ba(j)+Ba2Sc2O6W112_Scatoms*mu_Sc(j)+...
                Ba2Sc2O6W112_Oatoms*mu_O(j)))/Ba2Sc2O6W112_area)+(T(i)*Ba2Sc2O6W112_config);
            
           
%%%%comparing surface energies%%%%
            
            %{
            gammaratioW112_Ba2Sc2O6_bareW110(j,i) = W112_Ba2Sc2O6_gamma(j)/gamma_W110(i);
            gammaratioBa4O8W001_gamma_bareW110(j,i) = Ba4O8W001_gamma(j)/gamma_W110(i);
           
           
           if (gammaratioW112_Ba2Sc2O6_bareW110(j,i)>0.85) && (gammaratioW112_Ba2Sc2O6_bareW110(j,i)<1.05) ...
                && (gammaratioBa4O8W001_gamma_bareW110(j,i)>0.85) && (gammaratioBa4O8W001_gamma_bareW110(j,i)<1.05) %...
                %&& (gammaratioW112_Ba2Sc2O6_bareW110(j,i)<gammaratioBa4O8W001_gamma_bareW110(j,i))
                
                %disp('Yes!')
                if (~(goodshape))
                    
                q(k,1,:) = [T(i), P_torr(j), gammaratioBa4O8W001_gamma_bareW110(j,i), ...
                        gammaratioW112_Ba2Sc2O6_bareW110(j,i)];
                goodshape=true;
                end
                fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\n',q);
           else
               if (goodshape)
                   
                q(k,2,:)=[T(i), P_torr(j), gammaratioBa4O8W001_gamma_bareW110(j-1,i), ...
                        gammaratioW112_Ba2Sc2O6_bareW110(j-1,i)];
                k=k+1;
                goodshape=false;
               end 
           end
           
           if (Sc_O)
                if (Sccutoff(i)<mu_O(j))
                    Sc_O = false;
                    if exist('r','var')
                    r=[r;[T(i),P_torr(j),Sccutoff(i)]];
                    else
                    r=[T(i),P_torr(j),Sccutoff(i)];
                    end
                end
           end
           
           if (Ba_O)
                if (Bacutoff(i)<mu_O(j))
                    Ba_O = false;
                    if exist('s','var')
                    s=[s;[T(i),P_torr(j),Bacutoff(i)]];
                    else
                    s=[T(i),P_torr(j),Bacutoff(i)];
                    end
                end
           end
           
           if (W_O)
                if (Wcutoff(i)<mu_O(j))
                    W_O = false;
                    if exist('z','var')
                    z=[z;[T(i),P_torr(j),Wcutoff(i)]];
                    else
                    z=[T(i),P_torr(j),Wcutoff(i)];
                    end
                end
           end
           %}
           
           
           
		end
        
        
        
        %P is plot, L is label
        figure(i);
        hold on
        P1 = plot(mu_O,Ba4O8W001_gamma,'-.r','LineWidth',4);
        L1 = 'Ba_{0.50}O-top/W(0 0 1)';
        hold on
        P2 = plot(mu_O,Ba2O8W110_gamma,'-.b','LineWidth',4);
        L2 = 'Ba_{0.25}O-tri/W(1 1 0)';
        hold on
        P3 = plot(mu_O,Ba2O4W112_gamma,'-.g','LineWidth',4);
        L3 = 'Ba_{0.50}O-top/W(1 1 2)';
        hold on
        P4 = plot(mu_O,O8W001_gamma,':r','LineWidth',4);
        L4 = 'O-top/W(0 0 1)';
        hold on
        P5 = plot(mu_O,O8W110_gamma,':b','LineWidth',4);
        L5 = 'O-tri/W(1 1 0)';
        hold on
        P6 = plot(mu_O,O4W112_gamma,':g','LineWidth',4);
        L6 = 'O-top/W(1 1 2)';
        hold on
        P10 = plot(mu_O,Ba1Sc1O4W112_gamma,':g','LineWidth',2);
        L10 = 'Ba_{0.25}Sc_{0.25}O-top/W(1 1 2)';
        %hold on
        %P11 = plot(mu_O,Ba2Sc2O4W112_gamma,'-+g','LineWidth',2);
        %L11 = 'Ba_{0.5}Sc_{0.5}O-top/W(1 1 2)';
        %hold on
        P13 = plot(mu_O,Ba2O8W001_gamma,'-+r','LineWidth',2);
        L13 = 'Ba_{0.25}O-top/W(0 0 1)';
        hold on
        P14 = plot(mu_O,Ba2Sc2O8W001_gamma,'-<r','LineWidth',2);
        L14 = 'Ba_{0.25}Sc_{0.25}O-top/W(0 0 1)';
        hold on
        P12 = plot(mu_O,W112_Ba2Sc2O6_gamma,'-<g','LineWidth',2);
        L12 = '(BaSc)_{1/3} O-top/W(1 1 2)';
        hold on
        P7 = yline(gamma_bareW001(i,:), 'r', 'LineWidth',4);
        L7 = 'Bare W(0 0 1)';
        P8 = yline(gamma_W110(i,:), 'b', 'LineWidth',4);
        L8 = 'Bare W(1 1 0)';
        P9 = yline(gamma_bareW112(i,:), 'g', 'LineWidth',4);
        L9 = 'Bare W(1 1 2)';
        hold on
        
        line([Wcutoff(i,:) Wcutoff(i,:)], [-10 0.7],'Color','c', ...
            'LineWidth', 3, 'LineStyle','-');
        text([Wcutoff(i,:) + 0.05], 0.05, 'WO_3','fontsize', 28, 'Color', 'c');
        hold on
        text([Wcutoff(i,:) - 0.2], 0.05, 'W','fontsize', 28, 'Color', 'c');
        hold on
        
        line([Bacutoff(i,:) Bacutoff(i,:)], [-10 0.7],'Color','m', ...
            'LineWidth', 3, 'LineStyle','-');
        text([Bacutoff(i,:) + 0.05], 0.05, 'BaO','fontsize', 28, 'Color', 'm');
        hold on
        text([Bacutoff(i,:) - 0.2], 0.05, 'Ba','fontsize', 28, 'Color', 'm');
        
        line([Sccutoff(i,:) Sccutoff(i,:)], [-10 0.7],'Color','k', ...
            'LineWidth', 3, 'LineStyle','-');
        text([Sccutoff(i,:) + 0.05], 0.05, 'Sc_2O_3','fontsize', 24, 'Color', 'k');
        hold on
        text([Sccutoff(i,:) - 0.2], 0.05, 'Sc','fontsize', 24, 'Color', 'k');
        
       
        hold off;
        set(gca,'FontSize',32);
        xlabel({'\mu_O (eV)'},'fontsize', 32);
        ylabel({'Surface Energy (eV/ang^2)'},'fontsize', 32);
        axis([-12 -7 0 0.6])
        ax = gca;
        ax.LineWidth = 4;
        %axes('YColor','none');
        box on;
        %grid on;
        legend([P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P12; P13; P14], L1, L2, L3, L4, ...
            L5, L6, L7, L8, L9, L10, L12, L13, L14, 'fontsize', 18,'Location','north','NumColumns',5);
        %legend boxoff;
        txt = {['T =',num2str(T(i,:)),'K']};
        text(-8.5,0.48,txt,'fontsize', 26);
        set(gcf, 'Position',  [100, 100, 1800, 1000]);
        saveas(gcf,['allstablesurfaces_t',num2str(T(i,:)),'.png']);
        gamma_bareW001(i)
        gamma_W110(i)
        gamma_bareW112(i)
        Wcutoff(i,:)
       
end

%{

semilogx(q(7:109,1,2),q(7:109,1,1),'LineWidth',6,'LineStyle',':','color','red')
hold on
semilogx(q(:,2,2),q(:,2,1),'LineWidth',6,'LineStyle',':','color','blue')
hold on
semilogx(r(110:aa,2),r(110:aa,1),'LineWidth',6,'color','black')
hold on
semilogx(s(88:aa,2),s(88:aa,1),'LineWidth',6,'color','magenta')
hold on
semilogx(z(57:aa,2),z(57:aa,1),'LineWidth',6,'color','cyan')
ax.LineWidth = 3;
box on
xlim([1E-50 1])
%}