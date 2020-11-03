clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% mu O to PO2 %%%%%%%%%%%%%%%%%%%%%%%%%%
aa = 201;
cc = 401;
Temp = linspace(0,2000,aa);
Temp = transpose(Temp);
convertunit = (1/6.022E23)*1000*(1/1.60218E-19);
%mu_O = linspace(-18,-4.5,cc);
%mu_O = transpose(mu_O);


P_torr = logspace(-50,0,cc);
P_torr = transpose(P_torr);
for bb = 1:aa
    H(bb) = -8E-12*Temp(bb)^3+5E-08*Temp(bb)^2+0.0003*Temp(bb)-0.0889;
    TxS(bb)=3E-07*Temp(bb)^2+0.0023*Temp(bb)-0.0596;
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
            +k_B*Temp(ii)*log(P_MPa(kk)/P0_MPa))/2;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% surface energies %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% bulk %%%%%%%

W333_numunits = 54;
W333_F = readmatrix('thermalpropsW333.txt');
W333_F = W333_F(:,2);
W333_F = W333_F*convertunit;
W333_E0 = -.69946699E+03;
mu_metW = (W333_F+W333_E0)/W333_numunits;

WO3222_numunits = 16;
WO3222_F = readmatrix('thermalpropsWO3222.txt');
WO3222_F = WO3222_F(:,2);
WO3222_F = WO3222_F*convertunit;
WO3222_E0 = -.58083580E+03;
mu_WO3 = (WO3222_F+WO3222_E0)/WO3222_numunits;

Wcutoff = (mu_WO3-mu_metW)/3;

Ba333_numunits = 54;
Ba333_F = readmatrix('thermalpropsBa333.txt');
Ba333_F = Ba333_F(:,2);
Ba333_F = Ba333_F*convertunit;
Ba333_E0 = -.10386067E+03;
mu_metBa = (Ba333_F+Ba333_E0)/Ba333_numunits;

BaO222_numunits = 32;
BaO222_F = readmatrix('thermalpropsBaO222.txt');
BaO222_F = BaO222_F(:,2);
BaO222_F = BaO222_F*convertunit;
BaO222_E0 = -.37202020E+03;
mu_BaO = (BaO222_F+BaO222_E0)/BaO222_numunits;

Bacutoff = (mu_BaO-mu_metBa);

Sc333_numunits = 54;
Sc333_F = readmatrix('thermalpropsSc333.txt');
Sc333_F = Sc333_F(:,2);
Sc333_F = Sc333_F*convertunit;
Sc333_E0 = -.34197851E+03;
mu_metSc = (Sc333_F+Sc333_E0)/Sc333_numunits;

Sc2O3111_numunits = 8;
Sc2O3111_F = readmatrix('thermalpropsSc2O3111.txt');
Sc2O3111_F = Sc2O3111_F(:,2);
Sc2O3111_F = Sc2O3111_F*convertunit;
Sc2O3111_E0 = -.36907745E+03;
mu_Sc2O3 = (Sc2O3111_F+Sc2O3111_E0)/Sc2O3111_numunits;

Sccutoff = (mu_Sc2O3-2*mu_metSc)/3;

%%%%%%% (001) %%%%%%%

W001surfunpert_area = 80.5;
W001surfunpert_Watoms = 64;
W001surfunpert_F = readmatrix('thermalpropsW001-surfunpert.txt');
W001surfunpert_F = W001surfunpert_F(:,2);
W001surfunpert_F = W001surfunpert_F*convertunit;
W001surfunpert_E0 = -.80941964E+03;
W001surfunpert_F = W001surfunpert_F+W001surfunpert_E0;
W001surfunpert_gamma = (W001surfunpert_F-W001surfunpert_Watoms*mu_metW)/...
    W001surfunpert_area;

W001surfpert_area = 80.5;
W001surfpert_Watoms = 64;
W001surfpert_F = readmatrix('thermalpropsW001-surfpert.txt');
W001surfpert_F = W001surfpert_F(:,2);
W001surfpert_F = W001surfpert_F*convertunit;
W001surfpert_E0 = -.81006004E+03;
W001surfpert_F = W001surfpert_F+W001surfpert_E0;
W001surfpert_gamma = (W001surfpert_F-W001surfpert_Watoms*mu_metW)/...
    W001surfpert_area;

O8W001_area = 80.5;
O8W001_Watoms = 52;
O8W001_Oatoms = 8;
O8W001_F = readmatrix('thermalpropsO-W001.txt');
O8W001_F = O8W001_F(:,2);
O8W001_F = O8W001_F*convertunit;
O8W001_E0 = -.72150756E+03;
O8W001_F = O8W001_F+O8W001_E0;

Ba4O8W001_area = 80.5; 
Ba4O8W001_Watoms = 52;
Ba4O8W001_Baatoms = 4;
Ba4O8W001_Oatoms = 8;
Ba4O8W001_config = 5.808E-6;
Ba4O8W001_F = readmatrix('thermalpropsBa4O8-W001.txt');
Ba4O8W001_F = Ba4O8W001_F(:,2);
Ba4O8W001_F = Ba4O8W001_F*convertunit;
Ba4O8W001_E0 = -.74413878E+03;
Ba4O8W001_F = Ba4O8W001_F+Ba4O8W001_E0;

Ba2Sc2O8W001_area = 80.5;
Ba2Sc2O8W001_Watoms = 52;
Ba2Sc2O8W001_Baatoms = 2;
Ba2Sc2O8W001_Scatoms = 2;
Ba2Sc2O8W001_Oatoms = 8;
Ba2Sc2O8W001_config = 8.712E-6;
Ba2Sc2O8W001_F = readmatrix('thermalpropsBa2Sc2O8-W001.txt');
Ba2Sc2O8W001_F = Ba2Sc2O8W001_F(:,2);
Ba2Sc2O8W001_F = Ba2Sc2O8W001_F*convertunit;
Ba2Sc2O8W001_E0 = -.75521094E+03;
Ba2Sc2O8W001_F = Ba2Sc2O8W001_F+Ba2Sc2O8W001_E0;

%%%%%%% (110) %%%%%%%

W110_area = 56.88; 
W110_Watoms = 36;
W110_F = readmatrix('thermalpropsW110.txt');
W110_F = W110_F(:,2);
W110_F = W110_F*convertunit;
W110_E0 = -.45508870E+03;
W110_F = W110_F+W110_E0;
W110_gamma = (W110_F-W110_Watoms*mu_metW)/...
    W110_area;

O8W110_area = 56.88;
O8W110_Watoms = 36;
O8W110_Oatoms = 8;
O8W110_F = readmatrix('thermalpropsO-W110.txt');
O8W110_F = O8W110_F(:,2);
O8W110_F = O8W110_F*convertunit;
O8W110_E0 = -.52503875E+03;
O8W110_F = O8W110_F+O8W110_E0;

Ba2O8W110_area = 56.88; 
Ba2O8W110_Watoms = 36;
Ba2O8W110_Baatoms = 2;
Ba2O8W110_Oatoms = 8;
Ba2O8W110_config = 6.67E-6;
Ba2O8W110_F = readmatrix('thermalpropsBa2O8-W110.txt');
Ba2O8W110_F = Ba2O8W110_F(:,2);
Ba2O8W110_F = Ba2O8W110_F*convertunit;
Ba2O8W110_E0 = -.53288599E+03;
Ba2O8W110_F = Ba2O8W110_F+Ba2O8W110_E0;

Ba2Sc2O8W110_area = 56.88; 
Ba2Sc2O8W110_Watoms = 36;
Ba2Sc2O8W110_Baatoms = 2;
Ba2Sc2O8W110_Scatoms = 2;
Ba2Sc2O8W110_Oatoms = 8;
Ba2Sc2O8W110_config = 1.23E-5;
Ba2Sc2O8W110_F = readmatrix('thermalpropsBa2Sc2O8-W110.txt');
Ba2Sc2O8W110_F = Ba2Sc2O8W110_F(:,2);
Ba2Sc2O8W110_F = Ba2Sc2O8W110_F*convertunit;
Ba2Sc2O8W110_E0 = -.54556299E+03;
Ba2Sc2O8W110_F = Ba2Sc2O8W110_F+Ba2Sc2O8W110_E0;


%%%%%%% (112) %%%%%%%

W112_area = 49.27; 
W112_Watoms = 34;
W112_F = readmatrix('thermalpropsW112.txt');
W112_F = W112_F(:,2);
W112_F = W112_F*convertunit;
W112_E0 = -.43013015E+03;
W112_F = W112_F+W112_E0;
W112_gamma = (W112_F-W112_Watoms*mu_metW)/...
    W112_area;

O4W112_area = 49.27;
O4W112_Watoms = 38;
O4W112_Oatoms = 4;
O4W112_F = readmatrix('thermalpropsO-W112.txt');
O4W112_F = O4W112_F(:,2);
O4W112_F = O4W112_F*convertunit;
O4W112_E0 = -.51392254E+03;
O4W112_F = O4W112_F+O4W112_E0;

Ba2O4W112_2x2_area = 98.533; 
Ba2O4W112_2x2_Watoms = 76;
Ba2O4W112_2x2_Baatoms = 4;
Ba2O4W112_2x2_Oatoms = 8;
Ba2O4W112_2x2_config = 4.745E-6;
Ba2O4W112_2x2_F = readmatrix('thermalpropsBa2O4W112_2x2.txt');
Ba2O4W112_2x2_F = Ba2O4W112_2x2_F(:,2);
Ba2O4W112_2x2_F = Ba2O4W112_2x2_F*convertunit;
Ba2O4W112_2x2_E0 = -.10547239E+04;
Ba2O4W112_2x2_F = Ba2O4W112_2x2_F+Ba2O4W112_2x2_E0;

Ba1Sc1O4W112_area = 98.533;
Ba1Sc1O4W112_Watoms = 76;
Ba1Sc1O4W112_Baatoms = 2;
Ba1Sc1O4W112_Oatoms = 8;
Ba1Sc1O4W112_Scatoms = 2;
Ba1Sc1O4W112_config = 7.118E-6;
Ba1Sc1O4W112_F = readmatrix('thermalpropsBa1Sc1O4-W112.txt');
Ba1Sc1O4W112_F = Ba1Sc1O4W112_F(:,2);
Ba1Sc1O4W112_F = Ba1Sc1O4W112_F*convertunit;
Ba1Sc1O4W112_E0 = -.10629497E+04;
Ba1Sc1O4W112_F = Ba1Sc1O4W112_F+Ba1Sc1O4W112_E0;

Ba2Sc2O4W112_area = 98.533;
Ba2Sc2O4W112_Watoms = 76;
Ba2Sc2O4W112_Baatoms = 4;
Ba2Sc2O4W112_Oatoms = 8;
Ba2Sc2O4W112_Scatoms = 4;
Ba2Sc2O4W112_config = 4.745E-6;
Ba2Sc2O4W112_F = readmatrix('thermalpropsBa2Sc2O4-W112.txt');
Ba2Sc2O4W112_F = Ba2Sc2O4W112_F(:,2);
Ba2Sc2O4W112_F = Ba2Sc2O4W112_F*convertunit;
Ba2Sc2O4W112_E0 = -.10804814E+04;
Ba2Sc2O4W112_F = Ba2Sc2O4W112_F+Ba2Sc2O4W112_E0;

Ba2Sc2O6W112_area = 147.7998;
Ba2Sc2O6W112_Watoms = 114;
Ba2Sc2O6W112_Baatoms = 4;
Ba2Sc2O6W112_Oatoms = 12;
Ba2Sc2O6W112_Scatoms = 4;
Ba2Sc2O6W112_config = 7.52E-6;
Ba2Sc2O6W112_F = readmatrix('thermalpropsBa2Sc2O6-W112.txt');
Ba2Sc2O6W112_F = Ba2Sc2O6W112_F(:,2);
Ba2Sc2O6W112_F = Ba2Sc2O6W112_F*convertunit;
Ba2Sc2O6W112_E0 = -.16071218E+04;
Ba2Sc2O6W112_F = Ba2Sc2O6W112_F+Ba2Sc2O6W112_E0;


k = 1;
q(k,1,:) = zeros(1,4);
q(k,2,:) = zeros(1,4);
goodshape=false;
Sc_O = true;
Ba_O= true;
W_O = true;

for i = 1:201
    
    mu_O = plotmu_O(:,1);
    mu_O = transpose(mu_O);
    goodshape = true;
    Sc_O = true;
    Ba_O= true;
    W_O = true;
    
    for j = 1:length(mu_O)
            
            if mu_O(j)< Bacutoff(i)
                mu_Ba(j) = mu_metBa(i);
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
                +Ba4O8W001_Baatoms*mu_Ba(j)+Ba4O8W001_Oatoms*mu_O(j)))...
                /Ba4O8W001_area)+(Temp(i)*Ba4O8W001_config);
            
            Ba2Sc2O8W001_gamma(j) = ((Ba2Sc2O8W001_F(i)-(Ba2Sc2O8W001_Watoms*mu_W(j)...
                +Ba2Sc2O8W001_Baatoms*mu_Ba(j)+Ba2Sc2O8W001_Oatoms*mu_O(j)+Ba2Sc2O8W001_Scatoms*mu_Sc(j)))/...
                Ba2Sc2O8W001_area+(Temp(i)*Ba2Sc2O8W001_config));
            
            O8W001_gamma(j) = (O8W001_F(i)-(O8W001_Watoms*mu_W(j)...
                +O8W001_Oatoms*mu_O(j)))/O8W001_area;
            
            W001set = [W001surfpert_gamma(i); Ba4O8W001_gamma(j); ...
                Ba2Sc2O8W001_gamma(j); O8W001_gamma(j)];
            lowestW001 = min(W001set);
            
           
                    
%%%%110%%%%

          O8W110_gamma(j) = (O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Oatoms*mu_O(j)))/Ba2O8W110_area;
            
          Ba2O8W110_gamma(j) = ((Ba2O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Baatoms*mu_Ba(j)+Ba2O8W110_Oatoms*mu_O(j)))/Ba2O8W110_area)+(Temp(i)*Ba2O8W110_config);
                      
          Ba2Sc2O8W110_gamma(j) = ((Ba2Sc2O8W110_F(i)-(Ba2Sc2O8W110_Watoms*mu_W(j)...
                +Ba2Sc2O8W110_Baatoms*mu_Ba(j)+Ba2Sc2O8W110_Oatoms*mu_O(j)...
                +Ba2Sc2O8W110_Scatoms*mu_Sc(j)))/Ba2Sc2O8W110_area)+(Temp(i)*Ba2Sc2O8W110_config);
       
          W110set = [W110_gamma(i); Ba2O8W110_gamma(j); ...
                Ba2Sc2O8W110_gamma(j); O8W110_gamma(j)];
          lowestW110 = min(W110set);

    
%%%%112%%%%
            
          O4W112_gamma(j) = (O4W112_F(i)-(O4W112_Watoms*mu_W(j)...
                +O4W112_Oatoms*mu_O(j)))/O4W112_area;
            
            Ba2O4W112_2x2_gamma(j) = ((Ba2O4W112_2x2_F(i)-(Ba2O4W112_2x2_Watoms*mu_W(j)...
                +Ba2O4W112_2x2_Baatoms*mu_Ba(j)+Ba2O4W112_2x2_Oatoms*mu_O(j)))/Ba2O4W112_2x2_area)+(Temp(i)*Ba2O4W112_2x2_config);
               
            Ba1Sc1O4W112_gamma(j) = ((Ba1Sc1O4W112_F(i)-(Ba1Sc1O4W112_Watoms*mu_W(j)...
                +Ba1Sc1O4W112_Baatoms*mu_Ba(j)+Ba1Sc1O4W112_Oatoms*mu_O(j)+...
                Ba1Sc1O4W112_Scatoms*mu_Sc(j)))/Ba1Sc1O4W112_area)+(Temp(i)*Ba1Sc1O4W112_config);
            
            Ba2Sc2O4W112_gamma(j) = ((Ba2Sc2O4W112_F(i)-(Ba2Sc2O4W112_Watoms*mu_W(j)...
                +Ba2Sc2O4W112_Baatoms*mu_Ba(j)+Ba2Sc2O4W112_Oatoms*mu_O(j)+...
                Ba2Sc2O4W112_Scatoms*mu_Sc(j)))/Ba2Sc2O4W112_area)+(Temp(i)*Ba2Sc2O4W112_config);
            
            Ba2Sc2O6W112_gamma(j) = ((Ba2Sc2O6W112_F(i)-(Ba2Sc2O6W112_Watoms*mu_W(j)...
                +Ba2Sc2O6W112_Baatoms*mu_Ba(j)+Ba2Sc2O6W112_Scatoms*mu_Sc(j)+...
                Ba2Sc2O6W112_Oatoms*mu_O(j)))/Ba2Sc2O6W112_area)+(Temp(i)*Ba2Sc2O6W112_config);
                 
           
            W112set = [W112_gamma(i); Ba2O4W112_2x2_gamma(j); ...
                Ba1Sc1O4W112_gamma(j); Ba2Sc2O4W112_gamma(j); ...
                Ba2Sc2O6W112_gamma(j); O4W112_gamma(j)];
            lowestW112 = min(W112set);
          
%%%%comparing surface energies%%%%
            
            
            gammaratioW112_W110(j,i) = lowestW112(j)/lowestW110(j);
            gammaratioW001_W110(j,i) = lowestW001(j)/W110_gamma(j);
           
           
           if (gammaratioW112_W110(j,i)>0.85) && (gammaratioW112_W110(j,i)<1.05) ...
                && (gammaratioW001_W110(j,i)>0.85) && (gammaratioW001_W110(j,i)<1.05) %...
                %&& (gammaratioW112_Ba2Sc2O6_bareW110(j,i)<gammaratioBa4O8W001_gamma_bareW110(j,i))
                
                disp('Yes!')
                if (~(goodshape))
                    
                q(k,1,:) = [Temp(i), P_torr(j), gammaratioW001_W110(j,i), ...
                        gammaratioW112_W110(j,i)];
                goodshape=true;
                end
                fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\n',q);
           else
               if (goodshape)
                   
                q(k,2,:)=[Temp(i), P_torr(j), gammaratioW001_W110(j-1,i), ...
                        gammaratioW112_W110(j-1,i)];
                k=k+1;
                goodshape=false;
               end 
           end
           %}
            
           if (Sc_O)
                if (Sccutoff(i)<mu_O(j))
                    Sc_O = false;
                    if exist('r','var')
                    r=[r;[Temp(i),P_torr(j),Sccutoff(i)]];
                    else
                    r=[Temp(i),P_torr(j),Sccutoff(i)];
                    end
                end
           end
           
           if (Ba_O)
                if (Bacutoff(i)<mu_O(j))
                    Ba_O = false;
                    if exist('s','var')
                    s=[s;[Temp(i),P_torr(j),Bacutoff(i)]];
                    else
                    s=[Temp(i),P_torr(j),Bacutoff(i)];
                    end
                end
           end
           
           if (W_O)
                if (Wcutoff(i)<mu_O(j))
                    W_O = false;
                    if exist('z','var')
                    z=[z;[Temp(i),P_torr(j),Wcutoff(i)]];
                    else
                    z=[Temp(i),P_torr(j),Wcutoff(i)];
                    end
                end
           end
           
           
           
           
    end
        
        
        %{
        %P is plot, L is label
        figure(i);
        hold on
        P1 = plot(mu_O,Ba4O8W001_gamma,'-.r','LineWidth',4);
        L1 = 'Ba_{0.50}O-top/W(0 0 1)';
        hold on
        P2 = plot(mu_O,Ba2O8W110_gamma,'-.b','LineWidth',4);
        L2 = 'Ba_{0.25}O-tri/W(1 1 0)';
        hold on
        P3 = plot(mu_O,Ba2O4W112_2x2_gamma,'-.g','LineWidth',4);
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
        P7 = yline(W001surfpert_gamma(i,:), 'r', 'LineWidth',4);
        L7 = 'Bare W(0 0 1)';
        hold on
        P8 = yline(W110_gamma(i,:), 'b', 'LineWidth',4);
        L8 = 'Bare W(1 1 0)';
        hold on
        P9 = yline(W112_gamma(i,:), 'g', 'LineWidth',4);
        L9 = 'Bare W(1 1 2)';
        hold on
        P10 = plot(mu_O,Ba2Sc2O8W001_gamma,'-<r','LineWidth',4);
        L10 = 'Ba_{0.25}Sc_{0.25}O-top/W(0 0 1)';
        hold on
        P11 = plot(mu_O,Ba2Sc2O8W110_gamma,'-<b','LineWidth',4);
        L11 = 'Ba_{0.25}Sc_{0.25}O-top/W(1 1 0)';
        hold on
        P12 = plot(mu_O,Ba1Sc1O4W112_gamma,'-<g','LineWidth',4);
        L12 = 'Ba_{0.25}Sc_{0.25}O-top/W(1 1 2)';
        hold on
%         P13 = plot(mu_O,Ba2Sc2O4W112_gamma,'-+g','LineWidth',2);
%         L13 = 'Ba_{0.5}Sc_{0.5}O-top/W(1 1 2)';
%         hold on
%         P14 = plot(mu_O,Ba2Sc2O6W112_gamma,'-sg','LineWidth',2);
%         L14 = '(BaSc)_{1/3} O-top/W(1 1 2)';
%         hold on

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
        ylabel({'Surface Energy (eV/Ang.^2)'},'fontsize', 32);
        axis([-12 -6.5 0 0.6])
        ax = gca;
        ax.LineWidth = 4;
        %axes('YColor','none');
        box on;
        %grid on;
        legend([P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12], L1, L2, L3, L4, L5, L6, ...
            L7, L8, L9, L10, L11, L12, 'fontsize', 18, 'Location','northeast','NumColumns',4);
        legend boxon;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-7.55,0.105,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['cathodesurfpaper_t',num2str(Temp(i,:)),'.png']);
       %}
    
end


%{
figure
contourf(P_torr,Temp,plotmu_O,100,'LineColor','none')
colorbar
xlabel({'P_{O2} (Torr)'},'fontsize', 32);
ylabel({'Temperature (K)'},'fontsize', 32);
set(gca,'XScale','log','FontSize',28);
ax = gca;
ax.LineWidth = 3;
box on
hold on

semilogx(q(:,1,2),q(:,1,1),'LineWidth',6,'LineStyle',':','color','red')
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
