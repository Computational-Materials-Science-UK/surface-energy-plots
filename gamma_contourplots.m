clear all;
close all;

load areafrac.mat

aa = 201;
cc = 166;
Temp = linspace(0,2000,aa);
Temp = transpose(Temp);
convertunit = (1/6.022E23)*1000*(1/1.60218E-19);
SE_mu_O = linspace(-18,-4.5,cc);
SE_mu_O = transpose(SE_mu_O);

%%%%%%%%%%%%%%%%%%%%%%%%%% PO2 to mu O %%%%%%%%%%%%%%%%%%%%%%%%%%
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
Sc2O3111_E0 = -.36207745E+03;
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

W001surfpert_area = 80.5;
W001surfpert_Watoms = 64;
W001surfpert_F = readmatrix('thermalpropsW001-surfpert.txt');
W001surfpert_F = W001surfpert_F(:,2);
W001surfpert_F = W001surfpert_F*convertunit;
W001surfpert_E0 = -.81006004E+03;
W001surfpert_F = W001surfpert_F+W001surfpert_E0;
W001surfpert_WF = 4.48;

O8W001_area = 80.5;
O8W001_Watoms = 52;
O8W001_Oatoms = 8;
O8W001_F = readmatrix('thermalpropsO-W001.txt');
O8W001_F = O8W001_F(:,2);
O8W001_F = O8W001_F*convertunit;
O8W001_E0 = -.72150756E+03;
O8W001_F = O8W001_F+O8W001_E0;
O8W001_WF = 7.63;

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
Ba4O8W001_WF = 2.47;

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
Ba2Sc2O8W001_WF = 0.82;

Ba2O8W001_area = 80.5; 
Ba2O8W001_Watoms = 52;
Ba2O8W001_Baatoms = 2;
Ba2O8W001_Oatoms = 8;
Ba2O8W001_config = 6.67E-6;
Ba2O8W001_F = readmatrix('thermalpropsBa2O8-W001.txt');
Ba2O8W001_F = Ba2O8W001_F(:,2);
Ba2O8W001_F = Ba2O8W001_F*convertunit;
Ba2O8W001_E0 = -.73520456E+03;
Ba2O8W001_F = Ba2O8W001_F+Ba2O8W001_E0;
Ba2O8W001_WF = 1.23;

%%%%%%% (110) %%%%%%%

W110_area = 56.88; 
W110_Watoms = 36;
W110_F = readmatrix('thermalpropsW110.txt');
W110_F = W110_F(:,2);
W110_F = W110_F*convertunit;
W110_E0 = -.45508870E+03;
W110_F = W110_F+W110_E0;
W110_WF = 4.81;

O8W110_area = 56.88;
O8W110_Watoms = 36;
O8W110_Oatoms = 8;
O8W110_F = readmatrix('thermalpropsO-W110.txt');
O8W110_F = O8W110_F(:,2);
O8W110_F = O8W110_F*convertunit;
O8W110_E0 = -.52503875E+03;
O8W110_F = O8W110_F+O8W110_E0;
O8W110_WF = 7.27;

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
Ba2O8W110_WF = 2.67;

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
Ba2Sc2O8W110_WF = 2.47;


%%%%%%% (112) %%%%%%%

W112_area = 49.27; 
W112_Watoms = 34;
W112_F = readmatrix('thermalpropsW112.txt');
W112_F = W112_F(:,2);
W112_F = W112_F*convertunit;
W112_E0 = -.43013015E+03;
W112_F = W112_F+W112_E0;
W112_WF = 4.33;

O4W112_area = 49.27;
O4W112_Watoms = 38;
O4W112_Oatoms = 4;
O4W112_F = readmatrix('thermalpropsO-W112.txt');
O4W112_F = O4W112_F(:,2);
O4W112_F = O4W112_F*convertunit;
O4W112_E0 = -.51392254E+03;
O4W112_F = O4W112_F+O4W112_E0;
O4W112_WF = 7.94;

Sc2W112top_area = 49.27;
Sc2W112top_Watoms = 38;
Sc2W112top_Scatoms = 2;
Sc2W112top_F = readmatrix('thermalpropsSc2-W112.txt');
Sc2W112top_F = Sc2W112top_F(:,2);
Sc2W112top_F = Sc2W112top_F*convertunit;
Sc2W112top_E0 = -.49288297E+03;
Sc2W112top_F = Sc2W112top_F+Sc2W112top_E0;
Sc2W112top_WF= 2.887924;

Sc4W112top_area = 49.27;
Sc4W112top_Watoms = 38;
Sc4W112top_Scatoms = 4;
Sc4W112top_F = readmatrix('thermalpropsSc4-W112.txt');
Sc4W112top_F = Sc4W112top_F(:,2);
Sc4W112top_F = Sc4W112top_F*convertunit;
Sc4W112top_E0 = -.50486328E+03;
Sc4W112top_F = Sc4W112top_F+Sc4W112top_E0;
Sc4W112top_WF = 3.4339847;

Sc4triW112_area = 49.27;
Sc4triW112_Watoms = 38;
Sc4triW112_Scatoms = 4;
Sc4triW112_F = readmatrix('thermalpropsSc4triW112.txt');
Sc4triW112_F = Sc4triW112_F(:,2);
Sc4triW112_F = Sc4triW112_F*convertunit;
Sc4triW112_E0 = -.51116579E+03;
Sc4triW112_F = Sc4triW112_F+Sc4triW112_E0;
Sc4triW112_WF = 2.6720228;

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
Ba2O4W112_2x2_WF = 1.20;

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
Ba1Sc1O4W112_WF = 2.2;

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
Ba2Sc2O4W112_WF = 2.45;

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
Ba2Sc2O6W112_WF = 1.45;

Ba2Sc2O6halfdecW112_area = 147.7998/2;
Ba2Sc2O6halfdecW112_Watoms = 72;
Ba2Sc2O6halfdecW112_Baatoms = 2;
Ba2Sc2O6halfdecW112_Oatoms = 6;
Ba2Sc2O6halfdecW112_Scatoms = 2;
Ba2Sc2O6halfdecW112_config = 7.52E-6/2;
Ba2Sc2O6halfdecW112_F = readmatrix('thermalpropsBa2Sc2O6_dec-W112.txt');
Ba2Sc2O6halfdecW112_F = Ba2Sc2O6halfdecW112_F(:,2);
Ba2Sc2O6halfdecW112_F = Ba2Sc2O6halfdecW112_F*convertunit;
Ba2Sc2O6halfdecW112_E0 = -.98208329E+03;
Ba2Sc2O6halfdecW112_F = Ba2Sc2O6halfdecW112_F+Ba2Sc2O6halfdecW112_E0;
Ba2Sc2O6halfdecW112_WF = 1.45;

Ba2Sc2O6halfbareW112_area = 147.7998; 
Ba2Sc2O6halfbareW112_Watoms = 66;
Ba2Sc2O6halfbareW112_F = readmatrix('thermalpropsBa2Sc2O6_bare-W112.txt');
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F(:,2);
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F*convertunit;
Ba2Sc2O6halfbareW112_E0 = -.82364685E+03;
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F+Ba2Sc2O6halfbareW112_E0;
Ba2Sc2O6halfbareW112_gamma = (Ba2Sc2O6halfbareW112_F-Ba2Sc2O6halfbareW112_Watoms*mu_metW)/...
    Ba2Sc2O6halfbareW112_area;

Ba2O4Sc4W112_area = 49.27;
Ba2O4Sc4W112_Watoms = 38;
Ba2O4Sc4W112_Scatoms = 4;
Ba2O4Sc4W112_Oatoms = 4;
Ba2O4Sc4W112_Baatoms = 2;
Ba2O4Sc4W112_F = readmatrix('thermalpropsBa2O4Sc4-W112.txt');
Ba2O4Sc4W112_F = Ba2O4Sc4W112_F(:,2);
Ba2O4Sc4W112_F = Ba2O4Sc4W112_F*convertunit;
Ba2O4Sc4W112_E0 = -.55767079E+03;
Ba2O4Sc4W112_F = Ba2O4Sc4W112_F+Ba2O4Sc4W112_E0;
Ba2O4Sc4W112_WF = 1.92933196;

%ii=1;
%lowerboundshapegammas = zeros(ii,4);
pp=1;
allSEs=zeros(pp,9);


for i = 1:aa
     
    lowestgammas = zeros(cc,3);
    SE_mu_O = plotmu_O(i,:);
    Ba_O = false;
    W_O = false;
    Sc_O = false;
    for j = 1:length(SE_mu_O)
            
            mu_Ba(j) = mu_metBa(i);
            mu_W(j) = mu_metW(i);
            mu_Sc(j) = mu_metSc(i);
            
            %{
            if SE_mu_O(j)>= Bacutoff(i)
                if(~Ba_O) 
                    if exist('r','var')
                    r=[r;[Temp(i),P_torr(j),mu_metBa(i)]];
                    else
                    r=[Temp(i),P_torr(j),mu_metBa(i)];
                    end
                end
                mu_Ba(j) = mu_metBa(i);
                Ba_O = true;
            end
            %}
            
            if SE_mu_O(j)>= Bacutoff(i)
                if(~Ba_O) 
                    if exist('r','var')
                    r=[r;[Temp(i),P_torr(j),Bacutoff(i)]];
                    else
                    r=[Temp(i),P_torr(j),Bacutoff(i)];
                    end
                end
                mu_Ba(j) = mu_BaO(i)- SE_mu_O(j);
                Ba_O = true;
            end
            %}
            if SE_mu_O(j) >= Wcutoff(i)
                if(~W_O) 
                    if exist('s','var')
                    s=[s;[Temp(i),P_torr(j),Wcutoff(i)]];
                    else
                    s=[Temp(i),P_torr(j),Wcutoff(i)];
                    end
                end
                mu_W(j) = mu_WO3(i)- 3*SE_mu_O(j);
                W_O = true;
            end
            %{
            if SE_mu_O(j) >= Sccutoff(i)
                if(~Sc_O) 
                    if exist('z','var')
                    z=[z;[Temp(i),P_torr(j),mu_metSc(i)]];
                    else
                    z=[Temp(i),P_torr(j),mu_metSc(i)];
                    end
                end
                mu_Sc(j) = mu_metSc(i);
                Sc_O = true;
            end
            %}
            
            if SE_mu_O(j) >= Sccutoff(i)
                if(~Sc_O) 
                    if exist('z','var')
                    z=[z;[Temp(i),P_torr(j),Sccutoff(i)]];
                    else
                    z=[Temp(i),P_torr(j),Sccutoff(i)];
                    end
                end
                mu_Sc(j) = (mu_Sc2O3(i)- 3*SE_mu_O(j))/2;
                Sc_O = true;
            end
            %}
 %%%% 001 %%%%
 
            W001surfpert_gamma(j) = (W001surfpert_F(i)-...
                (W001surfpert_Watoms*mu_W(j)))/W001surfpert_area; 
            
            Ba4O8W001_gamma(j) = ((Ba4O8W001_F(i)-(Ba4O8W001_Watoms*mu_W(j)...
                +Ba4O8W001_Baatoms*mu_Ba(j)+Ba4O8W001_Oatoms*SE_mu_O(j)))/Ba4O8W001_area)+(Temp(i)*Ba4O8W001_config);
            
            Ba2Sc2O8W001_gamma(j) = ((Ba2Sc2O8W001_F(i)-(Ba2Sc2O8W001_Watoms*mu_W(j)...
                +Ba2Sc2O8W001_Baatoms*mu_Ba(j)+Ba2Sc2O8W001_Oatoms*SE_mu_O(j)+Ba2Sc2O8W001_Scatoms*mu_Sc(j)))/...
                Ba2Sc2O8W001_area+(Temp(i)*Ba2Sc2O8W001_config));
            
            O8W001_gamma(j) = (O8W001_F(i)-(O8W001_Watoms*mu_W(j)...
                +O8W001_Oatoms*SE_mu_O(j)))/O8W001_area;  
            
            Ba2O8W001_gamma(j) = ((Ba2O8W001_F(i)-(Ba2O8W001_Watoms*mu_W(j)...
                +Ba2O8W001_Baatoms*mu_Ba(j)+Ba2O8W001_Oatoms*SE_mu_O(j)))/Ba2O8W001_area)+(Temp(i)*Ba2O8W001_config);
            
            
            
            W001set = [W001surfpert_gamma(j); Ba4O8W001_gamma(j); ...
                Ba2Sc2O8W001_gamma(j); O8W001_gamma(j); Ba2O8W001_gamma(j)];
            lowestW001 = min(W001set);
            lowestW001pos = find(W001set==min(W001set));
            
            if lowestW001pos == 1
                W001WF = W001surfpert_WF;
            elseif lowestW001pos == 2
                W001WF = Ba4O8W001_WF;
            elseif lowestW001pos == 3
                W001WF = Ba2Sc2O8W001_WF;
            elseif lowestW001pos == 4
                W001WF = O8W001_WF;
            else
                W001WF = Ba2O8W001_WF;
            end
           
                    
%%%% 110 %%%%

        	W110_gamma(j) = (W110_F(i)-(W110_Watoms*mu_W(j)))/W110_area;
            
            Ba2O8W110_gamma(j) = ((Ba2O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Baatoms*mu_Ba(j)+Ba2O8W110_Oatoms*SE_mu_O(j)))/Ba2O8W110_area)+(Temp(i)*Ba2O8W110_config);
                      
            Ba2Sc2O8W110_gamma(j) = ((Ba2Sc2O8W110_F(i)-(Ba2Sc2O8W110_Watoms*mu_W(j)...
                +Ba2Sc2O8W110_Baatoms*mu_Ba(j)+Ba2Sc2O8W110_Oatoms*SE_mu_O(j)...
                +Ba2Sc2O8W110_Scatoms*mu_Sc(j)))/Ba2Sc2O8W110_area)+(Temp(i)*Ba2Sc2O8W110_config);
            
            O8W110_gamma(j) = (O8W110_F(i)-(Ba2O8W110_Watoms*mu_W(j)...
                +Ba2O8W110_Oatoms*SE_mu_O(j)))/Ba2O8W110_area;
       
          	W110set = [W110_gamma(j); Ba2O8W110_gamma(j); ...
                Ba2Sc2O8W110_gamma(j); O8W110_gamma(j)];
            lowestW110 = min(W110set);
            lowestW110pos = find(W110set==min(W110set));
            
            if lowestW110pos == 1
                W110WF = W110_WF;
            elseif lowestW110pos == 2
                W110WF = Ba2O8W110_WF;
            elseif lowestW110pos == 3
                W110WF = Ba2Sc2O8W110_WF;
            else
                W110WF = O8W110_WF;
            end
       
    
%%%% 112 %%%%
            
          W112_gamma(j) = (W112_F(i)-(W112_Watoms*mu_W(j)))/W112_area;
          
          O4W112_gamma(j) = (O4W112_F(i)-(O4W112_Watoms*mu_W(j)...
              +O4W112_Oatoms*SE_mu_O(j)))/O4W112_area;
          
          Sc2W112top_gamma(j) = (Sc2W112top_F(i)-(Sc2W112top_Watoms*mu_W(j)...
              +Sc2W112top_Scatoms*mu_Sc(j)))/Sc2W112top_area;
          
          Sc4W112top_gamma(j) = (Sc4W112top_F(i)-(Sc4W112top_Watoms*mu_W(j)...
              +Sc4W112top_Scatoms*mu_Sc(j)))/Sc4W112top_area;
          
          Sc4triW112_gamma(j) = (Sc4triW112_F(i)-(Sc4triW112_Watoms*mu_W(j)...
              +Sc4triW112_Scatoms*mu_Sc(j)))/Sc4triW112_area;
          
          Ba2O4W112_2x2_gamma(j) = ((Ba2O4W112_2x2_F(i)-(Ba2O4W112_2x2_Watoms*mu_W(j)...
              +Ba2O4W112_2x2_Baatoms*mu_Ba(j)+Ba2O4W112_2x2_Oatoms*SE_mu_O(j)))/Ba2O4W112_2x2_area)+(Temp(i)*Ba2O4W112_2x2_config);
             
          Ba1Sc1O4W112_gamma(j) = ((Ba1Sc1O4W112_F(i)-(Ba1Sc1O4W112_Watoms*mu_W(j)...
              +Ba1Sc1O4W112_Baatoms*mu_Ba(j)+Ba1Sc1O4W112_Oatoms*SE_mu_O(j)+...
              Ba1Sc1O4W112_Scatoms*mu_Sc(j)))/Ba1Sc1O4W112_area)+(Temp(i)*Ba1Sc1O4W112_config);
          
          Ba2Sc2O4W112_gamma(j) = ((Ba2Sc2O4W112_F(i)-(Ba2Sc2O4W112_Watoms*mu_W(j)...
              +Ba2Sc2O4W112_Baatoms*mu_Ba(j)+Ba2Sc2O4W112_Oatoms*SE_mu_O(j)+...
              Ba2Sc2O4W112_Scatoms*mu_Sc(j)))/Ba2Sc2O4W112_area)+(Temp(i)*Ba2Sc2O4W112_config);
          
          Ba2Sc2O6W112_gamma(j) = ((Ba2Sc2O6W112_F(i)-(Ba2Sc2O6W112_Watoms*mu_W(j)...
              +Ba2Sc2O6W112_Baatoms*mu_Ba(j)+Ba2Sc2O6W112_Scatoms*mu_Sc(j)+...
              Ba2Sc2O6W112_Oatoms*SE_mu_O(j)))/Ba2Sc2O6W112_area)+(Temp(i)*Ba2Sc2O6W112_config);
          
          Ba2Sc2O6halfdecW112_gamma(j) = (((Ba2Sc2O6halfdecW112_F(i)-(Ba2Sc2O6halfdecW112_Watoms*mu_W(j)...
              +Ba2Sc2O6halfdecW112_Baatoms*mu_Ba(j)+Ba2Sc2O6halfdecW112_Scatoms*mu_Sc(j)+...
              Ba2Sc2O6halfdecW112_Oatoms*SE_mu_O(j)))/Ba2Sc2O6halfdecW112_area)+(Temp(i)*Ba2Sc2O6halfdecW112_config))...
              -Ba2Sc2O6halfbareW112_gamma(i);
          
          Ba2O4Sc4W112_gamma(j) = ((Ba2O4Sc4W112_F(i)-(Ba2O4Sc4W112_Watoms*mu_W(j)...
              +Ba2O4Sc4W112_Baatoms*mu_Ba(j)+Ba2O4Sc4W112_Scatoms*mu_Sc(j)+...
              Ba2O4Sc4W112_Oatoms*SE_mu_O(j)))/Ba2O4Sc4W112_area);
                        
           
          W112set = [W112_gamma(j); O4W112_gamma(j); Sc2W112top_gamma(j);...
                Sc4W112top_gamma(j); Sc4triW112_gamma(j); Ba2O4W112_2x2_gamma(j); ...
                Ba1Sc1O4W112_gamma(j); Ba2Sc2O4W112_gamma(j); ...
                Ba2Sc2O6W112_gamma(j); Ba2Sc2O6halfdecW112_gamma(j); ...
                Ba2O4Sc4W112_gamma(j)];
          
          lowestW112 = min(W112set);
          lowestW112pos = find(W112set==min(W112set));
          
          
          if lowestW112pos == 1
              W112WF = W112_WF;
          elseif lowestW112pos == 2
              W112WF = O4W112_WF;
          elseif lowestW112pos == 3
              W112WF = Sc2W112top_WF;
          elseif lowestW112pos == 4
              W112WF = Sc4W112top_WF;
          elseif lowestW112pos == 5
              W112WF = Sc4triW112_WF;
          elseif lowestW112pos == 6
              W112WF = Ba2O4W112_2x2_WF;
          elseif lowestW112pos == 7
              W112WF = Ba1Sc1O4W112_WF;
          elseif lowestW112pos == 8
              W112WF = Ba2Sc2O4W112_WF;
          elseif lowestW112pos == 9
              W112WF = Ba2Sc2O6W112_WF;
          elseif lowestW112pos == 10
              W112WF = Ba2Sc2O6halfdecW112_WF;
          else
              W112WF = Ba2O4Sc4W112_WF;
          end
          %}
          
%%%% assemble data %%%%
          delta_W112_W110 = 0.5*((lowestW112*1.1/lowestW110)-lowestW112/(1.1*lowestW110));
          delta_W001_W110 = 0.5*((lowestW001*1.1/lowestW110)-lowestW001/(1.1*lowestW110));
          target001SE = 0.927;
          target112SE = 0.945;
          lowSEbound001 = target001SE-delta_W001_W110;
          highSEbound001 = target001SE+delta_W001_W110;
          lowSEbound112 = target112SE-delta_W112_W110;
          highSEbound112 = target112SE+delta_W112_W110;


          lowestgammas = [lowestW001; lowestW110; lowestW112]; 
          ratio_W001_W110 = lowestW001/lowestW110;
          ratio_W112_W110 = lowestW112/lowestW110;
            
          allSEs(pp,:) = [Temp(i); round(P_torr(j),8,'significant'); ...
              SE_mu_O(j); ratio_W001_W110; 1.0000; ratio_W112_W110; ...
              W001WF; W110WF; W112WF];
          pp = pp+1;
            
            
          if (ratio_W001_W110>lowSEbound001) && (ratio_W001_W110<highSEbound001) ...
                  && (ratio_W112_W110>lowSEbound112) && (ratio_W112_W110<highSEbound112)
                if exist('lowerboundshapegammas','var')
                    lowerboundshapegammas=[lowerboundshapegammas;...
                    [Temp(i), P_torr(j), ratio_W001_W110, ...
                    ratio_W112_W110]];
                else
                    lowerboundshapegammas=[Temp(i), P_torr(j), ratio_W001_W110, ...
                    ratio_W112_W110];
                end
          end
          
          %}
            
    end  
        
end

save('data.mat','allSEs');
save('allSEs.mat','allSEs');


allSEs(:,10) = areafrac(:,1);
allSEs(:,11) = areafrac(:,2);
allSEs(:,12) = areafrac(:,3);

points = [1500 1.00E-40 4.556;1500 5.72E-25 1.215;1500 1.63E-19 1.221;...
    1500 2.85E-15 1.266;1350 2.85E-15 2.653];

figure
contourf(P_torr,Temp,plotmu_O,100,'LineColor','none')
colorbar
h = colorbar;
set(get(h,'label'),'string','O Chemical Potential (eV)','fontsize', 28,'FontName','Tahoma');
xlabel({'P_{O_2} (Torr)'},'fontsize', 28);
ylabel({'Temperature (K)'},'fontsize', 28);
set(gca,'XScale','log','FontSize',28,'FontName','Tahoma');
ax = gca;
ax.LineWidth = 3;
box on
hold on
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
semilogx(lowerboundshapegammas(:,2),lowerboundshapegammas(:,1),'.',...
    'color',[128/255 128/255 128/255],'MarkerSize',20);
hold on
semilogx(r(89:aa,2),r(89:aa,1),'LineWidth',3,'color','magenta') %Ba cutoff
hold on
semilogx(s(57:aa,2),s(57:aa,1),'LineWidth',3,'color','cyan') %W cutoff
hold on
semilogx(z(101:aa,2),z(101:aa,1),'LineWidth',3,'color','black') %Sc cutoff
hold on
semilogx(points(:,2),points(:,1),'.',...
    'color','k','MarkerSize',12);
ax.LineWidth = 3;
box on
xlim([1E-50 1])
set(gcf, 'Position',  [0, 0, 1500, 800]);
saveas(gcf,['muocontour.png']);
%}