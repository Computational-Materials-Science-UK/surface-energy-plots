clear all;
close all;


Temp = linspace(0,2000,201);
Temp = transpose(Temp);
convertunit = (1/6.022E23)*1000*(1/1.60218E-19); %Unit conversion conversion from Kilojoules to eV
mu_O = linspace(-12,-6.5,111);
mu_O = transpose(mu_O);

%%%%%%% bulk %%%%%%%

W333_numunits = 54;  %Number of molecules per supercell
W333_F = readmatrix('thermalpropsW333.txt');  %Source File
W333_F = W333_F(:,2);  %Look at the second column of the input file
W333_F = W333_F*convertunit; %convert to Kilojoules
W333_E0 = -.69946699E+03; %Ground state energy of the box?
mu_metW = (W333_F+W333_E0)/W333_numunits;

[W333_F, mu_metW]=ReadFromFile(54, 'thermalpropsW333.txt', -.69946699E+03); %Custom Function test case.  If it works I'll replace all readins with it

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
Sc2O3111_E0 = -362.2070125; 
%Remember, this E0 value comes from the Sc2O3 2x2x2 structure.
%The F (phonon) values come from the 1x1x1 structure.
%We simply did not have the computational resources to do a 2x2x2 dfpt
%calculation.
mu_Sc2O3 = (Sc2O3111_F+Sc2O3111_E0)/Sc2O3111_numunits;

Sccutoff = (mu_Sc2O3-2*mu_metSc)/3;

<<<<<<< HEAD

=======
[Al222_F, mu_metAl]=ReadFromFile(32, 'thermalpropsAl222.txt', -.11970109E+03); %Bulk aluminum properties

[Al4O6_F, mu_Al4O6]=ReadFromFile(16, 'thermalpropsAl4O6.txt', -.59797974E+03); %Bulk alumina properties

Alcutoff = (mu_Al4O6-2*mu_metAl)/3; %Setting the Al/Alumina cutoff
>>>>>>> b22f8d397299c2d165925d680d0f2e685a601615

%%%%%%% (001) %%%%%%%

W001surfunpert_area = 80.5; %Units?
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
%W001surfpert_gamma = (W001surfpert_F-W001surfpert_Watoms*mu_metW)/...
   % W001surfpert_area;

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
Ba4O8W001_config = 5.808E-6; %What is this?  What does it do?
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

%%% Shankar's 1x1 001 Slabs %%%
%Bare 001 W slab
W001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
W001double_Watoms = 16;
W001double_F = readmatrix('thermalpropsO-W001.txt');
W001double_F = W001double_F(:,2);
W001double_F = W001double_F*convertunit;
W001double_E0 = -.20231144E+03;
W001double_F = W001double_F+W001double_E0;
%Sc Covered 001 W slab
ScW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
ScW001double_Watoms = 16;
ScW001double_Scatoms = 2;
ScW001double_F = readmatrix('thermalpropsO-W001.txt');
ScW001double_F = ScW001double_F(:,2);
ScW001double_F = ScW001double_F*convertunit;
ScW001double_E0 = -2.1860319E+02;
ScW001double_F = ScW001double_F+ScW001double_E0;
%O-Sc-W 001 slab
OScW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
OScW001double_Watoms = 16;
OScW001double_Scatoms = 2;
OScW001double_Oatoms = 2;
OScW001double_F = readmatrix('thermalpropsO-W001.txt');
OScW001double_F = OScW001double_F(:,2);
OScW001double_F = OScW001double_F*convertunit;
OScW001double_E0 = -2.3987622E+02;
OScW001double_F = OScW001double_F+OScW001double_E0;
%Ba-O-Sc-W 001 slab goes here

%%%%%%% (110) %%%%%%%

W110_area = 56.88; 
W110_Watoms = 36;
W110_F = readmatrix('thermalpropsW110.txt');
W110_F = W110_F(:,2);
W110_F = W110_F*convertunit;
W110_E0 = -.45508870E+03;
W110_F = W110_F+W110_E0;
%W110_gamma = (W110_F-W110_Watoms*mu_metW)/...
%    W110_area;

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

%%% Shankar's 1x1 110 Slabs %%%
%Bare 110 W slab
W110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
W110double_Watoms = 12;
W110double_F = readmatrix('thermalpropsO-W001.txt');
W110double_F = W110double_F(:,2);
W110double_F = W110double_F*convertunit;
W110double_E0 = -1.4977044E+02;
W110double_F = W110double_F+W110double_E0;
%Sc Covered 110 W slab
ScW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
ScW110double_Watoms = 12;
ScW110double_Scatoms = 4;
ScW110double_F = readmatrix('thermalpropsO-W001.txt');
ScW110double_F = ScW110double_F(:,2);
ScW110double_F = ScW110double_F*convertunit;
ScW110double_E0 = -1.7655339E+02;
ScW110double_F = ScW110double_F+ScW110double_E0;
%O-Sc-W 110 slab
OScW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
OScW110double_Watoms = 12;
OScW110double_Scatoms = 4;
OScW110double_Oatoms = 4;
OScW110double_F = readmatrix('thermalpropsO-W001.txt');
OScW110double_F = OScW110double_F(:,2);
OScW110double_F = OScW110double_F*convertunit;
OScW110double_E0 = -2.1577504E+02;
OScW110double_F = OScW110double_F+OScW110double_E0;
%Ba-O-Sc-W 110 slab goes here

%%%%%%% (112) %%%%%%%

W112_area = 49.27; 
W112_Watoms = 34;
W112_F = readmatrix('thermalpropsW112.txt');
W112_F = W112_F(:,2);
W112_F = W112_F*convertunit;
W112_E0 = -.43013015E+03;
W112_F = W112_F+W112_E0;
%W112_gamma = (W112_F-W112_Watoms*mu_metW)/...
%    W112_area;

O4W112_area = 49.27;
O4W112_Watoms = 38;
O4W112_Oatoms = 4;
O4W112_F = readmatrix('thermalpropsO-W112.txt');
O4W112_F = O4W112_F(:,2);
O4W112_F = O4W112_F*convertunit;
O4W112_E0 = -.51392254E+03;
O4W112_F = O4W112_F+O4W112_E0;

Sc2W112top_area = 49.27;
Sc2W112top_Watoms = 38;
Sc2W112top_Scatoms = 2;
Sc2W112top_F = readmatrix('thermalpropsSc2-W112.txt');
Sc2W112top_F = Sc2W112top_F(:,2);
Sc2W112top_F = Sc2W112top_F*convertunit;
Sc2W112top_E0 = -.49288297E+03;
Sc2W112top_F = Sc2W112top_F+Sc2W112top_E0;

Sc4W112top_area = 49.27;
Sc4W112top_Watoms = 38;
Sc4W112top_Scatoms = 4;
Sc4W112top_F = readmatrix('thermalpropsSc4-W112.txt');
Sc4W112top_F = Sc4W112top_F(:,2);
Sc4W112top_F = Sc4W112top_F*convertunit;
Sc4W112top_E0 = -.50486328E+03;
Sc4W112top_F = Sc4W112top_F+Sc4W112top_E0;

Sc4triW112_area = 49.27;
Sc4triW112_Watoms = 38;
Sc4triW112_Scatoms = 4;
Sc4triW112_F = readmatrix('thermalpropsSc4triW112.txt');
Sc4triW112_F = Sc4triW112_F(:,2);
Sc4triW112_F = Sc4triW112_F*convertunit;
Sc4triW112_E0 = -.51116579E+03;
Sc4triW112_F = Sc4triW112_F+Sc4triW112_E0;

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
Ba2Sc2O4W112_config = 0;
Ba2Sc2O4W112_F = readmatrix('thermalpropsBa2Sc2O4-W112.txt');
Ba2Sc2O4W112_F = Ba2Sc2O4W112_F(:,2);
Ba2Sc2O4W112_F = Ba2Sc2O4W112_F*convertunit;
Ba2Sc2O4W112_E0 = -.10807918E+04;
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

Ba2Sc2O6halfbareW112_area = 147.7998; 
Ba2Sc2O6halfbareW112_Watoms = 66;
Ba2Sc2O6halfbareW112_F = readmatrix('thermalpropsBa2Sc2O6_bare-W112.txt');
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F(:,2);
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F*convertunit;
Ba2Sc2O6halfbareW112_E0 = -.82364685E+03;
Ba2Sc2O6halfbareW112_F = Ba2Sc2O6halfbareW112_F+Ba2Sc2O6halfbareW112_E0;
Ba2Sc2O6halfbareW112_gamma = (Ba2Sc2O6halfbareW112_F-Ba2Sc2O6halfbareW112_Watoms*mu_metW)/...
    Ba2Sc2O6halfbareW112_area;

<<<<<<< HEAD
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
=======
%%% Shankar's 1x1 112 Slabs %%%
%Bare 112 W slab
W112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
W112double_Watoms = 22;
W112double_F = readmatrix('thermalpropsO-W001.txt');
W112double_F = W112double_F(:,2);
W112double_F = W112double_F*convertunit;
W112double_E0 = -2.7446276E+02;
W112double_F = W112double_F+W112double_E0;
%Sc Covered 112 W slab
ScW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
ScW112double_Watoms = 22;
ScW112double_Scatoms = 4;
ScW112double_F = readmatrix('thermalpropsO-W001.txt');
ScW112double_F = ScW112double_F(:,2);
ScW112double_F = ScW112double_F*convertunit;
ScW112double_E0 = -3.0382577E+02;
ScW112double_F = ScW112double_F+ScW112double_E0;
%O-Sc-W 112 slab
OScW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
OScW112double_Watoms = 22;
OScW112double_Scatoms = 4;
OScW112double_Oatoms = 4;
OScW112double_F = readmatrix('thermalpropsO-W001.txt');
OScW112double_F = OScW112double_F(:,2);
OScW112double_F = OScW112double_F*convertunit;
OScW112double_E0 = -3.4242413E+02;
OScW112double_F = OScW112double_F+OScW112double_E0;
%Ba-O-Sc-W 110 slab goes here
>>>>>>> b22f8d397299c2d165925d680d0f2e685a601615

%%%%%%% Plotting %%%%%%%


for i = 1
        %%% I think this is the 001 surface energy discovery place %%%
        for j = 1:length(mu_O)
            
            if mu_O(j,:)< Bacutoff(i,:)
                mu_Ba(j,:) = mu_metBa(i,:);
            else 
                %mu_Ba(j,:) = mu_metBa(i,:);
                mu_Ba(j,:) = mu_BaO(i,:)-mu_O(j,:);
            end
            
            if mu_O(j,:) < Wcutoff(i,:)
                mu_W(j,:) = mu_metW(i,:);
            else 
                %mu_W(j,:) = mu_metW(i,:);
                mu_W(j,:) = mu_WO3(i,:)-3*mu_O(j,:);
            end
            
            if mu_O(j,:) < Sccutoff(i,:)
                mu_Sc(j,:) = mu_metSc(i,:);
            else 
                %mu_Sc(j,:) = mu_metSc(i,:);
                mu_Sc(j,:) = (mu_Sc2O3(i,:)- 3*mu_O(j,:))/2;
            end
           
            W001surfpert_gamma(j,:) = (W001surfpert_F(i,:)-...
                (W001surfpert_Watoms*mu_W(j,:)))/W001surfpert_area; 
            
            Ba4O8W001_gamma(j,:) = ((Ba4O8W001_F(i,:)-(Ba4O8W001_Watoms*mu_W(j,:)...
                +Ba4O8W001_Baatoms*mu_Ba(j,:)+Ba4O8W001_Oatoms*mu_O(j,:)))/Ba4O8W001_area)+(Temp(i)*Ba4O8W001_config);
            
            Ba2Sc2O8W001_gamma(j,:) = ((Ba2Sc2O8W001_F(i,:)-(Ba2Sc2O8W001_Watoms*mu_W(j,:)...
                +Ba2Sc2O8W001_Baatoms*mu_Ba(j,:)+Ba2Sc2O8W001_Oatoms*mu_O(j,:)+Ba2Sc2O8W001_Scatoms*mu_Sc(j,:)))/...
                Ba2Sc2O8W001_area+(Temp(i)*Ba2Sc2O8W001_config));
            
            O8W001_gamma(j,:) = (O8W001_F(i,:)-(O8W001_Watoms*mu_W(j,:)...
                +O8W001_Oatoms*mu_O(j,:)))/O8W001_area;  
            
            Ba2O8W001_gamma(j,:) = ((Ba2O8W001_F(i,:)-(Ba2O8W001_Watoms*mu_W(j,:)...
                +Ba2O8W001_Baatoms*mu_Ba(j,:)+Ba2O8W001_Oatoms*mu_O(j,:)))/Ba2O8W001_area)+(Temp(i)*Ba2O8W001_config);

            %%% Shankar's 1x1 001 W slabs %%%
            W001double_gamma(j,:) = (W001double_F(i,:)-(W001double_Watoms*mu_W(j,:)))/W001double_area;

            ScW001double_gamma(j,:) = (ScW001double_F(i,:)-(ScW001double_Watoms*mu_W(j,:)...
                +ScW001double_Scatoms*mu_Sc(j,:)))/ScW001double_area;
            
            OScW001double_gamma(j,:) = ((OScW001double_F(i,:)-(OScW001double_Watoms*mu_W(j,:)...
                +OScW001double_Scatoms*mu_Ba(j,:)+OScW001double_Oatoms*mu_O(j,:)))/OScW001double_area);
           
        end
        
        %%% Here's the 110 Section %%%
        for k = 1:length(mu_O)
            
            if mu_O(k,:) < Bacutoff(i,:)
                mu_Ba(k,:) = mu_metBa(i,:);
            else 
                %mu_Ba(k,:) = mu_metBa(i,:);
                mu_Ba(k,:) = mu_BaO(i,:)- mu_O(k,:);
            end
            
            if mu_O(k,:) < Wcutoff(i,:)
                mu_W(k,:) = mu_metW(i,:);
            else 
                %mu_W(k,:) = mu_metW(i,:);
                mu_W(k,:) = mu_WO3(i,:)- 3*mu_O(k,:);
            end           
      
            if mu_O(j,:) < Sccutoff(i,:)
                mu_Sc(j,:) = mu_metSc(i,:);
            else 
                %mu_Sc(j,:) = mu_metSc(i,:);
                mu_Sc(j,:) = (mu_Sc2O3(i,:)- 3*mu_O(j,:))/2;
            end
            
            W110_gamma(k,:) = (W110_F(i,:)-(W110_Watoms*mu_W(k,:)))/W110_area;
            
            O8W110_gamma(k,:) = (O8W110_F(i,:)-(Ba2O8W110_Watoms*mu_W(k,:)...
                +Ba2O8W110_Oatoms*mu_O(k,:)))/Ba2O8W110_area;
            
            Ba2O8W110_gamma(k,:) = ((Ba2O8W110_F(i,:)-(Ba2O8W110_Watoms*mu_W(k,:)...
                +Ba2O8W110_Baatoms*mu_Ba(k,:)+Ba2O8W110_Oatoms*mu_O(k,:)))/Ba2O8W110_area)+(Temp(i)*Ba2O8W110_config);
                      
            Ba2Sc2O8W110_gamma(k,:) = ((Ba2Sc2O8W110_F(i,:)-(Ba2Sc2O8W110_Watoms*mu_W(k,:)...
                +Ba2Sc2O8W110_Baatoms*mu_Ba(k,:)+Ba2Sc2O8W110_Oatoms*mu_O(k,:)...
                +Ba2Sc2O8W110_Scatoms*mu_Sc(k,:)))/Ba2Sc2O8W110_area)+(Temp(i)*Ba2Sc2O8W110_config);

            %%% Shankar's 1x1 110 W slabs %%%
            W110double_gamma(j,:) = (W110double_F(i,:)-(W110double_Watoms*mu_W(j,:)))/W110double_area;

            ScW110double_gamma(j,:) = (ScW110double_F(i,:)-(ScW110double_Watoms*mu_W(j,:)...
                +ScW110double_Scatoms*mu_Sc(j,:)))/ScW110double_area;
            
            OScW110double_gamma(j,:) = ((OScW110double_F(i,:)-(OScW110double_Watoms*mu_W(j,:)...
                +OScW110double_Scatoms*mu_Ba(j,:)+OScW110double_Oatoms*mu_O(j,:)))/OScW110double_area);
        end
        
        for n = 1:length(mu_O)
            
            if mu_O(n,:) < Bacutoff(i,:)
                mu_Ba(n,:) = mu_metBa(i,:);
            else 
                %mu_Ba(n,:) = mu_metBa(i,:);
                mu_Ba(n,:) = mu_BaO(i,:)- mu_O(n,:);
            end
            
            if mu_O(n,:) < Wcutoff(i,:)
                mu_W(n,:) = mu_metW(i,:);
            else 
                %mu_W(n,:) = mu_metW(i,:);
                mu_W(n,:) = mu_WO3(i,:)- 3*mu_O(n,:);
            end
            
            if mu_O(n,:) < Sccutoff(i,:)
                mu_Sc(n,:) = mu_metSc(i,:);
            else 
                %mu_Sc(n,:) = mu_metSc(i,:);
                mu_Sc(n,:) = (mu_Sc2O3(i,:)- 3*mu_O(n,:))/2;
            end
            
            W112_gamma(n,:) = (W112_F(i,:)-(W112_Watoms*mu_W(n,:)))/W112_area;
            
            O4W112_gamma(n,:) = (O4W112_F(i,:)-(O4W112_Watoms*mu_W(n,:)...
                +O4W112_Oatoms*mu_O(n,:)))/O4W112_area;
            
            Sc2W112top_gamma(n,:) = (Sc2W112top_F(i,:)-(Sc2W112top_Watoms*mu_W(n,:)...
                +Sc2W112top_Scatoms*mu_Sc(n,:)))/Sc2W112top_area;
            
            Sc4W112top_gamma(n,:) = (Sc4W112top_F(i,:)-(Sc4W112top_Watoms*mu_W(n,:)...
                +Sc4W112top_Scatoms*mu_Sc(n,:)))/Sc4W112top_area;
            
            Sc4triW112_gamma(n,:) = (Sc4triW112_F(i,:)-(Sc4triW112_Watoms*mu_W(n,:)...
                +Sc4triW112_Scatoms*mu_Sc(n,:)))/Sc4triW112_area;
            
            Ba2O4W112_2x2_gamma(n,:) = ((Ba2O4W112_2x2_F(i,:)-(Ba2O4W112_2x2_Watoms*mu_W(n,:)...
                +Ba2O4W112_2x2_Baatoms*mu_Ba(n,:)+Ba2O4W112_2x2_Oatoms*mu_O(n,:)))/Ba2O4W112_2x2_area)+(Temp(i)*Ba2O4W112_2x2_config);
               
            Ba1Sc1O4W112_gamma(n,:) = ((Ba1Sc1O4W112_F(i,:)-(Ba1Sc1O4W112_Watoms*mu_W(n,:)...
                +Ba1Sc1O4W112_Baatoms*mu_Ba(n,:)+Ba1Sc1O4W112_Oatoms*mu_O(n,:)+...
                Ba1Sc1O4W112_Scatoms*mu_Sc(n,:)))/Ba1Sc1O4W112_area)+(Temp(i)*Ba1Sc1O4W112_config);
            
            Ba2Sc2O4W112_gamma(n,:) = ((Ba2Sc2O4W112_F(i,:)-(Ba2Sc2O4W112_Watoms*mu_W(n,:)...
                +Ba2Sc2O4W112_Baatoms*mu_Ba(n,:)+Ba2Sc2O4W112_Oatoms*mu_O(n,:)+...
                Ba2Sc2O4W112_Scatoms*mu_Sc(n,:)))/Ba2Sc2O4W112_area)+(Temp(i)*Ba2Sc2O4W112_config);
            
            Ba2Sc2O6W112_gamma(n,:) = ((Ba2Sc2O6W112_F(i,:)-(Ba2Sc2O6W112_Watoms*mu_W(n,:)...
                +Ba2Sc2O6W112_Baatoms*mu_Ba(n,:)+Ba2Sc2O6W112_Scatoms*mu_Sc(n,:)+...
                Ba2Sc2O6W112_Oatoms*mu_O(n,:)))/Ba2Sc2O6W112_area)+(Temp(i)*Ba2Sc2O6W112_config);
            
            Ba2Sc2O6halfdecW112_gamma(n,:) = (((Ba2Sc2O6halfdecW112_F(i,:)-(Ba2Sc2O6halfdecW112_Watoms*mu_W(n,:)...
                +Ba2Sc2O6halfdecW112_Baatoms*mu_Ba(n,:)+Ba2Sc2O6halfdecW112_Scatoms*mu_Sc(n,:)+...
                Ba2Sc2O6halfdecW112_Oatoms*mu_O(n,:)))/Ba2Sc2O6halfdecW112_area)+(Temp(i)*Ba2Sc2O6halfdecW112_config))...
                -Ba2Sc2O6halfbareW112_gamma(i);
            
            Ba2O4Sc4W112_gamma(n,:) = ((Ba2O4Sc4W112_F(i,:)-(Ba2O4Sc4W112_Watoms*mu_W(n,:)...
                +Ba2O4Sc4W112_Baatoms*mu_Ba(n,:)+Ba2O4Sc4W112_Scatoms*mu_Sc(n,:)+...
                Ba2O4Sc4W112_Oatoms*mu_O(n,:)))/Ba2O4Sc4W112_area);

            %%% Shankar's 1x1 112 W slabs %%%
            W112double_gamma(j,:) = (W112double_F(i,:)-(W112double_Watoms*mu_W(j,:)))/W112double_area;

            ScW112double_gamma(j,:) = (ScW112double_F(i,:)-(ScW112double_Watoms*mu_W(j,:)...
                +ScW112double_Scatoms*mu_Sc(j,:)))/ScW112double_area;
            
            OScW112double_gamma(j,:) = ((OScW112double_F(i,:)-(OScW112double_Watoms*mu_W(j,:)...
                +OScW112double_Scatoms*mu_Ba(j,:)+OScW112double_Oatoms*mu_O(j,:)))/OScW112double_area);
                  
        end
        
        %P is plot, L is label
        figure(i);
        hold on
        %Ba-O-covered
        P1 = plot(mu_O(1:111),Ba4O8W001_gamma(1:111),'-.r','LineWidth',4);
        L1 = 'Ba_{0.50}O-top/W(0 0 1)';
        hold on
        P2 = plot(mu_O(1:111),Ba2O8W110_gamma(1:111),'-.b','LineWidth',4);
        L2 = 'Ba_{0.25}O-tri/W(1 1 0)';
        hold on
        P3 = plot(mu_O(1:111),Ba2O4W112_2x2_gamma(1:111),'-.g','LineWidth',4);
        L3 = 'Ba_{0.50}O-top/W(1 1 2)';
        hold on
        
        %O-covered
        P4 = plot(mu_O,O8W001_gamma,':r','LineWidth',4);
        L4 = 'O-top/W(0 0 1)';
        hold on
        P5 = plot(mu_O,O8W110_gamma,':b','LineWidth',4);
        L5 = 'O-tri/W(1 1 0)';
        hold on
        P6 = plot(mu_O,O4W112_gamma,':g','LineWidth',4);
        L6 = 'O-top/W(1 1 2)';
        hold on    
        
        %Bare W
        P7 = plot(mu_O,W001surfpert_gamma, 'r', 'LineWidth',4);
        L7 = 'Bare W(0 0 1)';
        hold on
        P8 = plot(mu_O,W110_gamma, 'b', 'LineWidth',4);
        L8 = 'Bare W(1 1 0)';
        hold on
        P9 = plot(mu_O,W112_gamma, 'g', 'LineWidth',4);
        L9 = 'Bare W(1 1 2)';
        hold on
        
        %Sc-containing
        P10 = plot(mu_O,Ba2Sc2O8W001_gamma,'-^r','LineWidth',2);
        L10 = 'Ba_{0.25}Sc_{0.25}O-top/W(0 0 1)';
        hold on
        P11 = plot(mu_O,Ba2Sc2O8W110_gamma,'-^b','LineWidth',2);
        L11 = 'Ba_{0.25}Sc_{0.25}O-top/W(1 1 0)';
        hold on
        P12 = plot(mu_O,Ba1Sc1O4W112_gamma,'-^g','LineWidth',2);
        L12 = 'Ba_{0.25}Sc_{0.25}O-top/W(1 1 2)';
        hold on
        P13 = plot(mu_O,Ba2Sc2O4W112_gamma,':+g','LineWidth',2);
        L13 = 'Ba_{0.5}Sc_{0.5}O-top/W(1 1 2)';
        hold on
        P14 = plot(mu_O,Ba2Sc2O6W112_gamma,'-sg','LineWidth',2);
        L14 = '(BaSc)_{1/3} O-top/W(1 1 2) symmetric/wrong';
        hold on
        P15 = plot(mu_O,Ba2Sc2O6halfdecW112_gamma,'--g','LineWidth',2.85);
        L15 = '(BaSc)_{1/3} O-top/W(1 1 2)';
        hold on
        P16 = plot(mu_O,Ba2O8W001_gamma,'--r','LineWidth',2.85);
        L16 = 'Ba_{0.25}O-top/W(0 0 1)';
        hold on
        P17 = plot(mu_O,Sc4W112top_gamma,':g','LineWidth',2);
        L17 = 'Sc-top/W(1 1 2)';
        hold on
        P18 = plot(mu_O,Sc2W112top_gamma,'-.g','LineWidth',2);
        L18 = 'Sc_{0.5}-top/W(1 1 2)';
        hold on
        P19 = plot(mu_O,Sc4triW112_gamma,':^g','LineWidth',1.5);
        L19 = 'Sc-tri/W(1 1 2)';
        hold on
        P20 = plot(mu_O,Ba2O4Sc4W112_gamma,':*g','LineWidth',1.5);
        L20 = 'Ba_{0.5}-triOSc-top/W(1 1 2)';
        hold on
        %}
%         
        line([Wcutoff(i,:) Wcutoff(i,:)], [-10 0.7],'Color','c', ...
            'LineWidth', 4, 'LineStyle','-');
%         text([Wcutoff(i,:) + 0.05], 0.05, 'WO_3','fontsize', 28, 'Color', 'c');
%         hold on
%         text([Wcutoff(i,:) - 0.3], 0.05, 'W','fontsize', 28, 'Color', 'c');
%         hold on
%         
        line([Bacutoff(i,:) Bacutoff(i,:)], [-10 0.7],'Color','m', ...
            'LineWidth', 4, 'LineStyle','-');
%         text([Bacutoff(i,:) + 0.05], 0.05, 'BaO','fontsize', 28, 'Color', 'm');
%         hold on
%         text([Bacutoff(i,:) - 0.3], 0.05, 'Ba','fontsize', 28, 'Color', 'm');
%         
         line([Sccutoff(i,:) Sccutoff(i,:)], [-10 0.7],'Color','k', ...
             'LineWidth', 4, 'LineStyle','-');
%         text([Sccutoff(i,:) + 0.05], 0.05, 'Sc_2O_3','fontsize', 24, 'Color', 'k');
%         hold on
%         text([Sccutoff(i,:) - 0.2], 0.05, 'Sc','fontsize', 24, 'Color', 'k');
%         
       
        hold off;
        set(gca,'FontSize',28,'FontName','Tahoma');
        xlabel({'\mu_O (eV)'},'fontsize', 32);
        ylabel({'Surface Energy (eV/Ang.^2)'},'fontsize', 32);
        axis([-12 -7.5 0 0.6])
        ax = gca;
        ax.LineWidth = 4;
        %axes('YColor','none');
        box on;
        %grid on;
        legend([P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12; P13; ...
            P15; P16; P17; P18; P19; P20], L1, L2, L3, L4, L5, L6, L7, L8, L9, ...
            L10, L11, L12, L13, L15, L16, L17, L18, L19, L20, 'fontsize', 18, ...
            'Location','southoutside','NumColumns',5);
        legend off;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-11.75,0.05,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['allcalculatedsurfaces_t',num2str(Temp(i,:)),'.png']);
      
end

function [F, mu] = ReadFromFile(numunits, sourcefile, E0)
    %Inputs: numunits is the number of atoms in the supercell (int), sourcefile is the name of the input file (string), E0 is the base energy of the cell
    convertunit = (1/6.022E23)*1000*(1/1.60218E-19); %Unit conversion conversion from Kilojoules to eV
    F = readmatrix(sourcefile);  %Source File
    F = F(:,2);  %Look at the second column of the input file
    F = F*convertunit;  %Convert from kJ to eV
    mu = (F+E0)/numunits;  %find the chemical reactivity
end