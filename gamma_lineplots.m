clear all;
close all;


Temp = linspace(0,2000,201);
Temp = transpose(Temp);
convertunit = (1/6.022E23)*1000*(1/1.60218E-19);
mu_O = linspace(-12,-6.5,111);
mu_O = transpose(mu_O);

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

%%% Shankar's 1x1 001 Slabs %%%
%Bare 001 W slab
W001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
W001double_Watoms = 16;
W001double_F = readmatrix('thermalprops001-W-1x1_double.txt');
W001double_F = W001double_F(:,2);
W001double_F = W001double_F*convertunit;
W001double_E0 = -.20231144E+03;
W001double_F = W001double_F+W001double_E0;
%Sc Covered 001 W slab
ScW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
ScW001double_Watoms = 16;
ScW001double_Scatoms = 2;
ScW001double_F = readmatrix('thermalprops001-ScW-1x1_double.txt');
ScW001double_F = ScW001double_F(:,2);
ScW001double_F = ScW001double_F*convertunit;
ScW001double_E0 = -2.1860319E+02;
ScW001double_F = ScW001double_F+ScW001double_E0;
%O-Sc-W 001 slab
OScW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
OScW001double_Watoms = 16;
OScW001double_Scatoms = 2;
OScW001double_Oatoms = 2;
OScW001double_F = readmatrix('thermalprops001-OScW-1x1_double.txt');
OScW001double_F = OScW001double_F(:,2);
OScW001double_F = OScW001double_F*convertunit;
OScW001double_E0 = -2.3987622E+02;
OScW001double_F = OScW001double_F+OScW001double_E0;
%Ba-O-Sc-W 001 slab goes here
BaOScW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
BaOScW001double_Watoms = 16;
BaOScW001double_Scatoms = 2;
BaOScW001double_Oatoms = 2;
BaOScW001double_Baatoms = 2;
BaOScW001double_F = readmatrix('thermalprops001-BaOScW-1x1_double.txt');
BaOScW001double_F = BaOScW001double_F(:,2);
BaOScW001double_F = BaOScW001double_F*convertunit;
BaOScW001double_E0 = -2.4095725E+02;
BaOScW001double_F = BaOScW001double_F+BaOScW001double_E0;
%Sc-O-W 001 slab
ScOW001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
ScOW001double_Watoms = 16;
ScOW001double_Scatoms = 2;
ScOW001double_Oatoms = 2;
ScOW001double_F = readmatrix('thermalprops001-ScOW-1x1_double.txt');
ScOW001double_F = ScOW001double_F(:,2);
ScOW001double_F = ScOW001double_F*convertunit;
ScOW001double_E0 = -2.3987622E+02;
ScOW001double_F = ScOW001double_F+ScOW001double_E0;
%2-layer Sc Covered 001 W slab
Sc2W001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
Sc2W001double_Watoms = 16;
Sc2W001double_Scatoms = 4;
Sc2W001double_F = readmatrix('thermalprops001-Sc2W-1x1_double.txt');
Sc2W001double_F = Sc2W001double_F(:,2);
Sc2W001double_F = Sc2W001double_F*convertunit;
Sc2W001double_E0 = -2.3016138E+02;
Sc2W001double_F = Sc2W001double_F+Sc2W001double_E0;
%3-layer Sc Covered 001 W slab
Sc3W001double_area = 20.12316713; %Area in angstroms^2.  This may not be the right units 
Sc3W001double_Watoms = 16;
Sc3W001double_Scatoms = 6;
Sc3W001double_F = readmatrix('thermalprops001-Sc3W-1x1_double.txt');
Sc3W001double_F = Sc3W001double_F(:,2);
Sc3W001double_F = Sc3W001double_F*convertunit;
Sc3W001double_E0 = -2.4305146E+02;
Sc3W001double_F = Sc3W001double_F+Sc3W001double_E0;

%%% Shankar's 2x2 001 Slabs %%%
%Bare 001 W slab
W001double2_area = 4*20.12316713; %Area in angstroms^2.  This may not be the right units 
W001double2_Watoms = 4*16;
W001double2_F = readmatrix('thermalprops001-W-2x2_double.txt');
W001double2_F = W001double2_F(:,2);
W001double2_F = W001double2_F*convertunit;
W001double2_E0 = 4*-.20231144E+03;
W001double2_F = W001double2_F+W001double2_E0;
%Sc Covered 2x2 001 W slab
ScW001double2_area = 4*20.12316713; %Area in angstroms^2.  This may not be the right units 
ScW001double2_Watoms = 4*16;
ScW001double2_Scatoms = 4*2;
ScW001double2_F = readmatrix('thermalprops001-ScW-2x2_double.txt');
ScW001double2_F = ScW001double2_F(:,2);
ScW001double2_F = ScW001double2_F*convertunit;
ScW001double2_E0 = -8.6842424E+02;
ScW001double2_F = ScW001double2_F+ScW001double2_E0;
%O-Sc-W 001 slab
OScW001double2_area = 4*20.12316713; %Area in angstroms^2.  This may not be the right units 
OScW001double2_Watoms = 4*16;
OScW001double2_Scatoms = 4*2;
OScW001double2_Oatoms = 4*2;
OScW001double2_F = readmatrix('thermalprops001-OScW-2x2_double.txt');
OScW001double2_F = OScW001double2_F(:,2);
OScW001double2_F = OScW001double2_F*convertunit;
OScW001double2_E0 = 4*-2.3987622E+02;
OScW001double2_F = OScW001double2_F+OScW001double2_E0;
%Ba-O-Sc-W 001 slab goes here
BaOScW001double2_area = 4*20.12316713; %Area in angstroms^2.  This may not be the right units 
BaOScW001double2_Watoms = 4*16;
BaOScW001double2_Scatoms = 4*2;
BaOScW001double2_Oatoms = 4*2;
BaOScW001double2_Baatoms = 4*2;
BaOScW001double2_F = readmatrix('thermalprops001-BaOScW-2x2_double.txt');
BaOScW001double2_F = BaOScW001double2_F(:,2);
BaOScW001double2_F = BaOScW001double2_F*convertunit;
BaOScW001double2_E0 = 4*-2.4095725E+02;
BaOScW001double2_F = BaOScW001double2_F+BaOScW001double2_E0;
%Sc-O-W 001 slab
ScOW001double2_area = 4*20.12316713; %Area in angstroms^2.  This may not be the right units 
ScOW001double2_Watoms = 4*16;
ScOW001double2_Scatoms = 4*2;
ScOW001double2_Oatoms = 4*2;
ScOW001double2_F = readmatrix('thermalprops001-ScOW-2x2_double.txt');
ScOW001double2_F = ScOW001double2_F(:,2);
ScOW001double2_F = ScOW001double2_F*convertunit;
ScOW001double2_E0 = 4*-2.3987622E+02;
ScOW001double2_F = ScOW001double2_F+ScOW001double2_E0;

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

%%% Shankar's 1x1 110 Slabs %%%
%Bare 110 W slab
W110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
W110double_Watoms = 12;
W110double_F = readmatrix('thermalprops110-W-1x1_double.txt');
W110double_F = W110double_F(:,2);
W110double_F = W110double_F*convertunit;
W110double_E0 = -1.4977044E+02;
W110double_F = W110double_F+W110double_E0;
%Sc Covered 110 W slab
ScW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
ScW110double_Watoms = 12;
ScW110double_Scatoms = 4;
ScW110double_F = readmatrix('thermalprops110-ScW-1x1_double.txt');
ScW110double_F = ScW110double_F(:,2);
ScW110double_F = ScW110double_F*convertunit;
ScW110double_E0 = -1.7655339E+02;
ScW110double_F = ScW110double_F+ScW110double_E0;
%O-Sc-W 110 slab
OScW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
OScW110double_Watoms = 12;
OScW110double_Scatoms = 4;
OScW110double_Oatoms = 4;
OScW110double_F = readmatrix('thermalprops110-OScW-1x1_double.txt');
OScW110double_F = OScW110double_F(:,2);
OScW110double_F = OScW110double_F*convertunit;
OScW110double_E0 = -2.1577504E+02;
OScW110double_F = OScW110double_F+OScW110double_E0;
%Ba-O-Sc-W 110 slab goes here
BaOScW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
BaOScW110double_Watoms = 12;
BaOScW110double_Scatoms = 4;
BaOScW110double_Oatoms = 4;
BaOScW110double_Baatoms=4;
BaOScW110double_F = readmatrix('thermalprops110-BaOScW-1x1_double.txt');
BaOScW110double_F = BaOScW110double_F(:,2);
BaOScW110double_F = BaOScW110double_F*convertunit;
BaOScW110double_E0 = -2.2482369E+02;
BaOScW110double_F = BaOScW110double_F+BaOScW110double_E0;
%Sc-O-W 110 slab
ScOW110double_area = 28.61753414; %Area in angstroms^2.  This may not be the right units 
ScOW110double_Watoms = 12;
ScOW110double_Scatoms = 4;
ScOW110double_Oatoms = 4;
ScOW110double_F = readmatrix('thermalprops110-ScOW-1x1_double.txt');
ScOW110double_F = ScOW110double_F(:,2);
ScOW110double_F = ScOW110double_F*convertunit;
ScOW110double_E0 = -2.1577503E+02;
ScOW110double_F = ScOW110double_F+ScOW110double_E0;
%3O-2Sc-W 110 slab
O3Sc2W110double_area = 28.61753414; %Area in angstroms^2. 
O3Sc2W110double_Watoms = 12;
O3Sc2W110double_Scatoms = 4;
O3Sc2W110double_Oatoms = 6;
O3Sc2W110double_F = readmatrix('thermalprops110-O3Sc2W-1x1_double.txt');
O3Sc2W110double_F = O3Sc2W110double_F(:,2);
O3Sc2W110double_F = O3Sc2W110double_F*convertunit;
O3Sc2W110double_E0 = -.23253415E+03;
O3Sc2W110double_F = O3Sc2W110double_F+O3Sc2W110double_E0;
%Sc Covered 110 W slab
ScW110double2_area = 4*28.61753414; %Area in angstroms^2.  This may not be the right units 
ScW110double2_Watoms = 4*12;
ScW110double2_Scatoms = 4*4;
ScW110double2_F = readmatrix('thermalprops110-ScW-2x2_double.txt');
ScW110double2_F = ScW110double2_F(:,2);
ScW110double2_F = ScW110double2_F*convertunit;
ScW110double2_E0 = -7.0608688E+02;
ScW110double2_F = ScW110double2_F+ScW110double2_E0;

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

%%% Shankar's 1x1 112 Slabs %%%
%Bare 112 W slab
W112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
W112double_Watoms = 22;
W112double_F = readmatrix('thermalprops112-W-1x1_double.txt');
W112double_F = W112double_F(:,2);
W112double_F = W112double_F*convertunit;
W112double_E0 = -2.7446276E+02;
W112double_F = W112double_F+W112double_E0;
%Sc Covered 112 W slab
ScW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
ScW112double_Watoms = 22;
ScW112double_Scatoms = 4;
ScW112double_F = readmatrix('thermalprops112-ScW-1x1_double.txt');
ScW112double_F = ScW112double_F(:,2);
ScW112double_F = ScW112double_F*convertunit;
ScW112double_E0 = -3.0382577E+02;
ScW112double_F = ScW112double_F+ScW112double_E0;
%O-Sc-W 112 slab
OScW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
OScW112double_Watoms = 22;
OScW112double_Scatoms = 4;
OScW112double_Oatoms = 4;
OScW112double_F = readmatrix('thermalprops112-OScW-1x1_double.txt');
OScW112double_F = OScW112double_F(:,2);
OScW112double_F = OScW112double_F*convertunit;
OScW112double_E0 = -3.4242413E+02;
OScW112double_F = OScW112double_F+OScW112double_E0;
%Ba-O-Sc-W 110 slab goes here
BaOScW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
BaOScW112double_Watoms = 22;
BaOScW112double_Scatoms = 4;
BaOScW112double_Oatoms = 4;
BaOScW112double_Baatoms = 4;
BaOScW112double_F = readmatrix('thermalprops112-BaOScW-1x1_double.txt');
BaOScW112double_F = BaOScW112double_F(:,2);
BaOScW112double_F = BaOScW112double_F*convertunit;
BaOScW112double_E0 = -3.5455272E+02;
BaOScW112double_F = BaOScW112double_F+BaOScW112double_E0;
%Sc-O-W 112 slab
ScOW112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
ScOW112double_Watoms = 22;
ScOW112double_Scatoms = 4;
ScOW112double_Oatoms = 4;
ScOW112double_F = readmatrix('thermalprops112-ScOW-1x1_double.txt');
ScOW112double_F = ScOW112double_F(:,2);
ScOW112double_F = ScOW112double_F*convertunit;
ScOW112double_E0 = -3.4375551E+02;
ScOW112double_F = ScOW112double_F+ScOW112double_E0;
%O3-Sc2-W 112 slab
O3Sc2W112double_area = 50.11848623; %Area in angstroms^2.  This may not be the right units 
O3Sc2W112double_Watoms = 22;
O3Sc2W112double_Scatoms = 4;
O3Sc2W112double_Oatoms = 6;
O3Sc2W112double_F = readmatrix('thermalprops112-O3Sc2W-1x1_double.txt');
O3Sc2W112double_F = O3Sc2W112double_F(:,2);
O3Sc2W112double_F = O3Sc2W112double_F*convertunit;
O3Sc2W112double_E0 = -.36175427E+03;
O3Sc2W112double_F = O3Sc2W112double_F+O3Sc2W112double_E0;
%Sc Covered 112 W slab
ScW112double2_area = 4*50.11848623; %Area in angstroms^2.  This may not be the right units 
ScW112double2_Watoms = 4*22;
ScW112double2_Scatoms = 4*4;
ScW112double2_F = readmatrix('thermalprops112-ScW-2x2_double.txt');
ScW112double2_F = ScW112double2_F(:,2);
ScW112double2_F = ScW112double2_F*convertunit;
ScW112double2_E0 = -1.2154063E+03;
ScW112double2_F = ScW112double2_F+ScW112double2_E0;
Ba2O4Sc4W112_WF = 1.92933196;

%%%%%%% Plotting %%%%%%%


for i = 1:length(Temp)  % Loop through temperatures.  Turn off for debuging
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
                +OScW001double_Scatoms*mu_Sc(j,:)+OScW001double_Oatoms*mu_O(j,:)))/OScW001double_area);
            
            BaOScW001double_gamma(j,:) = ((BaOScW001double_F(i,:)-(BaOScW001double_Watoms*mu_W(j,:)...
                +BaOScW001double_Scatoms*mu_Sc(j,:)+BaOScW001double_Oatoms*mu_O(j,:)+BaOScW001double_Baatoms*mu_Ba(j,:)))/BaOScW001double_area);
            
            ScOW001double_gamma(j,:) = ((ScOW001double_F(i,:)-(ScOW001double_Watoms*mu_W(j,:)...
                +ScOW001double_Scatoms*mu_Sc(j,:)+ScOW001double_Oatoms*mu_O(j,:)))/ScOW001double_area);
            
            Sc2W001double_gamma(j,:) = (Sc2W001double_F(i,:)-(Sc2W001double_Watoms*mu_W(j,:)...
                +Sc2W001double_Scatoms*mu_Sc(j,:)))/Sc2W001double_area;
            
            Sc3W001double_gamma(j,:) = (Sc3W001double_F(i,:)-(Sc3W001double_Watoms*mu_W(j,:)...
                +Sc3W001double_Scatoms*mu_Sc(j,:)))/Sc3W001double_area;
            
            %Shankar's 2x2 001 W Slabs
            W001double2_gamma(j,:) = (W001double2_F(i,:)-(W001double2_Watoms*mu_W(j,:)))/W001double2_area;
            
            ScW001double2_gamma(j,:) = (ScW001double2_F(i,:)-(ScW001double2_Watoms*mu_W(j,:)...
                +ScW001double2_Scatoms*mu_Sc(j,:)))/ScW001double2_area;
            
            OScW001double2_gamma(j,:) = ((OScW001double2_F(i,:)-(OScW001double2_Watoms*mu_W(j,:)...
                +OScW001double2_Scatoms*mu_Sc(j,:)+OScW001double2_Oatoms*mu_O(j,:)))/OScW001double2_area);
            
            BaOScW001double2_gamma(j,:) = ((BaOScW001double2_F(i,:)-(BaOScW001double2_Watoms*mu_W(j,:)...
                +BaOScW001double2_Scatoms*mu_Sc(j,:)+BaOScW001double2_Oatoms*mu_O(j,:)+BaOScW001double2_Baatoms*mu_Ba(j,:)))/BaOScW001double2_area);
            
            ScOW001double2_gamma(j,:) = ((ScOW001double2_F(i,:)-(ScOW001double2_Watoms*mu_W(j,:)...
                +ScOW001double2_Scatoms*mu_Sc(j,:)+ScOW001double2_Oatoms*mu_O(j,:)))/ScOW001double2_area);
           
            %%%%%%%%%%%%%%%% Shankar's Other Calulations %%%%%%%%%%%%%%%%
            ScfromScW001double(j,:) = (ScW001double_F(i,:)-W001double_F(i,:))/ScW001double_Scatoms;
            ScfromSc2W001double(j,:) = (Sc2W001double_F(i,:)-ScW001double_F(i,:))/ScW001double_Scatoms;
            ScfromSc3W001double(j,:) = (Sc3W001double_F(i,:)-Sc2W001double_F(i,:))/ScW001double_Scatoms;
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
      
            if mu_O(k,:) < Sccutoff(i,:)
                mu_Sc(k,:) = mu_metSc(i,:);
            else 
                %mu_Sc(j,:) = mu_metSc(i,:);
                mu_Sc(k,:) = (mu_Sc2O3(i,:)- 3*mu_O(k,:))/2;
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
            W110double_gamma(k,:) = (W110double_F(i,:)-(W110double_Watoms*mu_W(k,:)))/W110double_area;

            ScW110double_gamma(k,:) = (ScW110double_F(i,:)-(ScW110double_Watoms*mu_W(k,:)...
                +ScW110double_Scatoms*mu_Sc(k,:)))/ScW110double_area;
            
            OScW110double_gamma(k,:) = ((OScW110double_F(i,:)-(OScW110double_Watoms*mu_W(k,:)...
                +OScW110double_Scatoms*mu_Sc(k,:)+OScW110double_Oatoms*mu_O(k,:)))/OScW110double_area);
            
            BaOScW110double_gamma(k,:) = ((BaOScW110double_F(i,:)-(BaOScW110double_Watoms*mu_W(k,:)...
                +BaOScW110double_Scatoms*mu_Sc(k,:)+BaOScW110double_Oatoms*mu_O(k,:)+BaOScW110double_Baatoms*mu_Ba(k,:)))/BaOScW110double_area);
            
            ScOW110double_gamma(k,:) = ((ScOW110double_F(i,:)-(ScOW110double_Watoms*mu_W(k,:)...
                +ScOW110double_Scatoms*mu_Sc(k,:)+ScOW110double_Oatoms*mu_O(k,:)))/ScOW110double_area);
            
            O3Sc2W110double_gamma(k,:) = ((O3Sc2W110double_F(i,:)-(O3Sc2W110double_Watoms*mu_W(k,:)...
                +O3Sc2W110double_Scatoms*mu_Sc(k,:)+O3Sc2W110double_Oatoms*mu_O(k,:)))/O3Sc2W110double_area);
            
            %Shankar's 2x2 110 W slabs
            ScW110double2_gamma(k,:) = (ScW110double2_F(i,:)-(ScW110double2_Watoms*mu_W(k,:)...
                +ScW110double2_Scatoms*mu_Sc(k,:)))/ScW110double2_area;
            
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
            W112double_gamma(n,:) = (W112double_F(i,:)-(W112double_Watoms*mu_W(n,:)))/W112double_area;

            ScW112double_gamma(n,:) = (ScW112double_F(i,:)-(ScW112double_Watoms*mu_W(n,:)...
                +ScW112double_Scatoms*mu_Sc(n,:)))/ScW112double_area;
            
            OScW112double_gamma(n,:) = ((OScW112double_F(i,:)-(OScW112double_Watoms*mu_W(n,:)...
                +OScW112double_Scatoms*mu_Sc(n,:)+OScW112double_Oatoms*mu_O(n,:)))/OScW112double_area);
            
            BaOScW112double_gamma(n,:) = ((BaOScW112double_F(i,:)-(BaOScW112double_Watoms*mu_W(n,:)...
                +BaOScW112double_Scatoms*mu_Sc(n,:)+BaOScW112double_Oatoms*mu_O(n,:)+BaOScW112double_Baatoms*mu_Ba(n,:)))/BaOScW112double_area);
            
            ScOW112double_gamma(n,:) = ((ScOW112double_F(i,:)-(ScOW112double_Watoms*mu_W(n,:)...
                +ScOW112double_Scatoms*mu_Sc(n,:)+ScOW112double_Oatoms*mu_O(n,:)))/ScOW112double_area);
            
            O3Sc2W112double_gamma(n,:) = ((O3Sc2W112double_F(i,:)-(O3Sc2W112double_Watoms*mu_W(n,:)...
                +O3Sc2W112double_Scatoms*mu_Sc(n,:)+O3Sc2W112double_Oatoms*mu_O(n,:)))/O3Sc2W112double_area);
            
            %Shanar's 2x2 112 W slabs
            ScW112double2_gamma(n,:) = (ScW112double2_F(i,:)-(ScW112double2_Watoms*mu_W(n,:)...
                +ScW112double2_Scatoms*mu_Sc(n,:)))/ScW112double2_area;
                  
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
        %movefile(['allcalculatedsurfaces_t',num2str(Temp(i,:)),'.png'], 'Mujan_Images');


        %Shankar's O-Sc-W series slabs
        figure(i+201);
        hold on
        P21 = plot(mu_O,W001double_gamma,'-r','LineWidth',2);
        L21 = 'Bare W(0 0 1)';
        hold on
        P22 = plot(mu_O,ScW001double_gamma,'--r','LineWidth',2);
        L22 = 'Sc-Covered W(0 0 1)';
        hold on
        P23 = plot(mu_O,OScW001double_gamma,':r','LineWidth',2);
        L23 = 'Layered O-Sc-W(0 0 1)';
        hold on
        P24 = plot(mu_O,BaOScW001double_gamma,'-or','LineWidth',2);
        L24 = 'Layered Ba-O-Sc-W(0 0 1)';
        hold on
        P25 = plot(mu_O,ScOW001double_gamma,'-.r','LineWidth',2);
        L25 = 'Layered Sc-O-W(0 0 1)';
        hold on
        P26 = plot(mu_O,W110double_gamma,'-g','LineWidth',2);
        L26 = 'Bare W(0 0 1)';
        hold on
        P27 = plot(mu_O,ScW110double_gamma,'--g','LineWidth',2);
        L27 = 'Sc-Covered W(1 1 0)';
        hold on
        P28 = plot(mu_O,OScW110double_gamma,':g','LineWidth',2);
        L28 = 'Layered O-Sc-W(1 1 0)';
        hold on
        P29 = plot(mu_O,BaOScW110double_gamma,'-og','LineWidth',2);
        L29 = 'Layered Ba-O-Sc-W(1 1 0)';
        hold on
        P30 = plot(mu_O,ScOW110double_gamma,'-.g','LineWidth',2);
        L30 = 'Layered Sc-O-W(1 1 0)';
        hold on
        P31 = plot(mu_O,O3Sc2W110double_gamma,'+g','LineWidth',2);
        L31 = 'Layered O3-Sc2-W(1 1 0)';
        hold on
        P32 = plot(mu_O,W112double_gamma,'-b','LineWidth',2);
        L32 = 'Bare W(1 1 2)';
        hold on
        P33 = plot(mu_O,ScW112double_gamma,'--b','LineWidth',2);
        L33 = 'Sc-Covered W(1 1 2)';
        hold on
        P34 = plot(mu_O,OScW112double_gamma,':b','LineWidth',2);
        L34 = 'Layered O-Sc-W(1 1 2)';
        hold on
        P35 = plot(mu_O,BaOScW112double_gamma,'-ob','LineWidth',2);
        L35 = 'Layered Ba-O-Sc-W(1 1 2)';
        hold on
        P36 = plot(mu_O,ScOW112double_gamma,'-.b','LineWidth',2);
        L36 = 'Layered Sc-O-W(1 1 2)';
        hold on
        P37 = plot(mu_O,O3Sc2W112double_gamma,'+b','LineWidth',2);
        L37 = 'Layered O3-Sc2-W(1 1 2)';
        hold on

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
        legend([P21; P22; P23; P24; P25; P26; P27; P28; P29; P30; P31; P32; P33; P34; P35; P36; P37], L21, L22, L23, L24, L25, L26, L27, L28, L29, L30, L31, L32, L33, L34, L35, L36, L37, 'fontsize', 18, ...
            'Location','southoutside');
        legend off;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-11.75,0.05,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['Shankarssurfaces_t',num2str(Temp(i,:)),'.png']);
        %movefile(['Shankarssurfaces_t',num2str(Temp(i,:)),'.png'], 'Shankar_Images');
        
        %%%%%%%%%%%Shankar's 1x1 vs. 2x2 Phonon Comparison%%%%%%%%%%%%%
        figure(i+402);
        hold on
        P40 = plot(mu_O,W001double_gamma,'--r','LineWidth',2);
        L40 = 'W(0 0 1)';
        hold on
        P41 = plot(mu_O,W001double2_gamma,':r','LineWidth',2);
        L41 = '2x2 W(0 0 1)';
        hold on
        P42 = plot(mu_O,ScW001double_gamma,'--g','LineWidth',2);
        L42 = 'Sc-Covered W(1 1 0)';
        hold on
        P43 = plot(mu_O,ScW001double2_gamma,':g','LineWidth',2);
        L43 = 'Sc-Covered 2x2 W(1 1 0)';
        hold on
        P44 = plot(mu_O,OScW001double_gamma,'--b','LineWidth',2);
        L44 = 'O-Sc-Covered W(1 1 2)';
        hold on
        P45 = plot(mu_O,OScW001double2_gamma,'ob','LineWidth',2);
        L45 = 'O-Sc-Covered 2x2 W(1 1 2)';
        P46 = plot(mu_O,BaOScW001double_gamma,'--c','LineWidth',2);
        L46 = 'Ba-O-Sc-Covered W(1 1 2)';
        hold on
        P47 = plot(mu_O,BaOScW001double2_gamma,':c','LineWidth',2);
        L47 = 'Ba-O-Sc-Covered 2x2 W(1 1 2)';
        P48 = plot(mu_O,ScOW001double_gamma,'--m','LineWidth',2);
        L48 = 'Sc-O-Covered W(1 1 2)';
        hold on
        P49 = plot(mu_O,ScOW001double2_gamma,':m','LineWidth',2);
        L49 = 'Sc-O-Covered 2x2 W(1 1 2)';
        hold on
        

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
        legend([P40; P41; P42; P43; P44; P45; P46; P47; P48; P49], L40, L41, L42, L43, L44, L45, L46, L47, L48, L49, 'fontsize', 18, ...
            'Location','southoutside');
        legend off;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-11.75,0.05,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['Shankars2x2Comparison_t',num2str(Temp(i,:)),'.png']);
        %movefile(['Shankars2x2Comparison_t',num2str(Temp(i,:)),'.png'], 'Shankar_Images');

      %%%%%%%%%%%Shankar's Layered Scandium Phonon Comparison%%%%%%%%%%%%%
        figure(i+603);
        hold on
        P50 = plot(mu_O,ScW001double_gamma,'--r','LineWidth',2);
        L50 = 'Single Layer';
        hold on
        P51 = plot(mu_O,Sc2W001double_gamma,'--g','LineWidth',2);
        L51 = 'Double Layer';
        hold on
        P52 = plot(mu_O,Sc3W001double_gamma,'--b','LineWidth',2);
        L52 = 'Triple Layer';
        hold on
        

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
        legend([P50; P51; P52], L50, L51, L52, 'fontsize', 18, ...
            'Location','southoutside');
        legend off;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-11.75,0.05,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['ShankarsScLayersComparison_t',num2str(Temp(i,:)),'.png']);
        %movefile(['Shankars2x2Comparison_t',num2str(Temp(i,:)),'.png'], 'Shankar_Images');
        
        %%%%%%%%%%%Shankar's Layered Scandium Per-Sc Atom Energy Comparison%%%%%%%%%%%%%
        figure(i+804);
        hold on
        P53 = plot(mu_O,ScfromScW001double,'--r','LineWidth',2);
        L53 = 'Single Layer';
        hold on
        P54 = plot(mu_O,ScfromSc2W001double,'--g','LineWidth',2);
        L54 = 'Double Layer';
        hold on
        P55 = plot(mu_O,ScfromSc3W001double,'--b','LineWidth',2);
        L55 = 'Triple Layer';
        hold on
        P56 = plot(mu_O, mu_Sc,'--k','LineWidth',2);
        L56 = 'Metal/Oxide Bulk';
        hold on

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
        ylabel({'Sc Energy (eV)'},'fontsize', 32);
        %axis([-12 -7.5 0 0.6])
        ax = gca;
        ax.LineWidth = 4;
        %axes('YColor','none');
        box on;
        %grid on;
        legend([P53; P54; P55; P56], L53, L54, L55, L56, 'fontsize', 18, ...
            'Location','southoutside');
        legend off;
        txt = {['T = ',num2str(Temp(i,:)),' K']};
        text(-11.75,0.05,txt,'fontsize', 35);
        set(gcf, 'Position',  [0, 0, 1500, 800]);
        saveas(gcf,['ShankarsScAtomEnergyComparison_t',num2str(Temp(i,:)),'.png']);
        %movefile(['ShankarsScAtomEnergyComparison_t',num2str(Temp(i,:)),'.png'], 'Shankar_Images');

end

function [F, mu] = ReadFromFile(numunits, sourcefile, E0)
    %Inputs: numunits is the number of atoms in the supercell (int), sourcefile is the name of the input file (string), E0 is the base energy of the cell
    convertunit = (1/6.022E23)*1000*(1/1.60218E-19); %Unit conversion conversion from Kilojoules to eV
    F = readmatrix(sourcefile);  %Source File
    F = F(:,2);  %Look at the second column of the input file
    F = F*convertunit;  %Convert from kJ to eV
    mu = (F+E0)/numunits;  %find the chemical reactivity
end