load allSEs.mat

aa = 201;
cc = 166;
Temp = linspace(0,2000,aa);
Temp = transpose(Temp);
P_torr = logspace(-50,0,cc);
P_torr = transpose(P_torr);

pp = 33367;
for i = 1:(pp-1)
    
    W001WF = allSEs(i,7);
    W110WF = allSEs(i,8);
    W112WF = allSEs(i,9);
    W001AF = allSEs(i,10);
    W110AF = allSEs(i,11);
    W112AF = allSEs(i,12);
    beta = 8.61733E-5*allSEs(i,1);
    
    effectiveWF = -1*beta*log((W001AF*exp(-1*...
        W001WF/beta))+(W110AF*exp(-1*W110WF/beta))...
        +(W112AF*exp(-1*W112WF/beta)));
    allSEs(i,13) = effectiveWF;
    
    
end
    
effWFlist = reshape(allSEs(:,13),[166, 201]);
effWFlist = transpose(effWFlist);

figure
contourf(P_torr,Temp,effWFlist,1000,'LineColor','none')
colorbar
h = colorbar;
set(get(h,'label'),'string','Effective Work Function (eV)','fontsize', 32,'FontName','Tahoma');
xlabel({'P_{O_2} (Torr)'},'fontsize', 32);
ylabel({'Temperature (K)'},'fontsize', 32);
set(gca,'XScale','log','FontSize',28,'FontName','Tahoma');
ax = gca;
ax.LineWidth = 3;
box on
hold on
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

