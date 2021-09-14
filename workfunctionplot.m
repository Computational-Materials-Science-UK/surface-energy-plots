%load allSEs.mat

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

points = [1500 1.00E-40 4.556;1500 5.72E-25 1.215;1500 1.63E-19 1.221;...
    1500 2.85E-15 1.266;1350 2.85E-15 2.653];

figure
contourf(P_torr,Temp,effWFlist,1000,'LineColor','none')
colorbar
h = colorbar;
set(get(h,'label'),'string','\phi_{eff} (eV)','FontName','Tahoma','FontSize',28);
xlabel({'P_{O_2} (Torr)'},'FontSize',28);
ylabel({'Temperature (K)'},'FontSize',28);
set(gca,'XScale','log','FontSize',28,'FontName','Tahoma');
ax = gca;
ax.LineWidth = 3;
box on
hold on
semilogx(lowerboundshapegammas(:,2),lowerboundshapegammas(:,1),'.',...
    'color',[128/255 128/255 128/255],'MarkerSize',14);
hold on
semilogx(r(89:aa,2),r(89:aa,1),'LineWidth',3,'color','magenta') %Ba cutoff
hold on
semilogx(s(57:aa,2),s(57:aa,1),'LineWidth',3,'color','cyan') %W cutoff
hold on
semilogx(z(101:aa,2),z(101:aa,1),'LineWidth',3,'color','black') %Sc cutoff
hold on
%semilogx(points(:,2),points(:,1),'.',...
%    'color','k','MarkerSize',12);
ax.LineWidth = 3;
box on
xlim([1E-50 1])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(gcf, 'Position',  [0, 0, 1500, 800]);
saveas(gcf,['workfunctioncontour.png']);
%}


figure
scatter3(allSEs(:,2),allSEs(:,1),allSEs(:,13),[],allSEs(:,13),'filled')
hold on
%scatter3(points(:,2),points(:,1),points(:,3))
%colorbar
xlabel({'P_{O_2} (Torr)'},'fontsize', 28);
ylabel({'Temperature (K)'},'fontsize', 28);
zlabel({'\phi_{eff} (eV)'},'fontsize', 28);
grid off
ax = gca;
ax.LineWidth = 3;
xlim([1E-50 1])
set(gca,'XScale','log','FontSize',28,'FontName','Tahoma');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
%}
