L1=csvread("./summaries/network_summaries_TR.csv")(2:end,:);
L2=csvread("./summaries/network_summaries_TR2.csv")(2:end,:);
L3=csvread("./summaries/network_summaries_TR3.csv")(2:end,:);
L4=csvread("./summaries/network_summaries_TR4.csv")(2:end,:);
L5=csvread("./summaries/network_summaries_TR5.csv")(2:end,:);
L6=csvread("./summaries/network_summaries_TR6.csv")(2:end,:);
Cn = [L1;L2;L3;L4;L5;L6];
#Cn = [L1;L2;L4;L5;L6];
#Cn = [L5];
Ns = unique(Cn(:,1));

for i = 1:length(Ns)
  idx = Cn(:,1)==Ns(i);
  Lb(i) = mean(Cn(idx,3));
endfor

pfit = polyfit(log2(Ns),log2(Lb'),1)
R = corrcoef(log2(Ns),log2(Lb'));
R2 = R(1,2).^2
plot(log2(Ns),log2(Lb),'*'); hold on;
plot(log2(Ns),log2(Ns)*pfit(1)+pfit(2),'-k');
hold off;
fs = 16
LegFit = ['Fit: m=',num2str(pfit(1)),', b=',num2str(pfit(2)),', R^2=',num2str(R2)]
hL=legend("Simulation Data",LegFit,'Location','northwest');
set(hL,'FontSize',fs);
xlabel("network size (log2 scale)",'FontSize',fs)
Nsr = 2.^(1:max(log2(Ns)));
xticks(log2(Nsr))
xticklabels(Nsr)
xlim([min(log2(Ns)),max(log2(Ns))])
ylabel("average maximal stable module number (log2 scale)",'FontSize',fs)
yticks([0,1,2])
yticklabels([1,2,4])
ylim([0 4])
