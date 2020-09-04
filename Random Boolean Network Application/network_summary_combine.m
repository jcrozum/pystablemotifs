cut = 50
replace = true


L1=csvread("./summaries/network_summaries_prev.csv")(2:end,:);
L2=csvread("./summaries/network_summaries_prev(1).csv")(2:end,:);
D1=csvread("./summaries/network_summaries_delproj.csv")(2:end,:);
D2=csvread("./summaries/network_summaries_delproj(1).csv")(2:end,:);
L3=csvread("./summaries/network_summaries_pass3.csv")(2:end,:);
L4=csvread("./summaries/network_summaries_pass3(1).csv")(2:end,:);

Cn = unique([L1;D1;L2;D2;L3;L4],'rows');

Ns = unique(Cn(:,1));

#           2     4     8    16    32    64    91   128  181  256  362  512 1024 2048 4096
Ncaps = [1000; 1000; 1000; 1000; 1000; 1000; 1000; 1000; 750; 500; 450; 400; 350; 300; 300];

k=0;
q=0;
for i = 1:length(Ns)
  idx = Cn(:,1)==Ns(i);
  #Lb(i) = mean(Cn(idx,3));
  #Ub(i) = mean(Cn(idx,4));
  mInd(i) = max(Cn(idx,2))+1;

  for j = 1:Ncaps(i)
    jdx = Cn(idx,2)==(j-1);
    if sum(jdx) > 0
      Cj = Cn(idx,:)(jdx,:);
      k = k + 1;
      F(k,:) = [Ns(i) (j-1) max(Cj(:,3)) min(Cj(:,4))];
      if min(Cj(:,4)) - max(Cj(:,3)) > cut
        q=q+1;
        Missed(q,:) = [Ns(i),j-1,min(Cj(:,4)) - max(Cj(:,3))];
      endif
    else
      q=q+1;
      Missed(q,:) = [Ns(i),j-1,NaN];
    endif
  endfor
  if q == 0
    MF(i) = 0;
  else
    MF(i) = sum(Missed(:,1)==Ns(i))/Ncaps(i);
  endif

  idx2 = (F(:,1)==Ns(i)) & ((F(:,4)-F(:,3))<=cut);
  idxp = (F(:,1)==Ns(i)) & (F(:,4)==F(:,3));
  Lb(i) = mean(F(idx2,3));
  Le(i) = std(F(idx2,3))/sqrt(sum(idx2));
  Ub(i) = mean(F(idx2,4));
  Ue(i) = std(F(idx2,4))/sqrt(sum(idx2));
  Mb(i) = mean(F(idx2,4) + F(idx2,3))/2;
  Me(i) = std((F(idx2,4) + F(idx2,3))/2)/sqrt(sum(idx2));
  P(i) = mean(F(idxp,3));
  Pe(i) = std(F(idxp,3))/sqrt(sum(idxp));
  SFL(i) = sum(F(idx2,3) == 1)/(sum(idx2));
  SFU(i) = sum(F(idx2,4) == 1)/(sum(idx2));

  ndx = (F(:,1)==Ns(i)) & ((F(:,4)-F(:,3))>cut);
  ndxL(i) = sum(ndx);
  if replace && sum(ndx) > 0
    top_replace = max(F(idx2,4)) * 1.1;
    bot_replace = 1;
    att_ratio = 1.25; # avg number of attractors per trap space (upper bound)

    if q == 0
      numSkip = 0;
    else
      numSkip = sum(isnan(Missed(Missed(:,1)==Ns(i),3)));
    endif
    Lb(i) = (Lb(i)*sum(idx2) + numSkip*bot_replace + sum(F(ndx,3)))/Ncaps(i);
    Ub(i) = (Lb(i)*sum(idx2) + numSkip*top_replace + sum(min([att_ratio*F(ndx,3),F(ndx,4)]')))/Ncaps(i);
  endif

endfor

fsummary = [Ns';ceil(Ncaps'.*MF);Ncaps';mInd;ndxL]

SFL = SFL';
SFU = SFU';
Lb = Lb';
Ub = Ub';
P = P';
Le = Le';
Ue = Ue';
Pe = Pe';
Mb = Mb';
Me = Me';
#Mb = 0.5 * (Ub + Lb);

Lfit = polyfit(log2(Ns),log2(Lb),1);
LR = corrcoef(log2(Ns),log2(Lb));
LR2 = LR(1,2).^2;
LegL = ['Fit: m=',num2str(Lfit(1)),', b=',num2str(Lfit(2)),', R^2=',num2str(LR2)]
Ufit = polyfit(log2(Ns),log2(Ub),1);
UR = corrcoef(log2(Ns),log2(Ub));
UR2 = UR(1,2).^2;
LegU = ['Fit: m=',num2str(Ufit(1)),', b=',num2str(Ufit(2)),', R^2=',num2str(UR2)]
Mfit = polyfit(log2(Ns),log2(Mb),1);
MR = corrcoef(log2(Ns),log2(Mb));
MR2 = MR(1,2).^2;
LegM = ['Fit: m=',num2str(Mfit(1)),', b=',num2str(Mfit(2)),', R^2=',num2str(MR2)]
Pfit = polyfit(log2(Ns),log2(P),1);
PR = corrcoef(log2(Ns),log2(P));
PR2 = PR(1,2).^2;
LegP = ['Fit: m=',num2str(Pfit(1)),', b=',num2str(Pfit(2)),', R^2=',num2str(PR2)]

fs=12;
errorbar(log2(Ns),log2(Lb),log2(Lb)-log2(Lb-Le),-log2(Lb)+log2(Lb+Le),'~*b'); hold on
errorbar(log2(Ns),log2(Ub),log2(Ub)-log2(Ub-Ue),-log2(Ub)+log2(Ub+Ue),'~*r')
errorbar(log2(Ns),log2(Mb),log2(Mb)-log2(Mb-Me),-log2(Mb)+log2(Mb+Me),'~*g')
errorbar(log2(Ns),log2(P),log2(P)-log2(P-Pe),-log2(P)+log2(P+Pe),'~*k')
plot(log2(Ns),log2(Ns)*Lfit(1)+Lfit(2),'-b')
plot(log2(Ns),log2(Ns)*Ufit(1)+Ufit(2),'-r')
plot(log2(Ns),log2(Ns)*Mfit(1)+Mfit(2),'-g')
plot(log2(Ns),log2(Ns)*Pfit(1)+Pfit(2),'-k')
hL=legend('Lower Bound','Upper Bound', 'Midpoint','Exact Count',LegL,LegU,LegM,LegP,'Location','northwest');
set(hL,'FontSize',fs);
xlabel("network size (log2 scale)",'FontSize',fs)
Nsr = 2.^(1:max(log2(Ns)));
xticks(log2(Nsr))
xticklabels(Nsr)
xlim([min(log2(Ns)),max(log2(Ns))])
ylabel("average attractor number (log2 scale)",'FontSize',fs)
yticks([0,1,2])
yticklabels([1,2,4])
ylim([0 4])
if replace
  title(["Attractor scaling for attractor bounds within ",num2str(cut),"\nwith timeout replacements"])
else
  title(["Attractor scaling for attractor bounds within ",num2str(cut)])
endif
set(gca,'FontSize',fs)
hold off
