load tor
load pol
load tem

plot(tor(:,1),tor(:,2)/sum(tor(:,2)),'-k',pol(:,1),pol(:,2)/sum(pol(:,2)),'-b',tem(:,1),tem(:,2)/sum(tem(:,2)),'-r')