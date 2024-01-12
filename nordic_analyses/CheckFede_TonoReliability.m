 fH  = squeeze(mean(Bn(:,[1:8],:),2));
 sH  = squeeze(mean(Bn(:,[1:8],:),2));

 selvvox1 = find(fH(1,:)>2);
 selvvox2 = find(fH(2,:)>2);
 ssvoxfH = intersect(selvvox1,selvvox2);
 ssvoxsH = intersect(selvvox1,selvvox2);


rC1 = (abs(fH(1,ssvoxfH)-fH(2,ssvoxfH)))./(abs(fH(1,ssvoxfH))+abs(fH(2,ssvoxfH)));

rC2 = (abs(sH(1,ssvoxsH)-sH(2,ssvoxsH)))./(abs(sH(1,ssvoxsH))+abs(sH(2,ssvoxsH)));

 fH  = squeeze(mean(An(:,[1:8],:),2));
 sH  = squeeze(mean(An(:,[1:8],:),2));

%  selvvox1 = find(fH(1,:)>0);
%  selvvox2 = find(fH(2,:)>0);
%  ssvoxfH = intersect(selvvox1,selvvox2);
% 
%  selvvox1 = find(sH(1,:)>0);
%  selvvox2 = find(sH(2,:)>0);
%  ssvoxsH = intersect(selvvox1,selvvox2);

rC1_An = (abs(fH(1,ssvoxfH)-fH(2,ssvoxfH)))./(abs(fH(1,ssvoxfH))+abs(fH(2,ssvoxfH)));

rC2_An = (abs(sH(1,ssvoxsH)-sH(2,ssvoxsH)))./(abs(sH(1,ssvoxsH))+abs(sH(2,ssvoxsH)));


indexv1 = find(rC1~=1);
indexv2 = find(rC1_An~=1);
indexv = intersect(indexv1,indexv2);





figure
subplot(1,2,1)
hist(rC1,200)
subplot(1,2,2)
hist(rC2,200)

gvox1 = union(find(rC1>0.8),find(rC1<-0.8));
ggvox1 = ssvoxfH(gvox1);

gvox2 = union(find(rC2>0.8),find(rC2<-0.8));
ggvox2 = ssvoxsH(gvox2);

length(intersect(ggvox1,ggvox2))




fH  = squeeze(mean(AnNS(:,[1:4],:),2));
sH  = squeeze(mean(AnNS(:,[5:8],:),2));

 selvvox1 = find(fH(1,:)>0);
 selvvox2 = find(fH(2,:)>0);
 ssvoxfH = intersect(selvvox1,selvvox2);

 selvvox1 = find(sH(1,:)>0);
 selvvox2 = find(sH(2,:)>0);
 ssvoxsH = intersect(selvvox1,selvvox2);

rC1 = (fH(1,ssvoxfH)-fH(2,ssvoxfH))./(fH(1,ssvoxfH)+fH(2,ssvoxfH));

rC2 = (sH(1,ssvoxsH)-sH(2,ssvoxsH))./(sH(1,ssvoxsH)+sH(2,ssvoxsH));
figure
subplot(1,2,1)
hist(rC1,200)
subplot(1,2,2)
hist(rC2,200)

gvox1 = union(find(rC1>0.8),find(rC1<-0.8));
ggvox1 = ssvoxfH(gvox1);

gvox2 = union(find(rC2>0.8),find(rC2<-0.8));
ggvox2 = ssvoxsH(gvox2);

length(intersect(ggvox1,ggvox2))