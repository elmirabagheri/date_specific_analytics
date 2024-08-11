clear all
close all
clc
%%load data
speed = xlsread('speed.xlsx');
time =  xlsread('time.xlsx');
%% make microtrips
it=[];
fp=1;
dd=0;
plot(time,speed);
xlabel('Time(s)','fontsize',20,'fontname','Times New Roman');
ylabel('Speed(m/s)','fontsize',20,'fontname','Times New Roman');
 for i=1:length(speed(:,1))   
 if speed(i,1)<=2
 it=[it i];
 end
 end
for j=2:length(it(1,:))
 if it(1,j)~=it(1,(j-1))+1
 fp=[fp it(j)];
 end
end
lp=fp-1;
for k=1:length(lp(1,:))-1
 %tedade microtripha=toole lp-1
 B= time([fp(k):lp(k+1)],[1]);
 D= speed([fp(k):lp(k+1)],[1]);
 d=lp(k+1)-fp(k)+1;
 dd=[dd d];
 %%tole microtrip ha
tmt([1:d],[k])=B;
vmt([1:d],[k])=D;
 for l=1:d
 if vmt(l,k)<=0
  c(l,k)=l;
 end
% max1=max(c);
% tidle(k)=tmt(max1(k),k)-tmt(1,k);
 tt(k)=tmt(d,k)-tmt(1,k);
 end
vmean(k,[1,2])=[mean(vmt([1:d],[k])),k];
% idleperc=tidle./tt*100;
 %scatter(vmean(1,k),idleperc(1,k),'k','*')
 % xlabel('average velocity'); ylabel('idle time percentage');
 % hold on
 
end
%figure
%  xlabel('average velocity'); ylabel('idle time percentage');
%  hold on
% figure
% plot (time,speed)
% xlabel('Time(s)','fontsize',20,'fontname','Times New Roman');
% ylabel('Speed(m/s)','fontsize',20,'fontname','Times New Roman');
% dd(1)=[];

plot (time,speed)
xlabel('Time(s)','fontsize',18,'fontname','Times New Roman');
ylabel('Speed(m/s)','fontsize',18,'fontname','Times New Roman');
dd(1)=[];

%% kmeans clustering
%Put PCA Results as A and B
vmean2=vmean(:,1);
A = xlsread('pcc_drive.xlsx');
B = xlsread('pcc_acceleration.xlsx');
kcluster=[A,B];
opts = statset('Display','final');
[idx,C] = kmeans(kcluster,4,'Distance','cityblock',...
    'Replicates',20,'Options',opts);
figure
plot(kcluster(:,1),kcluster(:,2),'r*','linewidth',1)
ylabel(' positive acceleration(m/s^2)','fontsize',18,'fontname','Times New Roman');
xlabel('Drive time percentage(%)','fontsize',18,'fontname','Times New Roman');
title('Stop time percentage-Deceleration time percentage','fontsize',18,'fontname','Times New Roman')
figure;
for i=1:4
    center(i,:)=[C(i,:),i];
end
center=sortrows(center,2);
plot(kcluster(idx==center(1,3),1),kcluster(idx==center(1,3),2),'yh');
xlabel('Stop time percentage(%)','fontsize',18,'fontname','Times New Roman');
ylabel('Cruise time percentage(%)','fontsize',18,'fontname','Times New Roman');
hold on
plot(kcluster(idx==center(2,3),1),kcluster(idx==center(2,3),2),'g+');
plot(kcluster(idx==center(3,3),1),kcluster(idx==center(3,3),2),'bo');
plot(kcluster(idx==center(4,3),1),kcluster(idx==center(4,3),2),'r*');
plot(C(:,1),C(:,2),'kx','linewidth',2,'markersize',25);
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
       'Location','NW')
  % title ('Cluster Assignments and Centroids','fontsize',18,'fontname','Times New Roman');  
hold off
%% speed acceleration matrix
acc(1)=0;
for i=1:length(speed)-1
 acc(i+1)=(speed(i+1)-speed(i))/(time(i+1)-time(i));
end
for i=1:length(acc)-1
 diffacc(i+1)=(acc(i+1)-acc(i))/(time(i+1)-time(i));
end
figure;
plot(speed,acc,'k*','markersize',5)
xlabel('Speed(m/s)','fontsize',18,'fontname','Times New Roman');
ylabel('Acceleration(m/s^2)','fontsize',18,'fontname','Times New Roman');
figure;
plot(time,acc)
xlabel('Time(s)','fontsize',18,'fontname','Times New Roman');
ylabel('Acceleration(m/s^2)','fontsize',18,'fontname','Times New Roman');
%% microtrip choosing 
aa=[kcluster,vmean(:,2),tt'];
clusterr1=find(idx==1);
clusterr2=find(idx==2);
clusterr3=find(idx==3);
clusterr4=find(idx==4);

 c1(:,:)=aa(idx==center(1,3),:); %% c1= tedade khoshe1
 c2(:,:)=aa(idx==center(2,3),:); %% c2= tedade khoshe2
 c3(:,:)=aa(idx==center(3,3),:); %% c3= tedade khoshe3
 c4(:,:)=aa(idx==center(4,3),:); %% c4= tedade khoshe4
 tt1=sum(c1(:,4));
 tt2=sum(c2(:,4));
 tt3=sum(c3(:,4));
 tt4=sum(c4(:,4));
 ttt=tt1+tt2+tt3+tt4;
 f1=tt1*1200/ttt;
 f2=tt2*1200/ttt;
 f3=tt3*1200/ttt;
 f4=tt4*1200/ttt;
 e1=tt1*100/ttt;
 e2=tt2*100/ttt;
 e3=tt3*100/ttt;
 e4=tt4*100/ttt;
%%
 for i=1:4
pcenter(i,1)=sqrt(center(i,1).^2+center(i,2).^2);
 pcenter(i,2)=center(i,3);
 end
 for i=1:length(c1(:,1))
 pc1(i,1)=sqrt(c1(i,1).^2+c1(i,2).^2);
 end
 pc11(:,1)=abs(pc1(:,1)-pcenter(1,1));
 for i=1:length(pc11(:,1))
 pc11(i,2)=i;
 end
 ct1=0;
 ct2=0;
 ct3=0;
 ct4=0;
%pc11=sortrows(pc11,1);
 i=0;
 while ct1<=f1 
i=i+1
xx1=pc11(i,2);
cycle1(i,:)=c1(xx1,:);
ct1=ct1+cycle1(i,4);
if i==100
    break
end
 end
 for iii1=1:dd(1,c1(xx1,3))
 if tmt(iii1,c1(xx1,3))-tmt(1,c1(xx1,3))>=f1-(ct1-cycle1(i,3)) 
 break
 end
 end
 

 for i =1:length(c2(:,1))
 pc2(i,1)=sqrt(c2(i,1).^2+c2(i,2).^2);
 end
 pc22(:,1)=abs(pc2(:,1)-pcenter(3,1));
 for i=1:length(pc22(:,1))
 pc22(i,2)=i;
 end
pc22=sortrows(pc22,1);
 
 i=0;

 while ct2<=f2
 i=i+1;
 xx2=pc22(i,2);
 cycle2(i,:)=c2(xx2,:);
 ct2=ct2+cycle2(i,4);
 
 end

 for iii2=1:dd(1,c2(xx2,3))
 if  tmt(iii2,c2(xx2,3))-tmt(1,c2(xx2,3))>=f2-(ct2-cycle2(i,4))  
 break
 end
 end
 
 for i =1:length(c3(:,1))
 pc3(i,1)=sqrt(c3(i,1).^2+c3(i,2).^2);
 end
 pc33(:,1)=abs(pc3(:,1)-pcenter(3,1));
 for i=1:length(pc33(:,1))
 pc33(i,2)=i;
 end
 pc33=sortrows(pc33,1);
 
 i=0;
 while ct3<=f3
 i=i+1;
 xx3=pc33(i,2);
 cycle3(i,:)=c3(xx3,:);
 ct3=ct3+cycle3(i,4);
 end
 for iii3=1:dd(1,c3(xx3,3))
 if tmt(iii3,c3(xx3,3))-tmt(1,c3(xx3,3))>=f3-(ct3-cycle3(i,4)) 
 break
 end
 end
 for i =1:length(c4(:,1))
 pc4(i,1)=sqrt(c4(i,1).^2+c4(i,2).^2);
 end
 pc44(:,1)=abs(pc4(:,1)-pcenter(4,1));
 for i=1:length(pc44(:,1))
 pc44(i,2)=i;
end
%pc44=sortrows(pc44,1);
 
 i=0;
 while ct4<=f4
 i=i+1;
 xx4=pc44(i,2);
 cycle4(i,:)=c4(xx4,:);
 ct4=ct4+cycle4(i,4);
 
 end
 for iii4=1:dd(1,c4(xx4,3))
 if tmt(iii4,c4(xx4,3))-tmt(1,c4(xx4,3))>=f4-(ct4-cycle4(i,4)) 
 break
 end
 end
 
cycle=[cycle1;cycle2;cycle3;cycle4];
sortcycle=sortrows(cycle,4);

AA=vmt(:,cycle(:,3));


vcycle=0;
for i=1:length(AA(1,:))
 clear vcycle1
 
 if cycle(i,3)==cycle1(end,3)
 vcycle1(:,1)=AA(1:iii1,i); 
 else if cycle(i,3)==cycle2(end,3)
 vcycle1(:,1)=AA(1:iii2,i);
 else if cycle(i,3)==cycle3(end,3)
 vcycle1(:,1)=AA(1:iii3,i);
 else if cycle(i,3)==cycle4(end,3)
 vcycle1(:,1)=AA(1:iii4,i);
 else
 vcycle1(:,1)=AA(1:dd(cycle(i,3)),i);
 end
 end
 end
 end
 vcycle=[vcycle;vcycle1];
end

timee=ct1+ct2+ct3+ct4;
tcycle=1:1:1200
figure
if length(vcycle(:,1))<=length(tcycle(1,:))
plot(tcycle(1,(1:length(vcycle))),vcycle)
xlabel('Time(s)','fontsize',18,'fontname','Times New Roman');ylabel('Speed(m/s)','fontsize',18,'fontname','Times New Roman');
else
plot(tcycle,vcycle((1:length(tcycle)),1))
xlabel('Time(s)','fontsize',18,'fontname','Times New Roman');ylabel('Speed(m/s)','fontsize',18,'fontname','Times New Roman');
end
