clc
clear all
close all
%%load data
speed= xlsread('speed.xlsx');
time= xlsread('time.xlsx');
%% make microtrip
DD = [1];
LL = [];
lp=[];
for ii=1:length(speed(:,1))
 if speed(ii,1)<=2
 DD=[DD,ii];
 end
end
h=[];
for jj=2:length(DD)
 if DD(1,jj)~=DD(1,(jj-1))+1
 h=[h DD(jj)];
 end
end

vmt =[];
time1=[];
REPORT=[];
%  LLl=[];
% for y=double(DD(1)):1:double(DD(2))
%         LLl = [LLl;speed(y)];        
% end
    DD=h; 
for kk = 2:length(DD)
    LL=[];
    time1=[];
    
    for jj=double(DD(kk-1)):1:double(DD(kk))
        LL = [LL;speed(jj)];
        time1=[time1;time(jj)];
        microtrip=[time1,LL];
        
    end


%% gain parameters
G = microtrip;
V_threshold=0.1; a_threshold=0.1;
%% main data parameters
%%%% Standard Deviation of Speed
i=G(:,2).^2;
Vsd=((1/(size(G,1)-2))*sum(i(:,1)))^(0.5); % Standard deviation of speed
%%%% Acceleration Matrix
A=zeros(size(G,1)-1,1); j=zeros(size(G,1)-1,1); k=zeros(size(G,1)-1,1);
for i=1:size(G,1)-1
    A(i,1)=(G(i+1,2)-G(i,2))/(G(i+1,1)-G(i,1)); % A: Matrix of accelerations
    if A(i,1)>0
        j(i,1)=A(i,1);j(1,1)=0;
    elseif A(i,1)<0
        k(i,1)=A(i,1);
    end
    if G(i,2)==0 && A(i)==0
        %         q=A==0;
    end
end
A(1:2,:)=0;
k(1:2,:)=[];
% stop_nr=sum(q); % stop_nr: Number of stops
j(j(:,1)==0)=[]; k(k(:,1)==0)=[];
A_mean_pos=sum(j)/size(j,1); % A_mean_pos: Average positive acceleration
A_mean_neg=sum(k)/size(k,1); % A_mean_neg: Average negative acceleration
%%%% Standard Deviation of Acceleration
i=A(:,1).^2;
Asd=((1/(size(A,1)-1))*(sum(i(:,1))))^(0.5); % Standard deviation of acceleration
%%%% Total Distance
Total_d=zeros(size(G,1)-1,1);
for i=1:size(G,1)-1
    Total_d(i,1)=((G(i+1,1)-G(i,1))*G(i,2))/3.6;
end
Total_d=sum(Total_d); % Total_d: Total distance
%%%% Total Time
Durations=zeros(size(G,1)-1,1);
for i=1:size(G,1)-1
    Durations(i,1)=G(i+1,1)-G(i,1);
end
Total_t=sum(Durations(:,1)); % Total_t: Total time
%%%% Mean Acceleration
% A_mean=sum(A(:,1))/Total_t; % A_mean: Mean acceleration
%%%% Standing Time
T_stop=zeros(size(G,1)-1,1);
for i=1:size(G,1)-1
    if G(i,2)<V_threshold && A(i,1)<a_threshold && A(i,1)>-a_threshold
        T_stop(i,1)=G(i+1,1)-G(i,1);
    else; T_stop(i,1)=0;
    end
end
T_stop=sum(T_stop); % T_stop: Standing time
%%%% Accelerating and Decelerating time
T_acceleration=zeros(size(A,1),1); T_Deceleration=zeros(size(A,1),1);
for i=1:size(G,1)-1
    if A(i,1)>a_threshold
        T_acceleration(i,1)=G(i+1,1)-G(i,1);
    else; T_acceleration(i,1)=0;
    end
    if A(i,1)<-a_threshold
        T_Deceleration(i,1)=G(i+1,1)-G(i,1);
    else; T_Deceleration(i,1)=0;
    end
end
T_acceleration=sum(T_acceleration); % T_acceleration: Accelerating time
T_Deceleration=sum(T_Deceleration); % T_Deceleration: Decelerating time
%%%% Driving and Cruise time
T_drive=Total_t-T_stop; % T_drive: Driving time
T_cruise=T_drive-(T_acceleration+T_Deceleration); % T_cruise: Cruise time
%%%% Time Percentage
prcnt_acc=(T_acceleration/Total_t)*100; % prcnt_acc: Accelerating time Percentage
prcnt_dec=(T_Deceleration/Total_t)*100; % prcnt_dec: Decelerating time Percentage
prcnt_c=(T_cruise/Total_t)*100; % prcnt_c: Cruise time Percentage
prcnt_D=(T_drive/Total_t)*100; % prcnt_D: Driving time Percentage
prcnt_S=(T_stop/Total_t)*100; % prcnt_S: Standing time Percentage
%%%% Other Parameters
Vm_trip=3.6*(Total_d/Total_t); % Vm_trip: Average Speed of Trip
Vm_drive=3.6*(Total_d/T_drive); % Vm_trip: Average Driving Speed
% stop_T_mean=T_stop/stop_nr; % stop_T_mean: Average Time of Standing
% Max_p_ac=max(A); % Max_p_ac: Maximum Positive Acceleration
% Max_n_ac=min(A); % Max_n_ac: Maximum Negative Acceleration
Max_V=max(G(:,2)); % Max_V: Maximum Speed
Main=[prcnt_D;prcnt_S;prcnt_c;prcnt_acc; prcnt_dec;Vm_drive;Vm_trip;Max_V;A_mean_pos;A_mean_neg;Vsd;Asd];
Parameters={'%Drive';'%Stop';'%Cruise';'%Accelerating';'%Decelerating';'Driving mean speed';'Trip mean speed';'Max Speed'
    'Mean Positive Acceleration';'Mean Negative Acceleration';'Standard deviation of Speed ';'Standard deviation of Acceleration'};
REPORT_TABLE=table(Parameters,Main);
REPORT=[REPORT,Main];
    if kk>1
        [a b]=size(vmt);
        while a<length(LL)
            vmt=[vmt;zeros(1,b)];
            a=a+1
        end
        while length(LL)<a
            LL=[LL;0];
        end
        vmt=[vmt,LL];
        
    else
        vmt=[vmt,LL];
        
    end
end

%% gain PCA
H = REPORT';
pcc_drive = H(:,1);
pcc_stop = H(:,2);
pcc_cruise = H(:,3);
pcc_acceleration = H(:,4);
pcc_deceleration = H(:,5);
driving_mean_spead = H(:,6);
trip_mean_speed = H(:,7);
max_speed = H(:,8);
mean_positive_acceleration = H(:,9);
mean_negative_acceleration = H(:,10);
standard_deviation_speed = H(:,11);
standard_deviation_acceleration = H(:,12);
%% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(H)
Xcentered = score*coeff'
figure(7)
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'PCC-D','PCC-S','PCC-C','PCC-ac','PCC-dec','dms','tms','ms','mpa','mna','sds','sda'});
xlabel('principal components','fontsize',28,'fontname','Times New Roman');
ylabel('explained','fontsize',28,'fontname','Times New Roman');
PCA_Percentages = [explained];
figure (6)
% PCA_name_tags = categorical({'Drive Percentage','Stop Percentage','Cruise Percentage','Acceleration Percentage','Decceleration Percentage','Driving mean speed','Trip mean speed','Max Speed','Mean Positive Acceleration','Mean Negative Acceleration','Standard deviation of Speed ','Standard deviation of Acceleration'});
% PCA_name_tags = reordercats(PCA_name_tags,{'Drive Percentage','Stop Percentage','Cruise Percentage','Acceleration Percentage','Decceleration Percentage','Driving mean speed','Trip mean speed','Max Speed','Mean Positive Acceleration','Mean Negative Acceleration','Standard deviation of Speed ','Standard deviation of Acceleration'});
bar(PCA_Percentages,0.5)
xlabel({'Principal Components'});
ylabel({'Explained'});

