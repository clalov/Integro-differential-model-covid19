% MAIN integro-differential

% YY = [I; Y; H; Q; C];
[YY,D,T]=integrodifferential();

figure;
plot(T,YY(:,1),'displayname','Asymp. infected'); hold on; 
plot(T,YY(:,2),'displayname','Symptomatic'); hold on;
plot(T,YY(:,3),'displayname','Hospitalized'); hold on;
plot(T,YY(:,4),'displayname','Quarantined'); hold on;
plot(T,YY(:,5),'displayname','Recovered'); hold on;
plot(T,D(1:end-1),'displayname','Deceased'); hold on;
legend('location','best');