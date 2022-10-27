omni_data = load('omni2_all_years.mat');omni_data = omni_data.omni_data;

F_IMF = omni_data(:,10);
Dst = omni_data(:,41);
V_bulk = omni_data(:,25);

dates = datetime(omni_data(:,1),1,omni_data(:,2),omni_data(:,3),0,0);

idx = find(F_IMF == 9.999000000000000e+02 | Dst == 9.999000000000000e+02 | dates < datetime(1986,1,1));

F_IMF(idx) = [];
Dst(idx) = [];
V_bulk(idx) = [];
dates(idx) = [];

idx_storm = find(Dst < -150);

%%
% Identify isolated peaks
diff_idx = find(diff(idx_storm) > 24);

% Identify maximum IMF within the geomagnetic storm
% Looks within +-24 hours window to account for time lag between IMF and Dst
Dst_storm = [];
IMF_storm = [];
V_storm = [];
date_storm = [];
for i=1:length(diff_idx)
    [Dst_storm(i), idx_min] = min(Dst(idx_storm(diff_idx(i))-48:idx_storm(diff_idx(i))+48));
    
    date_storm(i) = datenum(dates(idx_storm(diff_idx(i))-48+idx_min));
    
    IMF_storm(i) = max(F_IMF(idx_storm(diff_idx(i))-48:idx_storm(diff_idx(i))+48));
    V_storm(i) = max(V_bulk(idx_storm(diff_idx(i))-48:idx_storm(diff_idx(i))+48));
end

V_storm(V_storm == 9999) = NaN;

[p,S] = polyfit(IMF_storm, Dst_storm,1); 
[y_fit,delta] = polyval(p,sort(IMF_storm),S);

figure;
% subplot(211)
plot(IMF_storm, Dst_storm, 'o');hold all
plot(sort(IMF_storm),y_fit,'k-');
plot(sort(IMF_storm),y_fit+2*delta,'--','Color', [0.5,0.5,0.5]);
plot(sort(IMF_storm),y_fit-2*delta,'--','Color', [0.5,0.5,0.5]);
legend('Data','Linear Fit (r = -0.7545)','95% Prediction Interval')
xlabel('|B_{IMF}| (nT)');
ylabel('Dst (nT)')

% subplot(212)
% scatter(IMF_storm, Dst_storm./IMF_storm, 'o')