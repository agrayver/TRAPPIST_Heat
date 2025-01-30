%% OMNI data base for IMF
omni_data = load('data/omni2_all_years.mat');omni_data = omni_data.omni_data;

F_IMF = omni_data(:,10);
Dst = omni_data(:,41);
dates = datetime(omni_data(:,1),1,omni_data(:,2),omni_data(:,3),0,0);
V_bulk = omni_data(:,25);

idx = find(F_IMF == 9.999000000000000e+02 | Dst == 9.999000000000000e+02 | dates < datetime(1986,1,1));

F_IMF(idx) = [];
Dst(idx) = [];
dates(idx) = [];

idx_storm = find(Dst < -150);

% Identify isolated peaks
diff_idx = find(diff(idx_storm) > 24);

% Identify maximum IMF within the geomagnetic storm
% Look within +-48 hours window to account for the time lag between IMF and Dst
Dst_storm = [];
IMF_storm = [];
V_storm = [];
date_storm = [];
for i=1:length(diff_idx)
    [Dst_storm(i), idx_min] = min(Dst(idx_storm(diff_idx(i))-48:idx_storm(diff_idx(i))+48));
    date_storm(i) = datenum(dates(idx_storm(diff_idx(i))-48+idx_min));
    IMF_storm(i) = max(F_IMF(idx_storm(diff_idx(i))-48:idx_storm(diff_idx(i))+48));
end

IMF_storm(43) = [];
Dst_storm(43) = [];
date_storm(43) = [];

[p_IMF2Dst,S_IMF2Dst] = polyfit(IMF_storm, Dst_storm,1); 
[y_fit,p_IMF2Dst_se] = polyval(p_IMF2Dst,sort(IMF_storm),S_IMF2Dst);

% Check significance
x1 = ones(length(IMF_storm),1);
X = [x1 IMF_storm.'];
y = Dst_storm.';
[b,bint,r,rint,stats] = regress(y,X); % stats(3) is p-value and should be < 0.05 (alpha level) to be significant

figure;
plot(IMF_storm, Dst_storm, 'o');hold all
plot(sort(IMF_storm),y_fit,'k-');
plot(sort(IMF_storm),y_fit+2*p_IMF2Dst_se,'--','Color', [0.5,0.5,0.5]);
plot(sort(IMF_storm),y_fit-2*p_IMF2Dst_se,'--','Color', [0.5,0.5,0.5]);
legend('Data',sprintf('Linear Fit (r = %.2f, p-value=%.2e)', corr(IMF_storm', Dst_storm'), stats(3)),'2-\sigma Interval','location','best')
xlabel('|B_{IMF}| (nT)');
ylabel('Dst (nT)')

% print('IMF_vs_Dst.png','-dpng','-r300');