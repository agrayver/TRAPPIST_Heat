clear all
close all

colors = [0.0000    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];

Markers = {'+','o','*','x','v','d','^','s','>','<'};

planet_names = {'1b','1c','1d','1e','1f','1g'};

% Earth radius
a_Earth = 6371200;

% Reference value for the radiogenic heat of Earth
radiogenic_heat = 2 * 1e13;

% TRAPPIST-1 planets radii: b,c,d,e,f,g
% from https://arxiv.org/pdf/2010.01074.pdf
planetary_radii(1) = 1.116*a_Earth;
planetary_radii(2) = 1.097*a_Earth;
planetary_radii(3) = 0.788*a_Earth;
planetary_radii(4) = 0.920*a_Earth;
planetary_radii(5) = 1.045*a_Earth;
planetary_radii(6) = 1.129*a_Earth;

% switches on intrinsic magnetic field
with_intrinsic_field = 0;

% 1 -- 1-D pyrolitic mantle;
% 2 -- homogeneous sphere 0.01 S/m
% 3 -- homogeneous sphere 1 S/m
model_idx = 1;

%%
load('data/conductivity_models_1D_and_homogeneous.mat');
load(sprintf('data/TRAPPIST1_heating_intrinsic_field=%d_model_idx=%d.mat',with_intrinsic_field,model_idx));

n_planets = length(planet_data);

n_bins = 30;

% 'Normalization' = 'probability' is when
% the height of the histogram bar represents the proportion of the data in each class. 
figure('Position',[150 150 400 400]);
for pidx=1:n_planets

    histogram(log10(planet_data{pidx}.Q_total),n_bins,'HandleVisibility', 'off','Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth',1,'EdgeColor',colors(pidx,:));hold all;
    histogram(log10(planet_data{pidx}.Q_total),n_bins,'DisplayName',planet_names{pidx},'Normalization', 'probability', 'EdgeColor','none','FaceColor',colors(pidx,:));hold all;
    
    Q_mean(pidx) = mean(log10(planet_data{pidx}.Q_total));
    Q_std(pidx) = std(log10(planet_data{pidx}.Q_total));
    plot([Q_mean(pidx),Q_mean(pidx)],[0, 0.2],'--','LineWidth',1,'color',colors(pidx,:),'HandleVisibility','off');
    
    S_pl = 4*pi*planetary_radii(pidx)^2;
    Qs_dt(pidx) = 10^Q_mean(pidx) / S_pl;
end

plot(log10([radiogenic_heat,radiogenic_heat]),[0, 0.2],'--','LineWidth',2,'color','k','HandleVisibility','off');
ylim([0, 0.16]);
xlim([10, 15]);
% set(gca,'xscale','log');
ylabel('Relative probability');
xlabel('log_{10} Heat power (W)');
legend('show','location','best');

% print(sprintf('./figures/Energy_Histogram_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx),'-dpng','-r600');
% crop_png(sprintf('./figures/Energy_Histogram_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx)); 

% data_table = 10.^[Q_mean' (Q_mean - Q_std)' (Q_mean + Q_std)'];
% save(sprintf('./figures/Q_mean_std_intrinsic_field=%d_model_idx=%d.txt', with_intrinsic_field, model_idx), '-ascii', 'data_table');

%% Plot heat profiles
figure('Position',[50 50 400 400]);
% subplot(212);
for pidx=1:n_planets

    radii = r_TRAPPIST{pidx,1};
    
    dV = 4/3 * pi * (radii(1:end-1).^3 - radii(2:end).^3);
    
    Q_x = radii(1:end-1)/radii(1);
    Q_y = mean(planet_data{pidx}.Q_profile_mean,2);
    
    dQ_y = std(planet_data{pidx}.Q_profile_mean, 0, 2);
    
    stairs(Q_y./dV,Q_x,'DisplayName',planet_names{pidx}, 'Color', colors(pidx,:),'LineWidth',1);hold all

    Q_lower = (Q_y-dQ_y) ./ dV;
    Q_upper = (Q_y+dQ_y) ./ dV;
    
    y_l = [Q_lower(1);repelem(Q_lower(2:end),2)];
    y_u = [Q_upper(1);repelem(Q_upper(2:end),2)];
    x = [repelem(Q_x(1:end-1),2);Q_x(end)];
    fill([y_l;flipud(y_u)],[x;flipud(x)],colors(pidx,:),'FaceAlpha',0.1,'HandleVisibility','off','EdgeColor', 'none');
end

set(gca,'XScale','log');
xlabel('Heat power density (W/m^3)');
ylabel('Normalized radius');
legend('show','location','best');
ylim([0.48 1]);
xlim([1e-15 1e-5]);

% print(sprintf('./figures/Heat_Profile_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx),'-dpng','-r600');
% crop_png(sprintf('./figures/Heat_Profile_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx)); 

%% Plot conductivity models
figure('Position',[50 50 400 400]);
for pidx=1:n_planets

    sigma = sigma_TRAPPIST{pidx,model_idx};
    radii = r_TRAPPIST{pidx,1};
    
    stairs(sigma,radii/radii(1),...
        'DisplayName',planet_names{pidx},'LineWidth',1,'Color',colors(pidx,:));
    set(gca,'XScale','log');hold all

end

xlabel('Electrical Conductivity (S/m)');
ylabel('Normalized radius');
legend('show','location','best');
ylim([0.48 1]);
xlim([1e-3 1e1]);

% print(sprintf('./figures/Conductivity_Profile_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx),'-dpng','-r600');
% crop_png(sprintf('./figures/Conductivity_Profile_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx)); 

%% Plot production of heat per energy bands

n_simulation_years = length(planet_data{1}.n_collided_flares_per_year);

N = 30;
[E_bin_idx,E_bin_edges] = discretize(log10(planet_data{1}.Flare_energies),N);

figure('Position',[50 50 400 400]);
for pidx=1:n_planets
    
    Q_per_energy_bin = [];
    for i=1:N
        n = find(E_bin_idx == i);
        Q_per_energy_bin(i) = sum(planet_data{pidx}.Q_total_per_flare(E_bin_idx(n))) / n_simulation_years;
    end
    
    E_bin_mid = diff(E_bin_edges);
    scatter(10.^(E_bin_edges(1:end-1)+E_bin_mid/2), Q_per_energy_bin,Markers{pidx},'MarkerEdgeColor',colors(pidx,:),'DisplayName',planet_names{pidx});hold all;
    
end

set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('Flare Energy (erg)');
ylabel('Heat power (W)');
legend('show','location','best');
ylim([1e7 2e14]);
box on;

% print(sprintf('./figures/Energy_Heat_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx),'-dpng','-r600');
% crop_png(sprintf('./figures/Energy_Heat_intrinsic_field=%d_model_idx=%d.png', with_intrinsic_field, model_idx)); 

%%
figure;
for pidx=1:n_planets
    histogram(planet_data{pidx}.n_collided_flares_per_year, 'EdgeColor','none','DisplayName',planet_names{pidx},'FaceColor',colors(pidx,:));hold all
end