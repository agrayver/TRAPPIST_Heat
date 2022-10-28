Flare_energy = logspace(29,34);

% Radius of the sun
R_Sun = 696340*1e3;
% Radius of the TRAPPIST star
R_Star = R_Sun * 0.1192; % (https://arxiv.org/pdf/2010.01074.pdf)

% Astronomical unit
AU = 1.496e+11;
Gauss2nT = 1e5;

% From Samara et al. (2021)
mean_decay_rate = 1.6;
std_decay_rate = 0.1;

% Bolometric energy is roughly 30% of the kinetic energy
% See figure 3 in
% https://iopscience.iop.org/article/10.1088/0004-637X/759/1/71/meta
CME_bolometric2kinetic = 200 / 71;

% Semi-major orbital axes for TRAPPIST-1 planets (AU) (https://arxiv.org/pdf/2010.01074.pdf)
orbital_distances = [0.01154;
    0.01588;
    0.02227;
    0.02925;
    0.03849;
    0.04683];

B0_Sun = @(energy) 5.94 * 1e-9 * energy.^0.23; % eq.A6
B0_Star = B0_Sun(CME_bolometric2kinetic*Flare_energy) * (R_Sun / R_Star)^2; % eq. A5

colors = [0.0000    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
planet_names = {'1b','1c','1d','1e','1f','1g'};

n_planets = length(orbital_distances);

figure('Position',[150 150 400 400]);
for pidx=1:n_planets
    B_ICME_low = Gauss2nT * B0_Star * ((R_Star * 10) / (orbital_distances(pidx)*AU)) ^ (mean_decay_rate - 3*std_decay_rate); %eq. A4
    B_ICME_up = Gauss2nT * B0_Star * ((R_Star * 10) / (orbital_distances(pidx)*AU)) ^ (mean_decay_rate + 3*std_decay_rate); %eq. A4
    
    x2 = [Flare_energy, fliplr(Flare_energy)];
    inBetween = [B_ICME_low, fliplr(B_ICME_up)];
    fill(x2, inBetween, colors(pidx,:),'HandleVisibility','off','FaceAlpha',0.25,'EdgeColor','none');hold all
end
for pidx=1:n_planets
    B_ICME_mean = Gauss2nT * B0_Star * ((R_Star * 10) / (orbital_distances(pidx)*AU)) ^ mean_decay_rate; %eq. A4
    plot(Flare_energy,B_ICME_mean,'Color',colors(pidx,:),'DisplayName',planet_names{pidx},'LineWidth',2);hold all;
end
legend('show','location','northwest');
xlabel('E (erg)');
ylabel('|B| (nT)');
xlim([min(Flare_energy) max(Flare_energy)]);
set(gca,'xscale','log');

% print('./figures/Energy_vs_B_field.png','-dpng','-r600');
% crop_png('./figures/Energy_vs_B_field.png'); 