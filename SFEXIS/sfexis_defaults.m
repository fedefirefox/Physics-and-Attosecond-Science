function config = sfexis_defaults()
%SFEXIS_DEFAULTS Small constant-dipole, isolated-pulse streaking example.
% Original SFEXIS: Federico Vismarra. Cite Vismarra et al. (2022),
% Journal of Physics: Photonics 4, 034006, doi:10.1088/2515-7647/ac7991.
config.dt_fs = 0.005;
config.time_half_width_fs = 20;
config.delay_fs = -3:0.25:3;
config.energy_eV = 15:0.5:65;
config.ionization_potential_eV = 15.75962; % Ar threshold, constant dipole
config.xuv_photon_energy_eV = 50;
config.xuv_fwhm_fs = 0.25; % electric-field intensity envelope
config.xuv_cep_rad = 0;
config.ir_wavelength_nm = 800;
config.ir_fwhm_fs = 5; % vector-potential intensity envelope; see guide
config.ir_intensity_Wcm2 = 1e12; % actual peak of E(t)^2
config.ir_cep_rad = 0;
end
