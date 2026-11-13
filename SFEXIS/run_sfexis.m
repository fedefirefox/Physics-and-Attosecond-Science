function result = run_sfexis(overrides)
%RUN_SFEXIS Noninteractive 1D SFA calculation with a constant unit dipole.
% result = run_sfexis() or run_sfexis(struct('ir_intensity_Wcm2',0)).
% Positive delay places the XUV maximum later than the IR envelope maximum.
% No files, figures, external data, toolboxes, or workspace changes are made.
% Original SFEXIS: Federico Vismarra. Cite Vismarra et al. (2022),
% Journal of Physics: Photonics 4, 034006, doi:10.1088/2515-7647/ac7991.
if nargin == 0
    overrides = struct();
end
validateattributes(overrides, {'struct'}, {'scalar'});
config = sfexis_defaults();
names = fieldnames(overrides);
for k = 1:numel(names)
    if ~isfield(config, names{k})
        error('sfexis:UnknownOption', 'Unknown option: %s', names{k});
    end
    config.(names{k}) = overrides.(names{k});
end
names = fieldnames(config);
for k = 1:numel(names)
    value = config.(names{k});
    validateattributes(value, {'double'}, {'real','finite','nonempty'}, ...
        mfilename, names{k});
    if ismember(names{k}, {'delay_fs','energy_eV'})
        validateattributes(value, {'double'}, {'vector'});
    else
        validateattributes(value, {'double'}, {'scalar'});
    end
end
positive = {'dt_fs','time_half_width_fs','ionization_potential_eV', ...
    'xuv_photon_energy_eV','xuv_fwhm_fs','ir_wavelength_nm','ir_fwhm_fs'};
for k = 1:numel(positive)
    if config.(positive{k}) <= 0
        error('sfexis:InvalidValue', '%s must be positive.', positive{k});
    end
end
if config.ir_intensity_Wcm2 < 0 || any(config.energy_eV < 0) || ...
        any(diff(config.energy_eV) <= 0) || any(diff(config.delay_fs) <= 0)
    error('sfexis:InvalidGrid', ...
        'Intensity and energies must be nonnegative; axes must increase.');
end
if config.time_half_width_fs < max(abs(config.delay_fs)) + 4*config.xuv_fwhm_fs || ...
        config.time_half_width_fs < 4*config.ir_fwhm_fs
    error('sfexis:ShortWindow', 'Allow four pulse widths beyond the pulses.');
end
constants.time_au_fs = 0.02418884326505;
constants.energy_au_eV = 27.21138624;
constants.intensity_au_Wcm2 = 3.5094452e16;
constants.c_nm_fs = 299.792458;
hbar_eV_fs = constants.time_au_fs * constants.energy_au_eV;
omega_xuv_fs = config.xuv_photon_energy_eV / hbar_eV_fs;
omega_ir_fs = 2*pi*constants.c_nm_fs / config.ir_wavelength_nm;
% Resolve both the carrier and the largest kinetic/Volkov phase frequency.
a_bound = 2*sqrt(config.ir_intensity_Wcm2/constants.intensity_au_Wcm2) / ...
    (omega_ir_fs*constants.time_au_fs);
omega_action_fs = (((sqrt(2*max(config.energy_eV)/constants.energy_au_eV) + ...
    a_bound)^2/2)*constants.energy_au_eV + config.ionization_potential_eV) / hbar_eV_fs;
if config.dt_fs * max([omega_xuv_fs, omega_ir_fs, omega_action_fs]) >= pi/2 || ...
        config.dt_fs > min(config.xuv_fwhm_fs,config.ir_fwhm_fs)/10
    error('sfexis:Undersampled', 'Reduce dt_fs to resolve the fields and action.');
end
n = ceil(config.time_half_width_fs / config.dt_fs);
time_fs = (-n:n)' * config.dt_fs;
time_au = time_fs / constants.time_au_fs;
omega_ir_au = omega_ir_fs * constants.time_au_fs;
width_au = config.ir_fwhm_fs / constants.time_au_fs;
gaussian = exp(-2*log(2)*(time_au/width_au).^2);
gaussian_derivative = -4*log(2)*time_au/width_au^2 .* gaussian;
phase = omega_ir_au*time_au + config.ir_cep_rad;
potential = -gaussian .* sin(phase) / omega_ir_au;
electric = gaussian.*cos(phase) + gaussian_derivative.*sin(phase)/omega_ir_au;
scale = sqrt(config.ir_intensity_Wcm2/constants.intensity_au_Wcm2) / max(abs(electric));
potential = potential * scale;
electric = electric * scale;
shifted_fs = time_fs - config.delay_fs(:)';
xuv = exp(-2*log(2)*(shifted_fs/config.xuv_fwhm_fs).^2) .* ...
    cos(omega_xuv_fs*shifted_fs + config.xuv_cep_rad);
[signal, amplitude] = sfexis_spectrum(time_au, xuv, potential, ...
    config.energy_eV/constants.energy_au_eV, ...
    config.ionization_potential_eV/constants.energy_au_eV);
peak = max(signal(:));
if ~isfinite(peak) || peak <= 0
    error('sfexis:EmptySignal', 'No finite positive signal; check the energy range.');
end
result.config = config;
result.constants = constants;
result.time_fs = time_fs;
result.delay_fs = config.delay_fs(:);
result.energy_eV = config.energy_eV(:)';
result.ir_vector_potential_au = potential;
result.ir_electric_field_au = electric;
result.signal = signal; % delay x energy, p*|amplitude|^2, arbitrary scale
result.normalized_signal = signal / peak; % one global normalization
result.amplitude = amplitude;
result.citation = ['F. Vismarra, R. Borrego-Varillas, Y. Wu, D. Mocci, ' ...
    'M. Nisoli and M. Lucchini, Journal of Physics: Photonics 4, 034006 (2022). ' ...
    'https://doi.org/10.1088/2515-7647/ac7991'];
end
