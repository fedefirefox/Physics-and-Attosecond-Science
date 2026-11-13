function tests = test_sfexis
% Physics and regression checks, base MATLAB only.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'examples'));
testCase.TestData.root = root;
testCase.TestData.default = run_sfexis();
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testDefaultOutput(testCase)
r = testCase.TestData.default;
verifySize(testCase,r.signal,[25 101]);
verifyTrue(testCase,all(isfinite(r.signal(:))) && all(r.signal(:) >= 0));
verifyEqual(testCase,max(r.normalized_signal(:)),1,'AbsTol',eps);
verifyEqual(testCase,r.signal,abs(r.amplitude).^2.* ...
    sqrt(2*r.energy_eV/r.constants.energy_au_eV),'RelTol',1e-14);
verifyTrue(testCase,contains(r.citation,'10.1088/2515-7647/ac7991'));
end

function testZeroIRAnalyticSpectrum(testCase)
r = run_sfexis(struct('ir_intensity_Wcm2',0));
c = r.constants;
alpha = 2*log(2)/(r.config.xuv_fwhm_fs/c.time_au_fs)^2;
omega = r.config.xuv_photon_energy_eV/c.energy_au_eV;
q = (r.energy_eV+r.config.ionization_potential_eV)/c.energy_au_eV;
% Exact Fourier transform of exp(-alpha*t^2)*cos(omega*t).
fourier = sqrt(pi/alpha)/2 * (exp(-(q-omega).^2/(4*alpha)) + ...
    exp(-(q+omega).^2/(4*alpha)));
expected = fourier.^2 .* sqrt(2*r.energy_eV/c.energy_au_eV);
relative = max(abs(r.signal-expected),[],'all')/max(expected);
fprintf('Zero-IR analytic relative max error: %.3g\n',relative);
verifyLessThan(testCase,relative,1e-10);
end

function testDelayTranslationPhase(testCase)
r = run_sfexis(struct('ir_intensity_Wcm2',0,'delay_fs',[-0.237 0 0.411]));
q = (r.energy_eV+r.config.ionization_potential_eV)/r.constants.energy_au_eV;
expected = r.amplitude(2,:) .* exp(1i*r.delay_fs/r.constants.time_au_fs*q);
verifyLessThan(testCase,max(abs(r.amplitude-expected),[],'all') / ...
    max(abs(expected),[],'all'),1e-10);
end

function testConstantPotentialAnalyticAction(testCase)
t = (-100:0.1:100)';
alpha = 0.02;
omega = 1.8;
xuv = exp(-alpha*t.^2).*exp(-1i*omega*t);
energy = 0.4:0.1:2;
ip = 0.5;
potential = 0.12*ones(size(t));
signal = sfexis_spectrum(t,xuv,potential,energy,ip);
p = sqrt(2*energy);
q = (p+0.12).^2/2+ip;
expected = pi/alpha*exp(-(q-omega).^2/(2*alpha)).*p;
verifyLessThan(testCase,max(abs(signal-expected))/max(expected),1e-11);
end

function testIRFieldDerivativeAndIntensity(testCase)
r = testCase.TestData.default;
step = r.config.dt_fs/r.constants.time_au_fs;
from_potential = -gradient(r.ir_vector_potential_au,step);
relative = max(abs(from_potential-r.ir_electric_field_au))/ ...
    max(abs(r.ir_electric_field_au));
verifyLessThan(testCase,relative,5e-5);
verifyEqual(testCase,max(r.ir_electric_field_au.^2)*r.constants.intensity_au_Wcm2, ...
    r.config.ir_intensity_Wcm2,'RelTol',1e-12);
verifyLessThan(testCase,max(abs(r.ir_vector_potential_au([1 end]))),1e-9);
end

function testStreakingSign(testCase)
r = testCase.TestData.default;
free = run_sfexis(struct('ir_intensity_Wcm2',0));
centroid = sum(r.signal.*r.energy_eV,2)./sum(r.signal,2);
baseline = sum(free.signal.*free.energy_eV,2)./sum(free.signal,2);
a = interp1(r.time_fs,r.ir_vector_potential_au,r.delay_fs);
p0 = sqrt(2*(r.config.xuv_photon_energy_eV-r.config.ionization_potential_eV) / ...
    r.constants.energy_au_eV);
predicted = (-p0*a+a.^2/2)*r.constants.energy_au_eV;
correlation = corrcoef(centroid-baseline,predicted);
fprintf('Streaking-centroid/classical-shift correlation: %.6f\n',correlation(1,2));
verifyGreaterThan(testCase,correlation(1,2),0.98);
verifyGreaterThan(testCase,max(abs(centroid-baseline)),1);
end

function testTimeStepConvergence(testCase)
options = struct('dt_fs',0.01,'delay_fs',[-1.23 0 0.71],'energy_eV',20:1:55);
coarse = run_sfexis(options);
options.dt_fs = 0.005;
medium = run_sfexis(options);
options.dt_fs = 0.0025;
fine = run_sfexis(options);
error_coarse = norm(coarse.signal-fine.signal,'fro')/norm(fine.signal,'fro');
error_medium = norm(medium.signal-fine.signal,'fro')/norm(fine.signal,'fro');
fprintf('Time-step relative errors (10 as, 5 as vs 2.5 as): %.3g, %.3g\n', ...
    error_coarse,error_medium);
verifyLessThan(testCase,error_medium,1e-4);
verifyLessThan(testCase,error_medium,error_coarse/3);
end

function testTimeWindowConvergence(testCase)
r = testCase.TestData.default;
wide = run_sfexis(struct('time_half_width_fs',25));
verifyLessThan(testCase,norm(r.signal-wide.signal,'fro')/norm(wide.signal,'fro'),1e-9);
end

function testLegacyConstantDipoleQuadrature(testCase)
% Original SFEXIS.m velocity-scan formula, same fields and integration grid.
% Phase constants disappear in |a|^2. The original rectangular rule differs
% from the new trapezoidal rule, so compare converged spectra, not amplitudes.
r = run_sfexis(struct('dt_fs',0.0025,'delay_fs',0.7,'energy_eV',20:1:55));
t = r.time_fs/r.constants.time_au_fs;
dt = r.config.dt_fs/r.constants.time_au_fs;
s = r.time_fs-r.delay_fs;
xuv = exp(-2*log(2)*(s/r.config.xuv_fwhm_fs).^2).*cos( ...
    r.config.xuv_photon_energy_eV/(r.constants.energy_au_eV*r.constants.time_au_fs)*s);
p = sqrt(2*r.energy_eV/r.constants.energy_au_eV);
ip = r.config.ionization_potential_eV/r.constants.energy_au_eV;
original = zeros(size(p));
for k = 1:numel(p)
    shifted = r.ir_vector_potential_au-r.ir_vector_potential_au(1)+p(k);
    action = 0.5*flipud(cumsum(flipud(shifted.^2)))*dt;
    amplitude = -1i*sum(xuv.*exp(-1i*action).*exp(1i*ip*(t-t(1))))*dt;
    original(k) = abs(amplitude)^2*p(k);
end
relative = norm(original-r.signal)/norm(r.signal);
fprintf('Legacy rectangular/new trapezoidal relative error: %.3g\n',relative);
verifyLessThan(testCase,relative,0.005);
end

function testInvalidConfigurations(testCase)
verifyError(testCase,@()run_sfexis(struct('Ipeak',1)),'sfexis:UnknownOption');
verifyError(testCase,@()run_sfexis(struct('dt_fs',0)),'sfexis:InvalidValue');
verifyError(testCase,@()run_sfexis(struct('dt_fs',0.1)),'sfexis:Undersampled');
verifyError(testCase,@()run_sfexis(struct('time_half_width_fs',3)),'sfexis:ShortWindow');
verifyError(testCase,@()run_sfexis(struct('ir_intensity_Wcm2',-1)),'sfexis:InvalidGrid');
verifyError(testCase,@()run_sfexis(struct('energy_eV',[1 1])),'sfexis:InvalidGrid');
verifyError(testCase,@()run_sfexis(struct('delay_fs',[1 0])),'sfexis:InvalidGrid');
verifyError(testCase,@()sfexis_spectrum([0 1 1],ones(3,1),[0 0 0],1,1), ...
    'sfexis:InvalidGrid');
end

function testSavedExampleAndOverwriteProtection(testCase)
folder = tempname(fullfile(testCase.TestData.root,'results'));
r = minimal_streaking(folder);
loaded = load(fullfile(folder,'streaking.mat'),'result');
verifyEqual(testCase,loaded.result,r);
image_info = imfinfo(fullfile(folder,'streaking.png'));
verifyGreaterThan(testCase,image_info.Width,500);
verifyError(testCase,@()minimal_streaking(folder),'sfexis:ExistingOutput');
end
