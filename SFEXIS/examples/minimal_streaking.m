function result = minimal_streaking(output_directory)
%MINIMAL_STREAKING Run and save a deterministic isolated-pulse example.
% Add SFEXIS and SFEXIS/examples to your MATLAB path before calling.
% Cite Federico Vismarra et al. (2022), doi:10.1088/2515-7647/ac7991.
if nargin == 0
    output_directory = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'results','minimal_streaking');
end
if ~isfolder(output_directory)
    mkdir(output_directory);
end
mat_file = fullfile(output_directory,'streaking.mat');
png_file = fullfile(output_directory,'streaking.png');
if isfile(mat_file) || isfile(png_file)
    error('sfexis:ExistingOutput', 'Choose a new directory to preserve prior results.');
end
result = run_sfexis();
result.environment = version;
save(mat_file,'result');
fig = figure('Visible','off','Color','w');
cleanup = onCleanup(@() close(fig));
imagesc(result.delay_fs, result.energy_eV, result.normalized_signal.');
axis xy;
xlabel('XUV delay relative to IR (fs)');
ylabel('Electron kinetic energy (eV)');
title('SFEXIS: constant-dipole, single-point example');
colormap(parula);
bar = colorbar;
bar.Label.String = 'Globally normalized signal';
text(0.01,0.02,'Vismarra et al. (2022), DOI: 10.1088/2515-7647/ac7991', ...
    'Units','normalized','FontSize',8,'Color','w','Interpreter','none');
print(fig,png_file,'-dpng','-r150');
fprintf('Saved %d delays x %d energies to %s\n%s\n', ...
    numel(result.delay_fs),numel(result.energy_eV),output_directory,result.citation);
end
