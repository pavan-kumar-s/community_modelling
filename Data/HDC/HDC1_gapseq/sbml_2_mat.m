% matlab code to convert all the sbml files in the folder to mat files

% get all the sbml files in the folder
sbml_files = dir('*.xml');
sbml_files = {sbml_files(1:end).name};
for file =1:numel(sbml_files)
    model = readCbModel(sbml_files{file});
    save(replace(sbml_files{file},'.xml','.mat'),'model')
end