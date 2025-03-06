% IBD PM-13
% Optimization Script - Vc, Vp, CLp, and Qcp are being estimated!


%% Clear Workspace
clear %Removes all variables from the workspace.
close all % Closes all open figure windows.
clc %Clears the command window.
warning('off') %Turns off warnings (use with caution, as this may suppress important messages).
disp('1. Workspace cleared, ready for a fresh run.✅');

%% LOAD PROJECT FILE

% Define the path and name of the model
model_path = "G:\Shared drives\Pioneering Medicine (PM-13) - IBD\05 Team Member Spaces\Jeevan working folder\";
model_name = "IBD PM-13 PK Profiles.sbproj";
% Load the SimBiology project
proj = sbioloadproject(strcat(model_path, model_name)); %sbioloadproject is a MATLAB toolbox for biological systems modeling.
disp('2. project is loaded successfully!✅');
disp(proj);

%% SELECT MODEL

% Select the model
specificModelIndex = 1;  % Assuming the first model is to be used
modelFields = fieldnames(proj);
model = proj.(modelFields{specificModelIndex});
disp(['3. Model of interest is selected. It is :  ', model.Name, '✅']);

%% Model Properties

disp('4. Model properties and their values - Response species and the corresponding compartment is selected for responseMap');
modelProperties = get(model);
%disp(modelProperties);

% compartment
response_species_compartment = modelProperties.Compartments(1).Name;
disp([' --4a. Response Species Compartment :',response_species_compartment]);

% species
%response_species = modelProperties.Species;
response_species = modelProperties.Species(2).Name;
disp([' --4b. Response Species :', response_species]);

%% LOAD DATA (NOTE: CONSIDERING THE DIRECTORY IS THE SAME AS THE PROJECT FILE DIRECTORY)

% Load experimental data
dataFile = 'IBD PM-13 PK data -3 Tables - Fitting.xlsx'; % Replace with the actual file name
dataTable = readtable(dataFile);
disp(['5. Excel file is loaded ✅'])
disp([' --5a. Column Names (for verification) :',dataTable.Properties.VariableNames]);

%% CHANGE MODEL UNITS 

% Change model time units to days
configset = getconfigset(model);
X_axis_Units = 'hour';
set(configset, 'TimeUnits', X_axis_Units); % Set model's time units to "day"
disp('6. Model configuration is done ! ✅');
disp(' --6a. X_axis Units are declared (hour).')

%% ENABLE UNIT CONVERSION

% Enable unit conversion
configset.CompileOptions.UnitConversion = true;
disp('7. Unit Conversion is set to True! ✅ ');


%% GROUPED DATA

% Convert to groupedData object
if ~all(ismember({'Mouse_', 'time_hrs_', 'Conc_ng_mL_'}, dataTable.Properties.VariableNames))
    error('The required columns "ID", "Time_Day_", or "DrugConc_ug_mL_" are missing in the data file.');
end

%Creating a groupedDataObj
groupedDataObj = groupedData(dataTable, 'Mouse_', 'time_hrs_'); %While assigning GroupVariableName and IndependentVariableName
disp('8. Data successfully converted to groupedData.✅');

disp(' -- 8a. Checking Variable Names (the column names)');
disp(groupedDataObj.Properties.VariableNames); %column names
groupedDataObj.Properties.VariableUnits = {'','ID','hour','nanogram/milliliter','','','','','milligram'};
disp(' -- 8b. Checking Assigned Units to the variables');
disp(groupedDataObj.Properties.VariableUnits); %column units
disp(' -- 8c. Checking GroupedDataObj properties');
disp(groupedDataObj.Properties);


dependent_variable_name = 'Conc_ng_mL_';
disp(' -- 8x. Dependent variable name is also assigned (for responseMap). ');

%% RESPONSE MAP (AUTOMATED) -  RESPONSE MAP (SPECIES = COLUMN_NAME)


% Make the relation
the_relation = strcat(response_species_compartment , '.' , response_species ,' ', ' = ', ' '  , dependent_variable_name); %----------------->
% Display it
%disp(the_relation); 
%the_relation = response_species_compartment + "." + response_species + " = " + independent_variable_name;
%responseMap = {'Cen.Drug_conc_cen = DrugConc_ug_mL_'};

% Map the response
responseMap = {the_relation};
disp(['9. Response Mapping is Done ✅. The relation is [',the_relation,'].']);

%% PARAMETERS TO BE ESTIMATED

% Define objects being estimated and their initial estimates.
estimatedInfoObj = estimatedInfo({'log(Vc)', 'log(Vp)', 'log(k_CLp)', 'log(k_Q)'});

estimatedInfoObj(1).InitialValue = 0.9;   % initial value of Vc
estimatedInfoObj(1).Bounds       = [0.09 9];   % Untransformed lower & upper bounds for Vc
estimatedInfoObj(1).CategoryVariableName = 'GroupID';
%estimatedInfoObj(1).GroupVariableName = 'Mouse_';

estimatedInfoObj(2).InitialValue = 1.5;       % initial value of Vp
estimatedInfoObj(2).Bounds       = [0.15 15]; % Untransformed lower & upper bounds for Vp
estimatedInfoObj(2).CategoryVariableName = 'GroupID';

estimatedInfoObj(3).InitialValue = 0.004;         % initial value of k_CLp
estimatedInfoObj(3).Bounds       = [0.00004 0.4]; % Untransformed lower & upper bounds for k_CLp
estimatedInfoObj(3).CategoryVariableName = 'GroupID';

estimatedInfoObj(4).InitialValue = 0.03;       % initial value of k_Q
estimatedInfoObj(4).Bounds       = [0.0003 3];  % Untransformed lower & upper bounds for k_Q
estimatedInfoObj(4).CategoryVariableName = 'GroupID';

disp('10. Displaying the Parameters to be estimated and their associated data: ');
disp(estimatedInfoObj);


%% DOSES

% DOCUMENTATION : if the input groupedData object has dosing information, 
% you can use the 'createDoses' method to construct doses.

% Create the data dose for Cen.Dose_Cen.
doses1 = sbiodose('DoseIV_mg_');
doses1.TargetName = 'Vc.Vc_Drug_mg';
doses1 = createDoses(groupedDataObj, 'DoseIV_mg_', '', doses1);


% Convert doses1 to a cell array (still needed if you plan to process doses as cells).
doses1 = num2cell(doses1);

% If no SC dosing, dosesForFit is directly built from doses1.
dosesForFit = cell(size(doses1, 1), 1);
for i = 1:numel(dosesForFit)
    dosesForFit{i} = doses1{i}; % Directly assign doses from doses1
end

%disp(dosesForFit);
disp('11. Doses: Defined using sbiodose ✅ (No SC dosing)');

%% DOSES INFO IN TABULAR FORMAT

%Displaying Doses in Tabular Format
% Assuming dosesForFit is already defined and has the properties

% Initialize variables to store the data
numRows = size(dosesForFit, 1);
Names = cell(numRows, 1);
Amounts = zeros(numRows, 1);
AmountUnits = cell(numRows, 1);
TargetNames = cell(numRows, 1);
TimeUnits = cell(numRows, 1);

% Extract data from dosesForFit
for i = 1:numRows
    Names{i} = dosesForFit{i,1}.Name;
    Amounts(i) = dosesForFit{i,1}.Amount;
    AmountUnits{i} = dosesForFit{i,1}.AmountUnits;
    TargetNames{i} = dosesForFit{i,1}.TargetName;
    TimeUnits{i} = dosesForFit{i,1}.TimeUnits;
end

% Create the table
doseTable = table(Names, Amounts, AmountUnits, TargetNames, TimeUnits, ...
    'VariableNames', {'Name', 'Amount', 'AmountUnits', 'TargetName', 'TimeUnits'});

% Display the table
disp(' --11a. Displaying the doses data in a tabular format (for cross-checking).');
disp(doseTable);


%% OPTIONS [NOT IMPORTANT RIGHT NOW. THE DEFAULT ONES ARE SUITABLE]

% Define Algorithm options.
options                   = optimoptions('particleswarm');
options.FunctionTolerance = 1e-08;
options.MaxIterations     = 200;

disp('12. No issues until Options')

%% FIT 

% Define fit problem.
f              = fitproblem('FitFunction', 'sbiofit');
f.Data         = groupedDataObj;
f.Model        = model;
f.Estimated    = estimatedInfoObj;
f.ResponseMap  = responseMap;
f.ErrorModel   = 'exponential';
f.Doses        = dosesForFit;
f.FunctionName = 'particleswarm';
f.Options      = options;
f.ProgressPlot = true;
f.UseParallel  = true;
f.Pooled       = false;

% Estimate parameter values.
[results, simdataI] = f.fit;

% Configure the names of the resulting simulation data.
for i = 1:length(simdataI)
    simdataI(i).Name = 'SimData Individual';
end

% Assign output arguments.
args.output.results  = results;
args.output.simdataI = simdataI;

%% CREATING CV COLUMN IN RESULTS

% Assuming 'results.ParameterEstimates' is the name of your table
ParameterEstimates = results.ParameterEstimates;

% Calculate CV as StandardError / Estimate
ParameterEstimates.CV = (ParameterEstimates.StandardError ./ ParameterEstimates.Estimate)*100;

% Display the updated table
disp(ParameterEstimates);