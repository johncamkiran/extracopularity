function [E, k, S] = extracopularity(varargin)
%EXTRACOPULARITY Extracopularity coefficients for a 3D particle packing.
%
%   extracopularity(filename) returns extracopularity coefficients for the
%   particles described by a LAMMPS dump file with unscaled atomic coord-
%   inates. A copy of this file with a column for extracopularity coeff-
%   icients labelled "c_extra" is created at the path of the original.
%
%   extracopularity(S) returns extracopularity coefficients for the part-
%   icles described by an N x 3 coordinate matrix (a matrix whose N rows 
%   are coordinate vectors relative to some ordered basis).
%
%   E = extracopularity(...) returns a cell array of extracopularity coeff-
%   icients. If the input argument is a LAMMPS dump file with multiple
%   timesteps, the nth cell corresponds to the nth timestep.
%
%   [E, k] = extracopularity(...) returns cell arrays for extracopularity
%   coefficients and coordination numbers.
%
%   [E, k, S] = extracopularity(...) returns cell arrays for extracopul-
%   arity coefficients, coordination numbers, and coordinate matrices.
%
%   EXAMPLE:
%   [X, Y, Z] = meshgrid(1:20); S = [X(:), Y(:), Z(:)]; % defines system
%   E = extracopularity(S); % computes extracopularity coefficients
%   scatter3(S(:,1),S(:,2),S(:,3),[],E{1},'filled'); colorbar; % plots

%{
Copyright (c) 2022 John CAMKIRAN

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}

% ============================= PARSE INPUTS ==============================

% Print empty line
disp(' ');

% Throw error if no input argument
if nargin < 1
    error('Not enough input arguments.');
end

% Throw error if too many input arguments
if nargin > 1
    error('Too many input arguments.');
end

% Determine input type and throw error if invalid
if isnumeric(varargin{1})
    inputType = 'numeric';
else
    try char(varargin{1});
        inputType = 'string';
    catch
        error('Invalid data type. Argument must be a string or numeric.');
    end
end

% Proceed based on input type
switch inputType
    
    case 'string'
        
        % Read data and search for extracopularity coefficients
        filename = char(varargin{1});
        [S,dumpHeaders,dumpTables,~,E] = readLammpsDump(filename);
        
        % Return control if coefficients are found
        if ~isempty(E{1})
            return
        end
        
    case 'numeric'
        
        % Assign input to S
        S = varargin;
        
        % Throw error if input argument is complex
        if ~isequal(S{1},real(S{1}))
            error('Numeric argumnts must be real.');
        end
        
        % Throw error if input argument does not have 3 columns
        if size(S{1},2) ~= 3
            error('Numeric arguments must have 3 columns.');
        end
        
        % Throw error if input argument does not have more than 3 columns
        if size(S{1},1) <= size(S{1},2)
            error('Numeric arguments must have more than 3 rows.');
        end
        
end

% Throw error if system is degenerate to within machine precision
if min(svd(S{1})) < eps
    error('System must have 3 linearly independent dimensions.');
end

% ============================ COMPUTE OUTPUTS ============================

% Set number of perturbation samplings
NUM_SAMPLES = 8;

% Determine number of timesteps in input data
numTimesteps = numel(S);

% Initialize output cell arrays
E = cell(1,numTimesteps);
k = cell(1,numTimesteps);

% Precompute pair indices up to 14
pairs = nchoosek(1:14,2);

tic % records current time

% Analyze timesteps
for t = 1:numTimesteps
    
    % Display progress message
    disp(strcat("Analyzing timestep ",string(t)," of ", ...
        string(numTimesteps),"."));
    
    % Compute the partially robustified Voronoi adjacency matrix
    % (robustification is implicitly completed by computeHoodParms)
    A = computePartRobustVoronoi(S{t}, NUM_SAMPLES);
    
    % Determine the number of particles
    numParticles = size(S{t},1);
    
    % Initialize neighborhood parameters
    numDiffAngs = zeros(numParticles,1);
    numBonds = zeros(numParticles,1);
    
    % Create temporary variable to avoid large broadcast to parallel pool
    S_t = S{t};
    
    % Compute neighborhood parameters for each particle
    parfor i = 1:numParticles
        S_t_i = S_t;
        bonds = S_t_i(A(:,i),:) - S_t_i(i,:);
        [numBonds(i),numDiffAngs(i)]= computeHoodParms(bonds,pairs);
    end
    
    % Calculate the number of bond pairs for each neighborhood
    numBondPair = (numBonds.^2-numBonds)/2;
    
    % Bound numDiffAng from above by numBondPair (just in case)
    numDiffAngs = min(numDiffAngs,numBondPair);
    
    % Assign first output argument
    E{t} = log2(numBondPair) - log2(numDiffAngs);
    
    % Account for the trivial case (i.e. E := 0 if k = 1)
    E{t}(numBondPair == 0) = 0;
    
    % Assign second output argument
    k{t} = numBonds;
    
end

elapsedTime = toc; % measures elapsed time

% Display progress message
disp(strcat("Analysis completed in ",string(elapsedTime)," seconds."));

% Write to file (if particle position data was taken from a file)
if exist('dumpTables','var')
    writeLammpsDump(dumpTables,dumpHeaders,E,filename)
end

% =========================== TABULATE RESULTS ============================

% Calculate E for commonly encountered geometries
k_ref = [12 11 12 14 12 12  8 11 10 10  9  8  4  8  7  8  5 NaN  2];
m_ref = [ 3  3  4  6  6  7  3  6  5  6  5  4  1  5  4  6  3 NaN  1];
E_ref = transpose(log2((k_ref.^2-k_ref)./(2*m_ref)));

% Set registration precision
prc = 100*eps;

% Calculate initial, average, and final geometry fractions (known)
strucFrac{1} = 100*transpose(mean(abs(vertcat(E{1})-E_ref')<=prc));
strucFrac{2} = 100*transpose(mean(abs(vertcat(E{:})-E_ref')<=prc));
strucFrac{3} = 100*transpose(mean(abs(vertcat(E{end})-E_ref')<=prc));

% " (other)
strucFrac{1}(end-1,:)=100*mean(sum(abs(vertcat(E{1})-E_ref')<=prc,2)==0);
strucFrac{2}(end-1,:)=100*mean(sum(abs(vertcat(E{:})-E_ref')<=prc,2)==0);
strucFrac{3}(end-1,:)=100*mean(sum(abs(vertcat(E{end})-E_ref')<=prc,2)==0);

% Set table row names
rowNames = {'ICO' 'CPA' 'FCC' 'BCC' 'HCP' 'BPP' 'HXD,SA' 'CPP' 'BSP' ...
    'BSA,SC' 'CSA,CSP,TTP' 'HBP' 'TET,(9,6)' 'BTP' 'CTP,PBP' 'SDS' 'TBP'...
    'OTHER' 'TRIVIAL'};

% Set table column names
colNames = {'E' 'Initial (%)' 'Average (%)' 'Final (%)'};

% Construct table
tbl = table(E_ref,strucFrac{1},strucFrac{2},strucFrac{3}, ...
    'RowNames',rowNames,'VariableNames',colNames);

% Store existing display format
fmt = format;

% Temporarily change format for display purposes
format bank

% Display table
disp(' ');
disp(tbl);

% Revert to old format
format(fmt)
clear fmt

end

%==========================================================================
function [S,dumpHeaders,dumpTables,timeIdx,E] = readLammpsDump(filename)

% Read LAMMPS dump file
inputDumpId = fopen(filename,'r' );
if inputDumpId <= 0
    error(['No file found at the path ' filename]);
else
    disp('Reading LAMMPS dump file ...');
end
inputDump = fscanf(inputDumpId,'%c');
fclose(inputDumpId);

% Determine timestep indices
timeIdx = strfind(inputDump,'ITEM: TIMESTEP');
if isempty(timeIdx)
    error('File does not contain any timesteps.');
end
numTimeSteps = numel(timeIdx);
keyIndices   = [timeIdx numel(inputDump)];

% Initialize dump file parts
dumpHeaders = cell(1,numTimeSteps);
dumpTables  = cell(1,numTimeSteps);
S           = cell(1,numTimeSteps);
E           = cell(1,numTimeSteps);

for t = 1:numTimeSteps % cannot be parfor
    
    % Isolate region of inputDump corresponding to current timestep
    currentRegion = inputDump(keyIndices(t):keyIndices(t+1)-1);
    
    % Extract its header and body
    dumpHeaderEnd = strfind(currentRegion,'ITEM: ATOMS') + 11;
    dumpHeaders{t} = currentRegion(1:dumpHeaderEnd);
    dumpBody = currentRegion(dumpHeaderEnd+1:end);
    
    % Write as text to temporary file
    tempFileId = fopen('tmp.dump.txt','wt');
    fprintf(tempFileId, dumpBody);
    fclose(tempFileId);
    
    % Read as table from temporary file
    dumpTables{t} = readtable('tmp.dump.txt','Delimiter',' ');
    
    % Get table variable names
    columnNames = dumpTables{t}.Properties.VariableNames;
    
    % Search column names for unscaled coordinates
    coordinateColumns = strcmp('x',columnNames) ...
        | strcmp('y',columnNames) ...
        | strcmp('z',columnNames);
 
    % Read unscaled coordinates into variable S if found
    if sum(coordinateColumns)==3
        S{t} = table2array(dumpTables{t}(:,coordinateColumns));
    else
        error('Unscaled coordinates could not be found.');
    end
    
    % Search table for extracopularity coefficients
    extraColumn = strcmp('c_extra',columnNames);
    
    % Read extracopularity column into variable E if found
    if sum(extraColumn)==1
        E{t} = table2array(dumpTables{t}(:,extraColumn));
    else
        E{t} = [];
    end
    
    % Delete temporary file
    delete tmp.dump.txt
    
end

% Display progress message
disp('File successfully read.');

end

%==========================================================================
function writeLammpsDump(dumpTables,dumpHeaders,E,filename)

% Display progress message
disp('Writing results to LAMMPS dump file ...');

% Attempt to determine directory from filename
[directory, name, ext] = fileparts(filename);

% If unsuccesful use pwd
if isempty(directory)
    directory = pwd;
end

% Open file
newDumpName = [directory, '/', 'extra.', name,  ext];
outputDumpId = fopen(newDumpName,'wt');

for i = 1:numel(dumpTables) % cannot be parfor
    
    % Insert extracopularity coefficients to table
    newColumn = table(E{i},'VariableNames',{'c_extra'});
    newTable = [dumpTables{i}, newColumn];
    
    % Write as table to temporary file
    writetable(newTable,'tmp.dump.txt','Delimiter',' ');
    
    % Read as text from temporary file
    tempFileId  = fopen('tmp.dump.txt','r');
    tempFile = fscanf(tempFileId,'%c');
    
    % Add header and write to open file in directory
    fprintf(outputDumpId,[dumpHeaders{i} tempFile]);
    
end

% Close file
fclose(outputDumpId);

% Delete temporary file
delete tmp.dump.txt

% Display progress message
disp('File successfully written.');

end

%==========================================================================
function [k,m] = computeHoodParms(bonds,pairs)
% obtainHoodParms(bonds,pairs) returns the number of significant bonds k
% and different bond angles m for the neighborhood specified by "bonds"

% Set RMSE cutoff level
RMSE_CUT = 8.5;

% Set neighborhood tolerances
TAU_PLUS_ONE_DOUBLE = 1.50;
TAU_PLUS_ONE_SINGLE = 1.36;
TAU_PLUS_ONE_HALF   = 1.20;

% Compute bond lengths
bondLengths = sqrt(sum(bonds.*bonds,2));

% Sort bonds by increasing length
[bondLengths,idx] = sort(bondLengths);
bonds = bonds(idx,:);

% Remove outliers
isOutlier = bondLengths>(2*bondLengths(1));
bonds(isOutlier,:) = [];
bondLengths(isOutlier,:) = [];

% Determine number of bonds
k_raw = size(bonds,1);

% Determine nearest neighbor distance
k_eff = min(2,k_raw);
nnDist = sum(bondLengths(1:k_eff))/k_eff;

% Compute half-shell coordination number
k_half_shell = sum(bondLengths<=TAU_PLUS_ONE_HALF*nnDist);

% Check for trivial case
if k_raw < 2
    k = 2;
    m = 1;
    return
end

% Initialize RMSE vector
RMSE = zeros(1,33);

% ================================== 14 ===================================

% Copy and truncate pair matrix
pairs14 = pairs;
pairs14(any(pairs14 > min(14,k_raw),2),:) = [];
angles = computeAngle( bonds(pairs14(:,1),:) , bonds(pairs14(:,2),:) );

% BCC (body-centered cubic) [M]
BCC = [54.7356; 70.5288; 90; 109.4712; 125.2644; 180];
[MAT,idx] = min(abs(angles-BCC'),[],2);

if numelunique(idx) == numel(BCC) && k_half_shell <= 14
    SQR = MAT.*MAT;
    RMSE(1) = sqrt(sum(SQR)/numel(SQR));
    hasRightAngle = true;
else
    RMSE(1) = inf;
    hasRightAngle = false;
end

% ================================== 12 ===================================

% Copy and truncate pair matrix
pairs12 = pairs14;
pairs12(any(pairs12 > min(12,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs12(:,1),:),bonds(pairs12(:,2),:) );

% FCC (face-centered cubic) [M]
FCC = [60; 90; 120; 180];
[MAT,idx] = min(abs(angles-FCC'),[],2);

if numelunique(idx) == numel(FCC) && k_half_shell > 8 && k_half_shell <= 12
    SQR = MAT.*MAT;
    RMSE(2) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(2) = inf;
end

% HCP (hexagonal close-packed) [M, M, M*]
HCP{1} = [60; 90; 109.4712; 120; 146.4427; 180]; % 1.633
HCP{2} = [59.5133; 61.4601; 90; 107.6796; 120.4867; 145.6822; 180]; % 1.580
HCP{3} = [60; 90; 111.187; 119.5322; 147.1774; 180]; % 1.686

for ii = 1:numel(HCP)
    [MAT,idx] = min(abs(angles-HCP{ii}'),[],2);
    
    if numelunique(idx) == numel(HCP{ii}) && k_half_shell <= 12
        SQR = MAT.*MAT;
        RMSE(3+ii-1) = sqrt(sum(SQR)/numel(SQR));
    else
        RMSE(3+ii-1) = inf;
    end
end

% ICO (regular icosahedral) [M]
ICO = [63.4349; 116.5651; 180];

[MAT,idx] = min(abs(angles-ICO'),[],2);

if numelunique(idx) == numel(ICO) && k_half_shell <= 12
    SQR = MAT.*MAT;
    RMSE(6) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(6) = inf;
end

% BPP (bicapped pentagonal prismatic) [M*]
BPP = [60.4798; 91.0450; 110.9015; 120; 148.9550; 180];
[MAT,idx] = min(abs(angles-BPP'),[],2);

if numelunique(idx) == numel(BPP) && k_half_shell <= 12
    SQR = MAT.*MAT;
    RMSE(7) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(7) = inf;
end

% ================================== 11 ===================================

% Copy and truncate pair matrix
pairs11 = pairs12;
pairs11(any(pairs11 > min(11,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs11(:,1),:),bonds(pairs11(:,2),:) );

% CPP (capped pentagonal prismatic) [M*]
CPP = [60.5997; 91.0450; 110.9015; 120; 148.9550];
[MAT,idx] = min(abs(angles-CPP'),[],2);

if numelunique(idx) == numel(CPP) && k_half_shell <= 11
    SQR = MAT.*MAT;
    RMSE(8) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(8) = inf;
end

% ================================== 10 ===================================

% Copy and truncate pair matrix
pairs10 = pairs11;
pairs10(any(pairs10 > min(10,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs10(:,1),:),bonds(pairs10(:,2),:) );

% BSA (bicapped square antiprismatic) [M, LT, LJ, KA]
BSA{1}=[59.2640; 74.8585; 118.5280; 120.7360; 141.5925; 180];
BSA{2}=[64.0024; 67.7238; 78.9224; 115.9989; 128.0017; 139.7621; 179.2616];
BSA{3}=[60.0213; 73.6873; 75.5415; 120.0000; 141.2822; 179.9993];
BSA{4}=[60.9462; 72.2714; 76.3590; 119.0538; 121.8924; 140.9123; 179.9917];

for ii = 1:numel(BSA)
    [MAT,idx] = min(abs(angles-BSA{ii}'),[],2);
    
    if numelunique(idx) == numel(BSA{ii}) && k_half_shell <= 10
        SQR = MAT.*MAT;
        RMSE(9+ii-1) = sqrt(sum(SQR)/numel(SQR));
        
    else
        RMSE(9+ii-1) = inf;
    end
end

if any(RMSE(9:12)<inf)
    hasStraightAngle = true;
else
    hasStraightAngle = false;
end

% BSP (bicapped square prismatic) [M]
BSP = [54.7356; 70.5288; 109.4712; 125.2644; 180];
[MAT,idx] = min(abs(angles-BSP'),[],2);

if numelunique(idx) == numel(BSP) && ~hasRightAngle && k_half_shell <= 10
    SQR = MAT.*MAT;
    RMSE(13) = sqrt(sum(SQR)/numel(SQR));
    hasTheseFiveAngles = true;
else
    RMSE(13) = inf;
    hasTheseFiveAngles = false;
end

% =================================== 9 ===================================

% Copy and truncate pair matrix
pairs9 = pairs10;
pairs9(any(pairs9 > min(9,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs9(:,1),:),bonds(pairs9(:,2),:) );

% CSA (capped square antiprismatic) [M, LJ, KA*]
CSA{1}=[59.2641; 74.8585; 118.5283; 120.7359; 141.5925];
CSA{2}=[61.6446; 73.9485; 76.9638; 115.3109; 122.6594; 141.2438];
CSA{3}=[61.8542; 72.0865; 76.4563; 120; 123.7084; 140.8359];

for ii = 1:numel(CSA)
    [MAT,idx] = min(abs(angles-CSA{ii}'),[],2);
    
    if numelunique(idx) == numel(CSA{ii}) ...
            && ~hasStraightAngle  && k_half_shell <= 9
        
        SQR = MAT.*MAT;
        RMSE(14+ii-1) = sqrt(sum(SQR)/numel(SQR));
        
    else
        RMSE(14+ii-1) = inf;
    end
end

if any(RMSE(14:16)<inf)
    has60DegAngle = true;
else
    has60DegAngle = false;
end

% TTP (tricapped trigonal prismatic) [M, LT, LJ]
TTP{1} = [67.7923; 81.7868; 120; 135.5847; 139.1066];
TTP{2} = [68.8826; 77.2230; 87.7962; 120; 136.1017; 137.7634];
TTP{3} = [70.4623; 96.8836; 120; 131.5582; 141.2572];

for ii = 1:numel(TTP)
    [MAT,idx] = min(abs(angles-TTP{ii}'),[],2);
    
    if numelunique(idx) == numel(TTP{ii}) && k_half_shell <= 9
        SQR = MAT.*MAT;
        RMSE(17+ii-1) = sqrt(sum(SQR)/numel(SQR));
    else
        RMSE(17+ii-1) = inf;
    end
end

% =================================== 8 ===================================

% Copy and truncate pair matrix
pairs8 = pairs9;
pairs8(any(pairs8 > min(8,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs8(:,1),:),bonds(pairs8(:,2),:) );

% HXD (regular hexahedral) [M]
HXD = [70.5288; 109.4712; 180];
[MAT,idx] = min(abs(angles-HXD'),[],2);

if numelunique(idx) == numel(HXD) ...
        && ~hasRightAngle && ~hasTheseFiveAngles && k_half_shell <= 8
    
    SQR = MAT.*MAT;
    RMSE(20) = sqrt(sum(SQR)/numel(SQR));
    
else
    RMSE(20) = inf;
end

% SA (square antiprismatic) [M, LJ]
SA{1} = [74.8585; 118.5280; 141.5925];
SA{2} = [73.4274; 75.6923; 120.3805; 141.2138];

for ii = 1:numel(SA)
    [MAT,idx] = min(abs(angles-SA{ii}'),[],2);
    
    if numelunique(idx) == numel(SA{ii}) ...
            && ~has60DegAngle && ~hasStraightAngle && k_half_shell <= 8
        
        SQR = MAT.*MAT;
        RMSE(21+ii-1) = sqrt(sum(SQR)/numel(SQR));
        
    else
        RMSE(21+ii-1) = inf;
    end
end

% BTP (bicapped trigonal prismatic) [M, LJ**]
BTP{1} = [67.7923; 81.7868; 120; 135.5847; 139.1066];
BTP{2} =[69.1846; 79.8549; 82.5099; 85.5736; 122.4174; 130.3709; 138.9367];

for ii = 1:numel(BTP)
    [MAT,idx] = min(abs(angles-BTP{ii}'),[],2);
    
    if numelunique(idx) == numel(BTP{ii}) ...
            && ~has60DegAngle && k_half_shell <= 8
        
        SQR = MAT.*MAT;
        RMSE(23+ii-1) = sqrt(sum(SQR)/numel(SQR));
        
    else
        RMSE(23+ii-1) = inf;
    end
    
end

% HBP (hexagonal bipyramidal) [M]
HBP = [60; 90; 120; 180];

[MAT,idx] = min(abs(angles-HBP'),[],2);

if numelunique(idx)==numel(HBP) && k_half_shell <= 8
    SQR = MAT.*MAT;
    RMSE(25) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(25) = inf;
end

% SDS (snub disphenoidal) [M, LJ]
SDS{1} = [65.0604; 75.1576; 95.2967; 135.3026; 140.2180; 144.6244];
SDS{2} = [56.9514; 63.9619; 80.3951; 87.3046; 91.7857; 118.8147; ...
    121.9861; 132.4627; 160.4156; 172.9793];

for ii = 1:numel(SDS)
    [MAT,idx] = min(abs(angles-SDS{ii}'),[],2);
    
    if numelunique(idx)==numel(SDS{ii}) && k_half_shell <= 8
        
        SQR = MAT.*MAT;
        RMSE(26+ii-1) = sqrt(sum(SQR)/numel(SQR));
        
    else
        RMSE(26+ii-1) = inf;
    end
    
end

% =================================== 7 ===================================

% Copy and truncate pair matrix
pairs7 = pairs8;
pairs7(any(pairs7 > min(7,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs7(:,1),:),bonds(pairs7(:,2),:) );

% PBP (pentagonal bipyramidal) [M]
PBP = [72; 90; 144; 180];

[MAT,idx] = min(abs(angles-PBP'),[],2);

if numelunique(idx) == numel(PBP) && k_half_shell <= 7
    SQR = MAT.*MAT;
    RMSE(28) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(28) = inf;
end

% CTP (capped trigonal prismatic) [M, LJ]
CTP{1}=[67.7923; 81.7868; 135.5847; 139.1066];
CTP{2}=[73.3975; 77.0703; 78.9191; 83.7121; 86.9184; 129.7326; 141.4649;...
    146.7951];

for ii = 1:numel(CTP)
    [MAT,idx] = min(abs(angles-CTP{ii}'),[],2);
    
    if numelunique(idx) == numel(CTP{ii}) && k_half_shell <= 7
        SQR = MAT.*MAT;
        RMSE(29+ii-1) = sqrt(sum(SQR)/numel(SQR));
    else
        RMSE(29+ii-1) = inf;
    end
    
end

% =================================== 6 ===================================

% Copy and truncate pair matrix
pairs6 = pairs7;
pairs6(any(pairs6 > min(6,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs6(:,1),:),bonds(pairs6(:,2),:) );

% SC (regular octahedral / simple cubic) [M]
SC = [90; 180];

[MAT,idx] = min(abs(angles-SC'),[],2);

if numelunique(idx) == numel(SC) && k_half_shell <= 6
    SQR = MAT.*MAT;
    RMSE(31) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(31) = inf;
end

% =================================== 5 ===================================

% Copy and truncate pair matrix
pairs5 = pairs6;
pairs5(any(pairs5 > min(5,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs5(:,1),:),bonds(pairs5(:,2),:) );

% TBP (trigonal bipyramidal) [M]
TBP = [90; 120; 180];
[MAT,idx] = min(abs(angles-TBP'),[],2);

if numelunique(idx) == numel(TBP) && k_half_shell <= 5
    SQR = MAT.*MAT;
    RMSE(32) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(32) = inf;
end

% =================================== 4 ===================================

% Copy and truncate pair matrix
pairs4 = pairs5;
pairs4(any(pairs4 > min(4,k_raw),2),:) = [];
angles = computeAngle(bonds(pairs4(:,1),:),bonds(pairs4(:,2),:) );

% TET (tetrahedral) [M]
TET = 109.4712;
[MAT,idx] = min(abs(angles-TET'),[],2);

if numelunique(idx) == numel(TET) && k_half_shell <= 4
    SQR = MAT.*MAT;
    RMSE(33) = sqrt(sum(SQR)/numel(SQR));
else
    RMSE(33) = inf;
end


% =========================================================================

% Apply naive neighbor restriction (for BCC)
isFar = bondLengths > TAU_PLUS_ONE_DOUBLE*nnDist;
bonds(isFar,:) = [];
bondLengths(isFar,:) = [];
k_double_shell = size(bonds,1);

% Apply naive neighbor restriction (for all others)
isFar = bondLengths > TAU_PLUS_ONE_SINGLE*nnDist;
bonds(isFar,:) = [];
k_single_shell = size(bonds,1);

% Set correction factors
CORR_FACT = [ 1.3155    1.0000    1.3830    1.4386    1.3729    0.7769 ...
    1.3804    1.3333    1.1576    1.1618    1.1256    1.1476 ...
    1.1346    1.0937    1.1362    1.0642    0.9654    0.9727 ...
    0.9791    0.7218    0.7837    0.8054    0.9655    0.9561 ...
    1.0032    1.0677    1.6948    0.7964    0.8164    0.8528 ...
    0.4405    0.5423    0.3337];

% Adjust RMSE accoding to bias correction factors
RMSE = CORR_FACT.*RMSE;

% Compute minimum error to commonly encountered geometries
[minRmse,minIdx] = min(RMSE);

% Assign outputs (coordination number and angle count)
if minRmse < RMSE_CUT
    
    % Assign directly for commonly encountered geometries
    switch minIdx
        case 1 % BCC (body-centered cubic)
            k = min(14,k_double_shell);
            m = 6;
        case 2 % FCC (face-centered cubic)
            k = min(12,k_single_shell);
            m = 4;
        case {3,4,5} % HCP (hexagonal close-packed)
            k = min(12,k_single_shell);
            m = 6;
        case 6 % ICO,CPA (regular icosahedral, capped pentagonal antiprism)
            kk = [12 11];
            [~,k_idx] = min(abs(kk-k_single_shell));
            k = kk(k_idx);
            m = 3;
        case 7 % BPP (bicapped pentagonal prismatic)
            k = 12;
            m = 7;
        case 8 % CPP (capped pentagonal prismatic)
            k = 11;
            m = 6;
        case {9,10,11,12} % BSA (bicapped square antiprismatic)
            k = 10;
            m = 6;
        case 13 % BSP,CSP [(bi)capped square prismatic)]
            kk = [10 9];
            [~,k_idx] = min(abs(kk-k_single_shell));
            k = kk(k_idx);
            m = 5;
        case {14,15,16} % CSA (capped square antiprismatic)
            k = 9;
            m = 5;
        case {17,18,19} % TTP (tricapped trigonal prismatic)
            k = 9;
            m = 5;
        case 20 % HXD (regular hexahedral)
            k = 8;
            m = 3;
        case {21,22} % SA (square antiprismatic)
            k = 8;
            m = 3;
        case {23,24} % BTP (bicapped trigonal antiprismatic)
            k = 8;
            m =  5;
        case 25 % HBP (hexagonal bipyramidal)
            k = min(8,k_single_shell);
            m =  4;
        case {26,27} % SDS (snub disphenoidal)
            k = 8;
            m =  6;
        case 28 % PBP (pentagonal bipyramidal)
            k = 7;
            m = 4;
        case {29,30} % CTP (capped trigonal prismatic)
            k = 7;
            m = 4;
        case 31 % SC (regular octahedral / simple cubic)
            k = min(6,k_single_shell);
            m =  2;
        case 32 % TBP (trigonal bipyramidal)
            k = min(5,k_single_shell);
            m =  3;
        case 33 % TET (regular tetrahedral)
            k = min(4,k_single_shell);
            m =  1;
    end
    
else
    
    % Apply bond angle count estimator if geometry unrecognized
    if k_half_shell < 2 % trivial case
        k = 2;
        m = 1;
    else % ordinary case
        a = 11.407186472007339;
        b =  7.334483343628186;
        gamma = 3.4;
        k = min(12,k_half_shell);
        m = b - a*exp(-k/gamma);
    end
    
end

% =========================================================================

% Fast function for getting number of unique elements in an array
    function N = numelunique(X)
        N = nnz(diff(sort(X)))+1;
    end

% Function for computing angle between the row vectors of A and B
    function angles = computeAngle(A,B)
        angles = acosd(sum(A.*B,2) ./ sqrt(sum(A.*A,2).*sum(B.*B,2)));
    end

end

%==========================================================================
function A = computePartRobustVoronoi(S, M)
% computePartRobustVoronoi(S, M) returns the adjacency matrix implied by
% the partially robustified Voronoi tessellation of the system S, computed
% numerically from M perturbation samplings

% Set perturbation scale multiplier
PERTURBATION_SCALE_MULTIPLIER = 0.04; % = 1/25

% Determine the number of particles
numParticles = size(S,1);

% Compute the Voronoi adjacency matrix implied by the unperturbed system
try
    A0 = del2adj(delaunay(S));
catch
    S = addGaussianNoise(S, eps);
    A0 = del2adj(delaunay(S));
end

% Initalize nearest-neighbor distances
nnDist = zeros(numParticles,1);

% Compute nearest-neighbor distances
parfor i = 1:numParticles
    nnDist(i) = min(sqrt(sum( ( S(i,:) - S(A0(:,i),:) ).^2, 2)));
end

% Compute the mean nearest-neighbor distance
fullVector = full(nnDist);
avgNnDist = sum(fullVector)/numel(fullVector);

if M > 1
    
    % Determine perturbation scale
    perturbationScale = PERTURBATION_SCALE_MULTIPLIER*avgNnDist;
    
    % Initialize cell array for sample adjacency matrices
    AA = cell(M-1,1);
    
    % Sample adjacency matrices
    parfor m = 1:(M-1)
        S_perturbed = addGaussianNoise(S, perturbationScale);
        thisTriangulation = delaunay(S_perturbed);
        AA{m} = del2adj(thisTriangulation); %#ok<PFOUS>
    end
    
    % Append initial estimate of A to end of the cell array
    AA{end+1} = A0;
    
    % Take majority vote
    A = takeMajorityVote(AA);
    
else
    
    % Assign initial estimate to output
    A = A0;
    
end


end

%==========================================================================
function A = del2adj(T)
% del2adj(T) returns the Voronoi adjacency matrix A implied by the
% Delaunay triangulation matrix T

numParticles = max(T(:));
numDimensions = size(T,2);

% Initialize a sparse logical matrix
A = logical(sparse([], [], [], numParticles, numParticles));

% Populate the matrix
for i = 1:numDimensions
    x = mod(i-1, numDimensions)+1;
    y = mod(i, numDimensions)+1;
    A = A | sparse(T(:,x), T(:,y),1, numParticles, numParticles);
end

% Symmetrize the matrix
A = A | A';

end

%==========================================================================
function Y = addGaussianNoise(X,noiseScale)
% addGaussianNoise(X,noiseScale) returns a noisy version Y of the coord-
% inate matrix X

Y = X + noiseScale*randn(size(X));

end

%==========================================================================
function A = takeMajorityVote(AA)
% takeMajorityVote(AA) returns the adjacency matrix A that represents the
% majority "opinion" of the adjacency matrices in the cell array AA

% Determine the number of particles and samplings
numParticles = length(AA{1});
numSamplings = length(AA);

% Initialize output
A = logical(sparse([], [], [], numParticles, numParticles));

% Sum adjacency matrices
parfor i = 1:numSamplings
    A = A + AA{i};
end

% Determine the number of votes needed for a majority
majorityVote = 0.5*numSamplings;

% Overwrite each element of A to "true" if its value exceeds majority vote
% and to "false" otherwise
A = A > majorityVote;

end

