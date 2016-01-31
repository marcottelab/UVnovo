function UVnovo_partition(Meta)
% UVNOVO_PARTITION imports, partitions, and saves spectra & psm data.
% 
%	Vars are saved in a serialized matlab format. The is much more efficient in
%	time and space than standard matlab files, but io.loadSer() is required for
%	opening such files. See io.saveSer for more info.
% 


%% Import spectra and psm files.

msData = import_spectra(Meta.paths.spectra.path, Meta.params.pre.import_spectra);
psmData = import_psms(Meta.paths.psms.path);

% Partition for cross-validation.
nPartitions = Meta.params.pre.crossVal.nPartitions;
partitions = partition_psms(psmData, nPartitions);

% Index between psmData.scans & msData.scans. Assign later to same partitions.
mapi = cross_index(psmData, msData);
for i = 1:nPartitions
	partitions(i).test.msData =  mapi.psmData.msData(partitions(i).test.psmData);
	partitions(i).train.msData = mapi.psmData.msData(partitions(i).train.psmData);
end

% Create and save training and test sets.
%	Initialize variables for saving. We replace 'msData' and 'psmData' with
%	partition-specific data so that the variables are named consistently upon
%	reloading the data.
if ~exist('msData_ORIG', 'var'), msData_ORIG =  msData;  end
if ~exist('psmData_ORIG','var'), psmData_ORIG = psmData; end
msData.scans =  [];
psmData.scans = [];


for n = 1:nPartitions
	partitionDir = fullfile(Meta.paths.exp.path, sprintf('part%g', n));
	mkdir(partitionDir)
	
	% Create and save test partition vars.
	% @TODO save test partition as MS2 file. (?)
	Meta.partition.n = n;
	Meta.partition.type = 'test';
	msData.scans = msData_ORIG.scans( partitions(n).test.msData );
	
	fn_test = fullfile(partitionDir, sprintf('p%d_Test.mat', n));
	
	io.saveSer(fn_test, 'msData', 'Meta')
	
	% Create and save training partition vars.
	Meta.partition.type = 'train';
	msData.scans =  msData_ORIG.scans(  partitions(n).train.msData  );
	psmData.scans = psmData_ORIG.scans( partitions(n).train.psmData );
	
	fn_train = fullfile(partitionDir, sprintf('p%d_Train.mat', n));
	io.saveSer(fn_train, 'msData', 'Meta')
	
	Meta = rmfield(Meta, 'partition');
end

end %MAIN

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function partitions = partition_psms(psmData, nPartitions)
	% PARTITION_PSMS partitions psms into independent sets for cross validation.
	%	This is based on amino acid sequence and ignores PTM differences.
	%	A random hash for each psm sequence is calculated in import_psms(). We
	%	divide those hashes, by sorted order, into same-sized sets. This method
	%	for partitioning is mostly stable, in that most sequences will be 
	%	assigned to the same partition regardless of random differences in the
	%	overall sequence set.
	
	% Cannot create independent training and test sets without > 1 partitions.
	assert( nPartitions > 1, 'UVnovo:Partition:Singularity', '> 1 partition required')
	
	[h, ~, ic_hash] = unique([psmData.scans.seqHash]');
	partition = ceil( (1:numel(h)) / (numel(h)/nPartitions) )';
	
	%	map partition to psmData.scans
	nPsmScans = numel(psmData.scans);
	partition2psmData = accumarray(partition(ic_hash), 1:nPsmScans, [], @(x){x});
	
	partitions = struct;
	for i = 1:nPartitions
		testPart = false(nPartitions,1);
		testPart(i) = true;
		
		psmTest = sort(cat(1, partition2psmData{testPart}));
		partitions(i,1).test = struct( ...
			'nscan', [psmData.scans(psmTest).nscan]', ...
			'psmData', psmTest ...
			);
		
		psmTrain = sort(cat(1, partition2psmData{~testPart}));
		partitions(i,1).train = struct( ...
			'nscan', [psmData.scans(psmTrain).nscan]', ...
			'psmData', psmTrain ...
			);
	end
end

