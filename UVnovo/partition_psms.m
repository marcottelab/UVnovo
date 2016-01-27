function partitions = partition_psms(psmData, nPartitions)
% PARTITION_PSMS partitions psms into independent sets for cross validation.
%	This is based on amino acid sequence and ignores PTM differences.
%	A random hash for each psm sequence is calculated in import_psms(). We
%	divide those hashes, by sorted order, into same-sized sets. This method for
%	partitioning is mostly stable in that most sequences will be assigned to the
%	same partition regardless of random differences in the total sequence set.


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


