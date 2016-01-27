function [seqs, filename] = fasta_in(filename)
% fasta_in(): Read fasta file into cell array. This prompts for a file if none
%	is provided.
% 
%   much fasta' now (1-9-13)
%   FASTA format specified here: http://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml

if ~exist('filename','var') || isempty(filename)
    [fname, pname] = uigetfile('*.fasta;*.afa;*.fna;*.fas;*.faa;*.ffn;*.frn',...
        'select fasta file');
    if ~fname, disp('no file chosen'), return, end
    filename = fullfile(pname, fname);
end

fid = fopen(filename, 'r');
z_cleanupObj = onCleanup(@() fclose(fid));

fcontent = fread(fid, '*char')';

expr = '(?<=(.*)?)(?<seqid>(?<=>)([^\f\n\r]*)?)(?<sequence>[^>]*)';
qq = regexp(fcontent, expr, 'tokens'); % parse string with regexp
seqs = cat(1, qq{:});
seqs(:,1) = regexprep(seqs(:,1), '^[\s]*','');
seqs(:,2) = regexprep(seqs(:,2), '\s','');