
classdef CONSTS
	properties (Constant)
	
	%%% PARTICLE MASSES, Daltons
		% proton
		mProton = 1.007276466879; % NIST 2014
		% electron
		mElectron = 0.00548579909070; % NIST 2014
		% neutron
		mNeutron = 1.00866491588; % NIST 2014
	
	%%% MONOISOTOPIC MASSES, Daltons
		% hydrogen
		mH = 1.0078250321;
		% carbon
		mC = 12.00000000;
		% nitrogen, http://www.nndc.bnl.gov/masses/mass.mas03round
		mN = 14.0030740048;
		% oxygen
		mO = 15.99491463;
		% sulphur, http://www.nndc.bnl.gov/masses/mass.mas03round
		mS = 31.97207100;
		% water
		mH2O = 18.0105646942;
		% mass CO (Unimod)
		mCO  = 27.994915;
		% mass NH (Unimod)
		mNH  = 15.010899;
		% mass OH
		mOH  = 17.002740;
	
	%%% AVERAGE MASSES, Daltons
		% hydrogen, estimated from IUPAC 2009 atomic weight report
		mH_avg = 1.00797;
		% oxygen (Unimod)
		mO_avg = 15.9994;
		% water
		mH2O_avg = 18.0152833;
		% avg mass CO (http://library.med.utah.edu/masspec/mole.htm)
		mCO_avg = 28.01039;
		% avg mass NH (http://library.med.utah.edu/masspec/mole.htm)
		mNH_avg = 15.01468; 
		% avg mass OH
		mOH_avg = 17.00734; 
	
	
	%%% MONOISOTOPIC AMINO ACID MASSES, Daltons
	% Selenocysteine isotopes are common and all over the place. Not included now. http://www.ncbi.nlm.nih.gov/pubmed/22779694
		mAA = struct( ...
			'G',   57.021463715, ...
			'A',   71.037113779, ...
			'S',   87.032028401, ...
			'P',   97.052763843, ...
			'V',   99.068413907, ...
			'T',  101.047678465, ...
			'C',  103.009184469, ...
			'L',  113.084063971, ...
			'I',  113.084063971, ...
			'N',  114.042927438, ...
			'D',  115.026943023, ...
			'Q',  128.058577502, ...
			'K',  128.094963008, ...
			'E',  129.042593087, ...
			'M',  131.040484597, ...
			'H',  137.058911853, ...
			'F',  147.068413907, ...
			'R',  156.101111018, ...
			'Y',  163.063328529, ...
			'W',  186.079312944  ...
			);
	
	%%% MONOISOTOPIC PTM MASSES, Daltons
	% 	Common post translational modification masses defined at www.unimod.org
		mPTM = struct( ...
			'AMCA',      215.058243, ... 7-amino-4-methylcoumarin 3-acetic acid (C12H11NO4 - H2O)
			'Acetyl',    42.010565,  ... Unimod # 1
			'Amidated',  -0.984016,  ... Unimod # 2
			'Carbamyl',  43.005814,  ... Unimod # 5
			'Oxidation', 15.994915,  ... Unimod # 35
			'Carbamidomethyl', 57.021464, ... Unimod # 4. Iodoacetamide derivative (Cys alkylation)
			'Ethanolyl', 44.026215   ... Unimod # 278. Iodoethanol derivative (Cys alkylation)
			);
	
	%%% Peptide and AA constants
		% unit_g: mean Daltons per nucleon
		unit_g = 1.000468;
		%	This is for normalization similar to Kendrick mass conversion.
		%	https://en.wikipedia.org/wiki/Kendrick_mass
		%	See function estUnitG.m to calculate a data dependent unit_g.
		%	With 2400 E coli PSMs, it estimated 1.000488. Pretty close!
		% This is similar to the 'averagine scale'.
		%	Bajrami B, et.al. J Am Soc Mass Spectrom 2011, 20:2124–2134.
		
		% Averagine (the average amino acid)
		%	C:4.9384, H:7.7583, N:1.3577, O:1.4773, S:0.0417
		mAveragine = 111.031924816757; % monoisotopic mass
	
	end
end

% REFERENCES:
%	NIST 2014: http://physics.nist.gov/cuu/Constants/index.html
%	molecular properties, AA masses: http://www.chemicalize.org/
%	elemental masses & isotopic abundances:
%		http://www.nndc.bnl.gov/masses/mass.mas03round
%		http://www.sisweb.com/referenc/source/exactmaa.htm

