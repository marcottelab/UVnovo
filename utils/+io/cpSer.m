function bytestream = cpSer(toCopy)
% bytestream = IO.CPSER(toCopy) copy serialized bytestream of input.
%	Useful for transfering into another Matlab session.
%	Use IO.PASTESER to paste serialized clipboard contents and translate into
%	original representation.
% 
%	I haven't tested this a lot. It works for at least simple matrices and
%	structs but may have weird results on other objects.
% 
%	See also IO.PASTESER

%%

bytestream = getByteStreamFromArray(toCopy);

com.mathworks.page.utils.ClipboardHandler().copy(num2str(bytestream))

if ~nargout
	clear bytestream
end
