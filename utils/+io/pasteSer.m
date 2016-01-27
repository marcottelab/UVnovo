function varOut = pasteSer()
% varOut = IO.PASTESER()
%	Import serialized clipboard data that was copied using IO.CPSER
%	Useful for transfering from another Matlab session.
% 
%	I haven't tested this a lot. It works for at least transfering simple
%	self-contained data between sessions -- matrices, basic structs and such --
%	but this may have weird results on other objects.
% 
%	See also IO.CPSER

%%
bytestream = char(com.mathworks.page.utils.ClipboardHandler().paste);
varOut = getArrayFromByteStream( uint8( str2doublez( strsplit(bytestream))));

%%% I had this previously, not sure why...
% varOut = getArrayFromByteStream( uint8( str2doubleq( strsplit( bytestream(2:end-1)))));


