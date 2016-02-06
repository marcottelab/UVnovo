function oobErrorVals = plotOOBError(ensTB)
% Calculate and plot out-of-bag error vs. number of trees grown.

fprintf(1,'Calculating OOB error ... ')
oobErrorVals = oobError(ensTB);
fprintf(1,'done.\n')

ntrees = ensTB.NTrees;

figure;
plot(oobErrorVals, 'color', rgb('hotpink'), 'LineWidth',2, 'LineSmoothing','on')
set(gcf,'Renderer','OpenGL')
xlim([ min(max(0,ntrees-30),50), ntrees ])
xlabel('Number of Grown Trees');
ylabel('OOB Classification Error');
text(.95,.95, sprintf( 'final OOB error: %f\nnumber of predictors: %g', ...
	oobErrorVals(end), numel(ensTB.VarNames)), ...
	'Units','normalized', 'fontname','Consolas', ...
	'VerticalAlignment','Cap', 'HorizontalAlignment','Right')
set(gca, 'LooseInset', [0,0,0.05,0]);

if nargout == 0
	clear oobErrorVals
end
