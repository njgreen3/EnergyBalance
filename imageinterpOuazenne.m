function [fittedX, fittedE] = imageinterpOuazenne()
% format longg;
% format compact;
% fontSize = 20;
nth_order = 7; %poly fit
derivArray = nth_order:-1:0;
rowShift = 425; %to remove axis from image
columnShift = 125; %to remove axis from image
I1 = 1*35;
Iend = 15*35; 
E1 = 180;
Eend = 22;

%===============================================================================
% Read in a image.
% folder = pwd;
folder = 'C:\Users\njgreen3\Documents\Thesis stuff\papers';
baseFileName = 'magnetizationCurveNoTrend.PNG';
% baseFileName = 'magCurveDandanMa.PNG';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileNameOnSearchPath = baseFileName; % No path this time.
	if ~exist(fullFileNameOnSearchPath, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
grayImage = imread(fullFileName);
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[~, ~, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
    
end
% % Display the original gray scale image.
% subplot(2, 2, 1);
% imshow(grayImage, []);
% title('Original Grayscale Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
% Binarize to get rid of horrible jpeg noise
binaryImage = grayImage < 100;
% % Display the original gray scale image.
% subplot(2, 2, 2);
% imshow(binaryImage, []);
% title('Binary Image', 'FontSize', fontSize);
% axis on;

% Get the rows (y) and columns (x).
[rows, columns] = find(binaryImage(1:rowShift,columnShift:end));

rows2Eslope = (Eend-E1)/(rows(1)-rows(end));
rows2Eoffset = (E1 - rows(end)*rows2Eslope);
columns2Islope = (Iend-I1)/(columns(end)-columns(1));
columns2Ioffset = (I1 - columns(1)*columns2Islope);

Im = columns*columns2Islope + columns2Ioffset;
E = rows*rows2Eslope + rows2Eoffset;

% Fit a nth_order
coefficients = polyfit(Im, E, nth_order); % Gets coefficients of the formula.
% Fit a curve to 500 points in the range that x has.
fittedI = linspace(min(Im), max(Im), 500);
% Now get the y values.
fittedE = polyval(coefficients, fittedI);
Xcoeffs = coefficients.*derivArray;
Xcoeffs(end) = [];
dEdI = polyval(Xcoeffs, fittedI);
fittedX = fittedE./fittedI;

figure(1)
plot(fittedX, fittedE, dEdI,fittedE)
xlabel('Xm')
ylabel('E')
legend('E/I', 'dE/dI')
% Plot the fitting:
% subplot(2,2,3:4);
figure(2)
plot(fittedI, fittedE,Im,E);
grid on;
xlabel('Im')%, 'FontSize', fontSize);
ylabel('E')%, 'FontSize', fontSize);
% % Overlay the original points in red.
% hold on;
% plot(Im, E, 'r.', 'LineWidth', .1, 'MarkerSize', .1);