function [fittedX, fittedE] = imageinterp()
% format longg;
% format compact;
% fontSize = 20;
nth_order = 7; %poly fit
derivArray = nth_order:-1:0;
rowShift = 615; %to remove axis from image
columnShift = 99; %to remove axis from image
I1 = .0333;
Iend = 6.5; 
E1 = 300;
Eend = .0714;

%===============================================================================
% Read in a image.
% folder = pwd;
folder = 'C:\Users\njgreen3\Documents\Thesis stuff\papers';
baseFileName = 'magCurveDandanMa.PNG';
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
    % I want the blue curve. In the blue channel both blue and white are
    % equal to 255. In all the channels, the black grid lines and labels
    % are the same (about 0-100 depending on the thickness). I can get only
    % the blue curve by grabbing the values that aren't equal across
    % channels.
	greenImage = grayImage(:, :, 2); % Take green channel.
    blueImage = grayImage(:,:,3); %Take blue channel.
    
    binaryImage = blueImage;
    binaryImage(blueImage == greenImage) = 0;
    
end
% % Display the original gray scale image.
% figure(1)
% %subplot(1, 2, 1);
% imshow(grayImage, []);
% title('Original Grayscale Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Give a name to the title bar.
% % set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

% % Display the original gray scale image.
% figure(2)
% % subplot(1, 2, 2);
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
% figure(3)
% plot(fittedX, fittedE, dEdI,fittedE)
% xlabel('Xm', 'FontSize', fontSize);
% ylabel('E', 'FontSize', fontSize);
% legend('E/I', 'dE/dI')
% % Plot the fitting:
% % subplot(2,2,3:4);
% figure(4)
% plot(fittedI, fittedE,fittedI,fittedI*.5, 'b-', 'linewidth', 4);
% grid on;
% xlabel('Im', 'FontSize', fontSize);
% ylabel('E', 'FontSize', fontSize);
% % Overlay the original points in red.
% hold on;
% plot(Im, E, 'r.', 'LineWidth', .1, 'MarkerSize', .1);
% figure(5)
% plot(fittedI,fittedX/(2*pi*50),fittedI,dEdI/(2*pi*50))