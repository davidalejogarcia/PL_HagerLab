function [RegMatrix, fixedpts, movingpts] = TwoColorImageRegistration(I1,I2,TranType)

%Generates the Transformation information that is required to align input
%image I2 to input image I1;

%specify some constants
lp = 1;
hp = 5;
windowSz = 7;
dx = 5;
dy = 5;

max_disp = 10;
%Filter image 
I_filt1 = bpass(I1, lp, hp);
I_filt2 = bpass(I2, lp, hp);
%find threshold values using Otsu method for each image;
th1 = graythresh(I_filt1);
th1 = th1*(max(I_filt1(:)))/6;

th2 = graythresh(I_filt2);
th2 = th2*(max(I_filt2(:)))/6;

%Locate centroids
peaks1 = pkfnd(I_filt1, th1, hp);
if ~isempty(peaks1)
    Centroids1 = cntrd(I_filt1,peaks1,windowSz);
    
else
    Centroids1 = [];
end

peaks2 = pkfnd(I_filt2, th2, hp);
if ~isempty(peaks2)
    Centroids2 = cntrd(I_filt2,peaks2,windowSz);
else
    Centroids2 = [];
end
%fit the points to a psf
Particles1 = Centroids1;

sizeFrame = size(I1);
[x,y] = meshgrid(1:sizeFrame(2), 1:sizeFrame(1));

if ~isempty(Centroids1)
    
    img = double(I1);
    for i = 1:size(Centroids1,1)
        x_center = double(Centroids1(i,1));
        y_center = double(Centroids1(i,2));
        
        
        if x_center ~= 0
            
            
            x_center_round = round(x_center);
            y_center_round = round(y_center);
            y_sub = (y_center_round-dy:y_center_round+dy);
            x_sub = (x_center_round-dx:x_center_round+dx);
            
            img_sub = img(y_sub,x_sub);
            coordinates = {x_sub, y_sub};
            disp(['Fitting Gaussian to peak ', num2str(i)]);
            
            [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
            
            
            % Add fit parameters to the output
            Particles1(i, 8:13) = [parFit, ssr];
            %                 % Test plot
            %                 figure;
            %                 subplot(2,1,1);
            %                 imagesc(img);
            %                 hold on;
            %                 colormap(gray);
            %                 plot(PeaksInCurrFrame(10, iPeak),PeaksInCurrFrame(11, iPeak), 'or');
            %
            % Subtract found particle from image;
            exponent = ((x - parFit(4))/parFit(2)).^2 +...
                ((y - parFit(5))/parFit(2)).^2;
            
            img_del = parFit(1) * exp(-exponent/2);
            
            img =img - img_del;
        end
    end
end
            
Particles2 = Centroids2;
if ~isempty(Centroids2)
    
    img = double(I2);
    for i = 1:size(Centroids2,1)
        x_center = double(Centroids2(i,1));
        y_center = double(Centroids2(i,2));
        
        
        if x_center ~= 0
            
            
            x_center_round = round(x_center);
            y_center_round = round(y_center);
            y_sub = (y_center_round-dy:y_center_round+dy);
            x_sub = (x_center_round-dx:x_center_round+dx);
            
            img_sub = img(y_sub,x_sub);
            coordinates = {x_sub, y_sub};
            disp(['Fitting Gaussian to peak ', num2str(i)]);
            
            [parFit, ssr] = gauss_fit2D_cov(img_sub,coordinates);
            
            
            % Add fit parameters to the output
            Particles2(i, 8:13) = [parFit, ssr];
            %                 % Test plot
            %                 figure;
            %                 subplot(2,1,1);
            %                 imagesc(img);
            %                 hold on;
            %                 colormap(gray);
            %                 plot(PeaksInCurrFrame(10, iPeak),PeaksInCurrFrame(11, iPeak), 'or');
            %
            % Subtract found particle from image;
            exponent = ((x - parFit(4))/parFit(2)).^2 +...
                ((y - parFit(5))/parFit(2)).^2;
            
            img_del = parFit(1) * exp(-exponent/2);
            
            img =img - img_del;
        end
    end
end
if ~isempty(Particles1) && ~isempty(Particles2)
    Particles1 = Particles1(:,11:12);
    Particles2 = Particles2(:,11:12);
    P_indx = 1;
    used = [];
    %find particles in each image that most likely are the same points
    for i = 1:size(Particles1,1)
        x_cent = Particles1(i,1);
        y_cent = Particles1(i,2);
        
        d_min = 100;
        for j = 1:size(Particles2,1)
            if isempty(find(used == j,1))
                d_cur = sqrt((Particles2(j,1) - x_cent).^2 + (Particles2(j,2) - y_cent).^2);
                if d_cur < d_min
                    d_min = d_cur;
                    indx_cur = j;
                end
            end
        end
        if d_min < max_disp
            used = [used;indx_cur];
            fixedpts(P_indx,:) = Particles1(i,:);
            movingpts(P_indx,:) = Particles2(indx_cur,:);
            P_indx = P_indx + 1;
            
        end
        
    end
    
    RegMatrix{1} = fitgeotrans(movingpts,fixedpts,TranType);
    compareMat = eye(3);
    diffMat = RegMatrix{1}.T - compareMat;
    ind = 1;
    while max(abs(diffMat(:))) > 10*eps
        ind = ind+1;
        movingpts = transformPointsForward(RegMatrix{ind-1},movingpts);
        RegMatrix{ind} = fitgeotrans(movingpts,fixedpts,TranType);
        diffMat = RegMatrix{ind}.T - compareMat;
    end
else
    RegMatrix = [];
end

blah = 10;

