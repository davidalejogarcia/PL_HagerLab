function plotVariables(whatX,whatY, allParticles, filterParticles)

if nargin == 3
    Filter = 0;
elseif nargin == 4
    Filter = 1;
end

% remove the rows filled with zeros
Ix = find(allParticles(:,1) ~=0);
allParticles = allParticles(Ix, :);


if Filter ==1
    Ix = find(filterParticles(:,1) ~=0);
    filterParticles = filterParticles(Ix, :);
 
end
disp(whatY);
    % Case Intensity vs Sigma;
if (whatX == 1 && whatY == 2) || (whatX == 2 && whatY == 1);
    figure;
    plot(allParticles(:,7),allParticles(:,8), '.k');
    
    if Filter == 1
        hold on;
        plot(filterParticles(:,7),filterParticles(:,8), '.r');
    end
end 
    

