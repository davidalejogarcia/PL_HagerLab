function Particles = filterPeaks(allParticles, filterPkPar)

I_Threshold = filterPkPar(1);
SigmaLow = filterPkPar(2);
SigmaHigh = filterPkPar(3);

Index = find (allParticles(:,8)>SigmaLow & allParticles(:,8)<SigmaHigh & allParticles(:,7)>I_Threshold);
Particles = allParticles(Index,:);

% Put a line of 0's for frames without particles

  for i = 1: max(Particles(:,6))
        Index = find(Particles(:,6) == i);
        
        if isempty(Index);
            Particles = [Particles; 0 0 0 0 0 i 0 0 0 0 0 0 0];
        end
   end


Particles = sortrows(Particles, [6 13]);