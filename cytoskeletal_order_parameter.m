% --------------------------------------------------------------------
% Released under MIT License
% Author: Hamsini Suresh
% --------------------------------------------------------------------

% Post-process .mat files from MCMC and generate distributions of
% stress-fibre(SF) orientations in the homeostatic ensemble for 
% given stripe width(s)
%

clc; clear

%% Assembly

% -------- Inputs -------- %  
sNote  = '';
beta   = [0.14; 0.238];    % homeostatic temperature for each stripe width of interest
setnum = [4; 5];         % no. of .mat files for each stripe width of interest


% -------- Construct SF orientation distributions -------- %
cskphi = -180:0.25:180;  % discretise SF orientations
cskphi = cskphi(:);
rndphi = zeros(length(cskphi),length(setnum));  

for ru = 1:length(setnum)
    for ty = 1:setnum(ru)
       sfile = sprintf('cell_stripes_%s_beta%.3f.mat', sNote, beta(ru)); 
       load(sfile);
       
       % assemble overall SF orientation distribution for each stripe width
       rndphi(:,ru) = rndphi(:,ru) +  rndopcsk;   
    end
end


% rearrange each distribution so that SF orientation angles lie in interval
% [-90, 90] degrees
lower1 = 1;
[~, lower2] = min(abs(cskphi + 90));
lower2 = lower2 - 1;

[~, upper1] = min(abs(cskphi - 90));
upper1 = upper1 + 1;
upper2 = length(cskphi);

cskphi_pm90 = -90:0.25:90;
cskphi_pm90 = cskphi_pm90(:);
zero_index = (length(cskphi_pm90)+1)/2;
rndphi_pm90 = zeros(length(cskphi_pm90), length(setnum));
opcsksim = zeros(length(setnum),1);
for yy = 1:length(setnum)
    rndphi_pm90(:, yy) = rndphi(lower2+1:upper1-1, yy);  % copy from -90 to 90
    rndphi_pm90(zero_index:length(cskphi_pm90)-1, yy) = rndphi_pm90(zero_index:length(cskphi_pm90)-1, yy) + rndphi(lower1:lower2, yy);  % from -180 to -90.25
    rndphi_pm90(2:zero_index, yy) = rndphi_pm90(2:zero_index, yy) + rndphi(upper1:upper2, yy);  % from 90.25 to 180
    
    opcsksim(yy) = sqrt( (sum((cosd(2*cskphi_pm90)).*rndphi_pm90(:,yy))/(sum(rndphi_pm90(:,yy))))^2 + (sum((sind(2*cskphi_pm90)).*rndphi_pm90(:,yy))/(sum(rndphi_pm90(:,yy))))^2 );
end

% values at +/- 90
rndphi_pm90(1,:) = rndphi_pm90(2,:); 
rndphi_pm90(end,:) = rndphi_pm90(end-1,:); 

disp('Cytoskeletal order parameters for stripes of interest:')
disp(opcsksim)


%% Plotting

figure(1); clf
set(gcf,'color','w');

phicsk1 = rndphi_pm90(:,1)/trapz(cskphi_pm90, rndphi_pm90(:,1));
phicsk2 = rndphi_pm90(:,2)/trapz(cskphi_pm90, rndphi_pm90(:,2));

plot(cskphi_pm90,phicsk1,'LineWidth',2)
hold on
plot(cskphi_pm90,phicsk2,'LineWidth',2)
xlim([-90 90])
set(gca,'FontSize',24);