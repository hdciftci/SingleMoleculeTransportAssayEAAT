function Fit_Efflux ()
%clear all;
data=Trajectories();
minfret=0.4;
total_fit=zeros(size(data.fret,1),9);
for trace=1:size(data.fret,1)   
[z,y,resolution]=fetch_fret(data.fret,data.time,trace,minfret);
%set of stable parameters 
S_initial=100;
diff=0.22;
Kd=20;
%monomer=3;
nonopt0=[S_initial diff resolution Kd];
% set of parameters I would like to optimize
unbound_fret=0.63;
Km=30;
r=100;
vmax=10;
% defining initial parameter set to start optimization
opt0=[unbound_fret Km r vmax];
% merge parameter set
tot_opt0=double([nonopt0 opt0]);
lb=double([S_initial, diff-0.02 resolution Kd 0.4 30 100 0.1]); % lower bounds for parameters
ub=double([S_initial, diff+0.02 resolution Kd 0.75 30 100 100]); % upper bounds for parameters

%optimizing parameters
fitfun = @(opt_all,z) efflux_function(z,opt_all);
fit_fret=lsqcurvefit(fitfun,tot_opt0,z,y,lb,ub);
%truncating trajectories end based on fitted unbound fret state
fitted_y=efflux_function(z,fit_fret);
empty_frame=find(abs(fitted_y-fit_fret(5))<0.05,1, 'first');
if size(empty_frame, 1)==0
    empty_frame=0;
end
full_frame=find(abs(fitted_y-mean(fitted_y(1:10))<0.05),1, 'last');
trun_frame=min(length(z), (empty_frame+full_frame)); 

z=z(1:trun_frame);

PEB1a=(1/6.02e23/(4/3*pi*fit_fret(7)^3*10^(-24))*1e6);
total_fit(trace,1:8)=fit_fret;
total_fit(trace,9)=fit_fret(8)/PEB1a;
%plotting experimental and fitted values
plot_fit(z,y,efflux_function(z,fit_fret),z(end),trace);

end

cHeader = {'Substrate(µM)' 'FRET amplitude' 'resolution' 'Kd of PEB1a (µM)' 'FRET initial' 'Km of EAAT1' 'liposome radius (nm)' 'vmax (uM/sec)' 'vmax (# molecules/sec)'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

fid = fopen('rate_constants.csv','w'); 
fprintf(fid,'%s\n',textHeader);
dlmwrite('rate_constants.csv',total_fit,'-append') 
fclose(fid);
close all;
end
