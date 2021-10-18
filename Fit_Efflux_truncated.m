function Fit_Efflux_truncated ()
clear all;
data=Trajectories();
minfret=0.4;
total_fit=zeros(size(data.fret,1),10);
total_fit_truncated=zeros(size(data.fret,1),10);
truncated=0;
for trace=1:size(data.fret,1)   
total_fit(trace,1)=trace;
total_fit_truncated(trace,1)=trace;
[z,y,resolution]=fetch_fret(data.fret,data.time,trace,minfret);
%set of stable parameters 
S_initial=100;
Kd=20;
% set of parameters I would like to optimize
diff=0.22;
unbound_fret=0.63;
Km=30;
r=100;
vmax=10;
% defining initial parameter set to start optimization
nonopt0=[S_initial diff resolution Kd];
opt0=[unbound_fret Km r vmax];
% merge parameter set
tot_opt0=double([nonopt0 opt0]);
lb=double([S_initial diff-0.02 resolution Kd 0.4 30 1 0.1]); % lower bounds for parameters
ub=double([S_initial diff+0.02 resolution Kd 0.75 30 200 100]); % upper bounds for parameters

%optimizing parameters
fitfun = @(opt_all,z) efflux_function(z,opt_all);
fit_fret=lsqcurvefit(fitfun,tot_opt0,z,y,lb,ub);
%truncating trajectories end based on fitted unbound fret state
fitted_y=efflux_function(z,fit_fret);
empty_frame=find(abs(fitted_y-fit_fret(5))<0.05,1, 'first');
full_frame=find(abs(fitted_y-mean(fitted_y(1:10)))<0.05,1, 'last');
if length(z)>(empty_frame+full_frame)
z_truncated=z(1:(empty_frame+full_frame));
y_truncated=y(1:(empty_frame+full_frame));

%re-optimizing parameters with truncated trajectory
fitfun_truncated = @(opt_all,z_truncated) efflux_function(z_truncated,opt_all);
fit_fret_truncated=lsqcurvefit(fitfun_truncated,tot_opt0,z_truncated,y_truncated,lb,ub);

PEB1a_truncated=(1/6.02e23/(4/3*pi*fit_fret_truncated(7)^3*10^(-24))*1e6);
total_fit_truncated(trace,2:9)=fit_fret_truncated;
total_fit_truncated(trace,10)=fit_fret_truncated(8)/PEB1a_truncated;

truncated=1;
end
PEB1a=(1/6.02e23/(4/3*pi*fit_fret(7)^3*10^(-24))*1e6);
total_fit(trace,2:9)=fit_fret;
total_fit(trace,10)=fit_fret(8)/PEB1a;


%plotting experimental and fitted values
if truncated==1
  plot_fit_truncated(z,y,efflux_function(z,fit_fret),z_truncated,efflux_function(z_truncated,fit_fret_truncated), z(end),trace);
else 
  plot_fit(z,y,efflux_function(z,fit_fret),z(end),trace);  
end

truncated=0;
end

cHeader = {'trace #' 'Substrate(µM)' 'FRET amplitude' 'resolution' 'Kd of PEB1a (µM)' 'FRET initial' 'Km of EAAT1' 'liposome radius (nm)' 'vmax (uM/sec)' 'vmax (# molecules/sec)'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

fid = fopen('rate_constants.csv','w'); 
fprintf(fid,'%s\n',textHeader);
dlmwrite('rate_constants.csv',total_fit,'-append') 
fclose(fid);


fid = fopen('rate_constants_truncated.csv','w'); 
fprintf(fid,'%s\n',textHeader);
dlmwrite('rate_constants_truncated.csv',total_fit_truncated,'-append') 
fclose(fid);
%close all;
end
