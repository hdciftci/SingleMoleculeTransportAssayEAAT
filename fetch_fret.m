function [z,y,resolution]= fetch_fret (data_fret,data_time, trace, min_fret)
%fetching fret values bigger than params.min_fret (you can change this value in the beginning of Trajectories function) 
nonzero= data_fret(trace,:) >= min_fret; 
 %forming time(sec) FRET matrix named fret_trace 
fret_trace=horzcat(transpose(data_time(:,nonzero)/1000),transpose(data_fret(trace,nonzero)));
fret_trace=fret_trace(1:end-1,:);
%calculating the temporal resolution of the trace
resolution=fret_trace(2,1)-fret_trace(1,1);
%x is the time component of the trace
z=fret_trace(:,1);
z=double(z);
size(z)
%y is the fret component of the trace
y=fret_trace(:,2);
y=double(y);
size(y)
end