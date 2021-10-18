function data = Trajectories (varargin)
% Default parameter values
params.truncateLen = 3000;  %frames to calculate over
params.min_fret = 0.125;  % minimum fret value, below which we assume there is no FRET.

%% Process input arguments
narginchk(0,3);
nargoutchk(0,1);
[varargout{1:nargout}] = deal([]);
[cax,args] = axescheck(varargin{:});

switch numel(args)
    case 0
        files = getFiles();
    case 1
        files = args{1};
    case 2
        [files,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end

if ~iscell(files), files={files}; end
nFiles = numel(files);
if nFiles==0,  return;  end

for i=1:numel(files)% this file should only have traces that respond 
    % Load FRET data and truncate to target length
    data = loadTraces( files{1} );
    data.nFrames = params.truncateLen;
end
end
