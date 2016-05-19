function [ outDir ] = getGldasForcing(qNames, qLat, qLon, qStart, qEnd, outDir)
% GETGLDASFORCING This script will download GLDAS rainfall data from
%       Nasa's servers.
% Created by Peter J. Shellito 5/19/16
% 
% qNames: cell array of site names
% qLat: latitude corresponding to each site
% qLon: longitude corresponding to each site
% qStart: a vector [yyyy, mm, dd] specifying the day to start downloading.
%       Start date must be [2000, 2, 24] or later.
% qEnd: a vector [yyyy, mm, dd] specifying the day to stop downloading. End
%       date must be 20 days before today or earlier. Neither qStart nor 
%       qEnd support a starting hour and minute.
% outDir: Directory to place the output files. If nothing is provided,
%       default is to create a directory in the present directory titled
%       './outFiles/'
% 
% For your own reference, or if something goes wrong, the directory where
% you can find all the GLDAS data documentation:
% http://disc.sci.gsfc.nasa.gov/uui/datasets/GLDAS_NOAH025SUBP_3H_V001/summary
% And the location of a sample file:
% fileDir = 'ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/GLDAS_SUBP/GLDAS_NOAH025SUBP_3H/2016/121/GLDAS_NOAH025SUBP_3H.A2016121.0000.001.2016137055857.grb';
% 
% Initial testing revealed that processing one month took xx minutes using
% a CU internet connection.
% 
% ========================================================================
% NB: This script requires the user to have downloaded and installed the
% nctoolbox: http://nctoolbox.github.io/nctoolbox/
% Extract nctoolbox-1.1.3, and then in the matlab command line, navigate 
% to the nctoolbox: cd('Path/to/toolbox/nctoolbox-1.1.3')
% Then type 'setup_nctoolbox' and you should be able to use the toolbox.
% It may be necessary to run setup_nctoolbox every time you start a new 
% matlab session. To automate this, add the '/Path/to/toolbox/' to your 
% matlab path. Then edit (or create) a startup.m file that is also in your
% matlab path. Add the following to startup.m: setup_nctoolbox;
% ========================================================================
% 
% Copyright (c) 2016, Peter J. Shellito, University of Colorado Boulder
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% -------------------------------------------------------------------------
% If no output directory was provided
if nargin<6
    outDir = './outFiles';
end

% -------------------------------------------------------------------------
% Make sure data are of the proper type and dimensions
if ~iscellstr(qNames)
    error('Input names must be a cell array of strings.')
end
dnStart = floor(datenum(qStart));
dnEnd = floor(datenum(qEnd));
if dnStart < datenum(2000, 2, 24)
    dnStart = datenum(2000, 2, 24);
end
if dnEnd > floor(datenum(now))-20
    dnEnd = floor(datenum(now))-20;
end
if dnEnd < dnStart
    error('Start date cannot be after end date')
end
% -------------------------------------------------------------------------
% Set up some variables
% Number of sites requested
nSites = length(qNames);
% Initialize vectors for rounded lat and lon and lat/lon idcs
nearestLat = nan(nSites,1);
nearestLon = nan(nSites,1);
latIdcs = nan(nSites,1);
lonIdcs = nan(nSites,1);
% The years and days of years of the query
qDatenums = dnStart:dnEnd;
% Initialize a day of year vector
qDoy = nan(1,length(qDatenums));
% The years, months, and days of the query
[qYears, qMonths, qDays] = datevec(qDatenums);
% The unique years in the query
qUniqueYears = unique(qYears);
% For each unique year
for yy = 1:length(qUniqueYears)
    % Find all the indices where qYears matches this unique year
    sameYrIdcs = find(qYears == qUniqueYears(yy));
    % Subtract the datenum of day zero of the year from the datenums in the
    % year to get day of year
    qDoy(sameYrIdcs) = datenum(qDatenums(sameYrIdcs)) - datenum(qUniqueYears(yy),0,0);
end % Loop through each unique year
% Convert qYears to string
qYearStr = num2str(qYears');
% Convert qMonths to string
qMonthStr = num2str(qMonths', '%02d');
% Convert qDays to string
qDayStr = num2str(qDays', '%02d');
% Convert qDoy to string
qDoyStr = num2str(qDoy', '%03d');
% Create hour strings
qHourStr = num2str((0:300:2300)','%04d');
% The directory where nldas forcings are held
ftpBaseDir = '/data/s4pa/GLDAS_SUBP/GLDAS_NOAH025SUBP_3H/';
% The local directory where gldas forcings will be placed
localBaseDir = [pwd '/data'];
% If there is already a directory named data in the working direcotry, do not continue because
% it will be deleted at the end of this script and I don't want to delete
% anything that this script itself did not create.
if exist(localBaseDir, 'dir') == 7
    error('%s\n%s\n%s','This function requires use of the local directory named ''./data.'',', ...
        'and you already have a directory with that name. We are stopping here',...
        'because the end of this function will DELETE ./data and all files within it.')
end
% The beginning of the file name of forcings
ftpBaseFn = 'GLDAS_NOAH025SUBP_3H.A';
% The end of the file name of forcings
ftpEndFn = '*.grb';
% Strings of variables to read
latStr = 'lat'; % Center of 1/8 degree pixel
lonStr = 'lon'; % Center of 1/8 degree pixel
rainfallStr = 'VAR221-221-1-132_surface'; % Rainfall rate, past 3-hr average [kg/m^2/s]
snowfallStr = 'VAR221-221-1-131_surface'; % Snowfall rate, past 3-hr average [kg/m^2/s]
% Put all the strings into one cell
allStrings = {rainfallStr; snowfallStr};
% The number of variables I'm getting
nVars = length(allStrings);
% Names and units for the variables (for headers)
namesStrAll = 'Rainfall_Rate Snowfall_Rate';
unitsStrAll = '[kg/m^2/s]    [kg/m^2/s]';
% Variable formats
varFmt = '%13.4e %13.4e\n';
% Date formats
dateFmt = '%04d   %02d    %02d  %02d   %02d     ';
% -------------------------------------------------------------------------
% Create a directory to hold output data
if exist(outDir, 'dir') ~= 7
    disp(['Making an output directory here: ' outDir])
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Set up ftp connection and download the first file to get its lat/lon data
% The Nasa host that holds nldas forcings
nasaHost = 'hydro1.sci.gsfc.nasa.gov';
% Create an ftp object (open the connection)
ftpObj = ftp(nasaHost);
% The file name of the first date requested
qFileName = [ftpBaseDir qYearStr(1,:) '/' qDoyStr(1,:) '/' ftpBaseFn qYearStr(1,:) qDoyStr(1,:) '.' qHourStr(1,:) ftpEndFn];
% Get that first file
disp(['Getting ' qFileName '...'])
localFileName = mget(ftpObj,qFileName);
% Create ncdataset object
geo = ncdataset(localFileName{1});
% Extract lat and lon from the grib file
lat = geo.data(latStr);
lon = geo.data(lonStr);
% Size of these vectors
nLat = length(lat);
nLon = length(lon);

% -------------------------------------------------------------------------
% Set up output files. Write headers with lat/lon data in them. Get lat/lon
% idcs.
% Initialize a vector of file IDs
fid = nan(1,nSites);
% Loop through each site
for ss = 1:nSites
    % Get lat/lon idcs
    [latDiff latIdcs(ss)] = min(abs(lat-qLat(ss)));
    [lonDiff lonIdcs(ss)] = min(abs(lon-qLon(ss)));
    % Get the rounded lat and lon. This is the center point of the NLDAS
    % cell.
    nearestLat = lat(latIdcs(ss));
    nearestLon = lon(lonIdcs(ss));
    % If the difference between the queried and actual lat or lon is too
    % big, display a warning
    if latDiff>0.125 || lonDiff>0.125
        warning('%s site does not have a nearby nldas cell. \nQueried lat/lon: %f %f\nClosest NLDAS lat/lon: %f %f',...
            qNames{ss}, qLat(ss), qLon(ss), nearestLat, nearestLon)
    end
    % The name for this output file
    outFile = [outDir '/' qNames{ss} '.txt'];
    % Open the output file for the first time
    fid(ss) = fopen(outFile,'w');
    % Print header lines to the file
    fprintf(fid(ss),['%% Site: ' qNames{ss} '\n']);
    fprintf(fid(ss),['%% Site lat/lon: ' num2str([qLat(ss) qLon(ss)]) '\n']);
    fprintf(fid(ss),['%% Closest GLDAS pixel (1/8 degree) center: ' num2str([nearestLat nearestLon]) '\n']);
    fprintf(fid(ss),['%% File created on ' date '.\n']);
    fprintf(fid(ss),['%% Date                       ' namesStrAll '\n']);
    fprintf(fid(ss),['%% year month day hour minute ' unitsStrAll '\n']);
end % Loop through each site

% -------------------------------------------------------------------------
% Loop through each day in the record
for dd = 1:length(qDatenums)
    % Location of this day's data on the local machine
    localDir = [pwd ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:)];
    % Loop through 3 hours at a time
    for hh = 1:(24/3)
        % Create strings and such to define where the files are located on the host
        qFileName = [ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:) '/' ftpBaseFn qYearStr(dd,:) qDoyStr(dd,:) '.' qHourStr(hh,:) ftpEndFn];
        % Get the file from Nasa's server
        disp(['Getting ' qFileName '...'])
        localFileName = mget(ftpObj,qFileName);
        % Create ncgeodataset object
        geo = ncdataset(localFileName{:});
        % Initialize a vector to hold the timestep's data (all variables) for entire domain
        domainData = nan(nLat, nLon, nVars);
        % Loop through each variable
        for vv = 1:nVars
            % Get the data for each variable
            domainData(:,:,vv) = squeeze(geo.data(allStrings{vv}));
        end
        % Loop through each site 
        for ss = 1:nSites
            % Extract the point data at each site from the domain data
            siteMetData = squeeze(domainData(latIdcs(ss), lonIdcs(ss), :));
            % Print the data. Hour is hh-1 because hours are listed from
            % 00:00 to 23:00.
            fprintf(fid(ss), [dateFmt varFmt], [qYears(dd) qMonths(dd) qDays(dd) (hh-1) 0 siteMetData']);
        end % loop through each site
    end % Loop through each hour of the day
    % Delete this day's directory and all data within
    rmdir(localDir, 's')
end % Loop through each day in the record
% -------------------------------------------------------------------------
% Clean up and close files
% Delete all data stored in this session
disp(['Cleaning up...'])
rmdir(localBaseDir, 's')
% Loop through each site and close the files
for ss = 1:nSites
    outFile = qNames{ss};
    fclose(fid(ss));
end
% Close the ftp connection
close(ftpObj)

end
