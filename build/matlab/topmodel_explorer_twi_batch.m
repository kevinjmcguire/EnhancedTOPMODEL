function topmodel_explorer_twi_batch(varargin)
%TOPMODEL_EXPLORER_TWI_BATCH  30 m TWI for the 11 TOPMODEL Explorer watersheds.
%
%   topmodel_explorer_twi_batch                 % all 11 gauges
%   topmodel_explorer_twi_batch('02017500')     % one or more gauges by ID
%   topmodel_explorer_twi_batch(..., 'res', 60) % other cell size (m)
%
% Runs twi_watershed_usgs.m (unchanged; it must be on the MATLAB path) for
% each CAMELS gauge with the outlet snapped by drainage area, then writes a
% slim, compressed file per watershed:
%
%   twi_slim_<gauge>.mat   (-v7, zlib-compressed)
%     twiWS    single  TWI = ln(a/tan b), NaN outside the watershed (0.01 precision)
%     zWS      single  filled DEM, NaN outside (0.1 m)    - for a hillshade backdrop
%     distCh   single  D8 channel flow distance to the outlet (m), NaN outside
%     strahler uint8   Strahler order on stream cells, 0 elsewhere
%     xc, yc   double  cell-centre UTM coordinates of the cropped grid (m)
%     dx, epsg, zone, hemi, bx, by (boundary), outlet, stats, gauge, name
%
% All grids are cropped to the watershed bounding box plus one cell.
% At the end everything is zipped into topmodel_twi_30m.zip with a summary
% CSV that compares the delineated area with the CAMELS (GAGES-II) area.
%
% Gauges that fail are reported and skipped; re-running skips gauges whose
% slim file already exists (delete it to redo one).
%
% Memory: the MD-infinity step stores 8 single-precision weights per cell
% of the DEM box; the largest basin here (Fish River, 2,250 km^2) needs a
% few GB at 30 m. If a gauge runs out of memory, run it alone at 60 m:
%   topmodel_explorer_twi_batch('01013500', 'res', 60)

% ------------------------------------------------------------------ gauges
G = { ...  % id,       lat,       lon,      CAMELS area km^2, name
    '02017500', 37.50624,  -80.10671,   276.15, 'Johns Creek at New Castle, VA'
    '02028500', 37.86959,  -78.82335,   245.45, 'Rockfish River near Greenfield, VA'
    '02038850', 37.41542,  -78.63584,    22.00, 'Holiday Creek near Andersonville, VA'
    '01669000', 37.87708,  -76.90052,    71.84, 'Piscataway Creek near Tappahannock, VA'
    '01013500', 47.23739,  -68.58264,  2252.70, 'Fish River near Fort Kent, ME'
    '12040500', 47.53786, -124.31575,  1153.44, 'Queets River near Clearwater, WA'
    '11143000', 36.24580, -121.77329,   120.61, 'Big Sur River near Big Sur, CA'
    '11266500', 37.71687, -119.66628,   833.08, 'Merced River at Pohono Bridge near Yosemite, CA'
    '09352900', 37.47750, -107.54367,   188.15, 'Vallecito Creek near Bayfield, CO'
    '09512280', 33.88726, -111.95404,   188.71, 'Cave Creek below Cottonwood Creek near Cave Creek, AZ'
    '08202700', 29.37079,  -99.28782,   434.75, 'Seco Creek at Rowe Ranch near D''Hanis, TX'};

% ------------------------------------------------------------------ options
res = 30;
ids = {};
k = 1;
while k <= numel(varargin)
    a = varargin{k};
    if (ischar(a) || isstring(a)) && strcmpi(a, 'res')
        res = varargin{k+1};  k = k + 2;
    else
        ids = [ids, cellstr(a)]; %#ok<AGROW>
        k = k + 1;
    end
end
if isempty(ids), ids = G(:, 1)'; end

script = which('twi_watershed_usgs.m');
if isempty(script)
    error('twi_watershed_usgs.m is not on the MATLAB path.');
end
root   = pwd;
work   = fullfile(root, sprintf('twi_work_%gm', res));      % full outputs per gauge
slimD  = fullfile(root, sprintf('twi_slim_%gm', res));      % files to upload
if ~exist(work, 'dir'), mkdir(work); end
if ~exist(slimD, 'dir'), mkdir(slimD); end

% The TWI script reads the USGS key from the base workspace or the
% environment; inside this function only the environment is visible.
try
    key = evalin('base', 'USGS_API_KEY');
    if ~isempty(key), setenv('USGS_API_KEY', char(key)); end
catch
end

summary = cell(0, 9);
for i = 1:numel(ids)
    row = find(strcmp(G(:, 1), ids{i}));
    if isempty(row)
        warning('Unknown gauge %s - skipped.', ids{i});  continue
    end
    g = G(row, :);
    fslim = fullfile(slimD, sprintf('twi_slim_%s.mat', g{1}));
    if exist(fslim, 'file')
        fprintf('\n=== %s %s: slim file exists, skipping ===\n', g{1}, g{5});
        S = load(fslim, 'stats');
        summary(end+1, :) = summary_row(g, S.stats, res, 'existing', dir(fslim)); %#ok<AGROW>
        continue
    end
    fprintf('\n=== %s %s (CAMELS %.1f km^2) at %g m ===\n', g{1}, g{5}, g{4}, res);
    t0 = tic;
    try
        S = run_one(script, g, res, fullfile(work, g{1}));
        save(fslim, '-struct', 'S', '-v7');
        d = dir(fslim);
        fprintf('Saved %s (%.1f MB) in %.1f min\n', d.name, d.bytes/1e6, toc(t0)/60);
        summary(end+1, :) = summary_row(g, S.stats, res, 'ok', d); %#ok<AGROW>
    catch ME
        fprintf(2, 'FAILED %s: %s\n', g{1}, ME.message);
        summary(end+1, :) = {g{1}, g{5}, g{4}, NaN, NaN, res, NaN, NaN, ['failed: ' ME.message]}; %#ok<AGROW>
    end
    close all force
end

% ------------------------------------------------------------------ summary + zip
T = cell2table(summary, 'VariableNames', {'gauge', 'name', 'camels_area_km2', ...
    'delineated_area_km2', 'area_diff_pct', 'cell_m', 'twi_mean', 'file_MB', 'status'});
csvF = fullfile(slimD, 'twi_summary.csv');
writetable(T, csvF);
disp(T(:, [1 3 4 5 7 8 9]));
bad = ~strcmp(T.status, 'ok') & ~strcmp(T.status, 'existing');
big = abs(T.area_diff_pct) > 10;
if any(big & ~bad)
    fprintf(2, ['Check these outlets - delineated area differs from CAMELS by > 10%%: %s\n' ...
                'Re-run one with cfg.interactive = true in twi_watershed_usgs.m to pick it by hand.\n'], ...
            strjoin(T.gauge(big & ~bad), ', '));
end
dd = dir(fullfile(slimD, 'twi_slim_*.mat'));
files = [fullfile(slimD, {dd.name}), {csvF}];
zipF = fullfile(root, sprintf('topmodel_twi_%gm.zip', res));
zip(zipF, files);
d = dir(zipF);
fprintf('\nUpload this file: %s (%.1f MB)\n', zipF, d.bytes/1e6);
end

% =========================================================================
function S = run_one(script, g, res, outDir)
% Run the unchanged TWI script for one gauge in this function's workspace.
cfgOverride = struct( ...
    'lat', g{2}, 'lon', g{3}, ...
    'targetArea_km2', g{4}, ...           % snap to the cell matching the gauge area
    'forceDemRes_m', res, ...             % fixed 30 m DEM, no lidar tiers
    'interactive', false, ...
    'snapRadius_m', max(500, 4*res), ...
    'streamThresh_km2', max(0.1, 0.02*sqrt(g{4})), ...  % coarser network in big basins
    'outDir', outDir); %#ok<NASGU>
% Declare the script's outputs so MATLAB treats these names as variables
[twiWS, Zws, ws, distCh, strahler, isS, xc, yc, dx, epsg, zone, hemi, bx, by, stats] = deal([]);
run(script);                              % fills twiWS, ws, xc, yc, dx, ... here

% crop to the watershed bounding box + 1 cell
[r, c] = find(ws);
r1 = max(min(r) - 1, 1);  r2 = min(max(r) + 1, size(ws, 1));
c1 = max(min(c) - 1, 1);  c2 = min(max(c) + 1, size(ws, 2));
R = r1:r2;  C = c1:c2;
w = ws(R, C);

S = struct();
S.gauge = g{1};  S.name = g{5};
S.twiWS = single(round(twiWS(R, C) * 100) / 100);
z = Zws(R, C);  z(~w) = NaN;
S.zWS = single(round(z * 10) / 10);
dc = distCh(R, C);  dc(~w) = NaN;
S.distCh = single(round(dc));
so = strahler(R, C);  so(~(w & isS(R, C))) = 0;
S.strahler = uint8(so);
S.xc = xc(C);  S.yc = yc(R);
S.dx = dx;  S.epsg = epsg;  S.zone = zone;  S.hemi = hemi;
S.bx = single(bx);  S.by = single(by);
S.outlet = struct('lat', stats.lat_outlet, 'lon', stats.lon_outlet, ...
                  'snap_distance_m', stats.snap_distance_m);
S.stats = stats;
S.camels_area_km2 = g{4};
end

function row = summary_row(g, stats, res, status, d)
row = {g{1}, g{5}, g{4}, stats.area_km2, 100*(stats.area_km2/g{4} - 1), res, ...
       stats.twi_mean, d.bytes/1e6, status};
end
