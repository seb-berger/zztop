function edf = import_edf(filename, header_only)
%IMPORT_EDF Import physiological data from an EDF file.
%   EDF = IMPORT_EDF(FILENAME) returns a struct containing the data in
%   the EDF file with path FILENAME.
%
%   EDF = IMPORT_EDF(FILENAME, HEADER_ONLY) does neither read nor return
%   any signal samples if HEADER_ONLY is set to TRUE.

%   Copyright (c) 2017, Sebastian Berger.
%
%   Klinikum rechts der Isar der
%   Technischen Universitaet Muenchen
%   Munich, Germany
%
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met:
%       * Redistributions of source code must retain the above copyright
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright
%         notice, this list of conditions and the following disclaimer in
%         the documentation and/or other materials provided with the
%         distribution.
%       * Neither the names of the copyright holders nor the names of its
%         contributors may be used to endorse or promote products derived
%         from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
%   THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
%   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
%   THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.

narginchk(1, 2);

if nargin < 2
    header_only = false;
end

% Header format:<name>,    <len>, <is_num>
glob_hdr_fmt = {'ver',         8     false
                'pat',        80,    false
                'rec',        80,    false
                'date',        8,    false
                'time',        8,    false
                'len',         8,     true
                'reserved',   44,    false
                'num_rcd',     8,     true
                'dur_rcd',     8,     true
                'num_ch',      4,     true};

chan_hdr_fmt = {'label',      16,    false
                'transducer', 80,    false
                'unit',        8,    false
                'phys_min',    8,     true
                'phys_max',    8,     true
                'dig_min',     8,     true
                'dig_max',     8,     true
                'filter',     80,    false
                'num_smp',     8,     true
                'reserved',   32,    false};

glob_field_lens = [glob_hdr_fmt{:, 2}];
chan_field_lens = [chan_hdr_fmt{:, 2}];

glob_hdr_len = sum(glob_field_lens);
chan_hdr_len = sum(chan_field_lens);

% Open EDF file for reading.
[fid, msg] = fopen(filename, 'r');
if fid == -1
    error('Could not open file "%s": %s.', filename, msg);
end

% Make sure to close the file at function return.
fclose_callback = onCleanup(@() fclose(fid));

% Read global file header as one chunk of bytes.
[buf, n] = fread(fid, glob_hdr_len, '*uint8');
if n ~= glob_hdr_len
    error('Unexpected EOF.');
end

% Assert that header contains ASCII characters only.
if any(buf > 126)
    error('Corrupted file header.');
end

% Convert raw bytes to characters.
buf = char(buf);

% Split chunk into fields.
buf = mat2cell(buf, glob_field_lens);

% Trim any trailing whitespace of each field.
buf = cellfun(@(x) strtrim(x.'), buf, 'UniformOutput', false);

% Convert ASCII numbers to double values.
is_numeric = [glob_hdr_fmt{:, 3}];
buf(is_numeric) = num2cell(str2double(buf(is_numeric)));

% Interleave field names and data to create a struct.
buf = [glob_hdr_fmt(:, 1), buf].';
edf = struct(buf{:});

% Check header for consistency.
block_len = edf.num_ch * chan_hdr_len;

if ~strcmp(edf.ver, '0')                 ...
|| (edf.len ~= glob_hdr_len + block_len) ...
|| edf.dur_rcd <= 0                      ...
|| edf.num_ch < 0                        ...
|| edf.num_ch ~= floor(edf.num_ch)
    error('Corrupted file header.');
end

% Reject file formats other than plain EDF and EDF+C.
if ~isempty(edf.reserved)
    if numel(edf.reserved) < 5 || ~isequal(edf.reserved(1:5), 'EDF+C')
        error('File format not supported.');
    end
end

% Read channel headers as one chunk of bytes.
[buf, n] = fread(fid, block_len, '*uint8');
if n ~= block_len
    error('Unexpected EOF.');
end

% Assert that header contains ASCII characters only.
if any(buf < 32) || any(buf > 126)
    error('Corrupted file header.');
end

% Convert raw bytes to ASCII string.
buf = char(buf);

% Split chunk, one array per type of field.
buf = mat2cell(buf, chan_field_lens * edf.num_ch);

% Reshape arrays, one row per channel.
buf = cellfun(@(x) reshape(x, [], edf.num_ch).', buf, ...
    'UniformOutput', false);

% Convert multi-row strings to cell arrays of strings.
% Also trim any trailing whitespace
buf = cellfun(@cellstr, buf, 'UniformOutput', false);

% Convert ASCII numbers to double values
is_numeric = [chan_hdr_fmt{:, 3}];
fcn = @(x) num2cell(str2double(x));
buf(is_numeric) = cellfun(fcn, buf(is_numeric), 'UniformOutput', false);

% Interleave field names and data to create structs.
buf = [chan_hdr_fmt(:, 1), buf].';
edf.chan = struct(buf{:});

% Check channel headers for consistency.
num_smp = [edf.chan.num_smp];

if any(num_smp == 0) || any(num_smp ~= floor(num_smp))
    error('Corrupted file header.');
end

% Get number of samples available
stat = dir(filename);
pos  = ftell(fid);
num_bytes = stat.bytes - pos;

% Make sure file header and body are consistent.
block_len = sum([edf.chan.num_smp]);
num_bytes_expected = 2 * edf.num_rcd * block_len;

if num_bytes > num_bytes_expected
    error('Excess samples.');
elseif num_bytes < num_bytes_expected
    error('Missing samples.');
end

% Leave early if no signal data requested.
if header_only
    return
end

% The remainder of the file are the signal samples:
% read these little-endian 16 bit integers.
samples = fread(fid, Inf, '*int16', 0, 'l');

% Check if reading was successfull.
if 2 * numel(samples) ~= num_bytes
    error('Error reading from file.');
end

% Split raw chunk into per-channel data arrays.
samples = reshape(samples, block_len, edf.num_rcd);
samples = mat2cell(samples, [edf.chan.num_smp]);

% Calculate slopes and offsets.
slope  = ([edf.chan.phys_max] - [edf.chan.phys_min]) ./ ...
         ([edf.chan.dig_max]  - [edf.chan.dig_min]);

offset = [edf.chan.phys_max] - slope .* [edf.chan.dig_max];

% Convert 'digital' to 'physiological' data.
fcn     = @(x, m, t) m * double(x(:)) + t;
slope   = num2cell(slope.');
offset  = num2cell(offset.');
samples = cellfun(fcn, samples, slope, offset, 'UniformOutput', false);

% Add sample data to channel struct array.
edf.chan = cell2struct([struct2cell(edf.chan); samples.'], ...
                       [fieldnames(edf.chan); {'samples'}], 1);
