function [res, idx] = percentiles(x, y, steps, threshold)
%PERCENTILES Calculate percentile bands.
%   [RES, IDX] = PERCENTILES(X, Y, STEPS, THRESHOLD) takes a multi-class
%   population Y, where class membership is annotated by X, and calculates
%   its percentiles at a granularity of STEPS. Classes with less than
%   THRESHOLD samples are rejected. The function returns a matrix of
%   percentiles RES and a vector IDX of class membership annotations.

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

% Sort x and y with regard to the values in x.
[x, i] = sort(x);
y = y(i);

% Calculate percentile bounds.
bounds = linspace(0, 100, steps);

% Find blocks of equal values in x.
stops  = [find(diff(x) > 0), numel(x)];
starts = [1, stops(1:end-1) + 1];

% Reject blocks with less than 'threshold' number of samples.
reject = (stops - starts + 1) < threshold;
starts(reject) = [];
stops(reject) = [];

% Initialise results vector with NaNs
res = nan(numel(starts), steps);

% For each block of identical x, calculate percentiles.
for n = 1:numel(starts)
    res(n, :) = prctile(y(starts(n):stops(n)), bounds);
end

% Create index for each bin.
idx = x(starts).';
