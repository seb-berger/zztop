function y = zztop(x, wnd)
%ZZTOP Estimate peak pattern probability.
%   Y = ZZTOP(X) returns the peak pattern probability for each channel of
%   the multi-channel signal in X. Each column in X is treated as a
%   channel.
%
%   Y = ZZTOP(X, WND) uses a sliding window of length WND for the analysis.
%   The window is moved one sample per evaluation (maximum overlap).

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

if ~isreal(x)
    error('Input data must be real-valued.');
end

if ~ismatrix(x)
    error('Input data must be a vector or matrix.');
end

if ~all(isfinite(x(:)))
    error('Input data contains infinite values.');
end

% Transpose row vector input data.
if isrow(x)
    x = x.';
end

% By default, the window size equals the signal length.
if nargin == 1
    wnd = size(x, 1);
end

% Compensate for first and last sample.
wnd = wnd - 2;

if wnd < 1
    error('Window size must be greater than 2.');
elseif wnd > size(x, 1)
    error('Window size exceeds signal length.');
end

% Detect peaks;
peak = abs(diff(diff(x) >= 0));

% Prepare sliding window signals
head = peak(wnd + 1:end, :);
tail = peak(1:end - wnd, :);

if ~isempty(head)
    % Calculate peak probability using a sliding window.
    temp = [sum(peak(1:wnd, :)); head - tail];
    y = cumsum(temp) / wnd;
else
    % Calculate peak probability for the whole signal.
    y = sum(peak(1:end, :)) / wnd;
end
