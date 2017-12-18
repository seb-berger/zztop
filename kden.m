function [y, y_i] = kden(x, sigma, x_i)
%KDEN Calculate Gaussian kernel density estimate.
%   [Y, Y_I] = KDEN(X, SIGMA) calculates a probability density estimate for
%   the data in the vector X. The function uses a Gaussian kernel with
%   bandwidth SIGMA and evaluates the density at 100 points.
%   It returns a vector Y of densities and a vector Y_I of the respective
%   positions.
%
%   [Y, Y_I] = KDEN(X, SIGMA, X_I) evaluates the densities at the positions
%   provided in the vector X_I.

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

narginchk(2, 3);

if nargin == 2
    y_i = linspace(min(x) - 3*sigma, max(x) + 3*sigma, 100);
else
    y_i = x_i;
end

kernel = @(x) exp(-0.5 * x.^2) / sqrt(2 * pi);
y = zeros(size(y_i));

for n = 1:numel(y_i)
    y(n) = mean(kernel((x - y_i(n)) / sigma)) / sigma;
end
