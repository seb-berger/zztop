function y = symbolise(x, dim, dly, ties)
%SYMBOLISE Encode Ordinal Patterns.
%   Y = SYMBOLISE(X, DIM) encodes the samples in vector X into a series of
%   ordinal patterns of order DIM using time delay 1. Patterns are
%   represented by integers [1, 2, ..., factorial(DIM)]. If X is a matrix,
%   the function operates on each column individually.
%
%   Y = SYMBOLISE(X, DIM, DLY) uses the time delay DLY for encoding.
%
%   Y = SYMBOLISE(X, DIM, DLY, TIES) replaces ties by the value -1 if
%   TIES is set to TRUE. Otherwise (the default), ties are resolved by
%   sample positions: if X(i) == X(j), indices i and j are considered.
%
%   The encoding scheme used (Lehmer code) is thoroughly explained in:
%   Keller, K.; Sinn, M.; Emonds, J. Time series from the ordinal viewpoint.
%   Stochastics Dyn. 2007, 7, 247--272.

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

narginchk(2, 4);

if nargin < 4
    ties = false;
end

if nargin < 3
    dly = 1;
end

if ~isfloat(x) || ~isreal(x)
    error('X must contain real floating point values.');
end

if floor(dim) ~= dim || dim < 2
    error('Dim must be an integer > 1.');
end

if floor(dly) ~= dly || dly < 1
    error('Delay must be an integer > 0.');
end

if ~isscalar(ties) || ~islogical(ties)
    error('Ties must be a scalar logical value.');
end

if isrow(x)
    x = x.';
end

switch dim
    case 3
       y = encode_dim_3(x, dly, ties);
    case 4
       y = encode_dim_4(x, dly, ties);
    case 5
       y = encode_dim_5(x, dly, ties);
    otherwise
       error('Support for patterns of order %d is not implemented.', dim);
end
end


function y = encode_dim_3(x, dly, ties)
    % The Lehmer code is used for permutation encoding.
    x1 = x(1 + 0 * dly : end - 2 * dly, :);
    x2 = x(1 + 1 * dly : end - 1 * dly, :);
    x3 = x(1 + 2 * dly : end - 0 * dly, :);

    y  = 2 * ((x1 > x2) + (x1 > x3)) + ...
         1 * (            (x2 > x3));

    % Matlabian indexing.
    y = y + 1;

    % Mask ties if requested.
    if ties
        tie = (x1 == x2) | (x1 == x3) ...
                         | (x2 == x3);
        y(tie) = -1;
    end
end


function y = encode_dim_4(x, dly, ties)
    % The Lehmer code is used for permutation encoding.
    x1 = x(1 + 0 * dly : end - 3 * dly, :);
    x2 = x(1 + 1 * dly : end - 2 * dly, :);
    x3 = x(1 + 2 * dly : end - 1 * dly, :);
    x4 = x(1 + 3 * dly : end - 0 * dly, :);

    y  = 6 * ((x1 > x2) + (x1 > x3) + (x1 > x4)) + ...
         2 * (            (x2 > x3) + (x2 > x4)) + ...
         1 * (                        (x3 > x4));

    % Matlabian indexing.
    y = y + 1;

    % Mask ties if requested.
    if ties
        tie = (x1 == x2) | (x1 == x3) | (x1 == x4) ...
                         | (x2 == x3) | (x2 == x4) ...
                                      | (x3 == x4);
        y(tie) = -1;
    end
end


function y = encode_dim_5(x, dly, ties)
    % The Lehmer code is used for permutation encoding.
    x1 = x(1 + 0 * dly : end - 4 * dly, :);
    x2 = x(1 + 1 * dly : end - 3 * dly, :);
    x3 = x(1 + 2 * dly : end - 2 * dly, :);
    x4 = x(1 + 3 * dly : end - 1 * dly, :);
    x5 = x(1 + 4 * dly : end - 0 * dly, :);

    y  = 24 * ((x1 > x2) + (x1 > x3) + (x1 > x4) + (x1 > x5)) + ...
          6 * (            (x2 > x3) + (x2 > x4) + (x2 > x5)) + ...
          2 * (                        (x3 > x4) + (x3 > x5)) + ...
          1 * (                                    (x4 > x5));

    % Matlabian indexing.
    y = y + 1;

    % Mask ties if requested.
    if ties
        tie = (x1 == x2) | (x1 == x3) | (x1 == x4) | (x1 == x5) ...
                         | (x2 == x3) | (x2 == x4) | (x2 == x5) ...
                                      | (x3 == x4) | (x3 == x5) ...
                                                   | (x4 == x5);
        y(tie) = -1;
    end
end
