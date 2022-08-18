function f = prior_xi (xi, K)
  %% usage:  f = prior_xi (xi, K)
  %%
  %% Calculate the value of prior(xi) for a vector 
  %% of xi.
  %% In:
  %%    xi   -- matrix of any size
  %%    K    -- scalar, cardinality of the outcome space
  %% Out:
  %%    f  -- the value of the prior, same size as xi
  %%
  %%  Uses:  none.
  %% 
  %% (c) Ilya Nemenman, 2002--2006
  %% Distributed under GPL, version 2
  %% Copyright (2006). The Regents of the University of California.

  %% This material was produced under U.S. Government contract
  %% W-7405-ENG-36 for Los Alamos National Laboratory, which is
  %% operated by the University of California for the U.S.
  %% Department of Energy. The U.S. Government has rights to use,
  %% reproduce, and distribute this software.  NEITHER THE
  %% GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS
  %% OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
  %% SOFTWARE.  If software is modified to produce derivative works,
  %% such modified software should be clearly marked, so as not to
  %% confuse it with the version available from LANL.
  %%
  %% Additionally, this program is free software; you can redistribute
  %% it and/or modify it under the terms of the GNU General Public
  %% License as published by the Free Software Foundation; either
  %% version 2 of the License, or (at your option) any later version.
  %% Accordingly, this program is distributed in the hope that it will
  %% be useful, but WITHOUT ANY WARRANTY; without even the implied
  %% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  %% See the GNU General Public License for more details.

  if(length(K)~=1) error('Wrong dimension of K'); end;

  %% put any prior you wish instead of this line
  f = ones(size(xi));

