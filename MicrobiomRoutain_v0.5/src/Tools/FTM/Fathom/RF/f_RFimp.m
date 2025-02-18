function idx = f_RFimp(model,n,verb)
% - display table of variable importance from a Random Forest
% 
% USAGE: idx = f_RFimp(model,n);
% 
% model = Random Forest model generated by f_RFclassTrain
% n     = number of variables to list                          (default = all)
% verb  = list variable importance on display                  (default = 1)
% 
% idx   = index to variables, sorted decreasing by z-score (or meanAcc)
% 
% SEE ALSO: f_RFimp

% -----References:-----
% http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#micro

% -----Author:-----
% by David L. Jones, Sep-2009
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Feb-2014: transpose 'idx' when it's created


% Extract variable labels & importance measures as column vectors:
X_txt   = model.X_txt;
meanAcc = model.meanAcc';
SD      = model.meanAccSD';
zP      = model.zP';

% Create index to variable number:
idx = (1:size(model.X,2))';

% -----Check input & set defaults:-----
if (nargin < 2), n    = size(model.X,2); end % default list all variables
if (nargin < 3), verb = 1;  end              % default list output in display

if isnan(meanAcc(1))
   error('MODEL was generated with opt.importance=0!');
end

if (n>idx(end)), error('Size mismatch b/n N and model.X!'); end
% -------------------------------------


% List the 'n' most important variables, sorted by their z-scores:
% 
var = [idx meanAcc SD zP];
var = sortrows(var,-2); % sort descending by meanAcc
var = var(1:n,:);       % subset of most important variables
idx = var(:,1);         % extract index to these variables

if (verb>0)
   % Print table:
   fprintf('==================================================\n');
   fprintf('=====      Variable Importance Measures:     =====\n');
   fprintf('==================================================\n');
   fprintf('col 1: variable label\n');
   fprintf('col 2: meanAcc\n');
   fprintf('col 3: SD\n');
   fprintf('col 4: zP\n');
   fprintf('---------------------\n');
   for i=1:n
      impTxt = sprintf('%7.3f %7.3f %7.3f\n',var(i,2:end)'); % all fields = 7 digits
      fprintf(['%s ' impTxt],X_txt{idx(i)});
   end
   fprintf('==================================================\n');
end
