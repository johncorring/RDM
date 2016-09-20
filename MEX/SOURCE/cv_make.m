
function [objfiles, timestamp_out] = cv_make (f)
 
% Copyright 2006-2012, Timothy A. Davis, http://www.suitesparse.com
 
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end
 
if (~isempty (strfind (computer, '64')))
    fprintf ('Compiling (64-bit)\n') ;
    mexcmd = 'mex -largeArrayDims' ;
else
    fprintf ('Compiling (32-bit)\n') ;
    mexcmd = 'mex' ;
end
 
%create makefile: what needs to be compiled?
! sh build_make.sh
 
cvm = { };
cv = { };
files = textread('makefile', '%s', 'delimiter', ' ', ...
                'whitespace', '');
%parse makefile:
%(c)                
i = 1;
while(char(files{i}) ~= '#')
    cv = [cv; files(i)];
    i = i+1;
end
i = i+1;
%(mex)
while((files{i}) ~= '#')
    cvm = [cvm; files(i)];
    i = i+1;
end
 
 
details = 0 ;
kk = 0 ;
if (nargin == 0)
    force = 0 ;
elseif (ischar (f))
    fprintf ('cv_make: compiling ../../Source files and %s_mex.c\n', f) ;
    force = 0 ;
    cvm = {f} ;
else
    force = f ;
    details = details | (force > 1) ;                                       %#ok
    if (force & details)                                                    %#ok
        fprintf ('cv_make: re-compiling everything\n') ;
    end
end
if (force)
    fprintf ('(compiling all of CSparse from scratch)\n') ;
end
 
if (isempty (cvm))
    % mexFunctions, of the form cs_add_mex.c, etc, in this directory
    cvm = { } ;
        % add cs_mynewfunc to the above list
end
 
if (pc)
    obj = '.obj' ;
else
    obj = '.o' ;
end
 
srcdir = 'AUX/' ;
hfile = '../../Include/cs.h' ;
 
%-------------------------------------------------------------------------------
% With a current Microsoft compiler installed in the UF CISE lab, when CSparse
% is located in a shared network folder, the -I../../Include option fails.  The
% compiler -I or /I switch cannot use a relative path, even when using
% /I..\..\Include.  Relative paths work fine in the same setup for the source
% code files, which are in ../../Source/*.c.  This is very odd.  As a
% work-around, absolute paths are now used in this version.
 
% Change the following to 0 if relative paths work fine with /I on Windows:
relative_paths_do_not_work = 1 ;
 
if (pc & relative_paths_do_not_work)
    % begin pain
    here = pwd ;
    cd ../../Include
    % use quotes in case the path has spaces
    mexcmd = [mexcmd ' -I"' pwd '"'] ;
    cd (here)
    % end pain
else
    % Linux, Unix, and Mac, are just fine.
    mexcmd = [mexcmd ' -I../../Include'] ;
end
%-------------------------------------------------------------------------------
 
% compile each CSparse source file
if (nargout > 0)
    objfiles = [ obj] ;
end
timestamp = 0;
anysrc = 0;
CV = [];
for i = 1:length (cv)
    [s t kk] = compile_source (srcdir, cv {i}, obj, hfile, force, ...
        kk, details, mexcmd) ;
    timestamp = max (timestamp, t) ;
    anysrc = anysrc | s ;                                                   %#ok
    CV = [CV ' ' cv{i} obj] ;                                               %#ok
    if (nargout > 0)
        objfiles = [objfiles ' ../MEX/' cv{i} obj] ;   %#ok
    end
end
 
% compile each CSparse mexFunction
obj = ['.' mexext] ;
for i = 1:length (cvm)
    [s t] = cv_must_compile ('', cvm{i}, '.c', obj, hfile, force) ;
    timestamp = max (timestamp, t) ;
    if (anysrc | s)                                                         %#ok
        cmd = sprintf ('%s -O %s.c %s -output %s', ...
            mexcmd, cvm{i}, CV, cvm{i}) ;
        kk = do_cmd (cmd, kk, details) ;
    end
end
 
fprintf ('\n') ;
if (nargout > 1)
    timestamp_out = timestamp ;
end
 
if (force)
    fprintf ('Successfully compiled.\n') ;
end
end
 
%-------------------------------------------------------------------------------
function [s,t,kk] = compile_source (srcdir, f, obj, hfile, force, ...
    kk, details, mexcmd)
% compile a source code file in ../../Source, leaving object file in
% this directory.
[s t] = cv_must_compile (srcdir, f, '', obj, hfile, force) ;
if (s)
    cmd = sprintf ('%s -O -c %s%s.c', mexcmd, srcdir, f) ;
    kk = do_cmd (cmd, kk, details) ;
end
end
 
%-------------------------------------------------------------------------------
function kk = do_cmd (s, kk, details)
%DO_CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
        fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;
end