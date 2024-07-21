function [res1,res2,res3]=mexsvd(func,A,in2,in3)
% func : 'compile'... Only compiling mex
%        'gesdd', 'gesvd', 'gesvdx'
% A    : single or double matrix
% in2  : 'matrix', 'vector', 'econ', 0
% in3  : 'matrix', 'vector', 'econ', 0
% 'econ' and 0 are the same.
% The default is 'econ' for 'gesvdx'
% Only compiling mex: mexsvd('compile');
% The order of the arguments after in2 (including in2) is arbitrary.
% The arguments after in2 (including in2) are optional.
%
% written ... 2024-07-21 ... Yuki UCHINO

res1=0;
res2=0;
res3=0;

%% compile
if strcmp(func,'compile')
    disp('Compiling ... ')
    try
        mex -R2018a ./private/dgesdd.c -lmwlapack -silent -outdir ./private
        mex -R2018a ./private/dgesvd.c -lmwlapack -silent -outdir ./private
        mex -R2018a ./private/dgesvdx.c -lmwlapack -silent -outdir ./private
        mex -R2018a ./private/sgesdd.c -lmwlapack -silent -outdir ./private
        mex -R2018a ./private/sgesvd.c -lmwlapack -silent -outdir ./private
        mex -R2018a ./private/sgesvdx.c -lmwlapack -silent -outdir ./private
    catch
        error('error');
    end
    disp('Compilation Successful.')
    return
end

%% adjusting the arguments
if nargin == 3
    if in2 == 0
        in2 = 'econ';
    end
end
if nargin == 4
    if in2 == 0
        in2 = 'econ';
    elseif in3 == 0
        in3 = 'econ';
    end
end

%% swapping arguments
if nargin == 4
    if ~strcmp(in3,'econ')
        tmp = in2;
        in2 = in3;
        in3 = tmp;
    end
end

%% setting flags
if nargin == 2
    ecoS = 0;
    if nargout == 1
        MAT = 0;
        JOBU = 'N';
        JOBV = 'N';
        JOBZ = 'N';
    elseif nargout == 2
        MAT = 1;
        JOBU = 'V';
        JOBV = 'N';
        JOBZ = 'A';
    else
        MAT = 1;
        JOBU = 'V';
        JOBV = 'V';
        JOBZ = 'A';
    end
elseif nargin == 3
    if strcmp(in2,'econ')
        ecoS = 1;
        if nargout == 1
            MAT = 0;
            JOBU = 'N';
            JOBV = 'N';
            JOBZ = 'N';
        elseif nargout == 2
            MAT = 1;
            JOBU = 'S';
            JOBV = 'N';
            JOBZ = 'S';
        else
            MAT = 1;
            JOBU = 'S';
            JOBV = 'S';
            JOBZ = 'S';
        end
    elseif strcmp(in2,'matrix')
        ecoS = 0;
        MAT = 1;
        if nargout == 1
            JOBU = 'N';
            JOBV = 'N';
            JOBZ = 'N';
        elseif nargout == 2
            JOBU = 'V';
            JOBV = 'N';
            JOBZ = 'A';
        else
            JOBU = 'V';
            JOBV = 'V';
            JOBZ = 'A';
        end
    else
        ecoS = 0;
        MAT = 0;
        if nargout == 1
            JOBU = 'N';
            JOBV = 'N';
            JOBZ = 'N';
        elseif nargout == 2
            JOBU = 'V';
            JOBV = 'N';
            JOBZ = 'A';
        else
            JOBU = 'V';
            JOBV = 'V';
            JOBZ = 'A';
        end
    end
else
    ecoS = 1;
    if strcmp(in2,'matrix')
        MAT = 1;
    else
        MAT = 0;
    end
    if nargout == 1
        JOBU = 'N';
        JOBV = 'N';
        JOBZ = 'N';
    elseif nargout == 2
        JOBU = 'S';
        JOBV = 'N';
        JOBZ = 'S';
    else
        JOBU = 'S';
        JOBV = 'S';
        JOBZ = 'S';
    end
end

%% running LAPACK subroutines
if isa(A,'double')
    switch func
        case 'gesdd'
            if nargout == 1
                res1=dgesdd(A,JOBZ,MAT,ecoS);
            else
                [res1,res2,res3]=dgesdd(A,JOBZ,MAT,ecoS);
            end
            
        case 'gesvd'
            if strcmp(JOBU,'V')
                JOBU = 'A';
            end
            if strcmp(JOBV,'V')
                JOBV = 'A';
            end
            if nargout == 1
                res1=dgesvd(A,JOBU,JOBV,MAT,ecoS);
            else
                [res1,res2,res3]=dgesvd(A,JOBU,JOBV,MAT,ecoS);
            end
            
        case 'gesvdx'
            if strcmp(JOBU,'S')
                JOBU = 'V';
            end
            if strcmp(JOBV,'S')
                JOBV = 'V';
            end
            if nargout == 1
                res1=dgesvdx(A,JOBU,JOBV,MAT);
            else
                [res1,res2,res3]=dgesvdx(A,JOBU,JOBV,MAT);
            end
            
        otherwise
            error('error: func');
    end
else
    switch func
        case 'gesdd'
            if nargout == 1
                res1=sgesdd(A,JOBZ,MAT,ecoS);
            else
                [res1,res2,res3]=sgesdd(A,JOBZ,MAT,ecoS);
            end
            
        case 'gesvd'
            if strcmp(JOBU,'V')
                JOBU = 'A';
            end
            if strcmp(JOBV,'V')
                JOBV = 'A';
            end
            if nargout == 1
                res1=sgesvd(A,JOBU,JOBV,MAT,ecoS);
            else
                [res1,res2,res3]=sgesvd(A,JOBU,JOBV,MAT,ecoS);
            end
            
        case 'gesvdx'
            if strcmp(JOBU,'S')
                JOBU = 'V';
            end
            if strcmp(JOBV,'S')
                JOBV = 'V';
            end
            if nargout == 1
                res1=sgesvdx(A,JOBU,JOBV,MAT);
            else
                [res1,res2,res3]=sgesvdx(A,JOBU,JOBV,MAT);
            end
            
        otherwise
            error('error: func');
    end
end

%% transpose VT
if nargout == 3
    res3 = res3';
end

end
