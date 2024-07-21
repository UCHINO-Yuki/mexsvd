# mexeig - LAPACK's Singular solvers Interface for MATLAB

## Installation

1. Download all files.
2. Add the path of "mexsvd-main" folder to MATLAB's search path.
3. Execute the following command to compile the mex files:
```
mexsvd('compile');
```

## Usage

```
[Outputs,...] = mexsvd(FunctionName, Matrix, Options...);

%
% FunctionName: 'gesdd', 'gesvd', or 'gesvdx'
% Matrix      : single or double matrix
% Option1     : 'matrix', 'vector', 'econ', or 0
% Option2     : 'matrix', 'vector', 'econ', or 0
%
% 'econ' and 0 are the same
% The default is 'econ' for 'gesvdx'
%
```
 
