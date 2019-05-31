function [t,idat] = readNRNbin_Vclamp(fname,dtype)

% READNRNBIN_Vclamp(fname)	    Read binary file fname.out, output from NEURON, 
%                       into appropriate variables.
%
%  dtype = 0 if 'native'; = 1 if 'ieee-be'
%
%   modified from readNRNbin_Vonly.m by Christina Weaver
%   (christina.weaver@fandm.edu) on 4/5/12.
%

if( dtype == 1 )
    dstrg = 'ieee-be';
else dstrg = 'native';
end;

fnameOK = 0;

if( ~isempty(fname) ) 
    fnameOK = 1;
    finname = sprintf('%s.Ibin',fname);
    if( ~exist(finname,'file') ) 
        fnameOK = 0; 
    end;
    fprintf(1,'Reading data from binary file %s\n',finname);
end;

if( ~fnameOK )
    fprintf(1,'Error reading %s.Vbin, exiting readNRNbin_Vclamp\n',fname);
    t=[];   v=[];   st_data=[];
    return;
end;
fin = fopen(finname,'r',dstrg);

% file format taken from proc synTweak() in
%   ~/Neuron/LuebkeAm_forCluster/main_PFC_simEPSC_all.hoc:
% 
% vecsz	//size of all vectors
% vectors written:
% 
% t, i (of size vecsz)
% 


[npts] = fread(fin,1,'double');
[t, csz] = fread(fin,npts,'double');
[idat, csz] = fread(fin,npts,'double');

fclose(fin);

