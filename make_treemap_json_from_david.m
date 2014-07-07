function [js,err]=make_treemap_json_from_david(x)
%function [js,err]=make_treemap_json_from_david(x)
%
%IN:
%   x is a struct contining information for a JAVA string encoding a JSON
%   object to be generated
%
%js is a string in json format containing relevant data for a Jit
%icicle plot
% generating the JSON file to hold the treemap data may error out if the
% heap is too small and the gene list is too big


%JSON fields: id, name, data{area,dim,color}, children

%set up the parent node x, to be a summary of the gene list

err=0;
if isdeployed, addpath(fullfile(ctfroot,'json'));
else, addpath(fullfile(pwd,'json'));end
if isdeployed,javaaddpath(fullfile(ctfroot,'JSON-java.jar'));
else, javaaddpath(fullfile(pwd,'JSON-java.jar'));end
%addpath('./json');
%javaaddpath('./JSON-java.jar');
try 
    js=JSON.dump(x);
catch me
    err=1; js=[];
    return;
end
if isdeployed,javarmpath(fullfile(ctfroot,'HiTSelect','JSON-java.jar'));
else, javarmpath(fullfile(pwd,'JSON-java.jar'));end
%javarmpath('./JSON-java.jar');
js=strrep(js,'area','$area');%the $ prefix is needed by the Jit treemap code
js=strrep(js,'dim','$dim');
js=strrep(js,'color','$color');
%js=strrep(js,'[]','{}');%this is unlikely to be necessary
js=strrep(js,'"','\"');%the treemap code sets the data as a java string
js=['"' js '"'];