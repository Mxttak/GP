% generates a cpp file with access to vriables

function genMexTemplate(fname,numInputs,numOutputs,varargin)

  nvarargin = numel(varargin);
  if(nvarargin~=(numInputs+numOutputs))
    error('wrong number of input arguments')
  end

  % construct header
    header = ['// automatically generated template\n', ...
              '\n',...
              '#include <mex.h>\n\n',...
              'void yourfunc_('];
    for i = 1:(numInputs+numOutputs)
      if(i>1)
        header = [header,', \n               '];%#ok
      end
      %header = [header,'Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>* ', varargin{i},', mwSize* dims',varargin{i},', mwSize numDims',varargin{i}]; %#ok
      header = [header,'double* ', varargin{i},', mwSize* dims',varargin{i},', mwSize numDims',varargin{i}]; %#ok
    end
    header = [header,'){\n\n\t// put your code here\n\n}\n\n','/* ---------------------------------------------------------------------- */\n'];
  % END: construct header
  
  % construct body
    body = ['void mexFunction(int numOutputs, mxArray *outputs[],\n',...
            '                 int numInputs, const mxArray *inputs[]){\n'];
    % extract inputs
    for i = 1:numInputs      
      tmp = [
              % get pointer to the input
              '\t// get pointer to the input ',varargin{i},'\n',...
              '\tdouble* ptr',varargin{i},' = mexGetPr(inputs[',num2str(i-1),']);\n',...
              '\tmwSize numDims',varargin{i},' = mxGetNumberOfDimensions(inputs[',num2str(i-1),']);\n',...
              '\tmwSize* dims',varargin{i},' = (mwSize*) mxGetDimensions(inputs[',num2str(i-1),']);\n\n',...
              ];
      body = [body,tmp]; %#ok
    end
    
    body = [body,'\n'];
    
    % extract outputs
    for i = 1:numOutputs
      tmp = '';
      for j=1:numel(varargin{i+numInputs})
        if(j>1)
          tmp = [tmp,', '];%#ok
        end
        tmp = [tmp,num2str(varargin{i+numInputs}(j))];%#ok
      end
      
      tbody = ['\t// create output ',varargin{i+numInputs},'\n',...
              '\tmwSize numDims',varargin{i+numInputs},' = // <<<<< TODO\n',...
              '\tmwSize dims',varargin{i-1+numInputs},' = // <<<<< TODO\n',...
              '\toutputs[',num2str(i-1),'] = mxCreateNumericArray(numDims',varargin{i+numInputs},', dims',varargin{i+numInputs},' mxDOUBLE_CLASS, mxREAL);\n',...
              '\tdouble* ptr',varargin{i+numInputs},' = mxGetPr(outputs[',num2str(i-1),']);\n'];
      body = [body,tbody]; %#ok
    end
    
    % close body
    body = [body, '\n\tyourfunc_('];
    for i = 1:(numInputs+numOutputs)
      if(i>1)
        body = [body,', ']; %#ok
      end
      body = [body, 'ptr',varargin{i},', dims',varargin{i},', numDims',varargin{i}]; %#ok
    end
    body = [body,');\n\t\n\treturn;\n}']; % TODO: free allocated memory
  % END: construct body
  
  f = fopen([fname,'.cpp'],'w');
  fprintf(f,[header,body]);
  fclose(f);
end