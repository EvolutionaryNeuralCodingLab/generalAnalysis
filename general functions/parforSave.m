function parforSave(varargin)
evalExpression='save(';
for i=1:numel(varargin)
    evalExpression=[evalExpression '''' varargin{i} '''' ',' ''];
end
evalExpression=[evalExpression(1:end-1) ');'];
evalin('caller',evalExpression);