function varargout = write2txt(varargin)

fname = varargin{nargin};
y = varargin{nargin-1};
x = zeros(nargin-2,length(varargin{1}));
for i = 1:nargin-2
    x(i,:) = varargin{i};
end
if exist(fname,'file') == 2
    delete(fname);
end
str = cell(1,length(x(1,:)));
for i = 1:length(x(1,:))
    s = ' ';
    for j = 1:nargin-2
        s = [s,num2str(x(j,i)),' '];
    end
    s = s(2:end);
    str{i} = [s,num2str(y(i))];
end
fid = fopen(fname,'w');
for i = 1:length(x(1,:))
    fprintf(fid,'%s\n',str{i});
end
fclose(fid);
varargout{1} = [];