function varargout = kvf(varargin)
% Read and write key-value files.
  [varargout{1:nargout}] = feval(varargin{:});
end

% ------------------------------------------------------------------------------
% Public

function Write(fn, c, allow_overwrite)
% c is a struct containing string and numeric fields. Write it to the file
% fn. This fill will be used by the program fdra2c. The file is not
% necessarily portable to other machines.
  if (nargin < 3) allow_overwrite = false; end
  if (~allow_overwrite && exist(fn, 'file'))
    error(sprintf('Write:Kvf: I don''t want to overwrite %s!\n', fn));
  end
  fptr = fopen(fn, 'wb');
  if (fptr == -1) error(sprintf('Could not open %s for writing.', fn)); end
  kvf_Write(fptr, c);
  fclose(fptr);
end

function c = Read(fn)
% Load a struct from the file fn.
  fptr = fopen(fn, 'rb');
  if (fptr == -1) error(sprintf('Could not open %s for reading.', fn)); end
  c = kvf_Read(fptr);
  fclose(fptr);
end

% ------------------------------------------------------------------------------
% Private

function kvf_Write(fptr, c)
flds = fieldnames(c);
for (i = 1:numel(flds))
    d = c.(flds{i});
    if ((isnumeric(d) || islogical(d)) && all(isreal(d)))
        kvf_WriteString(fptr, flds{i});
        kvf_WriteCode(fptr, 'da');
        kvf_WriteInts(fptr, size(d));
        kvf_WriteArrayData(fptr, d);
    elseif (ischar(d))
        kvf_WriteString(fptr, flds{i});
        kvf_WriteCode(fptr, 'st');
        kvf_WriteString(fptr, d);
    elseif (isobject(d))
        kvf_Write(fptr, d)
    elseif (iscell(d))
    else
        fprintf(1, 'Write: Skipping %s.\n', flds{i});
    end
end
end

function kvf_WriteCode(fptr, code)
  if (fwrite(fptr, code, 'char') ~= 2)
    error(sprintf('Failed to write code %s.', code));
  end
end

function kvf_WriteString(fptr, s)
  if (fwrite(fptr, numel(s), 'int64') ~= 1 ||...
      fwrite(fptr, s, 'char') ~= numel(s))
    error(sprintf('Failed to write string %s.', s));
  end
end

function kvf_WriteInts(fptr, sz)
  nsz = numel(sz);
  if (fwrite(fptr, [nsz sz(:).'], 'int64') ~= nsz + 1)
    error('Failed to write ints.');
  end
end

function kvf_WriteArrayData(fptr, d)
  if (fwrite(fptr, d(:), 'double') ~= numel(d))
    error('Failed to write data.');
  end
end

function c = kvf_Read(fptr)
  c = struct();
  while (true)
    try
      fld = kvf_ReadString(fptr);
    catch
      if (~feof(fptr)) error('Error reading fieldname.');
      else break; end
    end
    code = fread(fptr, 2, 'char');
    code = char(code(:).');
    switch (code)
     case 'st'
      c.(fld) = kvf_ReadString(fptr);
     case 'da'
      c.(fld) = kvf_ReadArray(fptr);
     otherwise
      error(sprintf('Wrong code %s.', code));
    end
  end
end

function s = kvf_ReadString(fptr)
  n = kvf_ReadInts(fptr, 1);
  s = fread(fptr, n, 'char');
  if (numel(s) ~= n)
    error(sprintf('Read only %s of what should have length %d.', s, n));
  end
  s = char(s(:).');
end

function a = kvf_ReadArray(fptr)
  n = kvf_ReadInts(fptr, 1);
  sz = kvf_ReadInts(fptr, n);
  n = prod(sz);
  a = fread(fptr, n, 'double');
  if (numel(a) ~= n)
    error(sprintf('Read %d of expected %d.', numel(a), n));
  end
  a = reshape(a, sz(:).');
end
  
function a = kvf_ReadInts(fptr, n)
  a = fread(fptr, n, 'int64');
  if (any(a ~= round(a)))
    error('Expected integers.')
  end
end
