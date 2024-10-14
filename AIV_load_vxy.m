function a=AIV_load_vxy(fn)

if ~isfile(fn)  error(sprintf('AIV_load_vxy: File "%s" does not exist',fn)); end

a    = load(fn);

if ~isfloat(a.vx)
  if ~isfield(a,'vx') error(sprintf('AIV_load_vxy: v data is not float but no scale factor')); end
  a.vx = a.sf*double(a.vx);
  a.vy = a.sf*double(a.vy);
end

end

