function a=AIV_names(expno)

if expno(1:2)=='CR'  % Coral reef experiments   CR-2024-07-10-Exp-06  
  DD=CR_dirs;
  dir=expno;
  a.expno=expno(end-1:end);
else
  DD=AIV_dir;
  dir = AIV_date(expno);
  a.expno=expno;
end

% Where to look for the raw data experimental files
a.outdir=DD;
a.work     = sprintf('%s/data-raw/%s',DD,dir);
a.slide    = sprintf('%s/slide',   a.work);
a.null     = sprintf('%s/null',    a.work);
a.movie    = sprintf('%s/images',  a.work);

a.results  = sprintf('%s/data/%s',DD,dir);
a.velocity = sprintf('%s/velocity',a.results);

if ~isdir(a.results)  mkdir(a.results);  end;
if ~isdir(a.velocity) mkdir(a.velocity); end;

% Where to save the output data files
a.mat.param  = sprintf('%s/param.mat',  a.results);
a.mat.slide  = sprintf('%s/slide.mat',  a.results);
a.mat.tank   = sprintf('%s/tank.mat',   a.results);
a.mat.null   = sprintf('%s/null.mat',   a.results);
a.mat.coord  = sprintf('%s/coord.mat',  a.results);
a.mat.vstats = sprintf('%s/vstats.mat', a.results);
a.mat.traj   = sprintf('%s/traj-%%s.mat',  a.results);
a.mat.cam    = sprintf('%s/data/2024-04-26-Calibration/cam4.mat', DD);
end


