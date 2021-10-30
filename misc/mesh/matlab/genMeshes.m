% fd=@(p) dsphere(p,0,0,0,1);
% [p,t]=distmeshsurface(fd,@huniform,0.2,1.1*[-1,-1,-1;1,1,1]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\sphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\sphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\sphere.ply",fliplr(t),p);
% 
% a = 1.01;
% area = @(c) 2 * (a^2 + c^2 / sqrt( 1- c^2/a^2) * atanh( sqrt( 1- c^2/a^2) )) - 4;
% [c,~] = fsolve(area, 1)
% fd=@(p) p(:,1).^2/a^2 +p(:,2).^2/a^2 +p(:,3).^2/c^2 -1;
% [p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-2.1,-1.5; 2.1,2.1,1.5]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\slightlyOblate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\slightlyOblate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\slightlyOblate.ply",fliplr(t),p);
% 
% a = 1.1;
% area = @(c) 2 * (a^2 + c^2 / sqrt( 1- c^2/a^2) * atanh( sqrt( 1- c^2/a^2) )) - 4;
% [c,~] = fsolve(area, 1)
% fd=@(p) p(:,1).^2/a^2 +p(:,2).^2/a^2 +p(:,3).^2/c^2 -1;
% [p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-2.1,-1.5; 2.1,2.1,1.5]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\oblate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\oblate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\oblate.ply",fliplr(t),p);
% 
% a = 0.85;
% area = @(c) 2 * a^2 * ( 1 + c / a / sqrt( 1- a^2/c^2) * asin( sqrt( 1- a^2/c^2) )) - 4;
% [c,~] = fsolve(area, 1)
% fd=@(p) p(:,1).^2/a^2 +p(:,2).^2/a^2 +p(:,3).^2/c^2 -1;
%     [p,t]=distmeshsurface(fd,@huniform,0.2,[-1.5,-1.5,-2; 1.5,1.5,2]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\prolate.ply",fliplr(t),p);
% 
% a = 0.6;
% area = @(c) 2 * a^2 * ( 1 + c / a / sqrt( 1- a^2/c^2) * asin( sqrt( 1- a^2/c^2) )) - 4;
% [c,~] = fsolve(area, 1)
% fd=@(p) p(:,1).^2/a^2 +p(:,2).^2/a^2 +p(:,3).^2/c^2 -1;
%     [p,t]=distmeshsurface(fd,@huniform,0.2,[-1.5,-1.5,-2; 1.5,1.5,2]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\longProlate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\longPblate.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\longProlate.ply",fliplr(t),p);
% 
% fd=@(p) dsphere(p,0,0,0,1);
%     fh=@(p) 0.05+0.5*dsphere(p,0,0,1,0);
%     [p,t]=distmeshsurface(fd,fh,0.15,1.1*[-1,-1,-1;1,1,1]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\forTubeSphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\forTubeSphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\forTubeSphere.ply",fliplr(t),p);
% 
% fd=@(p) dsphere(p,0,0,0,1);
%     fh=@(p) 0.035+0.55*dsphere(p,0,0,1,0);
%     [p,t]=distmeshsurface(fd,fh,0.15,1.1*[-1,-1,-1;1,1,1]);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\finerForTubeSphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\finerForTubeSphere.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\finerForTubeSphere.ply",fliplr(t),p);

fd=@(p) sqrt(sum(p.^2,2))-1;
fh=@(p) 0.01+0.03*dcircle(p,0,0,0);
[p,t]=distmesh2d(fd,fh,0.02,[-1,-1;1,1],[]);
p = [p,zeros(size(p,1),1)];
[~,I] = min(sum(abs(p)')'); % find the center point index -> remember -1 for c++
plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\run\input-file\patch_flip.ply",t,p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\vesicle\input-file\patch.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\patch\input-file\patch.ply",fliplr(t),p);
% plywrite("C:\Users\Kieran\Dev\2020-Mem3DG-Applications\examples\tube\input-file\patch.ply",fliplr(t),p);
 