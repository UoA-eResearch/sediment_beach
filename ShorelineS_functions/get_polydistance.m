function [dmin,xcr,ycr,dmax,xcm,ycm] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
%function [dmin,xcr,ycr,dmax,xcm,ycm] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
%
% INPUT: 
%    Xr          : x-coordinate of reference line [m]
%    Yr          : y-coordinate of reference line [m]
%    Xc          : x-coordinate of coastline [m]
%    Yc          : y-coordinate of coastline [m]
%
% OUTPUT:
%    dmin        : distance from reference line point to the specified coastline [m] (negative is at right-side of the line)
%    xcr         : x-coordinate of coastline point as projected on the grid [m] (closest value to the grid)
%    ycr         : y-coordinate of coastline point as projected on the grid [m] (closest value to the grid)
%    dmin        : distance from reference line point to the specified coastline [m] (negative is at right-side of the line)(farthest crossing from the grid)
%    xcr         : x-coordinate of coastline point as projected on the grid [m] (farthest crossing from the grid)
%    ycr         : y-coordinate of coastline point as projected on the grid [m] (farthest crossing from the grid)
%
% EXAMPLE:
%    Xr=[20.5:185.5]';
%    Yr=1-Xr/185;
%    Xc=[1:200]';
%    Yc=2+sin(Xc.^0.5).^2;
%    [dmin,xcr,ycr]=get_segmentdistance(Xr, Yr, Xc, Yc);
%
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 Deltares
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    % dx,dy values
    if nargin<5
        Lcrit=500;  % maximum cross-shore distance to reference line (at both sides)
    end
    Xr=Xr(:);
    Yr=Yr(:);
    dx=Xr(2:end)-Xr(1:end-1); %for line in-between each reference point
    dy=Yr(2:end)-Yr(1:end-1);
    dx2 = [dx(1);(dx(2:end)+dx(1:end-1))/2;dx(end)];
    dy2 = [dy(1);(dy(2:end)+dy(1:end-1))/2;dy(end)];
    xcr=nan(length(Xr),1);
    ycr=nan(length(Xr),1);
    dmin=nan(length(Xr),1);
    if nargout>3
    xcm=nan(length(Xr),1);
    ycm=nan(length(Xr),1);
    dmax=nan(length(Xr),1);
    end
    Lxy=(dx2.^2+dy2.^2).^0.5;
    Xrnormal = [Xr,Xr]+[-Lcrit./Lxy.*dy2,Lcrit./Lxy.*dy2];
    Yrnormal = [Yr,Yr]+[Lcrit./Lxy.*dx2,-Lcrit./Lxy.*dx2];
    
    % Intersect all cross-shore normals with the coastline at once (in chunks
    % of normals to limit memory). This uses the same crossing computation as
    % get_intersections, and picks per normal the crossing closest to the
    % seaward end of the normal (smallest ui), as the former loop did.
    eps=1e-5;                         % same tolerance as in get_intersections
    Xc=Xc(:)';
    Yc=Yc(:)';
    nseg=length(Xc)-1;
    if nseg>=1
        xj=Xc(1:nseg);
        yj=Yc(1:nseg);
        dxc=Xc(2:end)-xj;
        dyc=Yc(2:end)-yj;
        xjmin=min(xj,Xc(2:end));
        xjmax=max(xj,Xc(2:end));
        yjmin=min(yj,Yc(2:end));
        yjmax=max(yj,Yc(2:end));
        rc2=dyc./dxc;
        y2r=yj-xj.*rc2;
        nchunk=max(1,floor(2e6/nseg));
        for i0=1:nchunk:length(Xr)
            ii=(i0:min(i0+nchunk-1,length(Xr)))';
            xi=Xrnormal(ii,1);
            yi=Yrnormal(ii,1);
            dx1=Xrnormal(ii,2)-xi;
            dy1=Yrnormal(ii,2)-yi;
            rc1=dy1./dx1;
            y1r=yi-xi.*rc1;
            
            xc=(y2r-y1r)./(rc1-rc2);
            yc=rc1.*xc+y1r;
            both=(dx1~=0) & (dxc~=0);
            xc(~both)=nan;
            yc(~both)=nan;
            id2=(dx1==0) & (dxc~=0);
            if any(id2(:))
                xcr2=xi+zeros(1,nseg);
                ycr2=rc2.*xcr2+y2r;
                xc(id2)=xcr2(id2);
                yc(id2)=ycr2(id2);
            end
            id3=(dx1~=0) & (dxc==0);
            if any(id3(:))
                xcr3=xj+zeros(length(ii),1);
                ycr3=rc1.*xcr3+y1r;
                xc(id3)=xcr3(id3);
                yc(id3)=ycr3(id3);
            end
            idnan=xc<max(min(xi,Xrnormal(ii,2)),xjmin)-eps | xc>min(max(xi,Xrnormal(ii,2)),xjmax)+eps;
            xc(idnan)=nan;
            yc(idnan)=nan;
            idnan=yc<max(min(yi,Yrnormal(ii,2)),yjmin)-eps | yc>min(max(yi,Yrnormal(ii,2)),yjmax)+eps;
            xc(idnan)=nan;
            yc(idnan)=nan;
            valid=~isnan(xc);
            
            % fraction along the normal (0 = seaward end, 1 = landward end)
            ui=((xc-xi).*dx1+(yc-yi).*dy1)./(dx1.^2+dy1.^2);
            ui=min(max(ui,0),1);
            
            uimin=ui;
            uimin(~valid)=inf;
            [umin,jmin]=min(uimin,[],2);
            has=isfinite(umin);
            idx=sub2ind(size(xc),find(has),jmin(has));
            xcr(ii(has))=xc(idx);
            ycr(ii(has))=yc(idx);
            dmin(ii(has))=(0.5-umin(has))*Lcrit*2;
            if nargout>3
                uimax=ui;
                uimax(~valid)=-inf;
                [umax,jmax]=max(uimax,[],2);
                idx=sub2ind(size(xc),find(has),jmax(has));
                xcm(ii(has))=xc(idx);
                ycm(ii(has))=yc(idx);
                dmax(ii(has))=(0.5-umax(has))*Lcrit*2;
            end
        end
    end
end 
