function imagesc_row(b,cm,scale,rshp)
%IMAGESC_ROW  Multiple 2D imagesc in rows 
%               Dim1+2=imagesc; dim3=rows
% imagesc_row(b,cm,scale)
%         b  Real-valued images (n1,n2,n3,n4)
%        cm  colormap (min max)             (opt; default=min+max of image)
%     scale  individual scaling (n3,n4)     (opt; default=no scaling)
%            'row' scale rows
%            'col' scale columns
%            'ind' scale subimages individually
%      rshp  Reshape data to 4D for best display (opt; default=false)
%
% 10/2013 Rolf Schulte
%  7/2014 RFS: introduce reshaping
if (nargin<1), help(mfilename); return; end;
if ~isreal(b(:)), warning('imagesc_row:real','b must be real'); end
if ~exist('cm','var'),    cm = []; end
if ~exist('scale','var'), scale = []; end
if ~exist('rshp','var'),  rshp = []; end
if isempty(rshp),         rshp = false; end
cb = false;
if ~isempty(findobj(gcf,'tag','Colorbar')), cb = true; end
% clf;

if rshp,
    [n1,n2,n3] = size(b);
    nn = ceil(sqrt(n3));
    b = reshape(b,[n1 n2 n3]);
    b1 = zeros(n1,n2,nn,nn);
    for l4=1:nn,
        for l3=1:nn,
            ll = l3 + (l4-1)*nn;
            if ll<=n3, b1(:,:,l3,l4) = b(:,:,ll); end
        end
    end
    b = b1;
end

n = size(b);
if (length(n)~=2 && length(n)~=3 && length(n)~=4), error('b must be 2D, 3D or 4D'); end
[n1,n2,n3,n4] = size(b);
if strcmpi(scale,'row'), 
    scale = 1./(ones(n3,1)*squeeze(max(max(max(b,[],1),[],2),[],3)).'); 
end
if strcmpi(scale,'col'),
    scale = 1./(squeeze(max(max(max(b,[],1),[],2),[],4))*ones(1,n4)); 
end
if strcmpi(scale,'ind'), scale = 1./squeeze(max(max(b,[],1),[],2)); end

if ~isempty(scale)
    ns = size(scale);
    if ns(1)~=n3, error('size(scale,1)=%g ~= size(b,3)=%g',ns(1),n3); end
    if ns(2)~=n4, error('size(scale,2)=%g ~= size(b,4)=%g',ns(2),n4); end
end
% bb = zeros(n1*n3,n2*n4);
bb = zeros(n1*n4,n2*n3);
for l3=1:n3,
    for l4=1:n4,
%         ind1 = (1:n1) + (l3-1)*n1;
%         ind2 = (1:n2) + (l4-1)*n2;
%         bb(ind1,ind2) = b(:,:,l3,l4);
        ind1 = (1:n1) + (l4-1)*n1;
        ind2 = (1:n2) + (l3-1)*n2;
        if isempty(scale),
            bb(ind1,ind2) = b(:,:,l3,l4);
        else
            bb(ind1,ind2) = scale(l3,l4)*b(:,:,l3,l4);
        end
    end
end
imagesc(bb);
hold on
for l3=(0:n3)+1/n2/2, plot(l3*n2*[1 1],[0 1]*n1*n4+0.5,'w'); end;
for l4=(0:n4)+1/n1/2, plot([0 1]*n2*n3+0.5,l4*n1*[1 1],'w'); end; 
hold off
if ~isempty(cm), caxis(cm); end
if cb, colorbar; end
axis image
set(gca,'XTick',-1,'YTick',-1,'XColor',[1 1 1],'YColor',[1 1 1]);
xlabel('','Color',[0 0 0])
ylabel('','Color',[0 0 0])
box off
% drawnow;
