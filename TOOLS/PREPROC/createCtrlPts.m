function ctrl = createCtrlPts(x, varargin)
%Written by John Corring; email me at johncorring@gmail.com. RDM 2015 v1.  
%This function creates a set of control points that will be used to form 
%the basis for the TPS. x is the template set, varargin consists up to 3 
%options:
%           length(varargin) = 1 => regularly sample x by varargin{1}
%           length(varargin) = 2 => cluster the points to belond to one of
%           varargin{2} clusters, and for clustering use the method
%           described by varargin{2}.
%           TODO::FINISH THIS STUFF
%           'kmeans' calls MATLAB's kmeans
%           'em' calls my EM (Gaussian)
%           'mml-em' calls Figueredo and Jain's algorithm (Gaussian)

    k = length(varargin);
    [n d] = size(x);

    if (k == 0)
        %%just use the points
        ctrl = x;
    elseif (k == 1)
        %%regularly sample the points
        m = varargin{1};
        if(m > n)
           error('Not enough points for that rate of control points.');
        end
        ctrl = x(1:m:end,:);
    elseif (k == 2)
        %%cluster the points and use the centroids
        m = varargin{1};
        method = varargin{2};

        %%TODO: clustering code here
    elseif (k > 2)
        
        dxsh = varargin{1};
        
        dysh = varargin{2};
        freq = varargin{3};
        
        if(d == 3)
            dzsh = varargin{3};
            freq = varargin{4};
        end
        
        
        minx = min(x(:,1));
        maxx = max(x(:,1));
        
        miny = min(x(:,2));
        maxy = max(x(:,2));
        
        if(d == 3)
            minz = min(x(:,3));
            maxz = max(x(:,3));
        end
        
        ctrlx = linspace(minx-dxsh, maxx+dxsh, freq);
        ctrly = linspace(miny-dysh, maxy+dysh, freq);
        
        if(d == 3)
            ctrlz = linspace(minz - dzsh, maxz + dzsh, freq);
        end
        
        [ctrlx2 ctrly2] = meshgrid(ctrlx, ctrly);
        
        
        ctrl = [ctrlx2(:) ctrly2(:)];
        
        if(d==3)
            [ctrlx2 ctrly2 ctrlz2] = meshgrid(ctrlx, ctrly, ctrlz);
            ctrl = [ctrlx2(:) ctrly2(:) ctrlz2(:)];
        end
    else
        %%failure case
        error('createCtrlPts does not recognize that option set. Please read the functions.');
    end

end