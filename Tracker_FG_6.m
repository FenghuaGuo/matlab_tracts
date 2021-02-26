 classdef Tracker_FG_6 < handle
    %
    % Abstract Tractography class
    %
    % see also DTITracker, FODTracker
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (Access = public) %protected)
        f;
        v2w;
        stepSize;
        threshold;
        maxAngle;
        lengthRange;
        fun;
        iterations = 1;
        pb = false;
    end
    
    methods (Access = public)
        function this = Tracker_FG_6(v2w)
            this.v2w = single(v2w);
        end 
        
        function this = setData(this, f)
            this.f = single(f);
        end
        
        function this = setParameters(this, stepSize, threshold, maxAngle, lengthRange) 
            if stepSize > 0
                this.stepSize = single(stepSize);
            else
                error('stepSize must be > 0');
            end
            if threshold >= 0
                this.threshold = single(threshold);
            else
                error('threshold must be >= 0');
            end
            if maxAngle > 0 && maxAngle <= 90
                this.maxAngle = single(maxAngle);
            else
                error('maxAngle must be in [0 90]');
            end
            if lengthRange(1) >= stepSize
                if (lengthRange(1) <= lengthRange(2))
                    this.lengthRange = single(lengthRange);
                else
                    error('lengthRange(2) must be >= lengthRange(1)');
                end
            else
                error('lengthRange(1) must be >= stepSize');
            end
        end
        
        function this = setProgressbar(this, val)
            this.pb = val;
        end
        
        function [tract, tractVal,tractCSD_FOD,tractDir,averageDir,stopVal,stopAngle] = track(this, point)
            point = single(point);
            if (isempty(this.stepSize) || isempty(this.threshold) || isempty(this.maxAngle) || isempty (this.lengthRange))
                error('set tracking parameters first, using setParameters');
            end
            if isempty(this.f)
                error('set input data first, using setData');
            end
            
            % interpolate
            this.interpolate(point);
            
            % mask out NaN values
            mask_ful = ~isnan(this.fun(1,:));
            point = point(:,mask_ful);
            this.fun = this.fun(:,mask_ful);
            
            % process function
            this.process();
            old_point = point; %
            % determine initial track direction(s)
            [point, dir, val, idk] = this.getInitDir(point);
            
            % add points_neighbour to mask, FG
               % here insert nb function, fguo
%                 [nsPoint, dir, val,idk] = getInitDir(this, sPoint); 
                udk = unique(idk);
                 [nb,nbv,ap,ip] = include_nb(old_point,udk,point);                  

            point_n6 = this.Get_neighbour(point);
            ibb = find(nb(:,ip));
            [ia,ib] = ismember(nbv(:,(ibb-1)*2+1)',point','rows');
            db = dir(:,ib);            
%             % add part too extract CSD_FOD, this.fun histogram
%             nfac = this.fun(1,:);
%             voxel_th = 0.1*nfac;
                
            % mask out all small peaks
            mask = val > this.threshold;
%             pointb = this.Get_neighbour(ap);
            
            npoint = point;
%             % add, mask &voxel_th              
%             mask_th = val > voxel_th;
%             mask = mask & mask_th;
                
            point = point(:,mask);
            dir = dir(:,mask);
            point = [point ap];
            dir = [dir db];
            val = [val val(:,ib)];
            
            % repeat data for probabilistic tractography
            point = repmat(point,[1 this.iterations]);
            dir = repmat(dir,[1 this.iterations]);
            val = repmat(val,[1 this.iterations]);
            
            % Track in both directions
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in first direction'); end;
            [tract1, tractVal1,tractCSD_FOD1,tractDir1,averageDir1] = this.trackOneDir(point, dir); %add voxel_th
            if this.pb; progressbar('ready'); end;
            
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in second direction'); end;
            [tract2, tractVal2,tractCSD_FOD2,tractDir2,averageDir2] = this.trackOneDir(point, -dir); %add voxel_th
            if this.pb; progressbar('ready'); end;
            
            % join tracts from both directions
            tract2 = cellfun(@flipud,tract2,'UniformOutput',false);
            tractVal2 = cellfun(@flipud,tractVal2,'UniformOutput',false);
            tract = cell(size(tract2));
            tractVal = cell(size(tractVal2));
            for j = 1:size(tract2,2)
                if ~isempty(tract2{j})
                    tract{j} = [tract2{j}; point(:,j)'];
                    tractVal{j} = [tractVal2{j}; val(:,j)'];
                    if ~isempty(tract1{j})
                        tract{j} = [tract{j}; tract1{j}];
                        tractVal{j} = [tractVal{j}; tractVal1{j}];
                    end
                else
                    if ~isempty(tract1{j})
                        tract{j} = [point(:,j)'; tract1{j}];
                        tractVal{j} = [val(:,j)'; tractVal1{j}];
                    end
                end
            end
            
            % enforce length limitations
            %             maska = cellfun('size',tract,1)*this.stepSize >= this.lengthRange(1);
            %             maskb = cellfun('size',tract,1)*this.stepSize <= this.lengthRange(2);
            %             mask = maska & maskb;
            %             tract = tract(mask);
            %             tractVal = tractVal(mask);
            mask = cellfun('size',tract,1)>1;
            tract = tract(mask);
            tractVal = tractVal(mask);
        end
    end
    
    methods (Access = private)
        function [tract, tractVal,tractCSD_FOD,tractDir,averageDir] = trackOneDir(this, point, dir)
            tract = cell(1,size(point,2));
            tractVal = cell(1,size(point,2));
            flist = 1:size(point,2);
            %
            tractCSD_FOD = cell(1,size(point,2)); % 45*it, it = tractLength & iteration
            tractDir = cell(1,size(point,2)); % 2*it, theta and phi, c2s(dir)
            averageDir = single(zeros(1,size(point,2))); % theta and phi, the dynamic angular threshold, 
%             stopVal = single(zeros(1,size(point,2)));
%             stopAngle = single(zeros(1,size(point,2)));
            rz = [0 0 1]; % compare dir to z-direction
            sPoint = point; 
            
            for it = 1:(this.lengthRange(2)/this.stepSize)
                if this.pb; progressbar(size(point,2),int2str(size(point,2))); end;
                % advance streamline

                
                
                point = point + this.stepSize .* dir; 
                theta_dir = acosd(sqrt(rz * dir)); %  3*TractNr   
                theta_dir = real(theta_dir);
                % interpolate
                this.interpolate(point);
              
                % mask out NaN values
                mask = ~isnan(this.fun(1,:)); 
                temp = mask;
                point = point(:,mask);
                dir = dir(:,mask);
                this.fun = this.fun(:,mask);
                flist = flist(mask);
                theta_dir = theta_dir(:,mask);
            
                % process function
                this.process();
                 % get new direction
                [newDir, val, angle] = this.getDir(dir);
               
                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*newDir,1));               
                % update dir
                newDir = flipsign([1 1 1],:).*newDir;
                
                

                
                % mask out small peaks
                mask = val > this.threshold;
%                 mask = mask | mask_nb1 | mask_nb2;
                opoint = point; % here insert nb function, fguo,20200
                temp = mask;
                fod = this.f;
                th = this.threshold.*10;
                stepsize = this.stepSize;
                [mask,ndir,nval,newangle,msk] = include_nbtrack3 (point,mask,dir,fod,th,stepsize);
                
%                 mask(ip) = true;
                % continue
                point = point(:,mask);            

                %                
%                 temp = mask;  
%                 point_maskP = temp & ~mask; 
%                 temp = temp(mask);
                theta_dir = theta_dir(:,mask);
%                 stopVal(flist(point_maskP)) = val(point_maskP); % 
                
                dir = dir(:,mask);
                newDir(:,msk) = ndir(:,msk);newDir = newDir(:,mask); 
                flist = flist(mask);
                angle(msk) = newangle(msk); angle = angle(mask); 
                val(msk) = nval(msk);val = val(:,mask); 
                
%                 dir(:,ip) = newd(:,ip);

                % mask out large angles
                
                mask = angle < this.maxAngle; %angle, fguo20200                            
%                 point_maskA = temp & ~mask;
%                 stopAngle(flist(point_maskA)) = angle(point_maskA); %angle, fguo20200
                theta_dir = theta_dir(:,mask);
                
                point = point(:,mask);
                dir = dir(:,mask);
                newDir = newDir(:,mask);
                flist = flist(mask);
                val = val(:,mask);


                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*newDir,1));
                
                % update dir
                dir = flipsign([1 1 1],:).*newDir;
                
                % stop if we are out of points
                if isempty(point)
                    break
                end

                % add points to the tracts
                 for i=1:length(flist)
                    tract{flist(i)}(it,:) = point(:,i);
                    tractVal{flist(i)}(it,:) = val(:,i);
                    tractDir{flist(i)}(it,:) = theta_dir(:,i);
                    averageDir(:,i) = nanmean(tractDir{flist(i)}); % this can be used earlier as an additional angular threshold
                    tractCSD_FOD{flist(i)}(it,:) = this.fun(:,i);
                end

            end
            % what to add between every step and whole trackOneDir function
        end
        
        function interpolate(this, point)
            point(4,:) = 1; voxel = this.v2w\point; voxel = voxel(1:3,:);
            this.fun = mex_interp(this.f, voxel);
        end
        %
        function point_n6 = Get_neighbour (this,point)
                point_n6 = cell(1,6);
                temp = point;
                temp(1,:) = point(1,:)-1;
                point_n6{1,1} = temp;
                
                temp = point;
                temp(1,:) = point(1,:)+1;
                point_n6{1,2} = temp;
                
                temp = point;
                temp(2,:) = point(2,:)-1;
                point_n6{1,3} = temp;
                
                temp = point;
                temp(2,:) = point(2,:)+1;
                point_n6{1,4} = temp;
                
                temp = point;
                temp(3,:) = point(3,:)-1;
                point_n6{1,5} = temp;
                
                temp = point;
                temp(3,:) = point(3,:)+1;
                point_n6{1,6} = temp; clear temp;
        end
    end
    methods (Access = public, Abstract = true) %protected
        process(this);       
        [point, dir, val,idk] = getInitDir(this, point);
        [dir, val, angle] = getDir(this, prevDir);
    end
end