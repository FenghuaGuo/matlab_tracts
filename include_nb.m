function [nb,nbv,ap,ip] = include_nb (oldpoint,idk,point)    

% % Get seedpoint
%                  p = Get_HBT_par_CSD;
%                  mask_s = ~isnan(this.f(:,:,:,1));
%                  SR = p.SeedPointRes;
%                  VDims = [1 1 1];
%                  sPoint = E_DTI_Get_Seeds_WBT(mask_s,SR,VDims,p);
%                  smask = ~isnan(this.fun(1,:));
%                  sPoint = sPoint(:,smask);

%                  [nsPoint, dir, val,idk] = getInitDir(this, sPoint);
%                 idk = unique(idk);
                tem = ones(1,size(oldpoint,2));
                tem(:,idk) = 0;
                nPoint = oldpoint(:,tem>0);
                temp = nPoint;
%                 temp = sPoint(:,~mask);
                point_n6 = Get_neighbour(temp);
%                 nb_mat = cell2mat(point_n6);
%                 nb_mat = reshape(nb_mat,[3 size(temp,2) 6]);
%                 nb_mat = permute(nb_mat,[1 3 2]);
                nb_mat = zeros(3,6,size(temp,2));
                for ii =1:6
                    nb_mat(:,ii,:) = point_n6{1,ii};
                    
                    
                end
                
%                 nb_temp = nb_mat(:,:,~mask);
               ot = zeros(6,size(temp,2));
                for ii = 1:size(temp,2)
                    
                    v1 = squeeze(nb_mat(:,:,ii));
                    m1 = ismember(v1',point','rows');
                    ot(:,ii) = m1;
                end
                
                nb = zeros(3,size(temp,2));
                nb(1,:) = ot(1,:)&ot(2,:);
                nb(2,:) = ot(3,:)&ot(4,:);
                nb(3,:) = ot(5,:)&ot(6,:);
                ap = [];
                ip = [];
                for ii = 1:size(temp,2)
                    if sum(nb(:,ii))>0
                        ap = [ap temp(:,ii)];
                        ip = [ip ii];
                    end
                end
                nbv = nb_mat(:,:,ip);
end