function resp_resid = mTRFpartial(stim,resp_full,fs,map,tmin,tmax,lambda,ver)

resp_resid = cell(size(resp_full));
switch ver
    case 1
        %% Original version
        for tt=1:length(resp_full)
        
            r_max=-1;
            for ll=1:length(lambda)
        
                model=mTRFtrain(stim{tt},resp_full{tt},fs,map,tmin,tmax,lambda(ll),'verbose',0);
                [pred,stats]=mTRFpredict(stim{tt},resp_full{tt},model,'verbose',0);
        
                if mean(stats.r)>r_max
                    r_max=mean(stats.r);
                    best_pred=pred;
                end
            end
        
            resp_resid{tt}=resp_full{tt}-best_pred;
        
        end
        
    case 2
        
        %% New version, individual trial models
        stats = mTRFcrossval(stim,resp_full,fs,map,tmin,tmax,lambda,'verbose',0);
        [~,I] = max(mean(mean(stats.r(:,:,:),1),3));
        
        for tt = 1:size(resp_full,2)
            model = mTRFtrain(stim{tt},resp_full{tt},fs,map,tmin,tmax,lambda(I),'verbose',0);
            resp_resid{tt} = resp_full{tt}-mTRFpredict(stim{tt},resp_full{tt},model,'verbose',0);
        end
        
    case 3
         %% New version, average trials model
        stats = mTRFcrossval(stim,resp_full,fs,map,tmin,tmax,lambda,'verbose',0);
        [~,I] = max(mean(mean(stats.r(:,:,:),1),3));

        for tt = 1:size(resp_full,2)
            model(tt) = mTRFtrain(stim{tt},resp_full{tt},fs,map,tmin,tmax,lambda(I),'verbose',0);
        end
        model = mTRFmodelAvg(model);
%       model = averageModels(model);
%         
        for tt = 1:size(resp_full,2)
            resp_resid{tt} = resp_full{tt}-mTRFpredict(stim{tt},resp_full{tt},model,'verbose',0);
        end

end
% 
% function model = averageModels(models)
%     model = models(1);
%     if length(models)>1
%         for n = 2:length(models)
%             model.w = model.w + models(n).w;
%             model.b = model.b + models(n).b;
%         end
%         model.w = model.w/n;
%         model.b = model.b/n;
%     end
% end

end