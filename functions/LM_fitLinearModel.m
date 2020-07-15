function model = LM_fitLinearModel(X,y,opt)
%
% LM_fitLinearModel
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement linear regression y = X * B with regularisation
%
% ------ sanity checks
assert( ~isempty(y) && ~isempty(X) );

switch opt.method.name
    
    case 'ridge-eig-XtX'
        tBegin=tic;
        % !!! Assuming N > p here !
        if opt.printOut
            fprintf('%s started...',opt.method.name);
        end
        
        [~,p] = size(X);
        [~,q] = size(y);
        
        % eigenvalue decomposition of X' * X & dropping numerical zeros 
        [V,S] = LM_regEigen(X, opt.accumulate, opt.method.removeEig);
        nL = mean(S);
        
        durSVD = toc(tBegin);
        tRest = tic;
        
        %% ------ ridge ------    
        % normalising lambdas
        lambda = nL * opt.method.lambda;
        
        if opt.accumulate
            % in this case y is actually X'*y precomputed!
            z = V' * y;
        else
            z = V' * (X' * y); % KEEP the parenthesis!
        end
        
        nLambda = length(lambda);
        coeffs = zeros(p,q,nLambda);
        
        for iL = 1:nLambda
            coeffs(:,:,iL) = V * ( z ./ (S + lambda(iL)) );
        end
        
        durTot = toc(tBegin);
        durRest = toc(tRest);
        
        % storing results
        model = struct();        
        model.meanEigenval = nL;
        model.normaliseLambdaType = 'lambda x mean(eigenval(XtX))';
        model.eigenVals = S;
        
        
        %%
        if opt.printOut
            fprintf('done.\nElapsed: total: %s, eig %s (%.1f %%), lambda loop %s (%.1f %%)\n',...
                seconds2HMS(durTot,{'h','mn','s'}),...
                seconds2HMS(durSVD,{'h','mn','s'}),durSVD/durTot*100,...
                seconds2HMS(durRest,{'h','mn','s'}),durRest/durTot*100);
        end
end

model.coeffs = coeffs;

end
%
%