function model = LM_fitLinearModel(X,y,trainOpt)
%
%
%
% ------ sanity checks
assert( ~isempty(y) && ~isempty(X) );

switch trainOpt.method.name
    
    case 'ridge-eig-XtX'
        tBegin=tic;
        % !!! Assuming N > p here !
        if trainOpt.printOut
            fprintf('%s started...',trainOpt.method.name);
        end
        
        [~,p] = size(X);
        [~,q] = size(y);
        
        % This would yield the singular values of X
        % ( eigenvalues of X'*X = (singular values of X)^2 ),
        % but is very costly when X is large:
        % [~,Sdd,~] = svd(X,'econ');
        % Sdd = diag(Sdd);
               
        % ensuring X'*X = V * S * V'
        if trainOpt.accumulate
            % in this case X is actually X' * X precomputed!
            [V,S] = eig(X,'vector');
        else
            [V,S] = eig(X'*X,'vector');
        end
        % !!! eig does not sort the eigenvalues as svd !!!
        [S,iS] = sort(S,'descend');
        V = V(:,iS);
        
        durSVD = toc(tBegin);
        tRest = tic;
        
        
        %% ------ dropping degenerate dimenions ("numerical zeros") ------
        switch trainOpt.method.removeEig.type
            case 'tol' % based on tolerance
                tol = eps(S(1));
                r1 = sum(S > tol)+1;
            case 'var' % based on total variance
                r1 = find(cumsum(S)/sum(S)>=trainOpt.removeEig.var,1)+1;
            case 'cond' % base on log10( condition number )
                r1 = find(log10(S(1)./S) > trainOpt.removeEig.cond,1);
        end

        S(r1:end) = [];
        V(:,r1:end) = [];
        
        nL = mean(S);
        

        %% ------ ridge ------    
        % normalising lambdas
        lambda = nL * trainOpt.method.lambda;
        
        if trainOpt.accumulate
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
        if trainOpt.printOut
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
%