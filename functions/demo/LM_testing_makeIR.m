function [impResponse,tResponse,nMaxResponse] = LM_testing_makeIR(type,Fs)

tMaxResponse = 400e-3; % in s
nMaxResponse = ceil(tMaxResponse*Fs)+1;
tResponse = (-nMaxResponse:nMaxResponse)/Fs;


switch type
    case 1
        % simple Gaussian response
        delayResponse = 200e-3; % in s
        widthResponse = 50e-3; % in s
        
        impResponse = normpdf(tResponse,delayResponse,widthResponse);
        
    case 2
        % bi-modal response:
        A1 = 1;
        delayResponse_1 = 80e-3; % in s
        widthResponse_1 = 20e-3; % in s
        A2 = -1.5;
        delayResponse_2 = 180e-3; % in s
        widthResponse_2 = 40e-3; % in s
        
        n1 = normpdf(delayResponse_1,delayResponse_1,widthResponse_1);
        n2 = normpdf(delayResponse_2,delayResponse_2,widthResponse_2);
        
        impResponse = A1 * normpdf(tResponse,delayResponse_1,widthResponse_1) / n1 + ...
            A2 * normpdf(tResponse,delayResponse_2,widthResponse_2) / n2;
end

impResponse = impResponse(:);
impResponse = impResponse ./ max(abs(impResponse),[],1);

end