function c = LM_convfft(x,y,type,maxMem)
if nargin < 3 % 'valid' convolution by default
    type = 'valid';
end
if nargin < 4
    maxMem = 15;
end

nx = size(x,1);
ny = size(y,1);
nFFT = pow2(nextpow2(nx+ny-1));

sy = size(y);

k = 16 * nFFT / 1e9; % 8 bytes * 2 (complex) * nFFT in Gb
ndims_y = prod(sy(2:end));
sc = 8 *  nFFT * ndims_y / 1e9; % to store the result
nd = floor((maxMem-sc)/k);

if maxMem < sc
    error('Will exceed max memmory');
end

if 0 < maxMem && nd < ndims_y
    
    assert(isvector(x) || size(x,2) == size(y,2));    
  
    % we'll need to work by batch
    y = y(:,:,:);

    c = nan([nFFT,size(y,[2,3])],'double');
    X = fft(x,nFFT);
    
    if ismatrix(y) % ndims(y) == 2
        nBatch = ceil(ndims_y / nd);
        bSize = distribute(ndims_y,nBatch);
        
        ie = 0;
        
        for iBatch = 1:nBatch
            ib = ie + 1;
            ie = ib + bSize(iBatch) - 1;
            
            c(:,ib:ie) = ifft( X .* fft(y(:,ib:ie),nFFT),nFFT,'symmetric');
        end
    else
        nBatch = ceil(size(x,2)*size(y,3) / nd);
        bSize = distribute(size(y,3),nBatch);
        
        ie = 0;
        
        for iBatch = 1:nBatch
            ib = ie + 1;
            ie = ib + bSize(iBatch) - 1;
            
            c(:,:,ib:ie) = ifft( X .* fft(y(:,:,ib:ie),nFFT),nFFT,'symmetric');
        end
    end
else % we can do everything at once
    c = ifft( fft(x,nFFT) .* fft(y,nFFT),nFFT,'symmetric');
end

m = min(nx,ny);
M = max(nx,ny);

switch type
    case 'valid'
        c = c( (1:(M-m+1)) + m-1,: );
    case 'full'
        c = c( (1:(M+m-1)),: ); 
    case 'same'
        c = c( (1:nx)+floor(ny/2),: );
end

c = reshape( c,[size(c,1),sy(2:end)] );

end