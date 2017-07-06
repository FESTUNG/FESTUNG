function localSolveBlockSize = determineLocalSolveBlockSize( K, trueLocalSolveSize )
    assert( K > 0, 'The mesh must have at least one element');
    
    if nargin == 1
      trueLocalSolveSize = 16;
    end
    
    localSolveBlockSize = min(K, trueLocalSolveSize);
    if ( trueLocalSolveSize > K )
        warning('Block size does not fit the problem. Reset block size to K');      
    end
    
    if (mod(K, localSolveBlockSize) ~= 0 )
        localSolveBlockSize = 1;
        warning('Block size does not fit the problem. Reset block size to 1');
    end
    
    assert( mod(K, localSolveBlockSize) == 0, ...
        'Block size does not fit to problem size!');
end