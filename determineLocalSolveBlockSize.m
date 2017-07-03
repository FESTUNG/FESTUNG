function localSolveBlockSize = determineLocalSolveBlockSize( K )
    assert( K > 0, 'The mesh must have at least one element');
    localSolveBlockSize = min(K, 16);
    if (mod(K, localSolveBlockSize) ~= 0)
      
%       if ( K < 16 )
%         localSolveBlockSize = K;
%       else
%         localSolveBlockSize = K / floor(K/16-0.5);
%       end
        
        localSolveBlockSize = 1;
        warning('Block size does not fit the problem. Reset block size to 1');
    end
    
%     if (mod(K, localSolveBlockSize) ~= 0)
%         localSolveBlockSize = 1;
%         warning('Block size does not fit the problem. Reset block size to 1');    
%     end
    
    assert( mod(K, localSolveBlockSize) == 0, ...
        'Block size does not fit to problem size!');
end